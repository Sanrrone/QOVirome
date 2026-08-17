#!/usr/bin/env bash
set -euo pipefail

# Host-side launcher for the Gutbusters virome-mining container (Mars-integrated).
# Mirrors qo-analysis-bgc/run.sh: download the assembly from S3, run the SIF with
# the DBs bind-mounted, upload results+logs, post the analysis webhook.

# --- Required variables from SQS / Mars worker ---
: "${RUN_ID:?}" "${ANALYSIS_WEBHOOK_URL:?}" "${S3_OUTPUT_DIR:?}" \
  "${QO_WORKDIR:?}" "${QO_ANALYSIS_SIF:?}"
# Default rather than require: an unset token under `set -u` would abort inside
# post_analysis_webhook — including when the ERR trap calls it — turning a clean
# failure into a run that never reports at all. Empty just yields a logged 401.
: "${INTERNAL_WORKER_AUTH_TOKEN:=}"
# Seeded so the input validations below can report before SAMPLE_ID is derived.
SAMPLE_ID="${SAMPLE_ID:-}"
S3_CONTIGS_URI="${S3_CONTIGS_URI:-${S3_INPUT_URI:-}}"

THREADS="${THREADS:-8}"
GB_MIN_LEN="${GB_MIN_LEN:-2000}"
GB_MIN_SCORE="${GB_MIN_SCORE:-0.9}"
# geNomad MMseqs2 marker-search sensitivity. Its own help: higher annotates more
# proteins but is slower AND uses more memory. Kept at geNomad's default 4.2 rather
# than qo-analysis-mge's 3.5, because mge's note warns the recall loss "falls hardest
# on divergent gut phages" — which is what this analysis exists to find. Exposed as an
# env knob (no CLI flag on the entrypoint) so it can be retuned without a .sif rebuild.
GB_GENOMAD_SENS="${GB_GENOMAD_SENS:-4.2}"
# Stage-5 corroboration floor: the score at which a classifier arm counts as a
# supporting vote. Kept separate from GB_MIN_SCORE, which used to double as VirSorter2's
# --min-score (removed 2026-08-05). It fed the legacy fallback gate's confidence threshold too,
# and that gate was removed 2026-08-07; the floor now only tags corroboration counts in
# pbb_viruses.tsv, which nothing gates on.
GB_CORROB_FLOOR="${GB_CORROB_FLOOR:-0.5}"
# How strict the stage-5 admission gate is. This is the ONE knob meant to be user-facing:
# the platform passes a LEVEL, and 5_eval.sh resolves it to a cut on the 0-1 viral score.
# Calibrated on the within-library benchmark (6,024 WGS contigs, 8 libraries, 33.3% viral in
# every one), which is the production-shaped question -- one mixed metagenome, decide per
# contig:
#     sensitive  0.15    46.5% precision, 42.5% recall, keeps 30.5% of candidates
#     balanced   0.267   51.1% precision, 26.8% recall, keeps 17.5%     <- default
#     strict     0.55    60.1% precision, 14.5% recall, keeps  8.1%
# The three sets are strictly NESTED, so "stricter" always means "a subset", which is what a
# UI control implies. Do not surface the VLP-vs-WGS precisions (91/93/95%) to a user -- that
# benchmark has one library per class, so its precision partly measures which prep was run.
GB_STRICTNESS="${GB_STRICTNESS:-balanced}"
# Escape hatch for reproducibility: an explicit cut wins over the level, so a completed run
# can be replayed exactly even if the levels are later recalibrated. Empty = use the level.
GB_VSCORE_MIN="${GB_VSCORE_MIN:-}"
ANALYSIS_API_BASE_URL="${ANALYSIS_WEBHOOK_URL%/api/webhooks/analysis}"

# DBs were moved to /data/containers/dbs and are bind-mounted read-only at
# /opt/gutbusters/db (the path the pipeline expects). Override with QO_VIROME_DB_ROOT.
QO_VIROME_DB_ROOT="${QO_VIROME_DB_ROOT:-/data/containers/dbs}"

# geNomad's DB is SHARED with qo-analysis-mge rather than copied — /data/containers
# is ~95% full and a second copy would be 1.4 GB for no benefit. Mars already
# exports QO_GENOMAD_DB_PATH for the mge analysis and passes its whole env to every
# container, so this resolves with no worker change.
QO_GENOMAD_DB_PATH="${QO_GENOMAD_DB_PATH:-/data/containers/qo-analysis-mge/refs}"
# `genomad download-database <dir>` writes a nested <dir>/genomad_db/; accept either
# form, exactly as qo-analysis-mge/run.sh does.
GENOMAD_DB_DIR="$QO_GENOMAD_DB_PATH"
[ -d "${GENOMAD_DB_DIR}/genomad_db" ] && GENOMAD_DB_DIR="${GENOMAD_DB_DIR}/genomad_db"

# Short job-specific TMPDIR to stay under the 108-char AF_UNIX socket-path limit.
_JOB_SHORT="${SLURM_JOB_ID:-$(printf '%s' "$RUN_ID" | md5sum | cut -c1-8)}"
TMPDIR="/tmp/qo-${_JOB_SHORT}"
TMP="$TMPDIR"
TEMP="$TMPDIR"

WORKDIR="${QO_WORKDIR}"
INPUT_DIR="${WORKDIR}/input"
OUTDIR="${WORKDIR}/output"
RESULTS_DIR="${OUTDIR}/results"
HOME_DIR="${WORKDIR}/home"
CACHE_DIR="${WORKDIR}/cache"
CONFIG_DIR="${WORKDIR}/config"
# Each Mars retry reuses the same per-job WORKDIR; clear stale output first.
rm -rf "$OUTDIR"
mkdir -p "$INPUT_DIR" "$RESULTS_DIR" "$TMPDIR" "$HOME_DIR" "$CACHE_DIR" "$CONFIG_DIR"

export TMPDIR TMP TEMP
export HOME="$HOME_DIR"
# Forward writable HOME/caches into the container (--cleanenv strips the host env).
export APPTAINERENV_TMPDIR="$TMPDIR"
export APPTAINERENV_TMP="$TMP"
export APPTAINERENV_TEMP="$TEMP"
export APPTAINERENV_HOME="$HOME_DIR"
export APPTAINERENV_XDG_CACHE_HOME="$CACHE_DIR"
export APPTAINERENV_XDG_CONFIG_HOME="$CONFIG_DIR"
# Without this the pipelines' ${GB_GENOMAD_SENS:-4.2} could never see anything but the
# default — --cleanenv drops the host env.
export APPTAINERENV_GB_GENOMAD_SENS="$GB_GENOMAD_SENS"
# Same reason: without the prefix 5_eval.sh could only ever see its own hard-coded
# default. (APPTAINERENV_ works for this name where it would NOT for GB_PHS, which
# entrypoint/gutbusters unconditionally re-exports from --score.)
export APPTAINERENV_GB_CORROB_FLOOR="$GB_CORROB_FLOOR"
export APPTAINERENV_GB_STRICTNESS="$GB_STRICTNESS"
# Same reason again: --cleanenv drops the host env, so without the prefix GB_KRAKEN=0 would be
# accepted here (skipping the DB check) and then silently ignored INSIDE the container, which is
# the worst of both -- kraken2 would run against a database run.sh never verified.
# Defaulted HERE rather than at the DB check below, because the export must not run first and
# read an unset variable.
GB_KRAKEN="${GB_KRAKEN:-1}"
case "$GB_KRAKEN" in
  0|1) ;;
  *) echo "[run.sh] invalid GB_KRAKEN='$GB_KRAKEN' (want 0 or 1)" >&2; exit 1;;
esac
export APPTAINERENV_GB_KRAKEN="$GB_KRAKEN"
# Only export the explicit cut when one was actually given, so the level stays in charge by
# default. Exporting it empty would in fact be harmless -- 5_eval.sh reads it with
# ${GB_VSCORE_MIN:-$_vsmin_default}, and :- treats unset and empty alike -- but leaving the
# variable absent keeps `env | grep GB_` honest about which knob decided the run.
# Written as an if, not `[ -n .. ] && export ..`: under `set -e` a trailing AND-list whose
# test fails is a footgun that reads as if it should abort.
if [ -n "$GB_VSCORE_MIN" ]; then
    export APPTAINERENV_GB_VSCORE_MIN="$GB_VSCORE_MIN"
    echo "[run.sh] viral score cut PINNED to $GB_VSCORE_MIN (overrides strictness=$GB_STRICTNESS)"
else
    echo "[run.sh] viral strictness: $GB_STRICTNESS"
fi

post_analysis_webhook() {
    local status="$1"
    local error_message="${2:-}"
    local report_url="${3:-}"
    local metrics_json="${4:-}"
    local aggregates_json="${5:-}"

    local _payload
    _payload="$(python3 - "$RUN_ID" "$status" "$S3_CONTIGS_URI" \
        "$SAMPLE_ID" \
        "${S3_OUTPUT_DIR}/results/" "${S3_OUTPUT_DIR}/logs/" \
        "$report_url" "$error_message" "$metrics_json" "$aggregates_json" <<'PYEOF'
import json
import sys

run_id, status, s3_contigs, sample_id, s3_results, s3_logs, report_url, error_message, metrics_json, aggregates_json = sys.argv[1:]
payload = {
    "runId": run_id,
    "status": status,
    "step": "analysis",
    "input": {
        "s3Contigs": s3_contigs,
        "sampleId": sample_id,
    },
    "output": {
        "s3Results": s3_results,
        "s3Logs": s3_logs,
    },
}
if report_url:
    payload["output"]["reportUrl"] = report_url
if status == "failed" and error_message:
    payload["error"] = error_message
if metrics_json:
    try:
        payload["metrics"] = json.loads(metrics_json)
    except Exception:
        pass
if aggregates_json:
    try:
        payload["aggregates"] = json.loads(aggregates_json)
    except Exception:
        pass
print(json.dumps(payload))
PYEOF
)"

    local _body_file="/tmp/qo-webhook-${RUN_ID}.body"
    local _http
    _http=$(curl -sS -o "$_body_file" -w '%{http_code}' \
        -X POST "$ANALYSIS_WEBHOOK_URL" \
        -H "Content-Type: application/json" \
        -H "Authorization: Bearer $INTERNAL_WORKER_AUTH_TOKEN" \
        --data "$_payload" 2>&1) || _http="000"
    echo "[webhook] status=$status http=$_http body=$(head -c 500 "$_body_file" 2>/dev/null || echo '')"
    if [[ "$_http" != "200" && "$_http" != "204" ]]; then
        echo "[webhook] WARN: webhook did not return 2xx (got $_http)"
    fi
    rm -f "$_body_file"
}

# --- ERROR HANDLING TRAP ---
on_error() {
    echo "[ERROR] Virome pipeline failed! Executing safety trap..."
    if [ -d "${RESULTS_DIR}/logs" ]; then
        # `|| true` keeps the failed-status webhook reachable if log sync fails.
        aws s3 sync "${RESULTS_DIR}/logs/" "${S3_OUTPUT_DIR}/logs/" || true
    fi
    post_analysis_webhook "failed" "Virome pipeline crashed"
    exit 1
}
trap on_error ERR

# --- VALIDATE INPUT ---
# Deliberately below the trap/webhook definition: these used to `exit 1` before
# either existed, so a malformed input left the run stuck in 'running' until the
# zombie sweeper reaped it instead of failing cleanly.
if [ -z "$S3_CONTIGS_URI" ]; then
    echo "[ERROR] no assembly input: both S3_CONTIGS_URI and S3_INPUT_URI are unset" >&2
    post_analysis_webhook "failed" "No assembly input: both S3_CONTIGS_URI and S3_INPUT_URI are unset"
    exit 1
fi

_RAW_BASENAME="$(basename "${S3_CONTIGS_URI%%\?*}")"

# Derive a human-readable SAMPLE_ID from the filename (strip the user-prefix UUID
# and the assembly extension); fall back to RUN_ID.
if [ -z "$SAMPLE_ID" ]; then
    _SAMPLE_FROM_NAME="$(printf '%s' "$_RAW_BASENAME" \
        | sed 's/^[0-9a-f]\{8\}-[0-9a-f]\{4\}-[0-9a-f]\{4\}-[0-9a-f]\{4\}-[0-9a-f]\{12\}_//' \
        | sed -E 's/\.(fna|fasta|fa)(\.gz)?$//I')"
    SAMPLE_ID="${_SAMPLE_FROM_NAME:-$RUN_ID}"
fi

# Fail fast on input shape — assembly FASTA only (plain or gzipped). Matched
# case-insensitively: accepted_data_types includes user-uploaded genomes, and a
# ".FASTA" upload is a valid assembly.
case "$(printf '%s' "$_RAW_BASENAME" | tr '[:upper:]' '[:lower:]')" in
    *.fna|*.fna.gz|*.fasta|*.fasta.gz|*.fa|*.fa.gz) ;;
    *)
        echo "[ERROR] S3_CONTIGS_URI must end in .fna/.fasta/.fa (gzipped allowed). got: $_RAW_BASENAME" >&2
        post_analysis_webhook "failed" "Unsupported input file type: ${_RAW_BASENAME} (expected .fna/.fasta/.fa, gzipped allowed)"
        exit 1
        ;;
esac

# --- VERIFY BIND-MOUNTED DATABASES ---
# All DBs live under QO_VIROME_DB_ROOT and are bind-mounted read-only at
# /opt/gutbusters/db. Fail fast (with a webhook) before launching the SIF.
_missing=""
# virsorter2 dropped from this list 2026-08-05 with VirSorter2 itself; requiring a DB no stage
# reads would fail a perfectly good host.
# kraken2/hgdb_unphaged dropped 2026-08-07: HRGMv2 is now the default for BOTH kraken2 roles
# (host assignment and the stage-2 bacterial tag), so the image needs one 19 GB reference instead
# of two. Requiring a DB no stage reads would fail a perfectly good host — same reasoning that
# removed virsorter2 above. Setting GB_KRAKEN_DB back to hgdb_unphaged still works, it is just no
# longer checked here.
# phrog ADDED 2026-08-07: stack9 consumes 15 PHROGs columns. Without this check a missing DB is
# silent — phrog_features.py warns, every contig gets median-imputed, and the run still reports
# "scored=N" while quietly scoring a degraded model.
# kraken2 is CONDITIONAL: it feeds no admission feature (stack9 has zero kraken2-derived
# columns), so with GB_KRAKEN=0 the 19 GB reference is not needed at all and requiring it would
# block the very deployment that flag exists for. Required only when GB_KRAKEN=1.
# blast/uhgv (4.5 GB) and phabox_db_v2_1 (1.7 GB) dropped 2026-08-09 with the move to the 3-stage
# pipeline: the UHGV BLAST claimer and PhaMer both lived in the retired stages 2-4, and neither
# feeds a stack9 feature. Restoring GB_BLAST=1 requires putting blast/uhgv back on this line.
_dbs="checkv_db/checkv-db-v1.5 deepvirfinder/models metaphapred/model stack4 vcontrast phrog"
if [ "$GB_KRAKEN" = "1" ]; then
    _dbs="$_dbs kraken2/HRGMv2"
else
    echo "[run.sh] GB_KRAKEN=0 — kraken2 skipped; host prediction unavailable, HRGMv2 not required"
fi
for sub in $_dbs; do
    [ -d "${QO_VIROME_DB_ROOT}/${sub}" ] || _missing="${_missing} ${sub}"
done
if [ -n "$_missing" ]; then
    echo "[ERROR] virome DBs missing under QO_VIROME_DB_ROOT=${QO_VIROME_DB_ROOT}:${_missing}" >&2
    post_analysis_webhook "failed" "Virome databases missing under ${QO_VIROME_DB_ROOT}:${_missing}"
    exit 1
fi

# geNomad DB (shared with qo-analysis-mge) — checked separately since it lives
# outside the db root. Same fail-fast contract as the block above.
if [ ! -d "$GENOMAD_DB_DIR" ] || [ -z "$(ls -A "$GENOMAD_DB_DIR" 2>/dev/null)" ]; then
    echo "[ERROR] geNomad DB missing or empty at QO_GENOMAD_DB_PATH=${QO_GENOMAD_DB_PATH}" >&2
    post_analysis_webhook "failed" "geNomad database missing or empty at ${QO_GENOMAD_DB_PATH} — run qo-analysis-mge/setup_db.sh on the host first"
    exit 1
fi

# VirSorter2's Done_all_setup marker used to be touched here so `vs2 run` never diverted into
# its internal setup workflow. VirSorter2 was removed from the pipeline 2026-08-05 and no stage
# invokes it, so the marker has no consumer.

# --- DOWNLOAD INPUT ---
echo "[1/4] Downloading assembly from S3..."
aws s3 cp "$S3_CONTIGS_URI" "$INPUT_DIR/"
_CANONICAL_BASENAME="$(printf '%s' "$_RAW_BASENAME" | sed 's/^[0-9a-f]\{8\}-[0-9a-f]\{4\}-[0-9a-f]\{4\}-[0-9a-f]\{4\}-[0-9a-f]\{12\}_//')"
if [ "$_CANONICAL_BASENAME" != "$_RAW_BASENAME" ]; then
    mv "${INPUT_DIR}/${_RAW_BASENAME}" "${INPUT_DIR}/${_CANONICAL_BASENAME}"
fi
LOCAL_INPUT="${INPUT_DIR}/${_CANONICAL_BASENAME}"
# Case-insensitive to match the extension whitelist above, which now accepts ".GZ".
case "$(printf '%s' "$LOCAL_INPUT" | tr '[:upper:]' '[:lower:]')" in
    *.gz) gunzip -f "$LOCAL_INPUT"; LOCAL_INPUT="${LOCAL_INPUT%.*}" ;;
esac

# --- LIVE PIPELINE CODE ---
# The .sif bakes pipelines/, scripts/ and the entrypoint at build time (Apptainer.def %files), so
# editing a pipeline step normally costs a full 2.7 GB rebuild. Bind the host-side copies over
# those paths instead: the image keeps supplying the conda envs and tool binaries (the expensive,
# rarely-changing half) while the orchestration shell/awk is read live from disk and can be
# changed with an scp.
#
# Guarded, and read-only. Binding a directory SHADOWS the baked one entirely, so a partial or
# missing host copy would silently run a crippled pipeline rather than fail — hence the explicit
# sentinel check on each path before it is bound, and the fall-back to the baked copy (with a
# loud warning) if the check does not pass. Never bind a directory that is missing a step.
#
# NOTE: this makes pipelines/ and scripts/ HOST-SIDE state, exactly like run.sh. Mars and Sun both
# claim from the same queue, so any edit here must be deployed to BOTH hosts or the same assembly
# will be analysed differently depending on which host wins the job.
PIPE_BINDS=()
_pipe_src="${SCRIPT_DIR_EARLY:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# The sentinels must name the stages that ACTUALLY RUN, which since 2026-08-09 are 1_prophage,
# 4_clean and 5_eval. They used to name 2_bactfilter/3_pullphages/4_pullviruses; those are retired
# but were kept on disk non-executable, so the old check still passed by accident. Deleting them --
# the obvious next cleanup -- would have failed this check and silently fallen back to the 5-stage
# pipeline baked into the .sif, which is exactly the failure the comment above warns about.
if [ -f "${_pipe_src}/pipelines/1_prophage.sh" ] && [ -f "${_pipe_src}/pipelines/4_clean.sh" ] \
   && [ -f "${_pipe_src}/pipelines/5_eval.sh" ]; then
    PIPE_BINDS+=(--bind "${_pipe_src}/pipelines:/opt/gutbusters/pipelines:ro")
    echo "[run.sh] pipelines/: LIVE from ${_pipe_src}/pipelines (no rebuild needed)"
else
    echo "[run.sh] WARN: ${_pipe_src}/pipelines incomplete — falling back to the copy baked in the .sif" >&2
fi
if [ -f "${_pipe_src}/scripts/pbbcalc.awk" ] && [ -f "${_pipe_src}/scripts/parsekr.awk" ]; then
    PIPE_BINDS+=(--bind "${_pipe_src}/scripts:/opt/gutbusters/scripts:ro")
    echo "[run.sh] scripts/: LIVE from ${_pipe_src}/scripts (no rebuild needed)"
else
    echo "[run.sh] WARN: ${_pipe_src}/scripts incomplete — falling back to the copy baked in the .sif" >&2
fi
if [ -f "${_pipe_src}/entrypoint/gutbusters" ]; then
    PIPE_BINDS+=(--bind "${_pipe_src}/entrypoint/gutbusters:/opt/gutbusters/entrypoint/gutbusters:ro")
    echo "[run.sh] entrypoint: LIVE from ${_pipe_src}/entrypoint/gutbusters (no rebuild needed)"
fi

# --- RUN CONTAINER ---
echo "[2/4] Running Gutbusters virome mining (threads=${THREADS})..."
apptainer run \
    --cleanenv \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "${QO_VIROME_DB_ROOT}:/opt/gutbusters/db:ro" \
    --bind "${GENOMAD_DB_DIR}:/opt/gutbusters/genomad_db:ro" \
    "${PIPE_BINDS[@]}" \
    "${QO_ANALYSIS_SIF}" \
    --in "$LOCAL_INPUT" \
    --outdir "$RESULTS_DIR" \
    --threads "$THREADS" \
    --min-len "$GB_MIN_LEN" \
    --score "$GB_MIN_SCORE" \
    --prefix "$SAMPLE_ID" \
    --keeptmp

VIRUSES_FASTA="${RESULTS_DIR}/${SAMPLE_ID}_viruses.fna"
VIRUSES_TABLE="${RESULTS_DIR}/${SAMPLE_ID}_table.tsv"
if [ ! -f "$VIRUSES_FASTA" ] || [ ! -f "$VIRUSES_TABLE" ]; then
    echo "[ERROR] expected outputs not found: $VIRUSES_FASTA / $VIRUSES_TABLE" >&2
    aws s3 sync "${RESULTS_DIR}/logs/" "${S3_OUTPUT_DIR}/logs/" || true
    post_analysis_webhook "failed" "Contract not fulfilled: missing ${SAMPLE_ID}_viruses.fna or ${SAMPLE_ID}_table.tsv"
    exit 1
fi

# --- UPLOAD ---
echo "[3/4] Uploading results and logs to S3..."
aws s3 cp "$VIRUSES_FASTA" "${S3_OUTPUT_DIR}/results/${SAMPLE_ID}_viruses.fna"
aws s3 cp "$VIRUSES_TABLE" "${S3_OUTPUT_DIR}/results/${SAMPLE_ID}_table.tsv"
[ -d "${RESULTS_DIR}/logs" ] && aws s3 sync "${RESULTS_DIR}/logs/" "${S3_OUTPUT_DIR}/logs/" || true

# Per-tool scores for every CANDIDATE contig, not just the admitted ones. The deliverable
# table carries only rows that PASSED the stage-5 gate, so without this the scores of
# rejected contigs are destroyed at the end of every run and any statement about what the
# gate does can only ever be a behaviour claim measured on the survivors. That is exactly
# why the corroboration gate's effect had to be simulated from stored deliverables rather
# than observed. `--keeptmp` above suppresses the entrypoint's own `rm -rf` of
# results/evaluation so these survive long enough to be copied out; the workdir is still
# deleted below, and since the entrypoint only cleans up AFTER every stage has run, this
# costs no peak disk — it just defers the final free by the length of this upload.
# Non-fatal on purpose: losing a diagnostic must never fail a run that produced its
# deliverables, but it is logged rather than swallowed by a bare `|| true`.
if [ -d "${RESULTS_DIR}/evaluation" ]; then
    aws s3 cp "${RESULTS_DIR}/evaluation/" "${S3_OUTPUT_DIR}/scores/" \
        --recursive --exclude '*' \
        --include '*_vscores.tsv' --include 'pbb_viruses.tsv' \
      || echo "[WARN] score-table upload failed (non-fatal; deliverables already uploaded)"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- PREVIEW + COMPLETE ---
echo "[4/4] Generating analysis preview..."
_PREVIEW_ERR=$(mktemp)
PREVIEW_DATA=$(python3 "${SCRIPT_DIR}/extract_preview.py" "${RESULTS_DIR}" 2>"$_PREVIEW_ERR") || {
    _contract_err=$(cat "$_PREVIEW_ERR"); rm -f "$_PREVIEW_ERR"
    aws s3 sync "${RESULTS_DIR}/logs/" "${S3_OUTPUT_DIR}/logs/" || true
    post_analysis_webhook "failed" "Contract not fulfilled: ${_contract_err}"
    exit 1
}
rm -f "$_PREVIEW_ERR"

post_analysis_webhook "completed" "" "" "${PREVIEW_DATA}" ""

trap - ERR
rm -rf "$WORKDIR" "$TMPDIR"
echo "Virome pipeline finished successfully."
