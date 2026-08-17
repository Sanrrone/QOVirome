#!/usr/bin/env bash
 
set -euo pipefail
unset PYTHONPATH 
export PYTHONNOUSERSITE=1
              
# ---------------- Env & inputs ----------------
: "${GB_OUTDIR:?Need GB_OUTDIR}"
: "${GB_THREADS:=4}"
: "${GB_DB_ROOT:?Need GB_DB_ROOT}"
: "${GB_IN:?Need GB_IN}"
: "${GB_MLEN:?Need GB_MLEN}"
: "${GB_HOME:?Need GUTBUSTERS_HOME}"
: "${GB_PHS:?Need PHAGE_SCORE}"
: "${GB_SCRIPTS:?Need GB_SCRIPTS}"
: "${GB_PREFIX:?Need GB_PREFIX}"


export TMPDIR=$GB_TMP
export GUTBUSTERS_TMP="$TMPDIR"
export TRANSFORMERS_CACHE="$GUTBUSTERS_TMP/hf"
export HF_HOME="$GUTBUSTERS_TMP/hf"
export TORCH_HOME="$GUTBUSTERS_TMP/torch"
export XDG_CACHE_HOME="$GUTBUSTERS_TMP/xdg"
export MPLCONFIGDIR="$GUTBUSTERS_TMP/matplotlib"
export SNAKEMAKE_SCHEDULER=greedy
#export PATH="/opt/micromamba/envs/vs2/bin:$PATH"
# ONE kraken2 DB (2026-08-07): HRGMv2 replaces hgdb_unphaged here too, so the image needs a
# single 19 GB reference instead of two (22 GB freed per host). hgdb_unphaged existed to be a
# PHAGE-SUBTRACTED bacterial reference for the old stage-2 veto; that veto became TAG-ONLY on
# 2026-08-05, so the subtraction no longer buys anything -- and it never worked: incomplete UHGV
# subtraction left residual prophage, so hgdb_unphaged mislabels 85.63% of TRUE VIRAL CHGV
# contigs as bacterial (tag AUC 0.5698, barely above chance). MEASURED head-to-head 2026-08-07:
#   false bacterial tag on true viral   CHGV 85.63% -> 57.54%   curated 4.17% -> 1.43%
#   tag-flag AUC                        CHGV 0.5698 -> 0.6448   (production-relevant)
#   k2_bact_frac AUC                    CHGV 0.6357 -> 0.6199   curated 0.8848 -> 0.9046
# HRGMv2 tags fewer true bacteria (86.51% vs 99.59%) but mislabels far fewer viruses; net
# discrimination improves. k2_bact_frac is worth +0.0005/+0.0008 on stack9 either way, i.e. the
# DB choice barely touches the model -- this is an ANNOTATION decision.
# ⚠️ HRGMv2 is GTDB-derived: rank "domain" not "superkingdom", names prefixed "d__". The
# taxonomy walk in 2_bactfilter.sh accepts both; with only the superkingdom test HRGMv2 tags
# ZERO contigs bacterial while still classifying 80.85% -- a SILENT no-op.
KRAKEN_DB="${GB_KRAKEN_DB:-${GB_DB_ROOT}/kraken2/HRGMv2}"
# Host calls use HRGMv2, bacterial tagging stays on hgdb_unphaged — see the comment block in
# 1_prophage.sh. This file MUST agree with stage 1: the two host tables are concatenated below
# into one deliverable column, and the two DBs use different nomenclature (GTDB vs NCBI), so
# splitting them would put "t__Megamonas funiformis" and "Megamonas" in the same column.
KRAKEN_HOST_DB="${GB_KRAKEN_HOST_DB:-${GB_DB_ROOT}/kraken2/HRGMv2}"

# ---------------- KRAKEN2 IS OPTIONAL (2026-08-07) ----------------
# GB_KRAKEN=1 (default) runs it; GB_KRAKEN=0 skips every kraken2 pass in the pipeline.
#
# It can be skipped because it feeds NOTHING that decides admission. stack9 contains zero
# kraken2-derived features (verified against the artifact feature list), and the only consumer
# that ever used one -- the 4-term fallback formula, via k2_bact_frac -- was removed on the same
# day. What is genuinely lost with GB_KRAKEN=0 is PROPHAGE HOST PREDICTION: the h_rank / h_name /
# h_score columns of the deliverable table, and n_with_predicted_host / top_hosts in the preview.
# Everything else kraken2 writes (k2_scores.tsv, k2_cand_scores.tsv, bact_contigs.tsv) is already
# unread by any downstream step.
#
# This is what lets the published image ship WITHOUT the 19 GB HRGMv2 reference: host prediction
# becomes an opt-in add-on rather than a hard dependency. run.sh only requires the DB when
# GB_KRAKEN=1. Production sets it to 1.
GB_KRAKEN="${GB_KRAKEN:-1}"
case "$GB_KRAKEN" in
  0|1) ;;
  *) echo "[$(basename "$0" .sh)] invalid GB_KRAKEN='$GB_KRAKEN' (want 0 or 1)" >&2; exit 1;;
esac
[ "$GB_KRAKEN" = "1" ] || echo "[$(basename "$0" .sh)] kraken2 DISABLED (GB_KRAKEN=0) -- no host prediction this run"


work="$GB_OUTDIR/evaluation"
prevwork="$GB_OUTDIR/pullviruses"
mkdir -p "$work"

echo "[5_evaluation] IN         : $GB_IN"
echo "[5_evaluation] OUTDIR     : $GB_OUTDIR"
echo "[5_evaluation] THREADS    : $GB_THREADS"

# Per-tool progress + failure localisation. This stage holds 11 of the run's 15 tool
# invocations, so "Step failed: 5_eval.sh" on its own points at an 80 KB log and nothing else.
if [ -r "${GB_SCRIPTS}/gb_progress.sh" ]; then
  . "${GB_SCRIPTS}/gb_progress.sh"
else
  gb_step() { echo "[5_evaluation] >> $1"; }; gb_tool_result() { :; }
fi

inputF=$prevwork/pulled_viruses.fna #input
virfasta=$GB_OUTDIR/minedviruses.fna #output

# STOP condition: every discovery tool (BLAST/CheckV/PhageBoost in stage 1,
# PhaMer/MetaPhaPred in stage 3, DeepVirFinder/VirSorter2 in stage 4) came up
# empty, so the merged candidate set has no viral contigs. Any single tool
# finding nothing is fine, but with nothing to score we stop here and emit empty
# deliverables — a virus-free sample is a valid, non-failing result (the scoring
# tools below would otherwise abort on an empty FASTA). Test for actual FASTA
# headers, not just file size, so a stray newline can't slip past.
if [ ! -s "$inputF" ] || ! grep -q '^>' "$inputF"; then
  echo "[5_evaluation] No viral contigs from any tool — stopping with an empty result (0 viruses)."
  : > "$GB_OUTDIR/${GB_PREFIX}_viruses.fna"
  : > "$GB_OUTDIR/${GB_PREFIX}_table.tsv"
  echo "[5_evaluation] Done :D"
  exit 0
fi

# ==============================================
# DEDUP
# ==============================================

echo "[5_evaluation] removing duplicated sequences"
awk '
  /^>/ {
    name = $1
    sub(/^>/, "", name)

    if (seen[name]++) {
      skip = 1   # we’ve already seen this ID → skip this header and its sequence
      next
    } else {
      skip = 0   # first time we see this ID → keep it
    }
  }

  skip == 0 { print }
' "$inputF" > "$work/nodup.fna"
inputF="$work/nodup.fna"

# ==============================================
# KRAKEN2 CANDIDATE PASS -- REMOVED 2026-08-07
# ==============================================
# This block subsetted stage 2's tagging pass and wrote k2_cand_scores.tsv (six columns of
# per-candidate bacterial k-mer content). Its ONLY consumer was the 4-term fallback formula,
# which was removed the same day, so the table became dead output: computed, logged, read by
# nothing. Deleted rather than left behind a GB_KRAKEN gate, because gating dead code only
# hides it.
#
# Nothing here fed the model: stack9 contains zero kraken2-derived features. The host
# attribution further down does NOT depend on this block -- it subsets stage 1's own host pass
# (HRGMv2 at --confidence 0.2), a different database at a different confidence, deliberately so.
# The measurement history for k2_bact_frac is preserved in BENCHMARKS.md.
# ==============================================
# DEEPVIRFINDER BLOCK
# ==============================================

rm -rf "${GB_OUTDIR:-.}/.gutbusters_dvf_cache/theano" 2>/dev/null || true
gb_step "DeepVirFinder — scoring contigs for viral origin"
dvf -i "$inputF" -o $work -l $GB_MLEN -c $GB_THREADS

mv $work/dvfpred.txt $work/dvp_vscores.tsv


# ==============================================
# VIRSORTER2 BLOCK -- SCORES FROM THE STAGE-4 CACHE, NO RUN HERE
# ==============================================
# VirSorter2 no longer RUNS in stage 5. It still runs in stage 4, where it MINES:
# 4_pullviruses.sh appends vs2viruses.fna to the candidate collection, so THE CANDIDATE
# POOL THIS STAGE RECEIVES IS UNCHANGED. Stage 5's own invocation only ever SCORED the
# remainder stage 4 had not already scored, and vs2 is not a term in the admission score
# below -- so nothing it produced here decides anything any more.
#
# MEASURED on a 53.05 Mbp pool (timings.tsv): 5249 s of 7802 s of stage-5 tool time,
# 67.2%. That predates the chunking fix; current-equivalent ~3490 s, still the single
# largest cost in the stage. Note the 5249 s came from a standalone run with NO stage-4
# cache, so the realised saving on a production run is smaller by the cache-hit rate --
# the log line below reports that rate so the real figure can be read off a live run.
#
# THE TABLE IS STILL WRITTEN, populated from the cache alone -- but NOT because the layout
# depends on it. CORRECTION to an earlier note here: pbbcalc.awk prints a FIXED 9-column
# header in its END block and defaults any arm it never saw to 0, so an absent
# vs2_vscores.tsv would NOT shift $6 confidence / $8 n_comp / $9 n_hom. The layout is safe
# either way. The file is kept because the cached scores are REAL DATA for the slice stage 4
# covered, and because n_hom (read by the REMOVED legacy fallback) was carried by vs2 alone (see the
# geNomad block below).
#
# INFORMATION LOSS, stated plainly: contigs stage 4 never scored now carry vs2=0.00 in the
# deliverable rather than a real score. pbbcalc treats a missing row as a silent arm, which
# is what it already did for any tool that did not call a contig. This affects the REPORTED
# column only -- admission is decided by the score gate, which has no vs2 term.
#
# 2026-08-05: VirSorter2 was removed from stage 4 as well, so the cache below is now ALWAYS
# absent and vs2_vscores.tsv is always header-only -- every contig reports vs2=0.00. The block
# is kept rather than deleted precisely because the guard already handles that case, and
# because deleting it would mean touching pbbcalc.awk's fixed 9-column contract and the legacy
# fallback gate's column indices for no gain. See 4_pullviruses.sh for the rationale.
echo "[5_evaluation] VirSorter2: removed from the pipeline (vs2 column reports 0.00; not an admission signal)"
rm -rf $work/sub
# Delete the VS2 scratch explicitly. A blanket `vs2_*` glob also matches the score
# table this block produces (vs2_vscores.tsv), and `rm -f` exits non-zero on the
# vs2_chunk* directories under `set -e`.
rm -rf $work/chunk* $work/vs2_chunk* $work/vs2_input.fna $work/vs2_cached_ids.txt $work/vs2_vscores.tsv

VS2_CACHE="$GB_OUTDIR/pullviruses/vs2_scores_cache.tsv"

# VirSorter2's own column order (no --include-groups => the 2.2.4 default
# dsDNAphage,ssDNA). pbbcalc.awk reads max_score by position, so this header must
# describe the real table.
echo -e "seqname\tdsDNAphage\tssDNA\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular" > $work/vs2_vscores.tsv
[ -s "$VS2_CACHE" ] && cat "$VS2_CACHE" >> $work/vs2_vscores.tsv || true

# Report the cache-hit rate: this is the number that says how much of the candidate pool
# still carries a real VS2 score, and it is what sizes the true saving of not running here.
_vs2_cached=$(( $(wc -l < "$work/vs2_vscores.tsv") - 1 ))
_vs2_cand=$(grep -c '^>' "$inputF" || true)
echo "[5_evaluation] VS2 scores: $_vs2_cached cached / $_vs2_cand candidates"


# ==============================================
# METAPHAPRED BLOCK
# ==============================================

rm -f $work/chunk* $work/*_mpp.tsv
# MetaPhaPred one-hot-encodes the sequence and has no symbol for N (or R/Y/S/W...). It dies
# on the whole CHUNK, not just the offending contig, so a single ambiguous base silently
# annihilates the scores of up to MAXN contigs that share that chunk.
# Mask non-ACGT to A for SCORING ONLY: $work/oneline.fna is a throwaway that is deleted right
# after the chunks are scored, and every downstream step still reads $inputF,
# so the sequence we actually ship is untouched.
# toupper() first is load-bearing -- soft-masked (lowercase) bases miss the encoder dict the
# same way N does, and without it the gsub would rewrite every lowercase base to A.
# For all-uppercase ACGT input (what MEGAHIT emits) this is byte-identical to the old line.
_mpp_amb=$(awk '!/^>/{s=toupper($0); n+=gsub(/[^ACGT]/,"",s)} END{print n+0}' $inputF)
if [ "$_mpp_amb" -gt 0 ]; then
  echo "[5_evaluation] MetaPhaPred input: masked $_mpp_amb ambiguous base(s) to A (scoring only; originals kept)"
fi
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { s=toupper($0); gsub(/[^ACGT]/,"A",s); printf("%s",s);}  END {printf("\n");}' < $inputF > $work/oneline.fna
#awk -v w=$work 'BEGIN {n=0;} /^>/ {if(n%30==0){file=sprintf("chunk%d.fa",n);} print >> w"/"file; n++; next;} { print >> w"/"file; }' < $work/oneline.fna

MAXBP=800000     # 8 Kb per chunk (tune)
MAXN=50          # max contigs per chunk (tune)
PREFIX="$work/chunk"
       
awk -v maxbp="$MAXBP" -v maxn="$MAXN" -v prefix="$PREFIX" '
BEGIN { chunk=1; bp=0; n=0; file=sprintf("%s%04d.fa", prefix, chunk) }
       
# new contig header
/^>/ {
  # if we already started a chunk, consider splitting BEFORE writing this contig
  if (seen && (bp >= maxbp || n >= maxn)) {
    chunk++; bp=0; n=0;
    file=sprintf("%s%04d.fa", prefix, chunk)
  }
  seen=1
  n++
  print $0 >> file
  next
}      
       
# sequence line(s)
{      
  bp += length($0)
  print $0 >> file
}      
' "$work/oneline.fna"


pushd "$work" >/dev/null
gb_step "MetaPhaPred — predicting phage contigs ($(ls chunk*.fa 2>/dev/null | wc -l) shards, ${GB_THREADS}-way)"
# Worker pool: keep up to GB_THREADS predict.py instances running continuously
# (wait -n frees a slot the moment any one finishes) instead of a batch barrier
# that idles every slot until the slowest chunk in the batch finishes.
running=0
for ff in chunk*.fa; do

        /opt/micromamba/envs/mpp/bin/python /opt/gutbusters/tools/MetaPhaPred/predict.py  -i $ff -o ${ff}_mpp.tsv &
        running=$((running+1))
	if [[ $running -ge $GB_THREADS ]]; then
                wait -n || true
                running=$((running-1))
        fi
done
wait
popd >/dev/null

echo -e "name\tclen\tscore" > $work/mpp_vscores.tsv
cat $work/chunk*_mpp.tsv >> $work/mpp_vscores.tsv
rm -f $work/oneline.fna $work/chunk*

# ==============================================
# PHABOX2 BLOCK — REMOVED FROM STAGE 5 (2026-07-27)
# ==============================================
# PhaMer no longer re-scores the pooled candidate set. Measured on 8,372 labelled contigs
# (Hannigan VLP/WGS physical fractionation): its marginal contribution to admission was
# 256 contigs at 43.0% precision (110 TP / 146 FP) -- the only arm scoring BELOW the 50%
# base rate, i.e. worse than a coin flip. 225 of those 256 landed in the weakest evidence
# tier and ZERO were homology-blind, so it was diluting the discovery tier without
# contributing discoveries. Dropping it moves overall precision 75.7% -> 79.8% and lifts
# every tier. Cost is 110 real viruses, 9 of them from the confident tier.
#
# This removes it from STAGE 5 ONLY. PhaMer still runs in stage 3, where it nominates
# candidates that the remaining four arms then judge -- that role is unaffected, and the
# arm is cheap there because it only sees what earlier stages left behind.
#
# Consequence for the column contract: no phamer_vscores.tsv is written, so pbbcalc.awk
# emits 9 columns instead of 10 and every positional read below shifted left by one.
# See the virome_stage5_benchmark_verdict memory for the full measurement.

# ==============================================
# GENOMAD: SCORING-ONLY PASS REINSTATED — feeds nn_chromosome (stack4_interaction)
# ==============================================
# geNomad's mining role in stage 1 is unaffected (unchanged). This second, scoring-only
# pass over the pooled candidate set was removed 2026-07-27 because nothing read its
# score under the old 4-term formula; it is reinstated because nn_chromosome is one of
# the 6 base signals of stack4_interaction (see STACK4 FEATURE EXTRACTION below). The
# gnmd column in pbb_viruses.tsv / the removed legacy fallback's n_hom is STILL 0 for every
# contig — this pass does not feed pbbcalc.awk, only the stack4 join.
echo "[5_evaluation] geNomad: scoring-only pass reinstated (feeds nn_chromosome for stack4_interaction)"

# ==============================================
# CHECKV BLOCK (annotation AND, as of 2026-08-09, 8 columns of the admission score)
# ==============================================
# CheckV runs twice for the same reason geNomad does: stage 1 uses it to DISCOVER and excise
# proviruses, and here it ANNOTATES the pooled candidate set so every admitted contig carries
# gene-level evidence -- not only the ones stage 1 happened to look at. Stage 1 only ever saw
# after_blast_remaining.fna, so BLAST-claimed contigs never reached it, and stage-1 provirus
# excisions carry new ids that would not join anyway.
#
# ⚠️ THAT LAST CLAUSE IS WHY THIS BLOCK MOVED UP HERE (2026-08-09). stack10 reads 8 CheckV
# columns, and the obvious source -- stage 1's prophage/checkv/ -- CANNOT cover the contigs stage 1
# itself created by excision. Measured on B13: all 135 contigs with no stage-1 CheckV row were
# `*_ckv_*` excisions, and scoring them with 8 median-imputed columns COST them admissions
# (72->62 sensitive, 65->46 balanced, 48->29 strict) against stack9. Excised proviruses are prime
# viral candidates, so that is precisely the wrong population to degrade. THIS pass runs on
# $inputF -- the deduplicated candidate pool, excisions included -- so it covers everything the
# join scores. It must therefore stay UPSTREAM of the stack4 join below.
#
# CheckV is NOT a noisy-OR arm and its output must NEVER be named *_vscores.tsv (that glob
# feeds pbbcalc.awk): CheckV emits no virus probability. What it emits is viral_genes /
# host_genes, and MEASURED on 8,372 labelled contigs (Hannigan VLP/WGS physical
# fractionation) their FRACTION f = viral/(viral+host) is the joint-strongest single signal
# available anywhere in this pipeline -- ranking AUC 0.7620, second only to kraken2's
# 1-bact_frac at 0.7801 and far ahead of every caller. Fitted inside stack10 the same census
# scores 0.9774 on curated as a single column, above geNomad's 0.9689.
#
# The three-way tier below is retained for display only:
#     viral_genes >= 1                      1025 contigs   93.5% precise   -> strong
#     no genes annotated at all              494 contigs   78.5%           -> unannotated
#     viral_genes == 0 AND host_genes >= 1   570 contigs   56.1%           -> weak
# The middle tier is AGNOSTIC (CheckV annotated nothing, so it holds no opinion), which is
# why it sits near the average rather than at the bottom -- do not collapse it into "weak".
#
# ⚠️ FAILURE SEMANTICS CHANGED with stack10. This block used to be annotation-only, so a CheckV
# failure cost display and nothing else. It now feeds 8 scored columns. A failure is still
# fail-OPEN -- ckv.tsv ends up empty, the join contributes no columns, and score_stack4.py
# imputes all 8 from the frozen training medians -- which degrades to roughly stack9 rather than
# failing the run. That is deliberate: CheckV is one tool among nine and losing it should cost
# accuracy, not the analysis.
ckv_dir="$work/checkv"
rm -rf "$ckv_dir"
mkdir -p "$ckv_dir"
: > "$work/ckv_quality.tsv"
: > "$work/ckv_completeness.tsv"
if [ -s "$inputF" ]; then
  echo "[5_evaluation] CheckV: annotating the candidate set (evidence tier + 8 stack10 columns)"
  # Same invocation as stage 1. CHECKVDB is baked into the image (Apptainer.def) and
  # re-exported by entrypoint/gutbusters, so there is no -d flag to pass here.
  gb_step "CheckV — annotating candidates (evidence tier + 8 stack10 columns)"
  _t0=$SECONDS; _rc=0
  /opt/micromamba/envs/gutbusters/bin/checkv end_to_end "$inputF" "$ckv_dir" \
      -t "$GB_THREADS" --remove_tmp || _rc=$?
  [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: CheckV failed -- tier empty, stack10's 8 columns median-imputed" >&2
  gb_tool_result "CheckV annotation" "$_rc" "$((SECONDS - _t0))" "evidence tier empty, 8 ckv_* columns median-imputed"
  if [ -s "$ckv_dir/quality_summary.tsv" ]; then
    cp -f "$ckv_dir/quality_summary.tsv" "$work/ckv_quality.tsv"
  else
    echo "[5_evaluation] WARN: no quality_summary.tsv -- evidence tier empty" >&2
  fi
  # completeness.tsv carries the 4 aai_* columns and was previously discarded by the rm below.
  [ -s "$ckv_dir/completeness.tsv" ] && cp -f "$ckv_dir/completeness.tsv" "$work/ckv_completeness.tsv"
else
  echo "[5_evaluation] candidate set empty -- skipping CheckV"
fi
rm -rf "$ckv_dir"

# ==============================================
# STACK4 FEATURE EXTRACTION
# ==============================================
# The signals stack9 reads that are not already computed above: geNomad (6 columns), the
# DeepVirFinder embedding (dvf_pca_01/03), km4 (136), and features.tsv, which now serves only as
# the id driver for the join below. mpp is already computed (mpp_vscores.tsv). These are
# independent, read-only passes over the same $inputF and run concurrently.
#
# ⚠️ MEASURED 2026-08-09, and NOT what the previous comment here claimed. geNomad does not
# dominate this block: on a 20,031-contig pool geNomad is 726s and DeepVirFinder is 1,665s, so DVF
# sets the block's wall time at 2.3x geNomad. Anything cheaper than DVF is free in its shadow, and
# anything that would shorten this block has to shorten DVF or drop contigs before it.
# Downstream, contrast (1,030s) and PHROGs (1,483s) run SERIALLY after this block and depend on
# geNomad's proteins, not on each other.
#
# Any one of these failing PER CONTIG does not fail the run: score_stack4.py fills a missing
# value from the artifact's frozen training median when the column is imputable (PHROGs, km4
# and the contrast block all are), and excludes the contig otherwise.
# ⚠️ A WHOLESALE failure is different and now ERRORS (changed 2026-08-07). If every contig ends
# up unscored the stage exits non-zero rather than falling back to a weaker formula -- see the
# ADMISSION SCORE block below for why the two fallback tiers were removed.
STACK4_DB="${GB_DB_ROOT}/stack4"
stack4_dir="$work/stack4"
mkdir -p "$stack4_dir"
: > "$stack4_dir/vscore_stack4.tsv"

GENOMAD_DB="${GB_GENOMAD_DB:-/opt/gutbusters/genomad_db}"
gn5_dir="$stack4_dir/genomad"
mkdir -p "$gn5_dir"

if [ -s "$inputF" ]; then
  # ONE gb_step for the whole block: these four run CONCURRENTLY, so numbering them separately
  # would advance the percentage in an order unrelated to wall clock. Their individual outcomes
  # are reported by gb_tool_result after the wait -- which is also the only place a user ever
  # learns that one of them failed open and left its columns to be median-imputed.
  gb_step "feature block — geNomad · DVF-embed · composition · km4 (4 in parallel)"
  _t0_block=$SECONDS
  # `|| _rc=$?` keeps the set -e exemption the bare `||` gave, and the subshell still ENDS on a
  # successful command, so `wait` below still sees 0 and the fail-open behaviour is unchanged.
  # The rc file is the only new thing: it is what lets the outcome be reported at all.
  ( _rc=0
    gnmd end-to-end --cleanup --threads "$GB_THREADS" --sensitivity "${GB_GENOMAD_SENS:-4.2}" \
      "$inputF" "$gn5_dir" "$GENOMAD_DB" || _rc=$?
    [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: stack4 geNomad pass failed -- nn_chromosome unavailable" >&2
    echo "$_rc" > "$stack4_dir/.rc_genomad"
  ) &
  _pid_gn=$!

  ( _rc=0
    /opt/micromamba/envs/dvf/bin/python "$GB_SCRIPTS/dvf_embed.py" \
      --fasta "$inputF" \
      --model "${GB_DB_ROOT}/deepvirfinder/models/model_siamese_varlen_1k_fl10_fn1000_dn1000_ep30_acc0.96.h5" \
      --pca-basis "$STACK4_DB/dvf_pca_frozen_basis.json" \
      --out "$stack4_dir/dvf_pca.tsv" || _rc=$?
    [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: dvf_embed.py failed -- dvf_pca_* unavailable" >&2
    echo "$_rc" > "$stack4_dir/.rc_dvfembed"
  ) &
  _pid_dvfembed=$!

  # extract_gene_features.py REMOVED 2026-08-09 -- it was dead compute. It produced
  # gc3_content_max / skew_concordance / strand_bias_B, which are stack4-era columns that left the
  # feature list at stack6 and are not among stack9's 168. gene_features.tsv was written, touched
  # for existence, and never read: it is not one of the 8 files in the join below, and
  # score_stack4.py reads only the joined table. Measured 1,270s of CPU per run on a 20k-contig
  # pool. It ran concurrently with DeepVirFinder (1,665s), so removing it does not shorten the
  # block on its own -- it hands 8 threads back to the passes that do.

  ( _rc=0
    /opt/micromamba/envs/gutbusters/bin/python3 "$GB_SCRIPTS/extract_features.py" \
      --fasta "$inputF" --pca-basis "$STACK4_DB/stack4_final_model.json" \
      --out "$stack4_dir/features.tsv" --workers "$GB_THREADS" || _rc=$?
    [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: extract_features.py failed -- delta_g_per_kb/mean_propeller_twist/km4dist_pca_00 unavailable" >&2
    echo "$_rc" > "$stack4_dir/.rc_features"
  ) &
  _pid_features=$!

  # stack9: 136 RAW canonical tetranucleotide frequencies. RAW, not PCA -- that distinction is the
  # entire result (the 2026-08-02 k-mer rejection used a dead label AND PCA-15'd the block before
  # any model saw it). Measured 1.8s for 14,916 contigs / 87.6 Mbp, so it joins this block rather
  # than earning its own; it is noise next to geNomad's ~509s.
  ( _rc=0
    /opt/micromamba/envs/gutbusters/bin/python3 "$GB_SCRIPTS/extract_km4.py" \
      --fasta "$inputF" --out "$stack4_dir/km4.tsv" --workers "$GB_THREADS" || _rc=$?
    [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: extract_km4.py failed -- km4_* unavailable (median-imputed)" >&2
    echo "$_rc" > "$stack4_dir/.rc_km4"
  ) &
  _pid_km4=$!

  wait "$_pid_gn" "$_pid_dvfembed" "$_pid_features" "$_pid_km4"

  # Report each tool's real outcome. A missing rc file means the subshell died before it could
  # write one (OOM-kill, for instance) -- that is a failure too, so it defaults to non-zero.
  _rc_of() { cat "$stack4_dir/.rc_$1" 2>/dev/null || echo 1; }
  _blk=$((SECONDS - _t0_block))
  gb_tool_result "geNomad"        "$(_rc_of genomad)"  "$_blk" "nn_chromosome unavailable — stack4_interaction median-imputed"
  gb_tool_result "DVF embeddings" "$(_rc_of dvfembed)" "$_blk" "dvf_pca_01 / dvf_pca_03 median-imputed"
  gb_tool_result "composition"    "$(_rc_of features)" "$_blk" "delta_g_per_kb / mean_propeller_twist / km4dist_pca_00 median-imputed"
  gb_tool_result "km4 tetranuc"   "$(_rc_of km4)"      "$_blk" "all 136 km4_* columns median-imputed — 77% of the feature space"

  # geNomad writes TWO files matching *_nn_classification.tsv into the same directory:
  #   <prefix>_nn_classification.tsv            one row per input contig
  #   <prefix>_provirus_nn_classification.tsv   only the excised provirus regions (often ~1 row)
  # `-print -quit` returned whichever the filesystem yielded first, which is directory order and
  # NOT deterministic -- observed picking each of the two on the same host. When it picked the
  # provirus file, nn_chromosome held a couple of rows, the join below dropped EVERY contig on
  # its `!(id in nnc)` test, and score_stack4.py reported scored=0: a silently empty admission
  # score for the whole run, with no error, because geNomad itself had succeeded.
  # Exclude the provirus variant and take the largest remaining match.
  # geNomad end-to-end already produces far more than the neural-network scores this stage used
  # to take. The marker branch -- 200k marker-protein profiles -- is where geNomad's published
  # MCC 95.3% comes from, and it was being computed and discarded on every run. Harvest it.
  #
  # In each case geNomad writes BOTH <prefix>_X.tsv and <prefix>_provirus_X.tsv into the same
  # directory. The old `-print -quit` returned whichever the filesystem yielded first, which is
  # directory order and NOT deterministic; when it returned the provirus file (a handful of
  # rows) the join dropped every contig and the gate silently fell back to the legacy formula.
  # Exclude the provirus variant and take the largest remaining match.
  _gn_pick() {
    find "$gn5_dir" -maxdepth 3 -iname "$1" ! -iname '*provirus*' \
      -printf '%s\t%p\n' 2>/dev/null | sort -rn | head -1 | cut -f2
  }
  gn5_agg_tsv="$(_gn_pick '*_aggregated_classification.tsv' || true)"
  gn5_mrk_tsv="$(_gn_pick '*_features.tsv' || true)"
  : > "$stack4_dir/genomad_agg.tsv"
  : > "$stack4_dir/genomad_marker.tsv"
  [ -n "$gn5_agg_tsv" ] && [ -s "$gn5_agg_tsv" ] && cp -f "$gn5_agg_tsv" "$stack4_dir/genomad_agg.tsv"
  [ -n "$gn5_mrk_tsv" ] && [ -s "$gn5_mrk_tsv" ] && cp -f "$gn5_mrk_tsv" "$stack4_dir/genomad_marker.tsv"
  echo "[5_evaluation] stack5 geNomad aggregated: ${gn5_agg_tsv:-NONE} ($(wc -l < "$stack4_dir/genomad_agg.tsv") rows)"
  echo "[5_evaluation] stack5 geNomad markers:    ${gn5_mrk_tsv:-NONE} ($(wc -l < "$stack4_dir/genomad_marker.tsv") rows)"

  [ -s "$stack4_dir/dvf_pca.tsv" ] || : > "$stack4_dir/dvf_pca.tsv"
  [ -s "$stack4_dir/features.tsv" ] || : > "$stack4_dir/features.tsv"

  # ==============================================
  # VIRAL-VS-CELLULAR SEQUENCE CONTRAST
  # ==============================================
  # Per predicted protein: delta = best_bits(BFVD, viral) - best_bits(AlphaFold/Swiss-Prot,
  # cellular), aggregated per contig into 4 features. See scripts/contrast_features.py.
  #
  # THE CONTRAST IS THE SIGNAL, NOT THE VIRAL HITS. Hits to BFVD alone rank at 0.4943 AUC --
  # chance -- because BFVD is full of AMGs and replication machinery whose sequences are shared
  # with cellular life. So BOTH searches are required; dropping the cellular one to halve the
  # runtime does not halve the signal, it deletes it.
  #
  # MEASURED standalone (grouped CV, honest out-of-fold): 0.8776 on CHGV virality and 0.9880 on
  # curated novel phage, above geNomad (0.7494/0.9689) and MetaPhaPred (0.8436/0.9629) on the
  # same contigs. Added to the 9-signal stack it moves CHGV 0.8864 -> 0.9099, paired donor
  # bootstrap 95% CI [+0.0196,+0.0275].
  #
  # STRUCTURE WAS TESTED AND REJECTED as the mechanism: Foldseek 3Di + ProstT5 over the SAME two
  # databases adds only +0.0025 AUC over this sequence search, and needs a GPU that Mars does not
  # have. The sequence control was run at -s 7.5 (default 5.7) deliberately, so the comparison
  # could not be won by handicapping it -- and those are the settings kept here, because the
  # frozen coefficients were fitted on the values they produce.
  #
  # GENE CALLING IS FREE: this consumes geNomad's own Prodigal-called proteins, which the pass
  # above already produced for nn_chromosome. No second gene-calling pass.
  #
  # Runs AFTER the parallel block rather than inside it: it depends on geNomad's output.
  VC_DB="${GB_DB_ROOT}/vcontrast"
  vc_dir="$stack4_dir/contrast"
  rm -rf "$vc_dir"; mkdir -p "$vc_dir"
  : > "$stack4_dir/contrast.tsv"

  # Match the annotate directory EXPLICITLY rather than globbing '*_proteins.faa': geNomad also
  # writes <prefix>_summary/<prefix>_virus_proteins.faa and _plasmid_proteins.faa (subsets of
  # this file, so a "largest match" tie-break is not robust) and find_proviruses/
  # <prefix>_provirus_proteins.faa (different id space entirely). Picking any of those would
  # silently compute the contrast over the wrong gene set. Same class of trap as the
  # provirus/aggregated ambiguity handled by _gn_pick above.
  gn5_faa="$(find "$gn5_dir" -maxdepth 3 -path '*_annotate/*_proteins.faa' ! -name '*provirus*' \
              -printf '%s\t%p\n' 2>/dev/null | sort -rn | head -1 | cut -f2)"

  if [ -n "${gn5_faa:-}" ] && [ -s "$gn5_faa" ] && [ -s "$VC_DB/bfvd_seq.dbtype" ] \
     && [ -s "$VC_DB/afdb_sp_seq.dbtype" ]; then
    echo "[5_evaluation] contrast: $(grep -c '^>' "$gn5_faa") proteins from $gn5_faa"
    _vc_ok=1
    _mm=/opt/micromamba/envs/gnmd/bin/mmseqs
    if ! "$_mm" createdb "$gn5_faa" "$vc_dir/q" > "$vc_dir/createdb.log" 2>&1; then
      echo "[5_evaluation] WARN: contrast createdb failed -- vc_* unavailable" >&2
      _vc_ok=0
    fi
    for _ref in bfvd afdb_sp; do
      [ "$_vc_ok" -eq 1 ] || break
      # -s 7.5 / -e 1e-3 / --max-seqs 300 are a CONTRACT with the frozen coefficients: they are
      # the settings the training features were computed under. Retuning them for speed silently
      # shifts the feature distribution out from under the model.
      if ! "$_mm" search "$vc_dir/q" "$VC_DB/${_ref}_seq" "$vc_dir/res_$_ref" "$vc_dir/tmp_$_ref" \
             -s 7.5 -e 1e-3 --max-seqs 300 --threads "$GB_THREADS" \
             > "$vc_dir/search_$_ref.log" 2>&1; then
        echo "[5_evaluation] WARN: contrast search vs $_ref failed -- vc_* unavailable" >&2
        _vc_ok=0; break
      fi
      # No taxid column: these targets are FASTA-derived MMseqs DBs and carry no taxonomy, so
      # requesting taxid makes convertalis abort with "mapping does not exist". The 4 shipped
      # features do not use it.
      if ! "$_mm" convertalis "$vc_dir/q" "$VC_DB/${_ref}_seq" "$vc_dir/res_$_ref" \
             "$vc_dir/$_ref.m8" --threads "$GB_THREADS" \
             --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov" \
             > "$vc_dir/conv_$_ref.log" 2>&1; then
        echo "[5_evaluation] WARN: contrast convertalis vs $_ref failed -- vc_* unavailable" >&2
        _vc_ok=0; break
      fi
      rm -rf "$vc_dir/tmp_$_ref" "$vc_dir/res_$_ref"*
      echo "[5_evaluation] contrast $_ref: $(wc -l < "$vc_dir/$_ref.m8") hits"
    done
    if [ "$_vc_ok" -eq 1 ]; then
      # --fasta is REQUIRED for coding_density: the contig-length definition reproduces the
      # training column 100.00% of the time, the max-gene-end fallback only 27.16%. $inputF is
      # the same nodup.fna geNomad was run on, so its contig ids match the protein headers.
      gb_step "MMseqs2 contrast — viral-vs-cellular homology"
      _t0=$SECONDS; _rc=0
      /opt/micromamba/envs/gutbusters/bin/python3 "$GB_SCRIPTS/contrast_features.py" \
        --faa "$gn5_faa" --viral-m8 "$vc_dir/bfvd.m8" --cellular-m8 "$vc_dir/afdb_sp.m8" \
        --fasta "$inputF" \
        --out "$stack4_dir/contrast.tsv" || _rc=$?
      [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: contrast_features.py failed -- vc_* unavailable" >&2
      gb_tool_result "MMseqs2 contrast" "$_rc" "$((SECONDS - _t0))" "the 4 vc_* columns median-imputed"
    fi
  else
    echo "[5_evaluation] WARN: contrast skipped (proteins='${gn5_faa:-NONE}', db=$VC_DB) -- vc_* unavailable" >&2
  fi
  # Scratch only; the m8 tables are kept for debugging but the query DB is the bulky part.
  rm -rf "$vc_dir/q" "$vc_dir/q_h" "$vc_dir/q".* "$vc_dir/q_h".* 2>/dev/null || true
  [ -s "$stack4_dir/contrast.tsv" ] || : > "$stack4_dir/contrast.tsv"

  # ==============================================
  # PHROGs PROFILE FAMILIES (stack9)
  # ==============================================
  # 38,880 viral protein families built by HMM-HMM remote homology comparison -- a different
  # MODALITY from the pairwise contrast above, reaching further into the twilight zone, and the
  # only feature block carrying FUNCTIONAL CATEGORIES. Measured +0.0065 CHGV AUC on top of stack7,
  # and it GROWS as the stack strengthens rather than being absorbed.
  #
  # ⚠️ HMMER, NEVER MMseqs2 -- a correctness constraint, not a preference. The Pharokka bundle's
  # phrogs_profile_db was built with MMseqs2 v18.8cc5c; this image ships 14.7e284. The prefilter
  # appears to succeed and the align stage dies with "Alignment died" -- a silent format
  # incompatibility, reproduced at 15,000 targets with --max-seqs 3000 and 24G split memory. The
  # .h3m in the same bundle is HMMER3 binary, whose format is stable across versions, and it is
  # what the shipped training features were computed on. pyhmmer is pinned 0.10.15 in its own env.
  #
  # GENE CALLING IS FREE: consumes geNomad's own proteins, the same .faa the contrast used.
  #
  # RUNS AFTER THE CONTRAST, NOT BESIDE IT. Both depend on gn5_faa and neither depends on the
  # other, so the parallel pattern above would apply -- but both are heavy and thread-hungry, and
  # Mars is a shared VM with measured ~40% CPU steal. Sequential at full GB_THREADS is predictable;
  # concurrent would oversubscribe and make both slower. Cost is ~6 min on a 14,916-contig set.
  # RETIRED 2026-08-13, DEFAULT OFF. The stack158 admission model reads none of the 15 phrog_*/
  # cat_* columns, so this pass became dead compute the moment that model shipped -- ~6 min/run and
  # a 1.3 GB .h3m for columns nothing consumes. Removing PHROGs was measured, not assumed: on the
  # development collections its removal is marginally POSITIVE (CHGV two-way bracket, §37), and the
  # frozen-artifact transfer to Hannigan and Shkoporov is neutral either way.
  #
  # GATED, NOT DELETED, following the GB_PHAGEBOOST pattern in 1_prophage.sh:305. The columns are
  # still needed to reproduce the published benchmarks, where PHROGs remains a baseline arm -- so
  # `GB_PHROGS=1` restores the pass exactly. The empty phrog.tsv left below is a pre-existing,
  # already-exercised code path (the `else` branch this replaces), so the awk join's column
  # contract is untouched and the 15 columns simply come out empty for a model that never asks.
  PHROG_DB="${GB_DB_ROOT}/phrog"
  : > "$stack4_dir/phrog.tsv"
  if [ "${GB_PHROGS:-0}" != "1" ]; then
    echo "[5_evaluation] PHROGs retired -- stack158 reads none of its 15 columns (GB_PHROGS=1 re-enables)"
  elif [ -n "${gn5_faa:-}" ] && [ -s "$gn5_faa" ] && [ -s "$PHROG_DB/all_phrogs.h3m" ] \
     && [ -s "$PHROG_DB/phrog_annot_v4.tsv" ]; then
    gb_step "pyhmmer/PHROGs — scanning proteins against 38k profiles"
    _t0=$SECONDS; _rc=0
    /opt/micromamba/envs/pyhmmer/bin/python "$GB_SCRIPTS/phrog_features.py" \
      --faa "$gn5_faa" --hmm "$PHROG_DB/all_phrogs.h3m" \
      --annot "$PHROG_DB/phrog_annot_v4.tsv" \
      --out "$stack4_dir/phrog.tsv" --cpus "$GB_THREADS" || _rc=$?
    [ "$_rc" -eq 0 ] || echo "[5_evaluation] WARN: phrog_features.py failed -- phrog_*/cat_* median-imputed" >&2
    gb_tool_result "pyhmmer/PHROGs" "$_rc" "$((SECONDS - _t0))" "the 15 phrog_*/cat_* columns median-imputed"
  else
    echo "[5_evaluation] WARN: PHROGs REQUESTED (GB_PHROGS=1) but skipped (proteins='${gn5_faa:-NONE}', db=$PHROG_DB) -- phrog_*/cat_* median-imputed" >&2
  fi
  [ -s "$stack4_dir/phrog.tsv" ] || : > "$stack4_dir/phrog.tsv"

  # stack10's 8 CheckV columns, distilled from the two tables the hoisted CheckV block kept.
  # ZERO new compute: CheckV already ran over this exact pool. Columns are located BY HEADER NAME,
  # never by position, so a CheckV version that reorders or inserts columns cannot silently rebind
  # a feature. viral_frac uses the +1/+2 smoothing the model was fitted with -- it is a feature
  # definition, not a display choice, and changing it here would silently mis-score every contig.
  # A contig present in quality_summary but absent from completeness (or vice versa) gets EMPTY
  # fields for the missing half, which score_stack4.py fills from the frozen training median.
  : > "$stack4_dir/ckv.tsv"
  if [ -s "$work/ckv_quality.tsv" ]; then
    awk -F'\t' -v OFS='\t' '
      { sub(/\r$/, "") }
      FILENAME==ARGV[1] {                        # completeness.tsv: the 4 aai_* columns
        if (FNR==1) { for (i=1;i<=NF;i++) h[$i]=i; next }
        ac[$1]=("aai_completeness" in h ? $(h["aai_completeness"]) : "")
        an[$1]=("aai_num_hits"     in h ? $(h["aai_num_hits"])     : "")
        ai[$1]=("aai_id"           in h ? $(h["aai_id"])           : "")
        af[$1]=("aai_af"           in h ? $(h["aai_af"])           : "")
        next
      }
      FILENAME==ARGV[2] {                        # quality_summary.tsv: the gene census
        if (FNR==1) {
          for (i=1;i<=NF;i++) q[$i]=i
          print "contig_id","ckv_viral_genes","ckv_host_genes","ckv_gene_count","ckv_viral_frac", \
                "ckv_aai_completeness","ckv_aai_num_hits","ckv_aai_id","ckv_aai_af"
          next
        }
        vg=("viral_genes" in q ? $(q["viral_genes"]) : "")
        hg=("host_genes"  in q ? $(q["host_genes"])  : "")
        gc=("gene_count"  in q ? $(q["gene_count"])  : "")
        # CheckV writes NA for a contig it could not annotate; that is a gap, not a zero.
        vf = (vg ~ /^[0-9]+$/ && hg ~ /^[0-9]+$/) ? (vg + 1) / (vg + hg + 2) : ""
        print $1, vg, hg, gc, vf, \
              ($1 in ac ? ac[$1] : ""), ($1 in an ? an[$1] : ""), \
              ($1 in ai ? ai[$1] : ""), ($1 in af ? af[$1] : "")
      }
    ' "$work/ckv_completeness.tsv" "$work/ckv_quality.tsv" > "$stack4_dir/ckv.tsv" \
      || : > "$stack4_dir/ckv.tsv"
  fi

  # Join the 5 inputs (mpp already computed above) keyed by contig id. A contig
  # missing ANY of the 6 base signals or gc3_content_max is dropped here rather
  # than imputed -- score_stack4.py's own missing-value handling is the second,
  # redundant line of defense for per-contig (not wholesale) tool failures.
  #
  # skew_concordance/strand_bias_B are deliberately NOT in the drop condition:
  # they are legitimately undefined for ~2% of contigs (no gene in a window with
  # |GC skew| > 0.05) and travel through as EMPTY fields, which score_stack4.py
  # fills from the model's frozen training median -- the same imputation the fit
  # used. Dropping those contigs instead would discard candidates the model was
  # trained to score. Column order here is a CONTRACT with score_stack4.py's
  # FEATURE_ORDER and the frozen coefficients; changing it silently rebinds
  # every coefficient to the wrong feature.
  # Column order is a CONTRACT with the artifact's "features" list; changing it silently
  # rebinds every coefficient to the wrong feature.
  #
  # n_uscg / n_virus_hallmarks / marker_enrichment_v are legitimately undefined for a contig on
  # which geNomad called no genes. They travel through as EMPTY and score_stack4.py fills them
  # from the model's frozen training median -- the same imputation the fit used. Dropping those
  # contigs would discard candidates the model was trained to score. The DROP condition covers
  # only the signals that are always defined when their tool ran at all.
  awk -F'\t' -v OFS='\t' -v nf="$stack4_dir/join_dropped.count" '
    # ⚠️ STRIP CR FIRST, on every record of every input. Two of these tables were written by a
    # csv.writer whose DEFAULT lineterminator is "\r\n", and the PHROGs and km4 blocks below are
    # concatenated into the output line VERBATIM as raw substrings. A trailing CR therefore landed
    # MID-LINE in joined_features.tsv -- and csv.DictReader (which score_stack4.py uses) treats a
    # bare CR as a row terminator, so the header stopped at the PHROGs block and ALL 136 km4_*
    # columns were invisible to the scorer, silently median-imputed on every run since 2026-08-07.
    # The writers were fixed too (lineterminator="\n"); this stays as the belt-and-braces, because
    # the failure is completely silent and any future producer could reintroduce it.
    { sub(/\r$/, "") }
    FILENAME==ARGV[1] { if (FNR>1) mppv[$1]=$3; next }
    FILENAME==ARGV[2] { if (FNR>1) { d01[$1]=$2; d03[$1]=$3 }; next }
    FILENAME==ARGV[3] {                       # geNomad aggregated: chromosome, plasmid, virus
      if (FNR==1) { for (i=1;i<=NF;i++) h[$i]=i; next }
      gchr[$1]=$(h["chromosome_score"]); gpla[$1]=$(h["plasmid_score"]); gvir[$1]=$(h["virus_score"]); next
    }
    FILENAME==ARGV[4] {                       # geNomad marker features (27 columns)
      if (FNR==1) { for (i=1;i<=NF;i++) m[$i]=i; next }
      if ("n_uscg" in m)              uscg[$1]=$(m["n_uscg"])
      if ("n_virus_hallmarks" in m)   vhall[$1]=$(m["n_virus_hallmarks"])
      if ("marker_enrichment_v" in m) menv[$1]=$(m["marker_enrichment_v"])
      next
    }
    FILENAME==ARGV[5] {                       # viral-vs-cellular contrast (8 features)
      if (FNR==1) { for (i=1;i<=NF;i++) v[$i]=i; next }
      vmd[$1]=$(v["vc_mean_delta"]); vp90[$1]=$(v["vc_delta_p90"])
      vmax[$1]=$(v["vc_delta_max"]); vfrc[$1]=$(v["vc_frac"])
      # stack7: three more columns out of the SAME two searches, at no extra cost.
      vdrk[$1]=$(v["dark_frac"]); vmbv[$1]=$(v["mean_bits_v"])
      vcd[$1]=$(v["coding_density"]); vng[$1]=$(v["n_genes"]); next
    }
    FILENAME==ARGV[6] {                       # stack9 PHROGs (15 features)
      # Header and blank-width come from the FILE, not from a hardcoded count: if this table is
      # empty (the pass failed or was skipped) the header contributes nothing and rows contribute
      # nothing, so the two stay aligned and the columns simply never appear. score_stack4.py then
      # takes its KeyError path into impute_median -- graceful degradation, not a broken join.
      if (FNR==1) { ph_h=substr($0,index($0,"\t")); ph_b=""; for(i=1;i<NF;i++) ph_b=ph_b"\t"; next }
      ph[$1]=substr($0,index($0,"\t")); next
    }
    FILENAME==ARGV[7] {                       # stack9 km4 RAW (136 features)
      if (FNR==1) { km_h=substr($0,index($0,"\t")); km_b=""; for(i=1;i<NF;i++) km_b=km_b"\t"; next }
      km[$1]=substr($0,index($0,"\t")); next
    }
    FILENAME==ARGV[8] {                       # stack10 CheckV (8 features)
      if (FNR==1) { ck_h=substr($0,index($0,"\t")); ck_b=""; for(i=1;i<NF;i++) ck_b=ck_b"\t"; next }
      ck[$1]=substr($0,index($0,"\t")); next
    }
    FILENAME==ARGV[9] {
      if (FNR==1) { print "id","gn_virus","gn_chr","gn_plasmid","mpp","n_uscg","n_virus_hallmarks","marker_enrichment_v","dvf_pca_01","dvf_pca_03","vc_mean_delta","vc_delta_p90","vc_delta_max","vc_frac","dark_frac","mean_bits_v","coding_density","n_genes" ph_h km_h ck_h; next }
      id=$1
      if (!(id in mppv) || !(id in d01) || !(id in d03) || !(id in gvir)) { dropped++; next }
      # The vc_*/dark_frac/mean_bits_v/coding_density/n_genes fields are NOT in the drop
      # condition: they all come from the contrast block, are legitimately undefined for a contig
      # on which geNomad called no genes (the same population n_uscg is undefined for), and
      # travel through EMPTY for score_stack4.py to fill from the frozen training median.
      #
      # PHROGs is absent for exactly that same population (no genes -> no proteins to scan), and
      # is likewise imputed. A contig that HAS genes but no PHROG hit arrives as real 0.0 from the
      # extractor, which is a measurement, not a gap -- do not conflate the two.
      _ph = (id in ph ? ph[id] : ph_b)
      _km = (id in km ? km[id] : km_b)
      # A contig absent from CheckV is NOT dropped: the 8 columns are imputable and CheckV
      # legitimately annotates nothing on some contigs. See the hoisted CheckV block for why the
      # SOURCE of this table matters more than the imputation.
      _ck = (id in ck ? ck[id] : ck_b)
      print id, gvir[id], gchr[id], gpla[id], mppv[id], \
            (id in uscg ? uscg[id] : ""), (id in vhall ? vhall[id] : ""), \
            (id in menv ? menv[id] : ""), d01[id], d03[id], \
            (id in vmd ? vmd[id] : ""), (id in vp90 ? vp90[id] : ""), \
            (id in vmax ? vmax[id] : ""), (id in vfrc ? vfrc[id] : ""), \
            (id in vdrk ? vdrk[id] : ""), (id in vmbv ? vmbv[id] : ""), \
            (id in vcd ? vcd[id] : ""), (id in vng ? vng[id] : "") _ph _km _ck
    }
    END { print dropped+0 > nf }
  ' "$work/mpp_vscores.tsv" "$stack4_dir/dvf_pca.tsv" "$stack4_dir/genomad_agg.tsv" "$stack4_dir/genomad_marker.tsv" "$stack4_dir/contrast.tsv" "$stack4_dir/phrog.tsv" "$stack4_dir/km4.tsv" "$stack4_dir/ckv.tsv" "$stack4_dir/features.tsv" \
    > "$stack4_dir/joined_features.tsv"
  # A wholesale failure of either new pass is SILENT otherwise: every contig would be scored with
  # median-filled columns, i.e. stack7-with-flat-blocks, and nothing would look wrong. Report the
  # coverage so a broken pass is visible in the log.
  # ⚠️ This reported 0 for EVERY table from the day it was written until 2026-08-09. `print NR>0 ?
  # NR-1 : 0` does not parse as a ternary: awk reads `print NR > 0` as output REDIRECTION to a file
  # named "0", then hits `? :` and dies with a syntax error. The 2>/dev/null hid the error and the
  # `|| echo 0` supplied a plausible-looking answer, so the one instrument guarding against a
  # WHOLESALE PHROGs/km4 failure -- which would score every contig with two blocks median-imputed
  # and otherwise look completely normal -- could never fire. Keep the parentheses.
  # Counts non-empty DATA rows rather than NR-1, so a table written with the interleaved blank rows
  # seen in the vscore tables reports contigs, not lines.
  _nrows() { awk -F'\t' 'NR>1 && $1!=""{n++} END{print (n+0)}' "$1" 2>/dev/null || echo 0; }
  echo "[5_evaluation] stack10 blocks: PHROGs $(_nrows "$stack4_dir/phrog.tsv") contig(s), km4 $(_nrows "$stack4_dir/km4.tsv") contig(s), CheckV $(_nrows "$stack4_dir/ckv.tsv") contig(s)"
  echo "[5_evaluation] stack5 join: $(cat "$stack4_dir/join_dropped.count" 2>/dev/null || echo 0) contig(s) dropped (missing an input signal)"

  gb_step "admission scoring — 158-signal ExtraTrees model"
  /opt/micromamba/envs/gutbusters/bin/python3 "$GB_SCRIPTS/score_stack4.py" \
    --input "$stack4_dir/joined_features.tsv" --model "$STACK4_DB/stack158_final_model.json" \
    --out "$stack4_dir/vscore_stack4.tsv" \
  || echo "[5_evaluation] WARN: score_stack4.py failed -- stack5 admission score unavailable" >&2
else
  echo "[5_evaluation] candidate set empty -- skipping stack4 feature extraction"
fi

_stack4_scored=$(awk -F'\t' 'NR>1 && $2!="" {c++} END{print c+0}' "$stack4_dir/vscore_stack4.tsv" 2>/dev/null || echo 0)
echo "[5_evaluation] stack4: $_stack4_scored contig(s) scored"

# CheckV ran earlier -- it was hoisted above the STACK4 FEATURE EXTRACTION block on 2026-08-09
# because stack10 reads 8 of its columns and the join needs them. See that block for the full
# rationale. $work/ckv_quality.tsv and $work/ckv_completeness.tsv are already populated here.

# ==============================================
# Calculation of viral score
# ==============================================

# Score at which an arm counts as a supporting vote. Read from the environment rather
# than added as an entrypoint flag, following the GB_GENOMAD_SENS precedent: the reader
# is baked into the .sif, but the VALUE stays retunable from run.sh (scp, no rebuild) or
# from the Mars worker .env.
#
# Deliberately NOT folded into GB_PHS. That variable used to double as VirSorter2's
# --min-score (VS2 removed 2026-08-05) and still sets the legacy fallback gate's confidence
# threshold below, so the two remain different quantities that must not share a knob.
GB_CORROB_FLOOR="${GB_CORROB_FLOOR:-0.5}"
# Reject a non-numeric value here rather than letting awk turn it into a lexicographic
# comparison that passes everything. "0,5" is one keystroke from "0.5" and is the
# decimal separator in several locales, so this is not hypothetical.
case "$GB_CORROB_FLOOR" in
  ''|*[!0-9.]*|*.*.*)
    echo "[5_evaluation] invalid GB_CORROB_FLOOR='$GB_CORROB_FLOOR' (want a decimal like 0.5)" >&2
    exit 1;;
esac
# Admission threshold on the additive viral score, which is normalised to 0-1 (see the
# ADMISSION SCORE block below). Same rationale as GB_CORROB_FLOOR for living in the
# environment: the reader is baked into the .sif, the VALUE stays retunable without a rebuild.
#
# The caller passes a LEVEL, not a number. Levels are CALIBRATED against the within-library
# benchmark -- 6,024 WGS contigs from 8 libraries, 33.3% viral in every one, so the score
# cannot win by detecting which library a contig came from. That is the production-shaped
# question: one mixed metagenome, decide per contig.
#
#   level       thr     precision   recall   keeps
#   sensitive   0.15       46.5%     42.5%   30.5% of candidates
#   balanced    0.267      51.1%     26.8%   17.5%     <- default
#   strict      0.55       60.1%     14.5%    8.1%
#
# ⚠️ Those precisions are the HONEST ones. The same thresholds read 91.0% / 93.3% / 95.5% on
# the headline VLP-vs-WGS benchmark, but that set has one library per class, so its precision
# partly measures "which prep was this". Do not surface those numbers to a user.
#
# Levels are NOT percentiles and must not become percentiles: per-sample adaptive thresholds
# were measured and precision collapses to the base rate at every quantile, because forcing a
# fixed admission rate on a virus-poor sample admits junk. Fixed cuts are the mechanism.
#
# ============================== 2026-08-04 MODEL REPLACEMENT ==============================
# The stack4 48-column model was RETIRED. It had been fitted against CHGV VLP enrichment, and
# VLP enrichment is not a virality measurement -- it measures abundance and prep recovery:
#
#   geNomad virus_score vs the VLP label        AUC 0.4883 / 0.4916 / 0.5007  (chance)
#   geNomad virus_score vs a virality label     AUC 0.7276          (same contigs)
#   69.9% of gc3_content_max's apparent skill disappears under GC-matching
#
# Measured consequence of that fit, against CURATED phage genomes vs CURATED gut bacterial
# genomes -- a comparison involving no VLP, no UHGV and no CheckV:
#
#   frozen stack4 vs legacy/mid/novel phage   AUC 0.5776 / 0.5302 / 0.6160
#   geNomad 0.9689,  MetaPhaPred 0.9629,  DeepVirFinder 0.9557 on the SAME fragments
#
# It scored phage only 0.02-0.06 above bacterial chromosome fragments, and at this file's 0.60
# default admitted 65.5% of bacterial fragments. It was barely filtering.
#
# The replacement (stack5_final_model.json) is 9 signals, plain standardized logistic, no
# interaction block: geNomad aggregated virus/chromosome/plasmid, mpp, geNomad marker
# n_uscg / n_virus_hallmarks / marker_enrichment_v, and dvf_pca_01/03. The marker branch is
# where geNomad's published MCC 95.3% comes from and this pipeline was already computing it on
# every run and throwing it away. No extra compute.
#
#   model                        CHGV virality   curated novel
#   frozen stack4 (retired)          0.6019      0.5776/0.5302/0.6160
#   shipped 12 features, refit       0.8697          0.9794
#   THIS (9 signals)                 0.8864          0.9877
#
# ============ 2026-08-05: + VIRAL-VS-CELLULAR CONTRAST (stack6_final_model.json) ============
# The 9-signal model above is superseded by the SAME 9 signals plus the 4 contrast features (see
# the VIRAL-VS-CELLULAR SEQUENCE CONTRAST block). Still a plain standardized logistic, still no
# interaction block -- 13 features.
#
#   model                          CHGV virality   curated novel (held out)
#   frozen stack4 (retired)            0.6019       0.5776/0.5302/0.6160
#   stack5, 9 signals                  0.8864           0.9874
#   THIS (13, + contrast)              0.9099           0.9893
#
# The gain is +0.0235 CHGV AUC over stack5, paired donor bootstrap 95% CI [+0.0196,+0.0275].
# The 4-feature quartet beat vc_mean_delta alone by +0.0045, CI [+0.0030,+0.0058], which is what
# justifies carrying 4 columns instead of 1.
#
# Calibrated on honest out-of-fold scores, CHGV grouped by donor (n=14,916, 45 donors,
# prevalence 19.55%) -- positives = contigs >=95% identical over >=80% of their length to a
# UHGV gut virus, negatives = contigs whose k-mers concentrate on one phage-free gut bacterial
# genome:
#
#   level       thr     precision   recall   keeps      (stack5 at ITS cut, same folds+rows)
#   sensitive   0.45      77.98%     62.07%   15.56%     74.74% / 58.33%  @0.45
#   balanced    0.60      85.61%     52.02%   11.88%     83.05% / 46.71%  @0.60   <- default
#   strict      0.80      93.04%     36.66%    7.70%     91.78% / 32.92%  @0.75
#
# Same cuts on the CURATED benchmark (NCBI phage genomes released 2023-2026, i.e. after every
# tool's training cutoff, vs 396 RefSeq complete gut bacterial genomes; grouped by source
# genome so no genome straddles a fold):
#
#   level       thr     precision   recall
#   sensitive   0.45      96.50%     96.50%
#   balanced    0.60      97.28%     95.53%
#   strict      0.80      98.02%     93.05%
#
# STRICT MOVED 0.75 -> 0.80, the other two are unchanged in value. The cuts are always
# re-derived rather than carried over -- a different model puts a different probability
# distribution under the same numbers -- and at 0.75 the new model would have traded 1.2 points
# of precision for 7.8 of recall. At 0.80 every level beats the outgoing model on BOTH precision
# and recall, so no user sees a regression on either axis at any strictness.
#
# ⚠️ The two tables differ enormously and the curated one is an UPPER BOUND, not a promise.
# Curated reference genomes are cleaner than assembly output, and the CHGV table's own labels
# are imperfect in both directions: its positives come from UHGV, which was built partly with
# CheckV (one of this pipeline's own stage-1 miners), and its "bacterial" negatives certainly
# contain novel phage that no database recognises. Quote the CHGV column to users.
#
# ============ 2026-08-05: + FREE CPU COLUMNS, GRADIENT-BOOSTED (stack7_final_model.json) =======
# Same 13 signals plus 4 columns that fall out of the contrast's OWN .m8 files at zero extra cost
# (dark_frac, mean_bits_v, coding_density, n_genes), and a gradient-boosted model in place of the
# logistic one. No new database, no new search, no GPU. 17 features.
#
#   model                          CHGV virality   curated novel (held out)
#   stack5, 9 signals                  0.8864           0.9874
#   stack6, 13 + contrast, linear      0.9099           0.9893
#   THIS (17, + free4, gbm)            0.9294           0.9903
#
# +0.0195 CHGV over stack6. Measured separately, the two ingredients stack rather than overlap:
# free4 contributes +0.0075 (nested CV, CI [+0.0054,+0.0096]) and the GBM +0.0116 (CI
# [+0.0088,+0.0144]), and the GBM's margin was identical on 13 and on 17 features.
#
# coding_density REQUIRES real contig lengths (--fasta): that definition reproduces the training
# column 100.00% of the time, max-gene-end only 27.16%.
#
#   level       thr     precision   recall   keeps      (stack6 at ITS cut, same folds+rows)
#   sensitive   0.45      78.28%     65.74%   16.42%     77.73% / 61.90%  @0.45
#   balanced    0.60      85.86%     55.38%   12.61%     85.80% / 52.23%  @0.60   <- default
#   strict      0.80      94.13%     41.22%    8.56%     92.91% / 36.87%  @0.80
#
# CUTS ARE UNCHANGED at 0.45/0.60/0.80 -- re-derived from scratch as always, and this time the
# same numbers cleared the bar. Every level beats stack6 on BOTH precision and recall, including
# strict, where stack6's own upgrade had to move 0.75 -> 0.80 to avoid a precision regression.
#
# Same cuts on CURATED (held out entirely -- the model is fitted on CHGV):
#
#   level       thr     precision   recall      (stack6 at the same cut)
#   sensitive   0.45      88.07%     98.92%     96.50% / 96.50%
#   balanced    0.60      92.93%     97.95%     97.28% / 95.53%
#   strict      0.80      96.75%     94.87%     98.02% / 93.05%
#
# ⚠️ On CURATED, stack7 trades precision for recall at fixed thresholds (-8.4 points of precision
# at sensitive, -1.3 at strict) even though its curated AUC is BETTER (0.9903 vs 0.9893). The
# ranking improved; the threshold placement moved. Cause: the cuts are calibrated on CHGV at
# 19.55% prevalence and curated runs at 42.9%, and a GBM's probabilities do not transfer across
# that gap the way the logistic's did. Calibrating on CHGV is deliberate -- it is the
# production-shaped benchmark -- but anyone quoting curated numbers must use the table above and
# not assume stack6's.
#
# Production parity before shipping: the extractor's columns reproduce the training columns
# (n_genes and coding_density 100.00% exact; the alignment-dependent ones to within the known
# 0.27% MMseqs2 run difference), and scoring both feature paths through the shipped pure-Python
# tree walker gave ZERO admission flips across all 14,916 contigs at all three cuts.
#
# ============ 2026-08-07: + PHROGs + km4 RAW (stack9_final_model.json) =========================
# 168 features: stack7's 17, plus 15 PHROGs profile-family columns and 136 RAW canonical
# tetranucleotide frequencies. Both blocks were REJECTED twice before shipping, and both times the
# rejection was a THRESHOLD artifact rather than a merit one -- see the cut note below.
#
#   model                          CHGV virality   curated novel (held out)
#   stack7, 17, gbm                    0.9284           0.9945
#   stack8, +PHROGs (32), gbm          0.9349           0.9952
#   THIS (168, +PHROGs +km4, gbm)      0.9391           0.9965
#
# The two ingredients stack: PHROGs +0.0065, km4 RAW +0.0046 on top of it. km4 is RAW per-k-mer
# frequency, NOT PCA -- the 2026-08-02 k-mer rejection was wrong twice (it scored against the VLP
# label invalidated the next day, and PCA-15'd each block before any model saw it). It saturates at
# k=4: k=5 and k=6 buy +0.0003, inside CI overlap. Mechanism is short-motif avoidance under host
# defence (restriction-modification 4-6 bp palindromes, CRISPR PAMs 2-6 bp), which is why 13
# restriction-enzyme densities, 5 palindrome summaries, a 29-column extended-palindrome panel and a
# 66-column physicochemical panel all added NOTHING on top: km4 is a superset, not a competitor.
#
#   level       thr     precision   recall   keeps      (stack7 at ITS cut, same folds+rows)
#   sensitive   0.45      78.63%     68.76%   17.10%     78.13% / 65.78%  @0.45
#   balanced    0.61      86.35%     58.57%   13.26%     86.19% / 55.45%  @0.60   <- default
#   strict      0.80      94.30%     42.56%    8.82%     94.03% / 41.08%  @0.80
#
# ⚠️ BALANCED MOVED 0.60 -> 0.61; the other two are unchanged. This is the whole reason PHROGs and
# km4 looked like failures twice. At the OLD cuts stack9 REGRESSES at balanced (86.04% vs 86.19%
# precision) -- shipping the artifact without moving the cut would have been a straight production
# regression. One hundredth of a threshold converts that into +0.16pp precision AND +3.12pp recall.
# Every level now beats stack7 on BOTH axes, which is the standing admissibility rule.
#
# ⚠️ The table above selects the threshold on the same out-of-fold scores it reports, so its
# absolute recall is optimistic. Measured honestly by nested CV (thresh_rederive.py: threshold
# chosen inside each outer training fold, applied to untouched test rows), in-sample recall at
# strict is overstated by 3.6-4.7pp. The ADMISSIBILITY comparison is unaffected -- both models get
# identical folds and identical treatment -- but the number to expect in production is the nested
# one: at stack7's own precision, stack9 buys +6.41 / +5.62 / +1.99pp recall, all significant.
#
# COST: PHROGs adds ~6 min/run (38,880 HMM profiles over geNomad's own proteins, no new gene
# calling) and a 1.3 GB .h3m in the DB root. km4 adds 1.8s. Both degrade to the frozen training
# median if their pass fails, so neither can empty the admission score.
GB_STRICTNESS="${GB_STRICTNESS:-balanced}"
if [ "$_stack4_scored" -gt 0 ]; then
  case "$GB_STRICTNESS" in
    # stack158 cuts, PREVALENCE-CALIBRATED 2026-08-14. These hold each level's precision PROMISE on
    # a real gut assembly instead of on the enriched development set. Product decision, taken
    # explicitly: a user is better served by fewer, more trustworthy calls than by a longer list.
    #
    # WHY THEY MOVED. Thresholds are inherited across a model swap as PRECISION TARGETS (0.7969 /
    # 0.8710 / 0.9453, read off stack10et's OOF at its installed cuts). Matching those targets on
    # CHGV put stack158 at 0.4138 / 0.5200 / 0.6562 -- but CHGV is 19.6% viral and a real assembly
    # is ~3% (UHGV-HQ label). Precision is NOT prevalence-invariant, so the same cut on real data
    # delivers far less precision than the level advertises: on Shkoporov (10 subjects, bulk
    # assemblies, prevalence 0.0304) strict measured 0.4656 precision against a 0.9461 promise,
    # even though its AUC (0.9698) EXCEEDS CHGV's (0.9577). Discrimination did not degrade; the base
    # rate changed. Resampling CHGV to 0.0304 reproduces 52% of that drop with the model, the cuts
    # and the collection all held fixed (qov_paper §53b) -- so half the loss is arithmetic, not
    # difficulty. The values below are where stack158 reaches each target AT 0.0304 prevalence.
    #
    # MEASURED EFFECT on real bulk assemblies (Shkoporov, 10 subjects, 48,304 contigs, prev 0.0304,
    # UHGV-HQ label) -- these are the numbers to quote, not the CHGV-resampled projection they were
    # derived from:
    #   level       cut 0.4138/0.5200/0.6562  ->  cut 0.7200/0.7850/0.9363
    #   sensitive     0.2440 prec / 0.9094 rec  ->  0.5428 prec / 0.7221 rec
    #   balanced      0.3270 prec / 0.8535 rec  ->  0.6149 prec / 0.6669 rec   <- default
    #   strict        0.4656 prec / 0.7643 rec  ->  0.6749 prec / 0.4074 rec
    # Balanced nearly DOUBLES precision (+28.8pp) for -18.7pp recall. At balanced and strict no
    # threshold of geNomad, CheckV, MetaPhaPred, contrast or PHROGs reaches these precisions at all.
    #
    # ⚠️ THE CALIBRATION UNDERSHOOTS ITS OWN TARGET AND THAT IS EXPECTED. It was derived by
    # resampling CHGV to 0.0304, which corrects ONLY the base-rate half of the gap -- §53b measured
    # base rate at 52% of the drop and genuine collection difficulty at 48%. So the derivation
    # PROJECTED ~0.87 precision / ~0.41 recall at balanced; real data delivers 0.61 / 0.67, i.e.
    # short on precision and much better on recall than projected. No available cut reaches the
    # 0.8710 target on Shkoporov -- even strict tops out at 0.6749. Do not chase the target by
    # raising the cut further; the residual gap is collection difficulty, not miscalibration.
    #
    # ⚠️ THIS IS A CEILING-SIDE CALIBRATION, i.e. deliberately conservative. The Shkoporov precision
    # it is anchored to is measured against a UHGV-HOMOLOGY label, so a genuinely novel virus with
    # no UHGV match is scored as a false positive. True precision at the old cuts was therefore
    # HIGHER than 0.4656 by an unmeasured amount, and these cuts are correspondingly tighter than
    # strictly necessary. Erring toward precision is the intended direction; do not "correct" this
    # by loosening without a label that can actually see novel viruses.
    # ⚠️ GB_VSCORE_MIN still overrides, and the pre-2026-08-14 behaviour is exactly
    # GB_VSCORE_MIN=0.4138 / 0.5200 / 0.6562 -- use that to reproduce any run scored before today.
    sensitive) _vsmin_default=0.7200 ;;
    balanced)  _vsmin_default=0.7850 ;;
    strict)    _vsmin_default=0.9363 ;;
    *)
      echo "[5_evaluation] invalid GB_STRICTNESS='$GB_STRICTNESS' (want sensitive|balanced|strict)" >&2
      exit 1;;
  esac
else
  case "$GB_STRICTNESS" in
    sensitive) _vsmin_default=0.15  ;;
    balanced)  _vsmin_default=0.267 ;;
    strict)    _vsmin_default=0.55  ;;
    *)
      echo "[5_evaluation] invalid GB_STRICTNESS='$GB_STRICTNESS' (want sensitive|balanced|strict)" >&2
      exit 1;;
  esac
fi
# An explicit GB_VSCORE_MIN still wins, so a run can be pinned to an exact cut for
# reproducibility after the levels are ever recalibrated.
GB_VSCORE_MIN="${GB_VSCORE_MIN:-$_vsmin_default}"
case "$GB_VSCORE_MIN" in
  ''|*[!0-9.]*|*.*.*)
    echo "[5_evaluation] invalid GB_VSCORE_MIN='$GB_VSCORE_MIN' (want a decimal like 0.267)" >&2
    exit 1;;
esac

awk -F'\t' -v floor="$GB_CORROB_FLOOR" -f $GB_SCRIPTS/pbbcalc.awk $work/*_vscores.tsv > $work/pbb_viruses.tsv

# pbbcalc emits ID,dvf,vs2,mpp,gnmd,confidence,n_arms,n_comp,n_hom — 9 columns, not 10:
# the phamer column is gone (see the PHABOX2 block above), so confidence is $6 (was $7),
# n_comp is $8 (was $9) and n_hom is $9 (was $10).
#
# The gate is confidence AND corroboration. Confidence alone barely discriminates:
# 67.9% of the contigs the old gate admitted sat at exactly 1.00, because noisy-OR is
# 1-PROD(1-s) and an absent arm is its identity element, so one tool shouting 0.99 into
# a silent room scores identically to five tools agreeing.
#
# What this rejects is precisely one class: contigs whose ONLY support is a single
# composition net (sole-dvf or sole-mpp) -- the one case where that net's own near-twin
# actively dissented, since dvf and mpp agree 88.7% of the time. A lone geNomad or
# VirSorter2 hit still passes: those are separate modalities, and a long contig
# carrying only protein-homology support is the signature of a divergent phage with no
# close reference, which is what this pipeline exists to find. Measured over 11.7k
# previously-admitted contigs from 52 runs, this removes ~5%.
#
# NR==1 is handled EXPLICITLY, not by accident. The original confidence-only test let the
# header through only because awk compared the string "confidence" against 0.9
# lexicographically -- a coincidence, not a design. A
# boolean corroboration test is false on the header and would drop it -- after which the
# join below mistakes the first data row for the header, and extract_preview.py's
# header.index("h_name") stops resolving, so n_with_predicted_host quietly reports 0
# with no error anywhere.
# ==============================================
# ADMISSION SCORE
# ==============================================
# The admission score is `stack10_final_model.json`, applied by score_stack4.py -- see STACK4
# FEATURE EXTRACTION above. A gradient-boosted tree ensemble over 176 features (168 + CheckV 8):
#
#   geNomad aggregated   gn_virus, gn_chr, gn_plasmid
#   geNomad markers      n_uscg, n_virus_hallmarks, marker_enrichment_v
#   composition nets     mpp, dvf_pca_01, dvf_pca_03
#   contrast             vc_mean_delta, vc_delta_p90, vc_delta_max, vc_frac
#   contrast, free       dark_frac, mean_bits_v, coding_density, n_genes
#   PHROGs               15 profile-family columns
#   km4 RAW              136 canonical tetranucleotide frequencies
#   CheckV census        ckv_viral_genes, ckv_host_genes, ckv_gene_count, ckv_viral_frac
#   CheckV aai           ckv_aai_completeness, ckv_aai_num_hits, ckv_aai_id, ckv_aai_af
#
# Fitted on CHGV against a VIRALITY label and frozen. CHGV OOF 0.9524, curated held out 0.9920
# (stack9 was 0.9391 / 0.9888; paired donor bootstrap +0.0137, 95% CI [+0.0109,+0.0166]).
#
# ⚠️ THE CUTS DID NOT MOVE, and that is a measured choice rather than an omission. stack10 at
# 0.45/0.61/0.80 beats stack9 on BOTH precision and recall on CHGV -- the population that
# resembles production -- at all three levels. On curated it TRADES instead (strict: precision
# 96.49% -> 94.32%, recall 95.43% -> 97.77%), because a stronger model at a fixed threshold is
# simply more permissive. Curated is fragmented isolate genomes and a model fitted on it learns
# provenance rather than virality, so CHGV governs. Re-deriving the cuts to match curated's
# precision would give back most of the recall this model was shipped for.
#
# ⚠️ The script name, the DB directory (stack4/) and the "stack4" variable names throughout this
# file identify the deployed model LINEAGE, not the design -- four models have now shipped
# through this path. What actually loads is named on the --model line above; read that, not the
# variable names.
#
# ⚠️ REMOVED 2026-08-07: the two fallback tiers that used to live here (a hand-weighted 4-term
# formula `(f + 0.5*(1-k2_bact_frac) + 0.3*dvf + 0.3*mpp)/2.1`, and a noisy-OR corroboration
# gate). Their full measurement history -- why those four terms, why f=0.75 encodes CheckV
# silence, why assembly coverage was rejected -- is preserved in BENCHMARKS.md rather than here,
# because the code they justified is gone. See the block below for why they were removed.
# ONE MODEL OR AN ERROR (2026-08-07). The two fallback tiers -- the legacy 4-term formula and
# the noisy-OR gate -- were REMOVED. What they actually did was convert a plumbing bug into a
# silent model substitution.
#
# The trigger was never "a tool failed". Every tool DB is fail-fast checked in run.sh, so a
# missing reference errors before the container starts. The one time this path fired for real,
# geNomad had SUCCEEDED: a non-deterministic `find` picked the provirus classification file, the
# join dropped every contig on its `!(id in nnc)` test, and score_stack4.py returned scored=0 --
# "a silently empty admission score for the whole run, with no error". The fallback's role in
# that incident was to emit plausible-looking admissions from a different formula instead of
# surfacing the bug. That is the opposite of a safety net.
#
# The removed tiers were also never re-derived against the current benchmarks: their 0.15/0.267/
# 0.55 cuts come from the legacy Hannigan set, whose headline numbers this file elsewhere says
# "partly measure 'which prep was this' -- never surface those to a user", while the shipping
# rule for every stack model is CUTS ARE RE-DERIVED, NEVER CARRIED OVER.
#
# A wholesale geNomad/DVF/mpp failure now ERRORS. That is the honest outcome: the contrast
# quartet, dark_frac, mean_bits_v, coding_density, n_genes AND PHROGs all consume geNomad's
# proteins, so a run in that state has lost most of the model regardless. PER-CONTIG gaps are
# unaffected and still degrade gracefully through the artifact's impute_median.
if [ "$_stack4_scored" -le 0 ]; then
  echo "[5_evaluation] FATAL: the admission model scored 0 of $(awk 'END{print NR-1}' "$stack4_dir/joined_features.tsv" 2>/dev/null || echo 0) candidates." >&2
  echo "[5_evaluation]   This means a REQUIRED input failed for every contig -- geNomad, DeepVirFinder" >&2
  echo "[5_evaluation]   or MetaPhaPred -- or the stack4 join dropped everything. Check the WARN lines" >&2
  echo "[5_evaluation]   above and '$stack4_dir/join_dropped.count'. Not falling back to a weaker model." >&2
  exit 1
fi
echo "[5_evaluation] gate: strictness=$GB_STRICTNESS, score >= $GB_VSCORE_MIN on 0-1" \
     "(stack158: geNomad agg+markers, mpp, viral-vs-cellular contrast x4, dark_frac, mean_bits_v," \
     "coding_density, 136 km4 RAW, 8 CheckV -- extremely randomized trees, 158 features)"

# ARGV[1]/ARGV[2] are matched by FILENAME rather than NR==FNR because either may legitimately
# be empty, and NR==FNR would then eat the first record of pbb_viruses.tsv -- its header.
# vscore.tsv is the ONLY place the formula is evaluated. The deliverable's viral_score and
# bact_frac columns are joined from this file rather than recomputed downstream, so the
# numbers a user reads are byte-identical to the ones that decided admission.
: > "$work/vscore.tsv"
awk -F'\t' -v vmin="$GB_VSCORE_MIN" -v vsf="$work/vscore.tsv" '
FILENAME==ARGV[1] { if (FNR>1 && $2!="") s4[$1]=$2+0; next }  # vscore_stack4.tsv: id, viral_score
FNR==1 { print; next }
{
  # score_stack4.py already excludes contigs missing a non-imputable input, so a missing lookup
  # here means "not scored", not "score of 0" -- skip it, do not silently reject it at the gate
  # as if it had genuinely scored low.
  if (!($1 in s4)) { unscored++; next }
  s = s4[$1]
  # Third column is bact_frac, retained for the deliverable join and deliberately EMPTY: it was
  # only ever populated by the removed 4-term fallback, so it has been blank in every stack-scored
  # production run. Kept rather than dropped so the deliverable column set does not change.
  if (s >= vmin+0) { printf "%s\t%.4f\t%s\n", $1, s, "" > vsf; print }
}
END {
  if (unscored+0 > 0) print "[5_evaluation] stack9 gate: " unscored+0 " candidate(s) had no stack9 score (excluded, not rejected -- see the join count above)" > "/dev/stderr"
}
' "$stack4_dir/vscore_stack4.tsv" "$work/pbb_viruses.tsv" > $work/table.tsv

# ids.txt is CUT from the already-gated table instead of re-evaluating the gate a second
# time, so the id list and the row list cannot drift apart. Behaviour is unchanged,
# including the header line "ID" landing in ids.txt (harmless: seqkit grep -nf matches
# nothing on it, and the grep -v chain below passes it through).
cut -f1 $work/table.tsv > $work/ids.txt

gb_step "extracting admitted viral contigs + building the deliverable table"
seqkit grep -w 0 -nf $work/ids.txt $inputF > $work/busted_viruses.fna

awk 'BEGIN{OFS="\t"} /^>/{id=substr($0,2); next} {print id, length($0)}' "$work/busted_viruses.fna" > "$work/summary_contigs.tsv"

# FILENAME==ARGV[1], not NR==FNR: if every contig were rejected, summary_contigs.tsv is
# empty and NR==FNR would then be true for the FIRST record of table.tsv, eating its
# header and emitting a 0-byte, headerless deliverable that extract_preview.py cannot
# parse. Matching on the filename makes the file discriminator explicit.
# cname["ID"] is seeded in BEGIN, not in the file-1 rule: if every contig were rejected
# then summary_contigs.tsv is empty, that rule never fires, and column 2 of the
# deliverable header would come out as an empty string instead of "length".
awk -F'\t' 'BEGIN{OFS="\t";cname["ID"]="length"}
FILENAME==ARGV[1]{cname[$1]=$2;next}
# ONLY the columns that still mean something are carried into the deliverable: $2 (dvf) and
# $4 (mpp), the two per-tool terms of the admission score. Everything else pbbcalc emits is
# dropped here.
#
# WHAT WAS DROPPED AND WHY. vs2 ($3) is populated only for the slice stage 4 cached, so it
# reads 0.00 on contigs VirSorter2 never examined -- a value, not a blank, which a reader
# takes as a measured "not viral". gnmd ($5) is now 0.00 for EVERY contig since geNomad no
# longer runs in this stage. confidence ($6) is the noisy-OR over all four arms, so with two
# of them degraded it is a number computed from data no longer shown; it also no longer
# gates anything. n_arms/n_comp/n_hom ($7..$9) are the corroboration counts for a gate that
# does not exist any more.
#
# They are still COMPUTED -- pbb_viruses.tsv keeps all nine columns -- but as of 2026-08-07
# NOTHING reads $6/$8/$9 any more: the legacy fallback that did was removed. They are inert
# columns kept so pbbcalc.awk and this file stay version-matched.
{if(FNR==1){print $1,cname[$1],$2,$4;next};printf "%s\t%d\t%0.2f\t%0.2f\n",$1,cname[$1],$2,$4}
' "$work/summary_contigs.tsv" "$work/table.tsv" > "$work/final_table.tsv"

# =============================================
# Mapping Host
# =============================================

# Only free viral contigs are classified here. Prophage-derived contigs (_ckv_ /
# _phage / _gnmd_ suffixes) already carry a host from stage 1, called on the SOURCE
# bacterial contig — classifying the excised phage DNA instead would overwrite that
# with a weaker call, since the join below strips the suffix and both keys collapse
# to the same contig id.
grep -v "_ckv_" "$work/ids.txt" | grep -v "_phage" | grep -v "_gnmd_" > "$work/noproids.txt" || true
kraken_out="$work/k2.kreads"
: > "$kraken_out"

# SUBSET of stage 1's single HOST pass (HRGMv2, --confidence 0.2), which classified the whole
# length-filtered assembly. NOT a subset of the candidate pass above: that one is a different
# database at a different confidence, and mixing them would put NCBI names at 0.5 in the same
# deliverable column as stage 1's GTDB names at 0.2.
#
# This is also where a PRE-EXISTING inconsistency was fixed (2026-08-06): the old code
# subsetted the --confidence 0.5 candidate pass while stage 1's host call has always used 0.2,
# so the two halves of one column were called under different thresholds. Both are 0.2 now.
#
# ids.txt carries its "ID" header line through to noproids.txt — harmless, no kreads row has
# that sequence id.
_k2_host_all="$GB_OUTDIR/prophage/k2_host_all.kreads"
if [ -s "$work/noproids.txt" ] && [ -s "$_k2_host_all" ]; then
  awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($2 in keep)' \
    "$work/noproids.txt" "$_k2_host_all" > "$kraken_out"
fi

# parsekr.awk emits no header (the sibling call in 1_prophage.sh reads it directly),
# so write straight to the host table — a `tail -n +2` here silently dropped the
# first real host assignment of every run.
#
# With GB_KRAKEN=0 both halves are empty and h_rank/h_name/h_score come out blank for every
# contig. That is the ONE thing disabling kraken2 costs; the join below already tolerates an
# empty full_host_kr.tsv (see its NR==FNR note), and extract_preview.py then reports
# n_with_predicted_host = 0 rather than failing.
: > "$work/notpro_host_kr.tsv"
if [ "$GB_KRAKEN" = "1" ] && [ -s "$kraken_out" ]; then
  gb_step "Kraken2 — assigning species-level bacterial host"
  gawk -F'\t' -f ${GB_SCRIPTS}/parsekr.awk $KRAKEN_HOST_DB/nodes.dmp $KRAKEN_HOST_DB/names.dmp $kraken_out | grep -v "unclassified" > "$work/notpro_host_kr.tsv" || true
fi
cat "$GB_OUTDIR/prophage/full_host_kr.tsv" "$work/notpro_host_kr.tsv" > "$work/full_host_kr.tsv"

# ============================================
# Final steps
# ===========================================



awk -F'\t' '
BEGIN{OFS="\t"}

# file1: full_host_kr.tsv  (ID, taxid, rank, name, score)
# FILENAME==ARGV[1] rather than NR==FNR for the same reason as the join above: a sample
# whose every contig is prophage-derived leaves full_host_kr.tsv empty, and NR==FNR
# would then consume the header of final_table.tsv as a host record and strip it.
FILENAME==ARGV[1]{
  hr[$1]=$3
  hn[$1]=$4
  hs[$1]=$5
  next
}

# file2: final_table.tsv
# NB: NO APOSTROPHES in comments inside this block -- it is a SINGLE-QUOTED awk program, so
# a stray single-quote character closes the quote and bash then chokes on the next paren.
# The explicit $1..$10 must match the final_table width EXACTLY (10 since phamer left; it
# was 11). One too many emits an empty field here, which shifts h_rank/h_name/h_score one
# place right in the HEADER while the data rows below still write them at 11..13 --
# extract_preview.py resolves them by header.index(), so it would then read the wrong
# column and report empty hosts with no error raised anywhere.
FNR==1 {if($1=="ID"){print $1,$2,$3,$4,"h_rank","h_name","h_score"}else{print}; next }

{
  id=$1
  if (match(id, /^(c_[0-9]+)_/)) {
    id=substr(id, 1, RLENGTH-1)     # c_1234
  }

  # final_table columns: 1=ID, 2=length, 3=dvf, 4=mpp; 5..7 are the host fields appended
  # here (may be empty). The explicit index list above and these three assignments must
  # ALWAYS move together -- writing the host block one place right of the header leaves a
  # blank column and shifts h_rank/h_name/h_score, and extract_preview.py resolves h_name
  # by header.index(), so it would silently read the wrong column and report empty hosts.
  if (id in hr){$5=hr[id];$6=hn[id];$7=hs[id]}else{$5=$6=$7=""}


  print
}
' "$work/full_host_kr.tsv" "$work/final_table.tsv" > "$work/busted_table.tsv"

# ---- evidence tier ------------------------------------------------------------------
# Appended AFTER the host block so the fragile explicit $1..$10 list above stays untouched.
# The remap and dedup steps below are column-count agnostic (they `print` whole records), so
# columns 14..16 flow through to the deliverable unchanged.
# quality_summary.tsv is: 1=contig_id .. 6=viral_genes, 7=host_genes.
# FILENAME==ARGV[1] rather than NR==FNR: ckv_quality.tsv is legitimately EMPTY whenever
# CheckV was skipped or failed, and NR==FNR would then be true for the FIRST record of
# busted_table.tsv -- eating its header and shipping a headerless table.
awk -F'\t' '
BEGIN{OFS="\t"}
# ARGV[1] = vscore.tsv (id, viral_score, bact_frac) -- headerless, written by the gate awk.
# Empty whenever the legacy fallback ran, in which case both columns come out blank rather
# than fabricated.
FILENAME==ARGV[1] { vs[$1]=$2; bfr[$1]=$3; next }
FILENAME==ARGV[2] {
  if (FNR==1) next
  vg[$1]=$6+0; hg[$1]=$7+0; seen[$1]=1
  next
}
# NB: NO APOSTROPHES in comments inside this block -- it is a SINGLE-QUOTED awk program, so
# a stray single-quote character closes the quote and bash then chokes on the next paren.
#
# evidence_tier is GONE. It was a 3-way bucketing (strong/weak/unannotated) of exactly the
# data in viral_genes + host_genes, which are still columns here, so no information is lost.
# It was also actively worse than what replaced it: as a ranking of viralness it scores AUC
# 0.6199 against 0.7620 for viral_frac, and the two ORDER contigs differently -- a "weak"
# contig (viral_genes 0, host >= 1) gets f = 0.00 while an "unannotated" one gets 0.75, so
# the tier labels invert the ranking the score itself uses on the silent contigs. Shipping
# both invited a reader to trust the wrong one. Nothing consumed it either: GetViromePreview
# ignored the preview key, and the counter in extract_preview.py was the only reader.
FNR==1 { print $0, "viral_genes", "host_genes", "viral_frac", "bact_frac", "viral_score"; next }
{
  sc = ($1 in vs)  ? vs[$1]  : ""
  bfq = ($1 in bfr) ? bfr[$1] : ""
  if ($1 in seen) {
    v=vg[$1]; h=hg[$1]
    # viral_frac is the SAME f the admission score uses, surfaced so a reader can see why a
    # contig was kept. Silence (CheckV annotated neither) is 0.75, NOT 0 -- see the
    # ADMISSION SCORE block for why absence of evidence is not evidence of absence.
    printf "%s\t%d\t%d\t%.4f\t%s\t%s\n", $0, v, h, (v + h > 0 ? v / (v + h) : 0.75), bfq, sc
  } else {
    # CheckV emitted no row for this contig at all -- distinct from "annotated nothing".
    print $0, "", "", "", bfq, sc
  }
}
' "$work/vscore.tsv" "$work/ckv_quality.tsv" "$work/busted_table.tsv" > "$work/tiered_table.tsv"

# viral_frac is $10 in the slimmed layout: 1=ID 2=length 3=dvf 4=mpp 5..7=host 8=viral_genes
# 9=host_genes 10=viral_frac 11=bact_frac 12=viral_score.
_ckv_n=$(awk -F'\t' 'NR>1 && $10!="" {if($10+0>=0.5)a++; else if($10+0>0)b++; else c++} END{printf "frac>=0.5=%d 0<frac<0.5=%d frac=0=%d", a+0, b+0, c+0}' "$work/tiered_table.tsv")
echo "[5_evaluation] CheckV viral fraction: ${_ckv_n}"


awk -F'\t' '
BEGIN{OFS="\t"}

# map file: header_map.tsv (col1 old, col2 new)
NR==FNR{
  old=$1; new=$2
  sub(/^>/,"",old)
  sub(/^>/,"",new)
  cname[new]=old
  next
}

# table file: busted_table.tsv
FNR==1 { print; next }

{
  id=$1

  # split id into base + suffix: c_1234 or c_1234_something
  base=id
  suffix=""
  if (match(id, /^(c_[0-9]+)_/)) {
    base=substr(id, 1, RLENGTH-1)     # c_1234
    suffix=substr(id, RLENGTH+1)      # after "c_1234_"
  }

  # replace base if mapping exists
  if (base in cname) {
    new_id = cname[base]
    if (suffix != "") new_id = new_id "_" suffix
    $1 = new_id
  }

  print
}
' "$GB_OUTDIR/prophage/header_map.tsv" "$work/tiered_table.tsv" > "$work/remapped_table.tsv"

awk -F'\t' '
BEGIN { OFS="\t" }

# First file: map (col1 realname=k..., col2 tmp=c_#)
NR==FNR {
  real=$1
  tmp=$2
  sub(/^>/,"",real)
  sub(/^>/,"",tmp)
  map[tmp]=real          # c_# -> k...
  next
}

# Second file: fasta
/^>/ {
  hdr = substr($0, 2)    # drop >
  base = hdr
  suffix = ""

  # capture base c_<digits> and anything after it (including nothing)
  if (match(hdr, /^(c_[0-9]+)(.*)$/, m)) {
    base = m[1]          # c_123
    suffix = m[2]        # _phage... or _... or empty
  }

  if (base in map) {
    print ">" map[base] suffix
  } else {
    print $0
  }
  next
}

{ print }
' "$GB_OUTDIR/prophage/header_map.tsv" "$work/busted_viruses.fna" > "$work/remapped_viruses.fna"

# ==============================================
# DEDUP deliverables by contig name
# ==============================================
# Internal c_N ids are unique, but remapping them back to original contig names
# (above) can re-collide when the INPUT assembly itself contains repeated contig
# names (e.g. a co-assembly / concatenated assembly). Duplicate-named records make
# the FASTA unusable for read mapping (minimap2/samtools reject duplicate @SQ
# entries). Collapse by name, keeping the first occurrence; keep the table in step.
awk '
  /^>/ {
    name = $1
    sub(/^>/, "", name)
    if (seen[name]++) { skip = 1; next } else { skip = 0 }
  }
  skip == 0 { print }
' "$work/remapped_viruses.fna" > "$GB_OUTDIR/${GB_PREFIX}_viruses.fna"

awk -F'\t' 'NR==1 { print; next } !seen[$1]++' \
  "$work/remapped_table.tsv" > "$GB_OUTDIR/${GB_PREFIX}_table.tsv"

_nin=$(grep -c '^>' "$work/remapped_viruses.fna" 2>/dev/null || true)
_nout=$(grep -c '^>' "$GB_OUTDIR/${GB_PREFIX}_viruses.fna" 2>/dev/null || true)
echo "[5_evaluation] deliverable viral contigs: ${_nout:-0} (removed $(( ${_nin:-0} - ${_nout:-0} )) duplicate-named record(s))"

echo "[5_evaluation] Done :D"
