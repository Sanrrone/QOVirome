#!/usr/bin/env bash

set -euo pipefail
#export PATH=/opt/micromamba/envs/gutbusters/bin:$PATH
# ---------------- Env & inputs ----------------
: "${GB_OUTDIR:?Need GB_OUTDIR}"
: "${GB_THREADS:=4}"
: "${GB_DB_ROOT:?Need GB_DB_ROOT}"
: "${GB_IN:?Need GB_IN}"
: "${GB_MLEN:?Need GB_MLEN}" 
: "${GB_SCRIPTS:?Need GB_SCRIPTS}"

export TMPDIR=$GB_TMP
export GUTBUSTERS_TMP="$TMPDIR"
export TRANSFORMERS_CACHE="$GUTBUSTERS_TMP/hf"
export HF_HOME="$GUTBUSTERS_TMP/hf"
export TORCH_HOME="$GUTBUSTERS_TMP/torch"
export XDG_CACHE_HOME="$GUTBUSTERS_TMP/xdg"
export MPLCONFIGDIR="$GUTBUSTERS_TMP/matplotlib"
export XDG_CONFIG_HOME="$GB_TMP/xdg"

work="$GB_OUTDIR/prophage"
mkdir -p "$work"
mkdir -p "$GB_TMP"/{hf,torch,matplotlib,xdg}
chmod -R 775 "$GB_TMP" 2>/dev/null 

# DBs (override via env if paths differ)
UHGV_DB_PREFIX="${GB_UHGV_DB_PREFIX:-${GB_DB_ROOT}/blast/uhgv/uhgv}"
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
# HOST CALLS USE A DIFFERENT DATABASE FROM BACTERIAL TAGGING, deliberately (measured
# 2026-08-06, see BENCHMARKS.md "Migrating kraken2 hgdb_unphaged -> HRGMv2"):
#   hgdb_unphaged  31,149 strains /  1,593 species  -> better bacterial TAG (+50.7% vs +42.5%)
#   HRGMv2          4,824 reps    /  4,654 species  -> better HOST call (+67% species-level)
# hgdb_unphaged packs ~20 strains onto each species, which makes k-mers ambiguous across taxa
# and drives kraken2's LCA up the tree: 200 of its 654 calls landed at order or superkingdom,
# i.e. "it is a bacterium". That is a property of leaf granularity, not of taxonomy quality —
# hgdb_unphaged is built on genuine NCBI taxonomy and has 3 "Unknown" strings in 33,597 names.
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

echo "[1_prophage] IN         : $GB_IN"
echo "[1_prophage] OUTDIR     : $GB_OUTDIR"
echo "[1_prophage] THREADS    : $GB_THREADS"
echo "[1_prophage] UHGV DB    : $UHGV_DB_PREFIX"
echo "[1_prophage] Kraken DB  : $KRAKEN_DB"
echo "[1_prophage] MIN LENGTH : $GB_MLEN"

# Per-tool progress + failure localisation. Sourced, never executed. Absence must not break the
# run -- the stage's job is to mine viruses, not to report on itself -- so every entry point is
# stubbed out if the helper is missing.
if [ -r "${GB_SCRIPTS}/gb_progress.sh" ]; then
  . "${GB_SCRIPTS}/gb_progress.sh"
else
  gb_step() { echo "[1_prophage] >> $1"; }; gb_tool_result() { :; }
fi

# Helper functions
subset_ids() { # subset fasta by ID list
  local ids="$1" fa="$2" out="$3"
  [[ -s "$ids" ]] || { : > "$out"; return 0; }
  seqkit grep -f "$ids" "$fa" > "$out"
}
drop_ids() { # drop sequences whose IDs are in list
  local ids="$1" fa="$2" out="$3"
  [[ -s "$ids" ]] || { cp -f "$fa" "$out"; return 0; }
  seqkit grep -v -f "$ids" "$fa" > "$out"
}

#prefilter
# Filter by length
gb_step "length filter — dropping contigs below ${GB_MLEN} bp"
perl "${GB_SCRIPTS}/removesmalls.pl" "$GB_MLEN" "$GB_IN" > "${work}/filt_mag.fna"

# Map old headers -> new headers (c_1, c_2, ...)
# Keep only the contig id (first whitespace-delimited token): assembler headers
# like ">k141_51352 flag=1 multi=3.0000 len=2342" carry trailing metadata, but
# only ">k141_51352" is the valid sequence name. Capturing the full line here
# would propagate the spaces+metadata back into the final FASTA/table names.
grep -E '^>' "${work}/filt_mag.fna" | awk '{print $1}' > "${work}/ids_original.txt"

awk '
  BEGIN { i=0 }
  /^>/ { print ">c_" ++i; next }
  { print }
' "${work}/filt_mag.fna" > "${work}/filt_mag.fna.tmp" \
  && mv -f "${work}/filt_mag.fna.tmp" "${work}/filt_mag.fna"

grep -E '^>' "${work}/filt_mag.fna" > "${work}/ids_new.txt"

paste "${work}/ids_original.txt" "${work}/ids_new.txt" | sed "s/>//g" > "${work}/header_map.tsv"
rm -f "${work}/ids_original.txt" "${work}/ids_new.txt"

# Update inputs
GB_IN="${work}/filt_mag.fna"
NEXT_IN="${work}/filt_mag.fna"

# ==============================================
# 0) Kraken2 HOST pass — ONCE, for the whole run
# ==============================================
# kraken2 costs a ~19 GB database LOAD and almost nothing per contig: measured 2026-08-06 on
# Sun, load = 8-11 s while classifying 20,000 contigs on top of it adds 0-2 s. So the number
# of PASSES is the entire cost, and the input size is free.
#
# kraken2 classifies each sequence independently, so a row is identical whether the contig was
# submitted alone or inside a 20,000-contig file — verified here, not assumed: 200 contigs run
# as a subset came back BYTE-IDENTICAL to the same 200 rows from the full-set run. That makes
# one pass over every contig, subset later by id, exactly equivalent to per-stage passes.
#
# This pass serves BOTH host call sites (here on provirus parents, and 5_eval.sh on the
# non-prophage candidates); both are subsets of filt_mag.fna and both use --confidence 0.2.
#
# ⚠️ 5_eval.sh carries a "WHY HERE AND NOT BEFORE STAGE 1" note saying kraken2 must not run
# first. That objection was about kraken2 DISCARDING contigs (it cost geNomad 43.6% of its
# calls). Since 2026-08-05 kraken2 removes nothing, so classifying early is inert — it only
# writes a table. Do not read that note as forbidding this pass.
k2_host_all="$work/k2_host_all.kreads"
: > "$k2_host_all"
if [ "$GB_KRAKEN" = "1" ] && [ -s "$GB_IN" ]; then
  gb_step "Kraken2 — host classification pass"
  _t0=$SECONDS
  kraken2 \
    --db "$KRAKEN_HOST_DB" \
    --threads "$GB_THREADS" \
    --confidence 0.2 \
    --report /dev/null \
    --output "$k2_host_all" \
    "$GB_IN" \
    || echo "[1_prophage] WARN: host kraken2 failed — hosts unassigned this run" >&2
  # Fail-open by design (the `||` above): a dead kraken2 costs the h_* columns, not the run.
  gb_tool_result "Kraken2 host" "$( [ -s "$k2_host_all" ] && echo 0 || echo 1 )" \
    "$((SECONDS - _t0))" "no predicted host for any contig this run"
fi
echo "[1_prophage] Kraken2 host pass: $(wc -l < "$k2_host_all") contigs classified once for the whole run"


# ==============================================
# 1) BLAST vs UHGV (viral/proviral detection)
# ==============================================
echo "[1_prophage] Step 1: BLAST vs UHGV"

blast_tsv="$work/blast.tsv"
# OFF BY DEFAULT since 2026-08-09. The UHGV BLAST claimer identified contigs homologous to known
# gut viruses so they could bypass the later stages -- but claims never bypassed the stack9 gate
# (5_eval.sh scores every pooled contig against vmin), and every contig now reaches stage 5
# regardless, so the claim changes nothing about which contigs are admitted. Removing it drops a
# 4.5 GB database that feeds no stack9 feature.
# An empty blast.tsv leaves hit_ids empty; subset_ids then writes an empty after_blast_viral.fna
# and the remainder awk passes every contig through -- byte-identical in shape to a zero-hit BLAST
# run, so nothing downstream needs to know this happened.
# Re-enabling needs GB_BLAST=1 *and* blast/uhgv restored to run.sh's _dbs list.
if [ "${GB_BLAST:-0}" = "1" ]; then
  raw_blast_tsv="$work/blast_uhgv.tsv"
  blastn -query "$NEXT_IN" -db "$UHGV_DB_PREFIX" \
    -num_threads "$GB_THREADS" -evalue 0.001 \
    -perc_identity 70 -qcov_hsp_perc 70 \
    -outfmt "6 qseqid sseqid pident qlen slen length" \
    > "$raw_blast_tsv"

  echo "[1_prophage] BLAST table -> $raw_blast_tsv"

  awk -F"\t" '{if($3>=80 && $6/$4 >= 0.8)print $1"\t"$2"\t"$3"\t"($6/$4)*100}' ${raw_blast_tsv} |
          awk -F"\t" 'BEGIN{getline;n[$1]=$3*$4;v[$1]=$0}{if($1 in n){if(n[$1]<$3*$4){n[$1]=$3*$4;v[$1]=$0}}else{n[$1]=$3*$4;v[$1]=$0}}END{for(k in n)print v[k]}'> "$blast_tsv"
  rm -f ${raw_blast_tsv}
else
  echo "[1_prophage] BLAST claimer DISABLED (GB_BLAST=0) -- no UHGV claims"
  : > "$blast_tsv"
fi
hit_ids="$work/after_blast_hit.ids"
awk -F"\t" '{print $1}' $blast_tsv > $hit_ids
 
# Extract full hits
after_blast_viral="$work/after_blast_viral.fna"
subset_ids $hit_ids "$NEXT_IN" "$after_blast_viral"

# ---------------------------------------------------------

# Remaining for CheckV
after_blast_remaining="$work/after_blast_remaining.fna"
awk -v elist="$hit_ids" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $NEXT_IN > $after_blast_remaining

#drop_ids "$hit_ids" "$NEXT_IN" "$after_blast_remaining"
NEXT_IN="$after_blast_remaining"

echo "[1_prophage] after_blast_viral.fna     : $( [[ -s $after_blast_viral ]] && wc -l < $after_blast_viral || echo 0 ) lines"
echo "[1_prophage] after_blast_remaining.fna : $( [[ -s $NEXT_IN ]] && wc -l < $NEXT_IN || echo 0 ) lines"

# ==============================================
# 2) CheckV on remaining contigs
# ==============================================
echo "[1_prophage] Step 2: CheckV on after_blast_remaining.fna"

checkv_dir="$work/checkv"
mkdir -p "$checkv_dir"

# Edge case: BLAST may match every contig, leaving nothing for CheckV. CheckV
# aborts on an empty FASTA, so skip it and synthesize the empty outputs the
# downstream parsing expects.
if [[ -s "$NEXT_IN" ]]; then
  # NOT fail-open: CheckV aborting kills the stage, so no gb_tool_result here -- the entrypoint
  # reports it from the progress state instead.
  gb_step "CheckV — finding viruses, excising prophage regions"
  /opt/micromamba/envs/gutbusters/bin/checkv end_to_end "$NEXT_IN" "$checkv_dir" -t "$GB_THREADS" --remove_tmp
else
  echo "[1_prophage] after_blast_remaining.fna is empty — skipping CheckV (all contigs matched UHGV)"
  : > "$checkv_dir/viruses.fna"
  : > "$checkv_dir/proviruses.fna"
  : > "$checkv_dir/contamination.tsv"
fi

prev=$(pwd)
cd "$checkv_dir"
touch proviruses.fna
# grep exits 1 (no match) on an empty proviruses.fna; tolerate it so set -e/pipefail
# don't kill the run when CheckV found no proviruses.
grep ">" proviruses.fna | sed "s/>//g" | sort -u > proids.txt || true

awk -F'\t' -v lfilt=$GB_MLEN '{
if(NR==FNR){
        n[$1]=0
}else{
if($1 in n == 0){
	if($6=="Yes" && $2>lfilt){
       		 if($4==0 && $5==1){print;next}
       		 if($4>$5){
       		         print $1
       		 }else{
       		 	split($9,a,",");
      	 	 	n[a[1]]=$7;
       	 		n[a[2]]=$8;
       		 	if(n["viral"]>n["host"]){print $1}
       	 	}
	}
}
}}' proids.txt contamination.tsv > nocont.txt

# Contigs whose viral portion CheckV harvested — i.e. exactly what ends up in
# after_checkv_viral.fna, which is proviruses.fna. CheckV writes those headers as
# "<contig>_<region> <start>-<end>/<contig_len>", so strip the region index to get
# an id that exists in the assembly: this list both excludes them from the remaining
# pool below and drives the Kraken2 host lookup on their source contigs.
# (The previous field-shuffling awk assumed a 3-token header and dropped the id
# itself, so neither consumer ever matched anything.)
awk '{id=$1; sub(/_[0-9]+$/,"",id); print id}' proids.txt | sort -u > checkv_valids.txt
sed -i "/>/ s/ /_ckv_/g" proviruses.fna
sed -i "/>/ s/\//__/g" proviruses.fna
grep ">" proviruses.fna | sed "s/>//g" | sort -u > proids.txt || true

cd $prev
awk -v elist="$checkv_dir/checkv_valids.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $NEXT_IN > $work/after_checkv_remaining.fna

if [[ -s "$checkv_dir/nocont.txt" && -s "$checkv_dir/viruses.fna" ]]; then
  seqkit grep -w 0 -nf "$checkv_dir/nocont.txt" "$checkv_dir/viruses.fna" > "$checkv_dir/tmp.fna"
else
  : > "$checkv_dir/tmp.fna"
fi
cat $checkv_dir/proviruses.fna $checkv_dir/tmp.fna > $work/after_checkv_viral.fna

after_checkv_viral="$work/after_checkv_viral.fna"
after_checkv_remaining="$work/after_checkv_remaining.fna"

# --------------------------------------------------------------------

# Remaining for PhageBoost
NEXT_IN="$after_checkv_remaining"

# ==============================================
# 3) PhageBoost on remaining contigs
# ==============================================
echo "[1_prophage] Step 3: PhageBoost on after_checkv_remaining.fna"

pb_dir="$work/phageboost"
mkdir -p "$pb_dir"

# PhageBoost RETIRED 2026-08-10, default OFF. Two arms on curated 20,000, identical except this
# gate (BENCHMARKS.md, "DROP PhageBoost?"): precision 95.64% -> 95.65%, recall 96.10% -> 99.12%,
# wall -393s. It improved BOTH axes, which for a tool REMOVAL needed explaining -- and the
# explanation is the reason it is off: of its 405 delivered excisions, 389 (96%) cut up contigs
# that were ALREADY FULLY VIRAL, truncating correct answers into fragments, and only 16 (4%)
# rescued a prophage from a bacterial host. CheckV is the inverse (82% bacterial parents).
# Cause: -cs/--mincontigsize defaults to 10000 and was never overridden here, so PhageBoost only
# ever saw 3,045 of 19,825 contigs -- the long tail, which skews viral -- and its features (CAI,
# GC-skew, feature z-scores) are all computed file-wide, i.e. blended across every organism in a
# metagenome rather than the single isolate it was designed for.
# Nothing downstream needs to know this happened: pb_dir stays empty, so after_phageboost_viral.fna
# is empty, phageboost_exclusive.txt is empty, and the subtraction awk passes every contig through
# untouched -- the same shape the GB_BLAST=0 path already relies on.
# Re-enable with GB_PHAGEBOOST=1.
if [ "${GB_PHAGEBOOST:-0}" != "1" ]; then
  echo "[1_prophage] PhageBoost DISABLED (GB_PHAGEBOOST=0) -- CheckV is the only excision tool"
elif [[ -s "$NEXT_IN" ]]; then
  PhageBoost -f "$NEXT_IN" -j "$GB_THREADS" -o "$pb_dir" -meta 0 -t 0.95
else
  echo "[1_prophage] after_checkv_remaining.fna is empty — skipping PhageBoost"
fi

# Normalize PB summary
after_phageboost_viral="$work/after_phageboost_viral.fna"
after_phageboost_remaining="$work/after_phageboost_remaining.fna"
find $pb_dir -name "*.fasta" -exec cat {} \; > $after_phageboost_viral
perl ${GB_SCRIPTS}/removesmalls.pl $GB_MLEN $after_phageboost_viral > $work/tmp.fna
rm -f $after_phageboost_viral
mv $work/tmp.fna $after_phageboost_viral

# Empty when PhageBoost found no phages; grep exits 1, so tolerate it.
grep ">" $after_phageboost_viral | sed "s/>//g" > $work/phageboost_valids.txt || true
awk -F"_phage" '{print $1}' $work/phageboost_valids.txt | sort -u > $work/phageboost_exclusive.txt
awk -v elist="$work/phageboost_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $NEXT_IN > $after_phageboost_remaining

#drop_ids "$pb_ids" "$NEXT_IN" "$after_phageboost_remaining"
NEXT_IN="$after_phageboost_remaining"

# ==============================================
# 4) geNomad on what BLAST/CheckV/PhageBoost left — VIRAL BRANCH ONLY
# ==============================================
# Scope: this stage takes geNomad's virus output and NOTHING else. geNomad also
# classifies plasmids and annotates conjugation machinery, but those are mobile
# genetic elements owned by the qo-analysis-mge analysis; duplicating them here
# would double-report the same loci across two analyses. The split is symmetric:
# mge runs geNomad with --disable-find-proviruses and never parses
# *_virus_summary.tsv, so prophage/virus is ours and plasmid/MGE is theirs.
# Accordingly we leave provirus finding ENABLED (that is the point of this stage)
# and delete the plasmid outputs before anything downstream can glob them.
echo "[1_prophage] Step 4: geNomad (viral branch) on after_phageboost_remaining.fna"

gn_dir="$work/genomad"
mkdir -p "$gn_dir"
after_genomad_viral="$work/after_genomad_viral.fna"
after_genomad_remaining="$work/after_genomad_remaining.fna"
: > "$after_genomad_viral"

# MMseqs2 marker-search sensitivity. qo-analysis-mge runs 3.5 as a speed trade and
# its own note warns the recall loss "falls hardest on divergent gut phages" — the
# exact contigs this analysis exists to find — so we keep geNomad's default 4.2.
GENOMAD_SENS="${GB_GENOMAD_SENS:-4.2}"
GENOMAD_DB="${GB_GENOMAD_DB:-/opt/gutbusters/genomad_db}"

# geNomad aborts on an empty FASTA, like every other tool in this stage.
# CLAIMER OFF BY DEFAULT since 2026-08-09: geNomad runs ONCE, in stage 5, for the 6 features stack9
# reads. As a stage-1 claimer it excised almost nothing (0 of 20,025 contigs on curated, 1 on B13)
# and its whole-contig claims are redundant now that every contig reaches stage 5 -- claims do not
# bypass the gate. This removes a pass, not a signal.
if [ "${GB_GENOMAD_CLAIM:-0}" != "1" ]; then
  echo "[1_prophage] geNomad CLAIMER DISABLED (GB_GENOMAD_CLAIM=0) -- stage 5 still runs it for features"
elif [[ -s "$NEXT_IN" ]]; then
  gb_step "geNomad — claimer pass (GB_GENOMAD_CLAIM=1)"
  gnmd end-to-end --cleanup \
      --threads "$GB_THREADS" \
      --sensitivity "$GENOMAD_SENS" \
      "$NEXT_IN" "$gn_dir" "$GENOMAD_DB"

  # Enforce "viral only" in code, not just by convention: drop the plasmid arm so
  # no later glob can pick it up.
  find "$gn_dir" -name '*_plasmid*' -delete 2>/dev/null || true

  _gn_virus="$(find "$gn_dir" -maxdepth 3 -type f -name '*_virus.fna' -print -quit 2>/dev/null || true)"
  if [[ -n "$_gn_virus" && -s "$_gn_virus" ]]; then
    # geNomad names excised proviruses "<contig>|provirus_<start>_<end>". A pipe
    # breaks seqkit and VirSorter2's "||" suffix handling downstream, so fold it
    # into the same "_<tag>_" convention CheckV (_ckv_) and PhageBoost (_phage)
    # already use, keeping the c_N prefix intact for the stage-5 header remap.
    # The second sub() is defensive: no "|" may survive into the pipeline.
    awk '/^>/{sub(/\|provirus_/,"_gnmd_"); sub(/\|.*$/,"")} {print}' "$_gn_virus" \
      > "$work/genomad_viral_raw.fna"
    perl "${GB_SCRIPTS}/removesmalls.pl" "$GB_MLEN" "$work/genomad_viral_raw.fna" > "$after_genomad_viral"
    rm -f "$work/genomad_viral_raw.fna"
  else
    echo "[1_prophage] geNomad found no viral contigs"
  fi
else
  echo "[1_prophage] after_phageboost_remaining.fna is empty — skipping geNomad"
fi

: > "$work/genomad_hits.txt"
grep ">" "$after_genomad_viral" | sed "s/>//g" | awk '{print $1}' > "$work/genomad_hits.txt" || true

# Two lists, because the two consumers need different sets.
# (a) Subtraction from the remaining pool: EVERY contig geNomad claimed, whether it
#     was excised as a provirus or called viral whole — none may be re-offered
#     downstream.
sed -E 's/_gnmd_[0-9]+_[0-9]+$//' "$work/genomad_hits.txt" | sort -u > "$work/genomad_exclusive.txt"
# (b) Kraken2 host attribution: PROVIRUS PARENTS ONLY. A whole-contig viral call is
#     the virus itself, not a bacterial host, so classifying it here would attribute
#     a "host" from the lax --confidence 0.2 prophage pass. Every other free-virus
#     class (BLAST/UHGV, PhaMer, MetaPhaPred, DVF, VS2) gets its host from stage 5's
#     --confidence 0.5 pass instead, and geNomad's whole-contig calls must match them
#     — otherwise they alone would keep a 0.2-confidence host whenever the stricter
#     stage-5 pass returns unclassified. Mirrors checkv_valids.txt (proviruses only)
#     and phageboost_exclusive.txt (_phage regions only).
grep -E '_gnmd_[0-9]+_[0-9]+$' "$work/genomad_hits.txt" \
  | sed -E 's/_gnmd_[0-9]+_[0-9]+$//' | sort -u > "$work/genomad_provirus_parents.txt" || true
awk -v elist="$work/genomad_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' \
  "$NEXT_IN" > "$after_genomad_remaining"

NEXT_IN="$after_genomad_remaining"

# ===================================================== #
# 5) Kraken2 (last): getting host from extracted phages #
# ===================================================== #
echo "[1_prophage] Step 5: Kraken2 on prophages contig source to track bacterial host"

cat $checkv_dir/checkv_valids.txt $work/phageboost_exclusive.txt $work/genomad_provirus_parents.txt \
  | sort -u > $work/exclusive_ids.txt
subset_ids "$work/exclusive_ids.txt" "$GB_IN" "$work/prophage_source_tok2.fna"

kraken_out="$work/k2.kreads"

# SUBSET of the single host pass at the top of this stage, not a second kraken2 run (see the
# comment there). Field 2 of a kreads line is the sequence id.
: > "$kraken_out"
if [ -s "$work/exclusive_ids.txt" ] && [ -s "$k2_host_all" ]; then
  awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($2 in keep)' \
    "$work/exclusive_ids.txt" "$k2_host_all" > "$kraken_out"
fi

: > "$work/full_host_kr.tsv"
if [ "$GB_KRAKEN" = "1" ] && [ -s "$kraken_out" ]; then
  gb_step "Kraken2 — resolving host taxonomy names"
  gawk -F'\t' -f ${GB_SCRIPTS}/parsekr.awk $KRAKEN_HOST_DB/nodes.dmp $KRAKEN_HOST_DB/names.dmp $kraken_out | grep -v "unclassified" > "$work/full_host_kr.tsv" || true
fi

echo "[1_prophage] Kraken host output -> $kraken_out"

# ==============================================
# Merge viral pieces for pipeline output
# ==============================================
echo "[1_prophage] Merging viral pieces to pipeline hand-off"

# Collect all viral outputs that you want to hand over as "viruses.fasta"
viral_collect="$work/viruses_prophages.fna"
: > "$viral_collect"
[[ -s "$after_blast_viral"       ]] && cat "$after_blast_viral"       >> "$viral_collect" || echo "[1_prophage] Warning: $after_blast_viral is empty"
[[ -s "$after_checkv_viral"     ]] && cat "$after_checkv_viral"     >> "$viral_collect" || echo "[1_prophage] Warning: $after_checkv_viral is empty"
[[ -s "$after_phageboost_viral" ]] && cat "$after_phageboost_viral" >> "$viral_collect" || echo "[1_prophage] Warning: $after_phageboost_viral is empty"
[[ -s "$after_genomad_viral"    ]] && cat "$after_genomad_viral"    >> "$viral_collect" || echo "[1_prophage] Warning: $after_genomad_viral is empty"

echo "[1_prophage] Done."

