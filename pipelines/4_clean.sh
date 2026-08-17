#!/usr/bin/env bash
# Hand stage 5 the WHOLE assembly. Replaces stages 2, 3 and 4 (shipped 2026-08-09).
#
# Keeps the benchmark's name -- this is the `clean` arm of the five-arm comparison, shipped
# unchanged, so the numbers in BENCHMARKS.md describe exactly this file.
#
# WHY THE THREE STAGES WENT AWAY. Stages 2-4 existed to select which contigs deserved scoring:
# a bacterial filter, then two claimers that pulled out likely viruses. That made sense when the
# gate was a weak formula. It stopped making sense once stack9 became the decision-maker, because
# claims never bypassed the gate -- 5_eval.sh scores every pooled contig against vmin, so a claimed
# contig and an unclaimed one are treated identically. The stages were therefore doing PREVALENCE
# ENGINEERING (feeding stage 5 a pool that was 86.2% viral instead of 60.1%), not selection.
#
# MEASURED (curated 20,000; five arms at the shipped cuts):
#   precision -0.44 / -0.15 / -0.08 pp  and recall -0.38 / -0.53 / -0.72 pp
#   at sensitive / balanced / strict, against the 5-stage pipeline.
# Bought with: 5 stages -> 3, and -6.2 GB of container (blast/uhgv 4.5 GB + phabox_db_v2_1 1.7 GB).
# On a real donor assembly (B13, 6,032 contigs) it admitted 315 contigs against the 5-stage
# pipeline's 308, at +0.8% wall -- inside the noise of a contended host.
#
# Cuts were deliberately NOT re-derived. Holding precision constant instead costs -0.50 / -0.54 /
# -3.12 pp of recall, i.e. far worse at strict, and it would silently redefine what a user's
# chosen strictness level means. See BENCHMARKS.md, "clean to production".
set -euo pipefail

work="$GB_OUTDIR/pullviruses"
mkdir -p "$work"

# after_genomad_remaining.fna is stage 1's final remainder whether or not the geNomad claimer ran:
# with GB_GENOMAD_CLAIM=0 the subtraction awk sees an empty exclusive list and passes everything
# through, so this filename is stable across both settings.
remain="$GB_OUTDIR/prophage/after_genomad_remaining.fna"
# The CheckV / PhageBoost excisions. These are NOT whole contigs -- they are prophage regions cut
# out of a parent, carrying new ids, and their parents were removed from the remainder. Dropping
# them here would lose every prophage the pipeline exists to find.
claims="$GB_OUTDIR/prophage/viruses_prophages.fna"
out="$work/pulled_viruses.fna"

: > "$out"
[ -s "$remain" ] && cat "$remain" >> "$out"
[ -s "$claims" ] && cat "$claims" >> "$out"

n_r=$(grep -c '^>' "$remain" 2>/dev/null || echo 0)
n_c=$(grep -c '^>' "$claims" 2>/dev/null || echo 0)
n_o=$(grep -c '^>' "$out"    2>/dev/null || echo 0)
n_u=$(grep '^>' "$out" | sort -u | wc -l)

echo "[4_clean] routing the whole assembly to stage 5"
echo "[4_clean]   stage-1 remainder : $n_r"
echo "[4_clean]   stage-1 claims    : $n_c   (CheckV excisions; PhageBoost retired 2026-08-10)"
echo "[4_clean]   candidate pool    : $n_o  (unique ids: $n_u)"
[ "$n_o" -ne "$n_u" ] && echo "[4_clean] WARN: $((n_o-n_u)) duplicate ids -- stage 5 will dedup" >&2

exit 0
