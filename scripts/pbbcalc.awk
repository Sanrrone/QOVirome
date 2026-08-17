BEGIN{
    OFS="\t"
    # Score at which an arm counts as a supporting vote for the corroboration gate in
    # 5_eval.sh. Defaulted here as well as there so a hand-run of this script against a
    # kept workdir reproduces the pipeline's numbers instead of silently counting every
    # non-zero arm as support.
    #
    # Validate the SHAPE of the string, then coerce. Neither step alone is enough:
    # -v hands the value over as a string, so a malformed floor="0,5" (a comma is the
    # decimal separator in several locales and one keystroke from "0.5") would be
    # compared against each arm LEXICOGRAPHICALLY -- "." > "," so every arm clears it.
    # But `floor+0` alone does not rescue that either: it parses "0,5" as 0, and a floor
    # of 0 also counts every arm. Both failures land in the same place, a gate that
    # silently degrades to confidence-only with nothing logged. So reject anything that
    # is not a plain decimal and fall back loudly. 5_eval.sh validates this too and
    # exits; this is the second line of defence, for a hand-run against a kept workdir.
    if (floor == "") floor = 0.5
    if (floor !~ /^[0-9]*\.?[0-9]+$/) {
        print "pbbcalc: malformed floor=\"" floor "\" -- falling back to 0.5" > "/dev/stderr"
        floor = 0.5
    }
    floor = floor + 0
}
# --- DeepVirFinder ---
# Every *_vscores.tsv carries a header (dvp/phamer from the tool, mpp/vs2 synthesised
# by 5_eval.sh), so skip the first line of each — not just of the files after the
# first, which let dvp_vscores.tsv's header through as a bogus "name" record.
FNR==1 { next }
FILENAME ~ /dvp_vscores.tsv$/ {
    gsub(/dvp_vscores.tsv/,"",f)
}
FILENAME ~ /dvp_vscores.tsv$/ {
    id=$1
    dvp[id]=$3
}

# --- VirSorter2 ---
# final-viral-score.tsv is: seqname, one column per --include-groups member, then
# max_score, max_score_group, length, hallmark, viral, cellular. With the default
# groups (dsDNAphage,ssDNA) max_score is $4; $2 is the dsDNAphage group score, which
# scores an ssDNA phage near zero.
FILENAME ~ /vs2_vscores.tsv$/ {
    # VS2 appends ||full / ||partial / ||lt2gene to seqname. sub() with an escaped
    # regex, not split(...,"||") — a bare "||" is an empty-alternation regex that
    # matches nothing, so the suffix survived and VS2 rows never joined the other
    # tools' contig ids.
    id=$1; sub(/\|\|.*$/,"",id)
    vs2[id]=$4
}

# --- MetaPhaPred ---
FILENAME ~ /mpp_vscores.tsv$/ {
    id=$1
    mpp[id]=$3
}

# --- PhaMer: REMOVED from stage 5 (2026-07-27) ---
# PhaMer no longer runs in stage 5, so no phamer_vscores.tsv is produced and no reader
# belongs here. Measured on 8,372 labelled contigs (Hannigan VLP/WGS): PhaMer's marginal
# contribution to admission was 256 contigs at 43.0% precision (110 TP / 146 FP) -- the
# only arm below the 50% base rate. 225 of those 256 landed in the lowest CheckV tier at
# 41.8%, and ZERO were homology-blind, so it was diluting the discovery tier without
# contributing any discoveries. Dropping it lifts every tier (overall 75.7% -> 79.8%).
# It still runs in STAGE 3, where it nominates candidates that the remaining four arms
# then judge. See the virome_stage5_benchmark_verdict memory for the full measurement.

# --- geNomad (viral branch) ---
# virus_summary.tsv is: seq_name, length, topology, coordinates, n_genes,
# genetic_code, virus_score, fdr, n_hallmarks, marker_enrichment, taxonomy.
# Only contigs geNomad called viral appear, so an absent id legitimately means
# "geNomad did not call this viral" -> 0, which is the correct noisy-OR identity.
# The sub() is defensive: stage 5 runs with --disable-find-proviruses so no
# "|provirus_" suffix should appear, but a stray one must not fork the id.
# fdr ($8) is deliberately NOT read and must not be: this pass calls geNomad without
# --enable-score-calibration, so geNomad writes the literal string "NA" into every fdr
# cell (its own --max-fdr help says the option is ignored when scores are uncalibrated).
# The column carries no information at all, and BOTH ways of testing it are wrong:
# `$8<=0.05` compares "NA" against "0.05" as strings and is always FALSE, while
# `$8+0<=0.05` coerces "NA" to 0 and is always TRUE. Enabling calibration is not a fix
# either -- this stage scores a pre-enriched candidate set, which violates the
# assumption geNomad's calibration is derived under.
FILENAME ~ /gnmd_vscores.tsv$/ {
    id=$1; sub(/\|.*$/,"",id)
    gnmd[id]=$7
}

END{
	# COLUMN CONTRACT CHANGED 2026-07-27: the phamer column is gone, so this table is
	# 9 wide, not 10, and confidence moved $7 -> $6 / n_comp $9 -> $8 / n_hom $10 -> $9.
	# The gate awk in 5_eval.sh reads all three POSITIONALLY and must move in lockstep.
	print "ID","dvf","vs2","mpp","gnmd","confidence","n_arms","n_comp","n_hom"
	# build union of IDs
	for (id in dvp)     all[id]
	for (id in vs2)     all[id]
	for (id in mpp)     all[id]
	for (id in gnmd)    all[id]


	# iterate once over union
	for (id in all) {
	    d = (id in dvp     ? dvp[id]     : 0)
	    v = (id in vs2     ? vs2[id]     : 0)
	    m = (id in mpp     ? mpp[id]     : 0)
	    g = (id in gnmd    ? gnmd[id]    : 0)

	    # probability all 4 tools fail
	    prob_fail = (1-d)*(1-v)*(1-m)*(1-g)

	    # Corroboration, counted by DETECTION MODALITY rather than by tool. Noisy-OR
	    # cannot express this at all: an absent or zero arm is its identity element, so
	    # one tool at 0.99 with four tools at 0.00 yields confidence 1.00 and the
	    # dissent leaves no trace (67.9% of admitted contigs sat at exactly 1.00).
	    #
	    # dvf and mpp are pooled because they are the same detector twice: both encode
	    # the contig and its reverse complement as k-mers into a deep net, and they
	    # agree 88.7% of the time. Counting them as two independent votes would make the
	    # most redundant pair available the CHEAPEST way for a contig to look
	    # corroborated. Pooling them is what makes n_comp>=2 mean "the twins AGREED".
	    #
	    # HONEST CAVEAT (measured 2026-07-27, 8,372 labelled contigs): this modality
	    # pooling does NOT pay off. At matched true-positive counts the gate performs the
	    # same as simply raising GB_PHS, and the 237 contigs it rejects are 51.9%
	    # positive against a 50.0% base rate -- i.e. it rejects at random with respect to
	    # truth. It is retained ONLY because changing it would change which contigs are
	    # admitted, and the CheckV evidence tier (added downstream) does the discriminating
	    # work far better: within the top tier even single-arm contigs are 85.7% precise.
	    # Do not invest further in tuning this gate; tune the tier instead.
	    #
	    # n_hom is now geNomad + VirSorter2 only (PhaMer left stage 5). Note vs2's term is
	    # operationally ">= 0.9", not ">= floor": VirSorter2 runs --min-score $GB_PHS
	    # --high-confidence-only and never emits a row below it, so GB_CORROB_FLOOR cannot
	    # move that arm at all.
	    #
	    # +0 on every arm: these come from fields, and a field that never looked numeric
	    # would otherwise compare against floor as a string.
	    ncomp = (d+0 >= floor) + (m+0 >= floor)                 # composition / k-mer nets
	    nhom  = (g+0 >= floor) + (v+0 >= floor)                 # marker + protein homology

	    print id, d, v, m, g, 1-prob_fail, ncomp+nhom, ncomp, nhom
	}

}
