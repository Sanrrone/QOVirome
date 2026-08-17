#!/usr/bin/env python3
"""15 PHROGs profile-family features per contig -- the PHROGs block of stack9.

PHROGs is a different MODALITY from everything else shipped: 38,880 viral protein families built by
HMM-HMM remote homology comparison, so it reaches further into the twilight zone than the pairwise
sequence search the shipped contrast uses. It also carries FUNCTIONAL CATEGORIES, which no other
feature has -- "this contig has genes that look like capsid, portal and tail" is a different claim
from "this contig has genes that score high against viral proteins". Measured: +0.0065 CHGV AUC on
top of stack7, and it GROWS rather than shrinks as the stack gets stronger.

⚠️ HMMER, NOT MMseqs2 -- this is a correctness constraint, not a preference. The Pharokka bundle's
`phrogs_profile_db` was built with MMseqs2 v18.8cc5c; this container ships 14.7e284. The prefilter
appears to succeed and the align stage then dies with "Alignment died" -- a silent format
incompatibility, not a resource limit (it reproduced at 15,000 targets with --max-seqs 3000 and 24G
split memory). The .h3m in the SAME bundle is HMMER3 binary, whose format is stable across versions.
The shipped training features were computed on the HMMER path, so it is also the only path that
reproduces the distribution the model was fitted on.

pyhmmer is pinned to 0.10.15 in its own env (Apptainer.def) precisely so nothing can drift it.

GENE CALLING IS FREE: this consumes geNomad's own Prodigal-called proteins, the same .faa the
viral-vs-cellular contrast already uses. No second gene-calling pass.

MISSING VS ZERO -- the distinction matters and is easy to get wrong:
  - A contig WITH called genes but NO PHROG hit gets 0.0 in every column. That is a real
    measurement ("nothing viral-looking here"), and the model was fitted on exactly that encoding.
  - A contig with NO called genes at all never appears in the .faa, so it is absent from this
    table entirely, travels through the join as EMPTY, and score_stack4.py fills it from the
    artifact's frozen training median. That is the same population n_uscg is undefined for.
Emitting 0.0 for the second case would tell the model "no viral genes" about a contig nobody looked
at, which is a different and false claim.
"""
import argparse
import csv
import re
import sys
import time
from collections import defaultdict

HDR = re.compile(r"^>(\S+)")

# The 9 annotated PHROG categories, mapped to the column names the model was fitted with.
# Anything else (overwhelmingly "unknown function", 104,113 of 109,331 families) lands in cat_unk.
CATS = {"head and packaging": "cat_head", "connector": "cat_conn", "tail": "cat_tail",
        "DNA, RNA and nucleotide metabolism": "cat_dna", "lysis": "cat_lysis",
        "integration and excision": "cat_integ", "transcription regulation": "cat_trans",
        "moron, auxiliary metabolic gene and host takeover": "cat_moron", "other": "cat_other"}
STRUCT = ["cat_head", "cat_conn", "cat_tail"]
CAT_COLS = list(CATS.values()) + ["cat_unk"]
COLS = (["id", "phrog_hit_frac", "phrog_mean_bits", "phrog_best_bits", "phrog_known_frac"]
        + CAT_COLS + ["phrog_struct_frac"])


def contig_of(protein_id):
    """Prodigal/geNomad protein ids are <contig>_<n>."""
    return protein_id.rsplit("_", 1)[0]


def load_annot(path):
    """category values containing commas are QUOTED in the source TSV; csv strips them, which is
    what pandas did when the training table was built. Parsing on raw split would mismatch."""
    annot = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t", quotechar='"')
        for row in rd:
            annot[str(row["phrog"])] = row["category"]
    return annot


def scan(faa, hmm, cpus, evalue):
    """Best-scoring PHROG per protein. Mirrors phrog_hmm.py, which produced the training table."""
    import pyhmmer
    from pyhmmer.easel import SequenceFile
    from pyhmmer.plan7 import HMMFile

    def query_name(hits):
        """pyhmmer moved this between 0.10 and 0.11; support both rather than pin a version."""
        q = getattr(hits, "query", None)
        if q is not None and getattr(q, "name", None) is not None:
            return q.name.decode()
        return hits.query_name.decode()

    alphabet = pyhmmer.easel.Alphabet.amino()
    t0 = time.time()
    with SequenceFile(faa, digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    print(f"[phrog] {len(seqs)} proteins loaded in {time.time()-t0:.1f}s", flush=True)

    best = {}
    n_hmm = 0
    t0 = time.time()
    with HMMFile(hmm) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=cpus, E=evalue):
            q = query_name(hits)
            n_hmm += 1
            for hit in hits:
                if hit.evalue > evalue:
                    continue
                nm = hit.name.decode()
                if nm not in best or hit.score > best[nm][0]:
                    best[nm] = (hit.score, q)
            if n_hmm % 5000 == 0:
                print(f"[phrog] {n_hmm} profiles, {len(best)} proteins hit, "
                      f"{time.time()-t0:.0f}s", flush=True)
    print(f"[phrog] {n_hmm} profiles vs {len(seqs)} proteins in {time.time()-t0:.0f}s -> "
          f"{len(best)} proteins with a hit "
          f"({100*len(best)/max(len(seqs),1):.1f}%)", flush=True)
    # ROUND TO 2dp -- a contract with the frozen coefficients, not cosmetics. The training table was
    # produced by phrog_hmm.py writing the best hit as f"{bits:.2f}" and the feature builder reading
    # that back, so the model was fitted on 2-decimal bit scores. Keeping full precision here shifts
    # phrog_best_bits by up to 0.005 against the distribution the trees were split on. Rounding
    # AFTER the argmax (round is monotone, so max of rounded == rounded max) reproduces the training
    # column exactly; phrog_mean_bits then sums the same rounded values the fit summed.
    return {k: (round(v[0], 2), v[1]) for k, v in best.items()}


def gene_counts(faa):
    genes = defaultdict(int)
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                genes[contig_of(HDR.match(line).group(1))] += 1
    return genes


def build(genes, best, annot):
    blank = {"n_hit": 0, "bits_sum": 0.0, "bits_max": 0.0, "known": 0,
             **{v: 0 for v in CAT_COLS}}
    per = defaultdict(lambda: dict(blank))
    for prot, (bits, ph) in best.items():
        d = per[contig_of(prot)]
        d["n_hit"] += 1
        d["bits_sum"] += bits
        d["bits_max"] = max(d["bits_max"], bits)
        cat = annot.get(ph.replace("phrog_", ""), "unknown function")
        col = CATS.get(cat)
        if col:
            d[col] += 1
            d["known"] += 1
        else:
            d["cat_unk"] += 1

    rows = []
    for contig, n in genes.items():
        d = per.get(contig, blank)
        r = {"id": contig,
             "phrog_hit_frac": d["n_hit"] / n,
             "phrog_mean_bits": d["bits_sum"] / n,
             "phrog_best_bits": d["bits_max"],
             "phrog_known_frac": d["known"] / n}
        for col in CAT_COLS:
            r[col] = d[col] / n
        r["phrog_struct_frac"] = sum(d[c] for c in STRUCT) / n
        rows.append(r)
    return rows


def selftest():
    annot = {"1": "integration and excision", "2": "head and packaging", "9": "unknown function"}
    genes = {"ctg1": 4, "ctg2": 2}
    best = {"ctg1_1": (100.0, "phrog_2"), "ctg1_2": (50.0, "phrog_1"),
            "ctg1_3": (10.0, "phrog_9")}
    rows = {r["id"]: r for r in build(genes, best, annot)}
    a = rows["ctg1"]
    assert abs(a["phrog_hit_frac"] - 3 / 4) < 1e-12, a["phrog_hit_frac"]
    assert abs(a["phrog_mean_bits"] - 160.0 / 4) < 1e-12, a["phrog_mean_bits"]
    assert abs(a["phrog_best_bits"] - 100.0) < 1e-12
    assert abs(a["phrog_known_frac"] - 2 / 4) < 1e-12          # phrog_9 is unknown -> not "known"
    assert abs(a["cat_head"] - 1 / 4) < 1e-12
    assert abs(a["cat_integ"] - 1 / 4) < 1e-12
    assert abs(a["cat_unk"] - 1 / 4) < 1e-12
    assert abs(a["phrog_struct_frac"] - 1 / 4) < 1e-12         # head only; no connector/tail
    b = rows["ctg2"]                                            # genes but NO hits -> real zeros
    assert b["phrog_hit_frac"] == 0.0 and b["phrog_best_bits"] == 0.0
    assert set(rows) == {"ctg1", "ctg2"}                        # nothing invented
    assert len(COLS) == 16, len(COLS)
    assert COLS[1] == "phrog_hit_frac" and COLS[-1] == "phrog_struct_frac"
    # a contig with no called genes must be ABSENT, not zero-filled
    assert "ctg3" not in rows
    print(f"SELFTEST PASS — {len(COLS)-1} columns, {COLS[1]}..{COLS[-1]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa")
    ap.add_argument("--hmm")
    ap.add_argument("--annot")
    ap.add_argument("--out")
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    for req in ("faa", "hmm", "annot", "out"):
        if not getattr(a, req):
            ap.error(f"--{req} is required unless --selftest")

    annot = load_annot(a.annot)
    print(f"[phrog] annotation: {len(annot)} phrogs", flush=True)
    genes = gene_counts(a.faa)
    print(f"[phrog] {sum(genes.values())} proteins on {len(genes)} contigs", flush=True)
    best = scan(a.faa, a.hmm, a.cpus, a.evalue)
    rows = build(genes, best, annot)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            # .17g, NOT .8g. The PHROGs training table was written by pandas at full precision, so
            # 8 significant digits loses up to 1e-5 on phrog_mean_bits (values run to ~1e3) and
            # perturbs 26,163 of 14,916x15 cells against the distribution the trees were split on.
            # ⚠️ extract_km4.py is the OPPOSITE case and must KEEP .8g: its training table came from
            # km5_extract.py, which itself wrote .8g, so .8g is what reproduces it exactly (verified
            # max|diff| = 0.0). Match the formatting of the table the model was FITTED on, per file.
            w.writerow({k: (f"{v:.17g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"[phrog] wrote {len(rows)} rows -> {a.out}", file=sys.stderr)
    print("PHROG_FEATURES_DONE", flush=True)


if __name__ == "__main__":
    main()
