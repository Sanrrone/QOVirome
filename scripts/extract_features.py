#!/usr/bin/env python3
"""
Per-contig delta_g_per_kb, mean_propeller_twist, km4dist_pca_00 -- the 3
covariates of the stack4_interaction admission score that are NOT gene-based
(gc3_content_max is extract_gene_features.py's job) and NOT the champion6 base signals
(mpp/nn_chromosome/dvf_pca_* come from MetaPhaPred/geNomad/dvf_embed.py).

nn_free_energy and kmer_occurrence_distance are the same pure-stdlib logic
validated in research (tier12_features.py, kmer_occurrence_distance.py).
km4dist_pca_00 applies the frozen median-impute + StandardScaler + PCA basis
from stack4_final_model.json (fit once on the full training population) --
NOT refit here, same reasoning as dvf_embed.py's frozen PCA.

Pure stdlib except numpy/json for the PCA projection. Parallelized across
contigs via multiprocessing.Pool.

Run in any env with numpy (e.g. the mpp or gutbusters env):
  python extract_features.py --fasta IN.fna \
      --pca-basis /opt/gutbusters/db/stack4/stack4_final_model.json \
      --out OUT.tsv
"""
import argparse
import itertools
import json
import multiprocessing as mp
import sys
from collections import Counter, defaultdict

import numpy as np

_VALID_BASES = set("ACGT")
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

# ============================================================
# nn_free_energy (SantaLucia 1998 unified NN thermodynamics)
# ============================================================
NN_DELTA_G37 = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}
INIT_TERM_GC = 0.98
INIT_TERM_AT = 1.03


def nn_free_energy(seq):
    n = len(seq)
    total = 0.0
    for i in range(n - 1):
        b1, b2 = seq[i], seq[i + 1]
        if b1 in _VALID_BASES and b2 in _VALID_BASES:
            total += NN_DELTA_G37[b1 + b2]
    if n >= 1:
        first, last = seq[0], seq[-1]
        if first in _VALID_BASES:
            total += INIT_TERM_GC if first in "GC" else INIT_TERM_AT
        if last in _VALID_BASES:
            total += INIT_TERM_GC if last in "GC" else INIT_TERM_AT
    return (total / n * 1000.0) if n > 0 else 0.0


# ============================================================
# dna_shape_features (propeller twist, El Hassan & Calladine 1996)
# ============================================================
PROPELLER_TWIST = {
    "AA": -18.66, "AC": -13.10, "AG": -14.00, "AT": -15.01,
    "CA": -9.45, "CC": -8.11, "CG": -10.03, "CT": -14.00,
    "GA": -13.48, "GC": -11.08, "GG": -8.11, "GT": -13.10,
    "TA": -11.85, "TC": -13.48, "TG": -9.45, "TT": -18.66,
}


def mean_propeller_twist(seq):
    values = []
    for i in range(max(len(seq) - 1, 0)):
        a, b = seq[i], seq[i + 1]
        if a in _VALID_BASES and b in _VALID_BASES:
            values.append(PROPELLER_TWIST[a + b])
    return sum(values) / len(values) if values else 0.0


# ============================================================
# kmer_occurrence_distance (k=4 only -- km5dist not part of the winning stack)
# ============================================================
def _reverse_complement(kmer):
    return "".join(_COMPLEMENT[b] for b in reversed(kmer))


def kmer_occurrence_distance_k4(seq):
    k = 4
    n = len(seq)
    positions_by_kmer = defaultdict(list)
    for i in range(n - k + 1):
        window = seq[i:i + k]
        if not _VALID_BASES.issuperset(window):
            continue
        rc = _reverse_complement(window)
        canonical = window if window <= rc else rc
        positions_by_kmer[canonical].append(i)

    result = {}
    for kmer, positions in positions_by_kmer.items():
        runs = []
        run_first = run_last = positions[0]
        for p in positions[1:]:
            if p - run_last < k:
                run_last = p
                continue
            runs.append((run_first, run_last))
            run_first = run_last = p
        runs.append((run_first, run_last))
        if len(runs) < 2:
            continue
        gaps = [runs[i + 1][0] - runs[i][1] - k for i in range(len(runs) - 1)]
        result[kmer] = sum(gaps) / len(gaps)
    return result


def _all_canonical_4mers():
    seen = set()
    for tup in itertools.product("ACGT", repeat=4):
        kmer = "".join(tup)
        rc = _reverse_complement(kmer)
        seen.add(kmer if kmer <= rc else rc)
    return sorted(seen)


CANONICAL_4MERS = _all_canonical_4mers()


def read_fasta(path):
    header = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if header is not None:
            yield header, "".join(buf)


def extract_raw(record):
    contig_id, seq = record
    seq = seq.upper()
    dist = kmer_occurrence_distance_k4(seq)
    return {
        "id": contig_id,
        "delta_g_per_kb": nn_free_energy(seq),
        "mean_propeller_twist": mean_propeller_twist(seq),
        "km4_dist": [dist.get(k, float("nan")) for k in CANONICAL_4MERS],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--pca-basis", required=True, help="stack4_final_model.json (uses its km4dist_pca section)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 2))
    args = ap.parse_args()

    with open(args.pca_basis) as f:
        model = json.load(f)
    basis = model["km4dist_pca"]
    assert basis["km4d_cols"] == [f"km4_dist_{k}" for k in CANONICAL_4MERS], "km4dist column order mismatch"
    median_impute = np.array(basis["median_impute"], dtype=np.float64)
    scaler_mean = np.array(basis["scaler_mean_"], dtype=np.float64)
    scaler_scale = np.array(basis["scaler_scale_"], dtype=np.float64)
    component_0 = np.array(basis["component_0"], dtype=np.float64)

    records = list(read_fasta(args.fasta))
    print(f"loaded {len(records)} contigs, workers={args.workers}", file=sys.stderr)

    rows = []
    with mp.Pool(args.workers) as pool:
        for i, feats in enumerate(pool.imap(extract_raw, records, chunksize=100)):
            rows.append(feats)
            if (i + 1) % 10000 == 0:
                print(f"  {i + 1}/{len(records)} done", file=sys.stderr)

    with open(args.out, "w") as out:
        out.write("id\tdelta_g_per_kb\tmean_propeller_twist\tkm4dist_pca_00\n")
        for feats in rows:
            raw = np.array(feats["km4_dist"], dtype=np.float64)
            imputed = np.where(np.isfinite(raw), raw, median_impute)
            standardized = (imputed - scaler_mean) / scaler_scale
            km4dist_pca_00 = float(component_0.dot(standardized))
            out.write("{}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(
                feats["id"], feats["delta_g_per_kb"], feats["mean_propeller_twist"], km4dist_pca_00))

    print(f"wrote {len(rows)} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
