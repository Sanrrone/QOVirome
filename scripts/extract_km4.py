#!/usr/bin/env python3
"""136 canonical tetranucleotide FREQUENCIES per contig -- the km4 RAW block of stack9.

RAW per-k-mer frequencies, NOT PCA components. That distinction is the whole result: k-mers were
believed rejected since 2026-08-02, but that screen was wrong twice -- it scored against the VLP
label invalidated the next day, and it PCA-15'd each block before any model saw it (km5's PC1 holds
24.4% of variance, so most of the block never reached a classifier). Given raw, km4 adds +0.0073 to
stack7 on CHGV and replicates on curated. It saturates at k=4: k=5 (512 cols) and k=6 (2080) buy
+0.0003, inside CI overlap.

MECHANISM: short-motif avoidance under host defence. Restriction-modification targets 4-6 bp
palindromes and CRISPR PAMs are 2-6 bp, so both pressures act on motifs that are definitionally
INSIDE a tetranucleotide vector. That is why 13 explicit restriction-enzyme densities, 5 palindrome
summary columns, a 29-column extended-palindrome panel and a 66-column physicochemical panel all
added nothing on top: this block is a SUPERSET of those features, not a competitor.

CANONICAL counts collapse each k-mer with its reverse complement -- the right invariance for
double-stranded DNA on contigs of arbitrary orientation. 256 raw -> 136 canonical (16 of the 256 are
their own reverse complement). Both the count and the COLUMN NAMES AND ORDER are a CONTRACT with the
artifact's "features" list; they are reproduced here by the identical construction the training
table used (`km5_extract.py`), and `--selftest` re-derives them and refuses to run on a mismatch.

COST: one vectorized pass per contig -- LUT-encode, roll the code, one 256-wide bincount, fold to
136. No Python loop over positions. Seconds for a full candidate set, against geNomad's ~509s and
the PHROGs search's ~6 min in the same stage.

MISSING: a contig shorter than 4 bp, or with no ACGT-only window, emits EMPTY fields (not "nan"),
which score_stack4.py fills from the artifact's frozen training median.
"""
import argparse
import csv
import multiprocessing as mp
import sys

import numpy as np

K = 4
NCAN = 136
LUT = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate("ACGT"):
    LUT[ord(b)] = i
    LUT[ord(b.lower())] = i

_codes = np.arange(4 ** K, dtype=np.int64)
_pw = 4 ** np.arange(K - 1, -1, -1, dtype=np.int64)
_digits = (_codes[:, None] // _pw[None, :]) % 4
_rc = ((3 - _digits)[:, ::-1] * _pw[None, :]).sum(1)
_uniq, INV = np.unique(np.minimum(_codes, _rc), return_inverse=True)
INV = INV.astype(np.int64)
assert len(_uniq) == NCAN, len(_uniq)

_rep = {}
for _code in range(4 ** K):
    _j = int(INV[_code])
    if _j not in _rep:
        _rep[_j] = "".join("ACGT"[(_code // int(p)) % 4] for p in _pw)
NAMES = [f"km{K}_{_rep[j]}" for j in range(NCAN)]
COLS = ["id"] + NAMES


def read_fasta(path):
    header, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line)
        if header is not None:
            yield header, "".join(buf)


def process_one(rec):
    cid, seq = rec
    arr = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    valid = arr < 4
    safe = np.where(valid, arr, 0).astype(np.int64)
    n = len(arr)
    if n < K:
        return [cid] + [None] * NCAN
    m = n - K + 1
    code = np.zeros(m, dtype=np.int64)
    ok = np.ones(m, dtype=bool)
    for i in range(K):
        code = code * 4 + safe[i:i + m]
        ok &= valid[i:i + m]
    if not ok.any():
        return [cid] + [None] * NCAN
    raw = np.bincount(code[ok], minlength=4 ** K)
    can = np.bincount(INV, weights=raw, minlength=NCAN)
    tot = can.sum()
    if not tot:
        return [cid] + [None] * NCAN
    return [cid] + (can / tot).tolist()


def selftest():
    """The names and their ORDER are a contract with the artifact. Verify, do not assume."""
    assert len(COLS) == NCAN + 1, len(COLS)
    assert NAMES[0] == "km4_AAAA", NAMES[0]
    assert NAMES[-1] == "km4_TTAA", NAMES[-1]
    # every name must be the lexicographic min of itself and its reverse complement
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for nm in NAMES:
        w = nm.split("_")[1]
        rc = "".join(comp[b] for b in reversed(w))
        assert w <= rc, f"{w} is not the canonical representative (rc {rc})"
    assert len(set(NAMES)) == NCAN
    # frequencies must sum to 1 and a homopolymer must put all its mass on one column
    row = process_one(("t", "ACGT" * 500))
    s = sum(row[1:])
    assert abs(s - 1.0) < 1e-9, s
    row = process_one(("p", "A" * 100))
    assert abs(row[1 + NAMES.index("km4_AAAA")] - 1.0) < 1e-9
    assert process_one(("short", "AC"))[1] is None
    assert process_one(("ns", "NNNNNNNN"))[1] is None
    print(f"SELFTEST PASS — {NCAN} canonical columns, {NAMES[0]}..{NAMES[-1]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta")
    ap.add_argument("--out")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    if not a.fasta or not a.out:
        ap.error("--fasta and --out are required unless --selftest")
    recs = list(read_fasta(a.fasta))
    print(f"{len(recs)} contigs, {sum(len(s) for _, s in recs)/1e6:.1f} Mbp, "
          f"{NCAN} k-mer columns", file=sys.stderr)
    n = 0
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        with mp.Pool(a.workers) as pool:
            for row in pool.imap(process_one, recs, chunksize=50):
                w.writerow([row[0]] + ["" if v is None else f"{v:.8g}" for v in row[1:]])
                n += 1
    print(f"wrote {n} rows -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
