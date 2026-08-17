#!/usr/bin/env python3
"""Per-contig viral-vs-cellular sequence contrast.

For each predicted protein, MMseqs2 is searched against TWO reference sets and the discriminative
quantity is the DIFFERENCE:

    delta(gene) = best_bits(BFVD, viral) - best_bits(AlphaFold/Swiss-Prot, cellular)

Counting hits to the viral database ALONE is worthless -- measured at 0.4943 AUC, i.e. chance --
because BFVD is full of AMGs and replication machinery whose sequences are shared with cellular
life, so "hits a viral protein" mostly means "has any recognisable protein at all". The contrast
is what separates a phage capsid (viral hit, no cellular counterpart) from a housekeeping enzyme
(hits both). Neither search can be dropped to save time; dropping the cellular one deletes the
signal entirely rather than halving it.

Measured standalone (grouped CV): 0.8776 on CHGV virality, 0.9880 on curated novel phage -- above
geNomad (0.7494 / 0.9689) and MetaPhaPred (0.8436 / 0.9629) on the same contigs.

Structure was tested and REJECTED as the mechanism: Foldseek 3Di + ProstT5 over the same two
databases adds +0.0025 AUC over this sequence search for a GPU that Mars does not have. See
BENCHMARKS.md.

INPUT is geNomad's own Prodigal-called protein FASTA (<prefix>_annotate/<prefix>_proteins.faa),
so gene calling costs nothing extra -- stage 5 already runs geNomad for nn_chromosome. The header
carries the coordinates and the contig id is the header up to the last '_N':

    >c_4_1 # 205 # 324 # -1 # ID=1_1;partial=00;...

As of 2026-08-05 (stack7) three more columns come out of the SAME two searches at no extra cost:

    dark_frac        fraction of genes hitting NEITHER reference
    mean_bits_v      mean best viral bits over ALL genes (non-hitters contribute 0)
    coding_density   coding bases / contig length

They were previously written by struc_features.py under a `struc_` prefix, which was an artifact of
that file being shared between the structural and sequence arms and caused a false "we need a GPU
for this" objection. They are plain CPU columns; the prefix is dropped here.

⚠️ coding_density REQUIRES a real contig length. Reproduction against the training column:
length-file definition 100.00% exact, max-gene-end 27.16%. --fasta is therefore not optional when
the model asks for coding_density; max gene end is only a per-contig fallback for a contig absent
from the FASTA, matching the fit's own fallback.

⚠️ The output columns and their NUMBER FORMATTING are a contract with the frozen model
(stack7_final_model.json). The fit consumed values rounded exactly as written below, so changing
a format string shifts the feature distribution out from under the model.
"""
import argparse
import re
from collections import defaultdict

HDR = re.compile(r"^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s")

COLS = ["id", "n_genes", "vc_mean_delta", "vc_delta_p90", "vc_delta_max", "vc_frac",
        "dark_frac", "mean_bits_v", "coding_density"]


def read_genes(faa):
    """contig -> [(gene id, start, end)]. Coordinates feed coding_density."""
    genes = defaultdict(list)
    with open(faa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = HDR.match(line)
            if not m:
                continue
            gid = m.group(1)
            genes[gid.rsplit("_", 1)[0]].append((gid, int(m.group(2)), int(m.group(3))))
    return genes


def read_lengths(fasta):
    """contig -> length. Empty dict when no FASTA is supplied."""
    lengths = {}
    if not fasta:
        return lengths
    name, n = None, 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def read_best_bits(m8, max_evalue):
    """gene_id -> bit score of its single best hit (max bits, e-value filtered)."""
    best = {}
    if not m8:
        return best
    with open(m8) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                ev, bits = float(f[10]), float(f[11])
            except ValueError:
                continue
            if ev > max_evalue:
                continue
            q = f[0]
            if q not in best or bits > best[q]:
                best[q] = bits
    return best


def pct(vals, q):
    """Linear-interpolated percentile. Kept identical to the fit's implementation."""
    if not vals:
        return 0.0
    v = sorted(vals)
    k = (len(v) - 1) * q
    lo, hi = int(k), min(int(k) + 1, len(v) - 1)
    return float(v[lo] + (v[hi] - v[lo]) * (k - lo))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True, help="geNomad <prefix>_proteins.faa")
    ap.add_argument("--viral-m8", required=True, help="MMseqs2 hits vs BFVD")
    ap.add_argument("--cellular-m8", required=True, help="MMseqs2 hits vs AlphaFold/Swiss-Prot")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--fasta", default="", help="contigs FASTA -- REQUIRED for coding_density")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes = read_genes(a.faa)
    vh = read_best_bits(a.viral_m8, a.evalue)
    ch = read_best_bits(a.cellular_m8, a.evalue)
    lengths = read_lengths(a.fasta)

    n_fallback = 0
    with open(a.out, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for contig, glist in genes.items():
            n = len(glist)
            deltas, vwins = [], 0
            dark, vbits, coding_bp = 0, 0.0, 0
            for gid, start, end in glist:
                bv = vh.get(gid, 0.0)
                bc = ch.get(gid, 0.0)
                # A gene hitting NEITHER reference contributes a 0 delta, not a missing value:
                # "no recognisable protein" is genuinely neutral evidence, and phage genomes are
                # full of such dark genes. Dropping them would rescale the mean by a
                # virus-correlated quantity.
                if bv > 0 and bv > bc:
                    vwins += 1
                if bv == 0 and bc == 0:
                    dark += 1
                if bv > 0:
                    vbits += bv
                coding_bp += end - start + 1
                deltas.append(bv - bc)
            # Fall back to the last gene end only when the contig is absent from the FASTA --
            # the fit used the same fallback, so this cannot drift from the training column.
            span = lengths.get(contig)
            if not span:
                span = max(e for _, _, e in glist)
                n_fallback += 1
            out.write("\t".join([
                contig, str(n),
                f"{(sum(deltas) / n):.4f}",
                f"{pct(deltas, 0.90):.2f}",
                f"{max(deltas):.2f}",
                f"{vwins / n:.6f}",
                f"{dark / n:.6f}",
                f"{vbits / n:.4f}",
                f"{coding_bp / max(span, 1):.6f}",
            ]) + "\n")
    if lengths and n_fallback:
        print(f"[contrast] WARN: {n_fallback} contig(s) absent from --fasta; coding_density "
              f"used the max-gene-end fallback for those")
    elif not lengths:
        print("[contrast] WARN: no --fasta given; coding_density uses max gene end for EVERY "
              "contig, which reproduces the training column only 27% of the time")

    print(f"[contrast] contigs={len(genes)} genes={sum(len(v) for v in genes.values())} "
          f"viral_hits={len(vh)} cellular_hits={len(ch)} -> {a.out}")


if __name__ == "__main__":
    main()
