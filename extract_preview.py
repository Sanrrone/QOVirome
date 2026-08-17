#!/usr/bin/env python3
"""extract_preview.py — emit the virome preview JSON for run.sh capture.

Host-side: invoked by run.sh after the container exits. Computes summary metrics
directly from the two deliverables in <results_dir>:
  - <prefix>_viruses.fna  (mined viral / phage contigs)
  - <prefix>_table.tsv    (per-contig scores + predicted bacterial host)

Prints a single JSON object to stdout (consumed as run.resultsPreview.metrics).
Exits 1 + prints to stderr if a deliverable is missing (contract violation).
Zero mined viruses is a valid outcome (exit 0 with n_viral_contigs = 0).
"""
import glob
import json
import os
import sys
from collections import Counter

TOP_HOSTS = 5

# Ranks too coarse to be a useful "predicted bacterial host". A superkingdom call
# ("Bacteria") or a clade ("Terrabacteria group") just says "some bacterium" and,
# being the modal Kraken2 outcome on divergent viral contigs, otherwise dominates
# n_with_predicted_host and top_hosts. Phylum and finer name a real lineage and are
# kept. The full per-contig h_rank/h_name is still in the deliverable table.
_COARSE_HOST_RANKS = {"superkingdom", "domain", "kingdom", "clade", "no rank", ""}


def _one(results_dir, pattern):
    hits = sorted(glob.glob(os.path.join(results_dir, pattern)))
    return hits[0] if hits else None


def _fasta_stats(path):
    n = 0
    total = 0
    longest = 0
    cur = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    total += cur
                    longest = max(longest, cur)
                n += 1
                cur = 0
            else:
                cur += len(line.strip())
        if n:
            total += cur
            longest = max(longest, cur)
    return n, total, longest


def _host_stats(path):
    n_with_host = 0
    hosts = Counter()
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx = header.index("h_name")
        except ValueError:
            return 0, []  # no host column → no host info
        # h_rank sits alongside h_name; absent on older tables → no rank filter.
        try:
            rank_idx = header.index("h_rank")
        except ValueError:
            rank_idx = None
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= idx:
                continue
            name = cols[idx].strip()
            if not name or name.upper() == "NA":
                continue
            if rank_idx is not None and len(cols) > rank_idx:
                if cols[rank_idx].strip().lower() in _COARSE_HOST_RANKS:
                    continue  # too coarse to count as a host
            n_with_host += 1
            hosts[name] += 1
    top = [{"name": name, "count": count} for name, count in hosts.most_common(TOP_HOSTS)]
    return n_with_host, top


def main():
    if len(sys.argv) < 2:
        print("Usage: extract_preview.py <results_dir>", file=sys.stderr)
        sys.exit(1)
    results_dir = sys.argv[1]

    fasta = _one(results_dir, "*_viruses.fna")
    table = _one(results_dir, "*_table.tsv")
    if not fasta:
        print(f"Contract violation: no *_viruses.fna in {results_dir}", file=sys.stderr)
        sys.exit(1)
    if not table:
        print(f"Contract violation: no *_table.tsv in {results_dir}", file=sys.stderr)
        sys.exit(1)

    n_viral, total_bp, longest_bp = _fasta_stats(fasta)
    n_with_host, top_hosts = _host_stats(table)

    print(json.dumps({
        "preview_type": "virome",
        "n_viral_contigs": n_viral,
        "total_viral_bp": total_bp,
        "longest_contig_bp": longest_bp,
        "n_with_predicted_host": n_with_host,
        "top_hosts": top_hosts,
    }))


if __name__ == "__main__":
    main()
