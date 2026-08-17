#!/usr/bin/env python
"""
Per-contig DVF embedding -> frozen-PCA projection, producing the 4 dvf_pca_*
columns the stack4 admission score needs (mpp/nn_chromosome/dvf_pca_01/03/06/13
are its 6 base signals).

Same tapped-model approach validated in research (see extract_dvf_embedding_chgv.py):
taps global_max_pooling1d_1 at both call sites (forward input_1, reverse-complement
input_2) of DeepVirFinder's siamese model, giving a 1000-dim vector per direction.
ALWAYS uses the 1k model regardless of contig length -- that is what the frozen PCA
basis was fit against, so switching models by length here would be a train/inference
mismatch.

Run in the dvf conda env (matches DVF's own runtime requirements):
  /opt/micromamba/envs/dvf/bin/python dvf_embed.py --fasta IN.fna \
      --model /opt/gutbusters/db/deepvirfinder/models/model_siamese_varlen_1k_fl10_fn1000_dn1000_ep30_acc0.96.h5 \
      --pca-basis /opt/gutbusters/db/stack4/dvf_pca_frozen_basis.json \
      --out OUT.tsv
"""
import argparse
import json
import os
import sys
import traceback

os.environ["KERAS_BACKEND"] = "theano"
import theano  # noqa: E402
theano.config.gcc.cxxflags = "-Wno-c++11-narrowing"

import numpy as np  # noqa: E402
from keras.models import load_model, Model  # noqa: E402
from Bio.Seq import Seq  # noqa: E402


def encode_seq(seq):
    seq_code = list()
    for pos in range(len(seq)):
        letter = seq[pos]
        if letter in ["A", "a"]:
            code = [1, 0, 0, 0]
        elif letter in ["C", "c"]:
            code = [0, 1, 0, 0]
        elif letter in ["G", "g"]:
            code = [0, 0, 1, 0]
        elif letter in ["T", "t"]:
            code = [0, 0, 0, 1]
        else:
            code = [1 / 4, 1 / 4, 1 / 4, 1 / 4]
        seq_code.append(code)
    return seq_code


def iter_fasta(path):
    header = None
    seq_chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    ap = argparse.ArgumentParser(description="Extract dvf_pca_01/03/06/13 via DVF's tapped siamese embedding + frozen PCA.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--model", required=True, help="Path to the 1k DVF .h5 model (always used, regardless of contig length).")
    ap.add_argument("--pca-basis", required=True, help="Frozen PCA basis JSON (mean_ + the 4 needed component rows).")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.pca_basis) as f:
        basis = json.load(f)
    emb_cols = basis["emb_cols"]
    mean_ = np.array(basis["mean_"], dtype=np.float64)
    comp_names = list(basis["components"].keys())
    comp_mat = np.array([basis["components"][n] for n in comp_names], dtype=np.float64)  # (4, 2000)
    assert len(emb_cols) == mean_.shape[0] == comp_mat.shape[1], "PCA basis shape mismatch"
    n_fw = sum(1 for c in emb_cols if c.startswith("fw_"))
    assert emb_cols[:n_fw] == [f"fw_{i:03d}" for i in range(n_fw)]
    assert emb_cols[n_fw:] == [f"rc_{i:03d}" for i in range(len(emb_cols) - n_fw)]

    print("Loading model: {}".format(args.model))
    model = load_model(args.model)
    gmp = model.get_layer("global_max_pooling1d_1")
    fw_pool = gmp.get_output_at(0)
    rc_pool = gmp.get_output_at(1)
    tapped = Model(inputs=model.input, outputs=[fw_pool, rc_pool])

    n_scored = 0
    n_failed = 0
    with open(args.out, "w") as out:
        out.write("\t".join(["id"] + comp_names) + "\n")
        for cid, seq in iter_fasta(args.fasta):
            try:
                L = len(seq)
                if L == 0:
                    raise ValueError("zero-length sequence")
                n_count = seq.count("N") + seq.count("n")
                if n_count / L > 0.3:
                    raise ValueError(">30% N bases ({}/{})".format(n_count, L))

                code_fw = encode_seq(seq)
                seq_rc = str(Seq(seq).reverse_complement())
                code_rc = encode_seq(seq_rc)

                x_fw = np.array([code_fw], dtype="float32")
                x_rc = np.array([code_rc], dtype="float32")
                f0, f1 = tapped.predict([x_fw, x_rc], batch_size=1)

                emb = np.concatenate([f0[0], f1[0]]).astype(np.float64)
                proj = comp_mat.dot(emb - mean_)
                out.write("\t".join([cid] + ["{:.6f}".format(v) for v in proj]) + "\n")
                n_scored += 1
            except Exception as e:
                n_failed += 1
                sys.stderr.write("FAILED contig {}: {}\n".format(cid, e))
                traceback.print_exc(file=sys.stderr)

    print("DONE. scored={} failed={}".format(n_scored, n_failed))
    print("Output written to: {}".format(args.out))


if __name__ == "__main__":
    main()
