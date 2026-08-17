#!/usr/bin/env python3
"""
Applies the frozen stack4 logistic-regression model to a joined per-contig
feature table, producing a per-contig viral_score in (0, 1).

As of 2026-08-03 the model carries 6 covariates (was 4), so the design is 48
columns: 6 base signals + 6 covariates + their 36 pairwise interactions. The
"stack4" name is kept throughout (script name, DB directory, artifact name) --
it identifies the deployed model lineage, not the covariate count, and renaming
a bind-mounted DB path across two hosts buys nothing.

predict_proba is reproduced BY HAND (standardize -> build the interaction
design -> dot with the frozen coefficients -> sigmoid) rather than loaded via
sklearn, so a mismatched sklearn version in the runtime env can never silently
change scores.

MISSING VALUES, two distinct cases:
  - skew_concordance / strand_bias_B are legitimately undefined for some
    contigs (no gene sits in a window with |GC skew| > 0.05, ~2% of the
    training population). These are imputed from the model's frozen training
    median -- exactly what the fit did -- rather than failing the contig.
    The set of imputable columns comes from the artifact, not from a list
    hardcoded here, so it cannot drift from what the model was fit with.
  - Every other feature is produced by a single upstream tool call per contig,
    so a missing value means that tool failed on this specific contig. Those
    contigs fail OUT of scoring (empty score) and are excluded -- not
    imputed, and not silently rejected as if they had scored low (see the
    stage-5 gate's handling).
"""
import argparse
import csv
import json
import math

FEATURE_ORDER = [
    "mpp", "nn_chromosome", "dvf_pca_01", "dvf_pca_03", "dvf_pca_06", "dvf_pca_13",
    "delta_g_per_kb", "gc3_content_max", "mean_propeller_twist", "km4dist_pca_00",
    "skew_concordance", "strand_bias_B",
]
N_BASE = 6
N_COV = 6

# As of 2026-08-04 an artifact may instead declare "design": "linear", carrying its own
# "features" list, a single "scaler", and one coefficient per feature -- no interaction block.
#
# The 6+6+36 design above existed to squeeze signal out of a weak compositional stack. Measured
# against VIRALITY (curated phage genomes vs curated bacterial genomes) that stack scored
# 0.5776/0.5302/0.6160 while geNomad, mpp and DeepVirFinder each reached ~0.96 on the same
# fragments -- it had been fitted against CHGV VLP enrichment, which is an abundance/recovery
# measurement rather than a virality one. The replacement uses nine individually strong signals
# where interactions bought nothing.
#
# As of 2026-08-05 an artifact may declare "design": "gbm" instead, carrying its own "features"
# list and a "gbm" block of {baseline, trees, node_layout}. Trees are walked in pure Python.
#
# As of 2026-08-11 an artifact may declare "design": "forest" -- a bagged ensemble (ExtraTrees)
# carrying {trees, node_layout, aggregation, n_estimators}. Same routing rule as "gbm"; what
# differs is that the per-tree leaf values are PROBABILITIES that get averaged, with no baseline
# and no sigmoid. The "aggregation" tag is asserted so the two can never be silently confused.
#
# The legacy path is kept intact so an old artifact still scores identically.


def sigmoid(x):
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def score_row(values, base_mean, base_scale, cov_mean, cov_scale, coef, intercept):
    x6 = [(values[i] - base_mean[i]) / base_scale[i] for i in range(N_BASE)]
    c = [(values[N_BASE + j] - cov_mean[j]) / cov_scale[j] for j in range(N_COV)]
    design = x6 + c
    for j in range(N_COV):
        design.extend(x6[i] * c[j] for i in range(N_BASE))
    z = intercept + sum(k * d for k, d in zip(coef, design))
    return sigmoid(z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="TSV with columns: id + the FEATURE_ORDER columns")
    ap.add_argument("--model", required=True, help="frozen model JSON")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.model) as f:
        model = json.load(f)

    if model.get("design") == "gbm":
        # Gradient-boosted trees, walked BY HAND for the same reason the logistic paths are
        # reproduced by hand: a sklearn version bump must never silently move a score. The trees
        # are plain numbers in the artifact, so there is no sklearn (or numpy) at runtime.
        # Parity with predict_proba is verified at FIT time and recorded in _provenance.
        feature_order = model["features"]
        impute_median = model.get("impute_median", {})
        _trees = model["gbm"]["trees"]
        _baseline = model["gbm"]["baseline"]
        _LEAF, _VAL, _FIDX, _TH, _L, _R, _MGL = range(7)
        assert model["gbm"]["node_layout"] == [
            "is_leaf", "value", "feature_idx", "num_threshold",
            "left", "right", "missing_go_to_left"], "gbm node layout drifted from the frozen model"

        def score_values(values):
            raw = _baseline
            for t in _trees:
                i = 0
                while not t[i][_LEAF]:
                    nd = t[i]
                    v = values[nd[_FIDX]]
                    if v != v:                      # NaN -- follow the split's missing direction
                        i = nd[_L] if nd[_MGL] else nd[_R]
                    elif v <= nd[_TH]:
                        i = nd[_L]
                    else:
                        i = nd[_R]
                raw += t[i][_VAL]
            return sigmoid(raw)
    elif model.get("design") == "forest":
        # A bagged forest (ExtraTrees). ROUTING is identical to the gbm path -- sklearn also splits
        # `v <= threshold -> left` -- so the traversal below is the same walk. AGGREGATION is not:
        # a forest is the MEAN of per-tree P(class 1) at the reached leaf, with no baseline and no
        # sigmoid. Feeding that through the gbm path would sigmoid an already-probability mean and
        # squash every score toward 0.5, so the aggregation is tagged in the artifact and asserted
        # here: a mismatched artifact/scorer pair fails loudly rather than scoring wrongly.
        feature_order = model["features"]
        impute_median = model.get("impute_median", {})
        _trees = model["forest"]["trees"]
        _ntree = len(_trees)
        _LEAF, _VAL, _FIDX, _TH, _L, _R, _MGL = range(7)
        assert model["forest"]["node_layout"] == [
            "is_leaf", "value", "feature_idx", "num_threshold",
            "left", "right", "missing_go_to_left"], "forest node layout drifted from the frozen model"
        assert model["forest"]["aggregation"] == "mean_proba", \
            "unknown forest aggregation -- refusing to guess how to combine the trees"
        assert _ntree == model["forest"].get("n_estimators", _ntree), "tree count drifted"

        def score_values(values):
            acc = 0.0
            for t in _trees:
                i = 0
                while not t[i][_LEAF]:
                    nd = t[i]
                    v = values[nd[_FIDX]]
                    if v != v:                      # NaN -- follow the split's missing direction
                        i = nd[_L] if nd[_MGL] else nd[_R]
                    elif v <= nd[_TH]:
                        i = nd[_L]
                    else:
                        i = nd[_R]
                acc += t[i][_VAL]
            return acc / _ntree
    elif model.get("design") == "linear":
        feature_order = model["features"]
        mean = model["scaler"]["mean_"]
        scale = model["scaler"]["scale_"]
        coef = model["logreg"]["coef_"]
        intercept = model["logreg"]["intercept_"]
        impute_median = model.get("impute_median", {})
        assert len(coef) == len(feature_order), "coefficient count does not match feature count"
        assert len(mean) == len(feature_order) and len(scale) == len(feature_order), \
            "scaler width does not match feature count"

        def score_values(values):
            z = intercept + sum(
                coef[i] * ((values[i] - mean[i]) / scale[i]) for i in range(len(coef)))
            return sigmoid(z)
    else:
        feature_order = model["champion_features"] + model["covariates"]
        assert feature_order == FEATURE_ORDER, "feature order drifted from the frozen model"
        impute_median = model.get("impute_median", {})
        _bm = model["base_scaler"]["mean_"]
        _bs = model["base_scaler"]["scale_"]
        _cm = model["covariate_scaler"]["mean_"]
        _cs = model["covariate_scaler"]["scale_"]
        _coef = model["logreg"]["coef_"]
        _int = model["logreg"]["intercept_"]
        assert len(_coef) == N_BASE + N_COV + N_BASE * N_COV, \
            "coefficient count does not match the design width"

        def score_values(values):
            return score_row(values, _bm, _bs, _cm, _cs, _coef, _int)
    n_scored = 0
    n_missing = 0
    n_imputed = 0
    with open(args.input, newline="") as fin, open(args.out, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "viral_score"])
        for row in reader:
            values = []
            usable = True
            row_imputed = False
            for col in feature_order:
                try:
                    v = float(row[col])
                    if not math.isfinite(v):
                        raise ValueError("non-finite")
                except (KeyError, TypeError, ValueError):
                    if col in impute_median:
                        v = impute_median[col]
                        row_imputed = True
                    else:
                        usable = False
                        break
                values.append(v)
            if not usable:
                writer.writerow([row.get("id", ""), ""])
                n_missing += 1
                continue
            score = score_values(values)
            writer.writerow([row["id"], "{:.6f}".format(score)])
            n_scored += 1
            if row_imputed:
                n_imputed += 1

    print("scored={} (of which {} used an imputed covariate) missing={} -> {}".format(
        n_scored, n_imputed, n_missing, args.out))


if __name__ == "__main__":
    main()
