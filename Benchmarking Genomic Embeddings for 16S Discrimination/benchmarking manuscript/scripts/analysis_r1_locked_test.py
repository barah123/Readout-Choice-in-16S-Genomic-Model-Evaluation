"""POST HOC: the R1 -> R2 readout effect measured on the LOCKED TEST.

NOT part of the frozen confirmatory analysis.

The manuscript's title foregrounds readout choice, but the +2.89 pp R1 -> R2 effect is a
DEVELOPMENT result and was never confirmed. Embeddings are cached, so the confirmatory
measurement costs nothing but is missing. This supplies it.

Fidelity to the frozen run: the outer folds are rebuilt with the frozen splitter
(StratifiedKFold(5, shuffle=True, random_state=0) over the isolate-level table) and each
fold reuses the region the frozen run actually selected, read from `_chosen_regions` in
locked_test_scores.json rather than re-derived. R2 is therefore expected to reproduce the
frozen per-block values, which is checked and reported as a fidelity test before R1 is
interpreted.

R1 is the frozen definition: tie-aware nearest-reference on cosine distance over L2-normalised
embeddings, all training neighbours within 1e-9 of the minimum distance contributing
fractional credit, scored with soft balanced accuracy (run_readouts_full.nn_credit).

Writes results/r1_locked_test_results.json.
"""
import json
import os
import sys
import warnings
from collections import Counter, defaultdict

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
import run_readouts_full as rr

warnings.filterwarnings("ignore", category=ConvergenceWarning)
SEED, N_SPLITS = 0, 5
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "r1_locked_test_results.json"))


def probe():
    return make_pipeline(StandardScaler(), LogisticRegression(
        max_iter=2000, class_weight="balanced", random_state=SEED))


def main():
    rows = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_blocks.json")))
    scores = json.load(open(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json")))

    emb = {}
    for m in MODELS:
        f = os.path.join(ROOT, "phase4_locked_test", f"emb_test_{m}.npz")
        if os.path.exists(f):
            z = np.load(f, allow_pickle=True)
            for r in [f"V{i}" for i in range(1, 10)]:
                if f"{r}__keys" in z:
                    emb[(m, r)] = dict(zip(z[f"{r}__keys"], z[f"{r}__vecs"]))

    per = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        per[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    out = {}
    for b in sorted(scores):
        v = scores[b]
        if v.get("_status") != "ok" or not v.get("_block_eligible"):
            continue
        pr = per[b]
        labels = {i: lab for region, d in pr.items() for i, (_, lab) in d.items()}
        isolates = sorted(labels)
        y_all = np.array([labels[i] for i in isolates])
        chosen = v["_chosen_regions"]

        acc = defaultdict(list)
        for k, (tr_idx, te_idx) in enumerate(StratifiedKFold(
                N_SPLITS, shuffle=True, random_state=SEED).split(
                    np.zeros(len(isolates)), y_all)):
            if k >= len(chosen):
                break
            region = chosen[k]
            d = pr[region]
            tr = [isolates[i] for i in tr_idx if isolates[i] in d]
            te = [isolates[i] for i in te_idx if isolates[i] in d]
            ytr = np.array([d[i][1] for i in tr])
            yte = np.array([d[i][1] for i in te])
            if len(set(ytr)) < 2 or len(set(yte)) < 2:
                continue
            str_tr = [d[i][0] for i in tr]
            str_te = [d[i][0] for i in te]

            for m in MODELS:
                lut = emb.get((m, region))
                if lut is None or not all(s in lut for s in str_tr + str_te):
                    continue
                Utr = np.vstack([lut[s] for s in str_tr])
                Ute = np.vstack([lut[s] for s in str_te])
                # R1: tie-aware nearest reference on cosine distance
                a = Utr / (np.linalg.norm(Utr, axis=1, keepdims=True) + 1e-12)
                c = Ute / (np.linalg.norm(Ute, axis=1, keepdims=True) + 1e-12)
                D = 1 - a @ c.T                      # train x test
                cp = rr.nn_credit(D, ytr)
                acc[f"{m}:R1"].append(rr.soft_balanced_accuracy(
                    yte, np.where(yte == 1, cp, 1 - cp)))
                # R2: linear probe
                p = probe(); p.fit(Utr, ytr)
                acc[f"{m}:R2"].append(rr.balanced_accuracy(yte, p.predict(Ute)))

        out[b] = {k: float(np.mean(x)) for k, x in acc.items() if x}
        out[b]["_frozen_dnabert2_R2"] = v["dnabert2:R2"]

    names = sorted(out)

    def boot(x, seed=1729, n=20000):
        rng = np.random.default_rng(seed)
        dr = x[rng.integers(0, len(x), size=(n, len(x)))].mean(axis=1)
        return float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))

    # fidelity: does the rebuilt R2 reproduce the frozen per-block values?
    delta = np.array([out[b]["dnabert2:R2"] - out[b]["_frozen_dnabert2_R2"] for b in names])
    fidelity = {"max_abs_deviation": float(np.abs(delta).max()),
                "exact_reproduction": bool(np.abs(delta).max() < 1e-12)}

    summary = {"n_blocks": len(names), "fidelity_vs_frozen_R2": fidelity, "per_model": {}}
    for m in MODELS:
        r1 = np.array([out[b][f"{m}:R1"] for b in names])
        r2 = np.array([out[b][f"{m}:R2"] for b in names])
        g = r2 - r1
        lo, hi = boot(g)
        summary["per_model"][m] = {
            "R1": float(r1.mean()), "R2": float(r2.mean()),
            "R1_to_R2": float(g.mean()), "ci95": [lo, hi],
            "sd": float(g.std(ddof=1)),
            "blocks_improved": int((g > 0).sum())}
    allg = np.array([np.mean([out[b][f"{m}:R2"] - out[b][f"{m}:R1"] for m in MODELS])
                     for b in names])
    lo, hi = boot(allg)
    summary["mean_over_models"] = {"R1_to_R2": float(allg.mean()), "ci95": [lo, hi],
                                   "sd": float(allg.std(ddof=1))}
    summary["development_reference"] = {
        "R1_to_R2_mean_over_models": 0.0289,
        "note": "development-set value quoted in the manuscript, 25 species-rank blocks"}

    json.dump({"summary": summary, "blocks": out}, open(OUT, "w"), indent=2, default=float)

    pp = 100.0
    print(f"fidelity: rebuilt R2 vs frozen R2, max |deviation| = "
          f"{fidelity['max_abs_deviation']:.2e}  exact={fidelity['exact_reproduction']}")
    print(f"\nLOCKED TEST, n = {len(names)} blocks")
    print(f"{'model':12}{'R1':>9}{'R2':>9}{'R1->R2 pp':>12}{'95% CI pp':>22}{'blocks up':>11}")
    for m in MODELS:
        d = summary["per_model"][m]
        print(f"{m:12}{pp*d['R1']:9.2f}{pp*d['R2']:9.2f}{pp*d['R1_to_R2']:+12.2f}"
              f"  [{pp*d['ci95'][0]:+7.2f},{pp*d['ci95'][1]:+7.2f}]{d['blocks_improved']:11}")
    a = summary["mean_over_models"]
    print(f"{'MEAN':12}{'':9}{'':9}{pp*a['R1_to_R2']:+12.2f}"
          f"  [{pp*a['ci95'][0]:+7.2f},{pp*a['ci95'][1]:+7.2f}]")
    print(f"\ndevelopment reference (manuscript): +2.89 pp")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
