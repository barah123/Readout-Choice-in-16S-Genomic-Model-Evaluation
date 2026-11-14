"""POST HOC: per-neighbour decomposition and divergence stratification.

NOT part of the frozen confirmatory analysis.

Balanced accuracy pools every neighbour species into one class, which discards the structure
the study premise is about. A method that cleanly separates the target from two neighbours and
fails on two others scores identically to one that half-fails on all four. This decomposes the
pooled metric without changing it.

Method. The classifier is EXACTLY the one the paper reports: pooled target-vs-rest, trained on
the frozen isolate-grouped folds at each block's modal chosen region. Nothing is retrained
per neighbour. At test time the predictions are partitioned by neighbour species, and balanced
accuracy is recomputed over {test targets} union {test isolates of that species}. This
decomposes the reported metric's errors rather than substituting a different task.

Divergence is the project's canonical identity (run_readouts_full.identity_dist, the same
measure C_align and the similarity-disjoint analysis use) between the target's haplotypes and
that neighbour species' haplotypes at that region. Both the mean and the MINIMUM over
cross-pairs are reported; the minimum is the hardest pair and is the one that should drive
difficulty.

Writes results/divergence_stratified_results.json.
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
SEED, N_SPLITS, MIN_N = 0, 5, 5
OUT = os.path.abspath(os.path.join(HERE, "..", "results",
                                   "divergence_stratified_results.json"))


def probe():
    return make_pipeline(StandardScaler(), LogisticRegression(
        max_iter=2000, class_weight="balanced", random_state=SEED))


def bal(y, p):
    y, p = np.asarray(y), np.asarray(p)
    out = []
    for c in (0, 1):
        m = y == c
        if m.sum():
            out.append((p[m] == c).mean())
    return float(np.mean(out)) if len(out) == 2 else float("nan")


def main():
    rows = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_blocks.json")))
    scores = json.load(open(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json")))
    z = np.load(os.path.join(ROOT, "phase4_locked_test", "emb_test_dnabert2.npz"),
                allow_pickle=True)
    emb = {r: dict(zip(z[f"{r}__keys"], z[f"{r}__vecs"]))
           for r in [f"V{i}" for i in range(1, 10)] if f"{r}__keys" in z}

    per = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        per[r["block_id"]][r["region"]][r["isolate_group"]] = (
            r["slice"], r["is_target"], r["species"])

    pairs = []
    for b in sorted(scores):
        v = scores[b]
        if v.get("_status") != "ok" or not v.get("_block_eligible"):
            continue
        region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        d = per[b][region]
        iso = sorted(d)
        slices = [d[i][0] for i in iso]
        y = np.array([d[i][1] for i in iso])
        sp = [d[i][2] for i in iso]
        lut = emb.get(region)

        # divergence: identity between target and each neighbour species, at this region
        uniq = sorted(set(slices))
        idx = {s: k for k, s in enumerate(uniq)}
        I = 1.0 - rr.identity_dist(uniq)
        np.fill_diagonal(I, 1.0)
        t_hap = sorted({slices[i] for i in range(len(iso)) if y[i] == 1})
        ti = [idx[s] for s in t_hap]

        # pooled classifier, frozen folds, predictions collected per test isolate
        pred = {"C_kmer": np.empty(len(iso), object),
                "dnabert2": np.empty(len(iso), object)}
        for tr, te in StratifiedKFold(N_SPLITS, shuffle=True,
                                      random_state=SEED).split(np.zeros(len(y)), y):
            if len(set(y[tr])) < 2 or len(set(y[te])) < 2:
                continue
            str_tr = [slices[i] for i in tr]
            str_te = [slices[i] for i in te]
            Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
            c = probe(); c.fit(Ktr, y[tr])
            for j, i in enumerate(te):
                pred["C_kmer"][i] = int(c.predict(Kte[j:j + 1])[0])
            if lut is not None and all(s in lut for s in str_tr + str_te):
                p = probe(); p.fit(np.vstack([lut[s] for s in str_tr]), y[tr])
                P = p.predict(np.vstack([lut[s] for s in str_te]))
                for j, i in enumerate(te):
                    pred["dnabert2"][i] = int(P[j])

        tgt = [i for i in range(len(iso)) if y[i] == 1 and pred["C_kmer"][i] is not None]
        for species in sorted({sp[i] for i in range(len(iso)) if y[i] == 0}):
            nb = [i for i in range(len(iso))
                  if y[i] == 0 and sp[i] == species and pred["C_kmer"][i] is not None]
            if len(nb) < MIN_N or len(tgt) < MIN_N:
                continue
            n_hap = sorted({slices[i] for i in nb})
            ni = [idx[s] for s in n_hap]
            sub = I[np.ix_(ti, ni)]
            rec = {"block": b, "region": region, "neighbour_species": species,
                   "n_target": len(tgt), "n_neighbour": len(nb),
                   "mean_identity": float(sub.mean()), "max_identity": float(sub.max()),
                   "n_target_haplotypes": len(t_hap), "n_neighbour_haplotypes": len(n_hap),
                   "shares_haplotype_with_target": bool(set(t_hap) & set(n_hap))}
            sel = tgt + nb
            yy = np.array([1] * len(tgt) + [0] * len(nb))
            for arm in ("C_kmer", "dnabert2"):
                if any(pred[arm][i] is None for i in sel):
                    rec[arm] = float("nan")
                    continue
                rec[arm] = bal(yy, np.array([pred[arm][i] for i in sel]))
            pairs.append(rec)
        print(f"  {b:26} {region}  "
              f"{sum(1 for p in pairs if p['block'] == b)} neighbour species", flush=True)

    ok = [p for p in pairs if p["C_kmer"] == p["C_kmer"] and p["dnabert2"] == p["dnabert2"]]
    k = np.array([p["C_kmer"] for p in ok])
    d = np.array([p["dnabert2"] for p in ok])
    mi = np.array([p["max_identity"] for p in ok])

    # bins on MAXIMUM target-neighbour identity: the hardest pair drives difficulty
    edges = [0.0, 0.95, 0.98, 0.99, 0.999, 1.0001]
    lab = ["<95%", "95-98%", "98-99%", "99-99.9%", ">=99.9%"]
    bins = []
    for lo, hi, L in zip(edges[:-1], edges[1:], lab):
        m = (mi >= lo) & (mi < hi)
        if m.sum():
            bins.append({"bin": L, "n_pairs": int(m.sum()),
                         "C_kmer": float(k[m].mean()), "dnabert2": float(d[m].mean()),
                         "difference_pp": float(100 * (d[m] - k[m]).mean()),
                         "n_blocks": len({ok[i]["block"] for i in np.nonzero(m)[0]})})

    from scipy import stats as st
    rho_k = st.spearmanr(mi, k)
    rho_d = st.spearmanr(mi, d)
    summary = {
        "n_target_neighbour_pairs": len(ok),
        "n_blocks": len({p["block"] for p in ok}),
        "pooled_mean": {"C_kmer": float(k.mean()), "dnabert2": float(d.mean())},
        "bins_by_max_target_neighbour_identity": bins,
        "spearman_accuracy_vs_max_identity": {
            "C_kmer": {"rho": float(rho_k.statistic), "p": float(rho_k.pvalue)},
            "dnabert2": {"rho": float(rho_d.statistic), "p": float(rho_d.pvalue)}},
        "pairs_sharing_a_haplotype_with_target": int(
            sum(1 for p in ok if p["shares_haplotype_with_target"])),
    }
    json.dump({"summary": summary, "pairs": ok}, open(OUT, "w"), indent=2, default=float)

    print(f"\n{len(ok)} target-neighbour pairs across {summary['n_blocks']} blocks")
    print(f"pooled means: k-mer {100*k.mean():.2f}%  DNABERT-2 {100*d.mean():.2f}%")
    print(f"\n{'max target-neighbour identity':32}{'pairs':>7}{'blocks':>8}"
          f"{'k-mer':>9}{'DNABERT-2':>11}{'diff pp':>9}")
    for bn in bins:
        print(f"{bn['bin']:32}{bn['n_pairs']:7}{bn['n_blocks']:8}"
              f"{100*bn['C_kmer']:9.2f}{100*bn['dnabert2']:11.2f}{bn['difference_pp']:+9.2f}")
    print(f"\nSpearman accuracy vs max identity: k-mer rho={rho_k.statistic:+.3f} "
          f"(p={rho_k.pvalue:.2e}), DNABERT-2 rho={rho_d.statistic:+.3f} "
          f"(p={rho_d.pvalue:.2e})")
    print(f"pairs where the neighbour shares a haplotype with the target: "
          f"{summary['pairs_sharing_a_haplotype_with_target']} of {len(ok)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
