"""POST HOC: does the comparison change when each representation gets its own regularisation?

NOT part of the frozen confirmatory analysis.

The frozen plan fixes StandardScaler -> LogisticRegression(C=1.0) for every arm. That is a
legitimate FIXED-READOUT comparison and answers "how do these representations perform under
one specified readout". It does not answer "how well can each perform after comparable
optimisation", because L2 strength interacts with effective dimension and the arms differ by
an order of magnitude in width: canonical 4/5/6-mers give roughly 2,700 sparse dimensions,
against 768 (DNABERT-2, DNABERT-S), 1,024 (NTv2) and 256 (HyenaDNA) dense ones.

This runs the second comparison. Every arm gets the SAME grid, the SAME selection criterion
and the SAME budget, chosen by inner cross-validation on training folds only:

    C in {0.01, 0.1, 1, 10, 100}, inner 3-fold, criterion = balanced accuracy

Both outer splitters from the disjointness analyses are used, because the question "was the
fixed readout unfair" and the question "does tuning rescue the disjoint verdict" are
different and both matter. The fixed-C result is reported alongside, not replaced.

Writes results/tuned_readout_results.json.
"""
import json
import os
import sys
import warnings
from collections import Counter, defaultdict

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
import run_readouts_full as rr

warnings.filterwarnings("ignore", category=ConvergenceWarning)

SEED, N_SPLITS, N_INNER = 0, 5, 3
GRID = [0.01, 0.1, 1.0, 10.0, 100.0]
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "tuned_readout_results.json"))


def probe(C):
    return make_pipeline(StandardScaler(), LogisticRegression(
        C=C, max_iter=2000, class_weight="balanced", random_state=SEED))


def pick_C(X, y, groups):
    """Inner CV over the shared grid. Grouped when the outer split is grouped."""
    y = np.asarray(y)
    use_grp = groups is not None and len(set(groups)) >= N_INNER
    if use_grp:
        g = np.asarray(groups)
        ok = all(len(set(g[y == c])) >= N_INNER for c in set(y))
        if not ok:
            use_grp = False
    sp = (StratifiedGroupKFold(n_splits=N_INNER, shuffle=True, random_state=SEED)
          if use_grp else StratifiedKFold(n_splits=N_INNER, shuffle=True, random_state=SEED))
    it = list(sp.split(np.zeros(len(y)), y, np.asarray(groups)) if use_grp
              else sp.split(np.zeros(len(y)), y))
    best, bestC = -np.inf, 1.0
    for C in GRID:
        sc = []
        for tr, te in it:
            if len(set(y[tr])) < 2 or len(set(y[te])) < 2:
                continue
            m = probe(C); m.fit(X[tr], y[tr])
            sc.append(rr.balanced_accuracy(y[te], m.predict(X[te])))
        if sc and np.mean(sc) > best:
            best, bestC = float(np.mean(sc)), C
    return bestC


def run(slices, y, emb_at, splitter, groups=None):
    y = np.asarray(y)
    if splitter == "isolate":
        it = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED).split(
            np.zeros(len(y)), y)
    else:
        it = StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=SEED).split(
            np.zeros(len(y)), y, groups)

    fixed, tuned, chosen = defaultdict(list), defaultdict(list), defaultdict(list)
    for tr, te in it:
        ytr, yte = y[tr], y[te]
        if len(set(ytr)) < 2 or len(set(yte)) < 2:
            continue
        str_tr = [slices[i] for i in tr]
        str_te = [slices[i] for i in te]
        gtr = [groups[i] for i in tr] if groups is not None else None

        feats = {}
        Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)   # featurise once, reuse per C
        feats["C_kmer"] = (np.asarray(Ktr), np.asarray(Kte))
        for m in MODELS:
            lut = emb_at.get(m)
            if lut is not None and all(s in lut for s in str_tr + str_te):
                feats[f"{m}:R2"] = (np.vstack([lut[s] for s in str_tr]),
                                    np.vstack([lut[s] for s in str_te]))
        for arm, (Xtr, Xte) in feats.items():
            f = probe(1.0); f.fit(Xtr, ytr)
            fixed[arm].append(rr.balanced_accuracy(yte, f.predict(Xte)))
            C = pick_C(Xtr, ytr, gtr)
            chosen[arm].append(C)
            t = probe(C); t.fit(Xtr, ytr)
            tuned[arm].append(rr.balanced_accuracy(yte, t.predict(Xte)))

    return ({k: float(np.mean(v)) for k, v in fixed.items()},
            {k: float(np.mean(v)) for k, v in tuned.items()},
            {k: Counter(v).most_common() for k, v in chosen.items()})


def main():
    rows = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_blocks.json")))
    scores = json.load(open(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json")))
    sd = json.load(open(os.path.join(HERE, "..", "results", "sequence_disjoint_results.json")))
    evaluable = set(sd["blocks"])     # the 14 blocks that support haplotype-disjoint folds

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
    for b in sorted(evaluable):
        v = scores[b]
        region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        d = per[b][region]
        iso = sorted(d)
        slices = [d[i][0] for i in iso]
        y = [d[i][1] for i in iso]
        emb_at = {m: emb.get((m, region)) for m in MODELS}

        fi, ti, ci = run(slices, y, emb_at, "isolate")
        fs, ts, cs = run(slices, y, emb_at, "haplotype", groups=slices)
        out[b] = {"region": region,
                  "isolate_disjoint": {"fixed_C1": fi, "tuned": ti, "C_chosen": ci},
                  "sequence_disjoint": {"fixed_C1": fs, "tuned": ts, "C_chosen": cs}}
        print(f"  {b:26} iso fixed {fi['C_kmer']:.4f}/{fi['dnabert2:R2']:.4f} "
              f"tuned {ti['C_kmer']:.4f}/{ti['dnabert2:R2']:.4f} | "
              f"seq fixed {fs['C_kmer']:.4f}/{fs['dnabert2:R2']:.4f} "
              f"tuned {ts['C_kmer']:.4f}/{ts['dnabert2:R2']:.4f}", flush=True)

    names = sorted(out)

    def boot(x, seed=1729, n=20000):
        rng = np.random.default_rng(seed)
        dr = x[rng.integers(0, len(x), size=(n, len(x)))].mean(axis=1)
        return float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))

    summary = {"grid": GRID, "inner_folds": N_INNER, "n_blocks": len(names),
               "criterion": "balanced accuracy, inner CV on training folds only"}
    for split in ("isolate_disjoint", "sequence_disjoint"):
        s = {}
        for mode in ("fixed_C1", "tuned"):
            x = np.array([out[b][split][mode]["dnabert2:R2"] - out[b][split][mode]["C_kmer"]
                          for b in names])
            lo, hi = boot(x)
            s[mode] = {"mean": float(x.mean()), "sd": float(x.std(ddof=1)), "ci95": [lo, hi],
                       "verdict": ("EQUIVALENT" if (lo > -0.010 and hi < 0.010)
                                   else "INCONCLUSIVE at this margin"),
                       "mean_accuracy": {a: float(np.mean([out[b][split][mode][a]
                                                           for b in names]))
                                         for a in out[names[0]][split][mode]}}
        gain = {a: float(np.mean([out[b][split]["tuned"][a] - out[b][split]["fixed_C1"][a]
                                  for b in names]))
                for a in out[names[0]][split]["fixed_C1"]}
        s["tuning_gain_by_arm"] = gain
        picks = defaultdict(Counter)
        for b in names:
            for a, cc in out[b][split]["C_chosen"].items():
                for C, n in cc:
                    picks[a][C] += n
        s["C_selected_counts"] = {a: dict(sorted(c.items())) for a, c in picks.items()}
        summary[split] = s

    json.dump({"summary": summary, "blocks": out}, open(OUT, "w"), indent=2, default=float)

    pp = 100.0
    print("\nPRIMARY CONTRAST  DNABERT-2:R2 - C_kmer")
    print(f"{'split':20}{'readout':10}{'mean pp':>9}{'95% CI pp':>22}{'SD':>7}  verdict")
    for split, lab in (("isolate_disjoint", "isolate"), ("sequence_disjoint", "sequence")):
        for mode in ("fixed_C1", "tuned"):
            c = summary[split][mode]
            print(f"{lab:20}{mode:10}{pp*c['mean']:+9.2f}  "
                  f"[{pp*c['ci95'][0]:+7.2f},{pp*c['ci95'][1]:+7.2f}]"
                  f"{pp*c['sd']:7.2f}  {c['verdict']}")
    print("\nTUNING GAIN per arm (tuned minus fixed C=1.0, percentage points)")
    for split, lab in (("isolate_disjoint", "isolate-disjoint"),
                       ("sequence_disjoint", "sequence-disjoint")):
        g = summary[split]["tuning_gain_by_arm"]
        print(f"  {lab:20} " + "  ".join(f"{a}={pp*v:+.2f}" for a, v in sorted(g.items())))
    print("\nC SELECTED (counts across blocks x folds)")
    for split, lab in (("isolate_disjoint", "isolate-disjoint"),
                       ("sequence_disjoint", "sequence-disjoint")):
        print(f"  {lab}")
        for a, c in sorted(summary[split]["C_selected_counts"].items()):
            print(f"    {a:14} {c}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
