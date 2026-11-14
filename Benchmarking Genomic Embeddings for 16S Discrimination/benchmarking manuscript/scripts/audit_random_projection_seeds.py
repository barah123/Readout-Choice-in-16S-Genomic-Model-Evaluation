"""AQ24 / check 3: how much does C_random depend on the random projection that defines it?

`score_test.py` re-seeds the generator inside each outer fold, so the locked test evaluates
C_random at ONE projection, not an average over projections. This script re-runs exactly that
arm over many projection seeds, holding everything that does not depend on the projection
fixed: the frozen outer split, the inner-CV region choice (C_kmer only, training isolates
only) and the k-mer features themselves are computed once per fold and reused across seeds.

Reports the per-block and across-block spread of C_random, and the spread of the derived
quantity the manuscript interprets, best-of-four R2 minus C_random.

Run from this directory.  conda run -n glm-benchmark python audit_random_projection_seeds.py
"""
import json
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import StratifiedKFold

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
sys.path.insert(0, os.path.join(ROOT, "phase4_locked_test"))
import run_readouts_full as rr          # noqa: E402
import score_test as st                 # noqa: E402

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEEDS = list(range(20))
TEST = os.path.join(ROOT, "phase4_locked_test")


def block_c_random(per_region, meta):
    """C_random per seed for one block, replicating score_test.score_block exactly."""
    labels = {}
    for region, d in per_region.items():
        for iso, (_, lab) in d.items():
            labels[iso] = lab
    isolates = sorted(labels)
    y_all = np.array([labels[i] for i in isolates])
    if len(np.unique(y_all)) < 2 or min(np.bincount(y_all)) < st.N_SPLITS:
        return None
    eligible = [r for r in st.REGIONS
                if meta.get("_eligible_regions") and r in meta["_eligible_regions"]]
    if not eligible:
        return None

    skf = StratifiedKFold(n_splits=st.N_SPLITS, shuffle=True, random_state=st.SEED)
    per_seed = defaultdict(list)
    for tr_idx, te_idx in skf.split(np.zeros(len(isolates)), y_all):
        tr_iso = {isolates[i] for i in tr_idx}
        te_iso = {isolates[i] for i in te_idx}

        best_r, best_s = None, -np.inf
        for r in eligible:
            d = per_region[r]
            sl = [d[i][0] for i in isolates if i in tr_iso and i in d]
            yy = [d[i][1] for i in isolates if i in tr_iso and i in d]
            if len(set(yy)) < 2:
                continue
            s = st.kmer_score_inner(sl, yy)
            if s == s and s > best_s:
                best_s, best_r = s, r
        if best_r is None:
            continue

        d = per_region[best_r]
        tr = [i for i in isolates if i in tr_iso and i in d]
        te = [i for i in isolates if i in te_iso and i in d]
        ytr = np.array([d[i][1] for i in tr])
        yte = np.array([d[i][1] for i in te])
        if len(set(ytr)) < 2 or len(set(yte)) < 2:
            continue

        # k-mer features do not depend on the projection: compute once, reuse per seed
        Ktr, Kte = rr.kmer_fit_transform([d[i][0] for i in tr], [d[i][0] for i in te])
        p = Ktr.shape[1]
        for seed in SEEDS:
            rp = np.random.default_rng(seed).standard_normal((p, 256)) / np.sqrt(max(p, 1))
            c = st.probe()
            c.fit(Ktr @ rp, ytr)
            per_seed[seed].append(rr.balanced_accuracy(yte, c.predict(Kte @ rp)))

    if not per_seed:
        return None
    return {s: float(np.mean(v)) for s, v in per_seed.items()}


def main():
    rows = json.load(open(os.path.join(TEST, "test_blocks.json")))
    summ = json.load(open(os.path.join(TEST, "test_block_summary.json")))
    scores = json.load(open(os.path.join(TEST, "locked_test_scores.json")))

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    out = {}
    for bid in sorted(by_block):
        s = summ.get(bid, {})
        elig = [r for r, v in s.get("regions", {}).items() if v.get("eligible")]
        meta = {"_eligible_regions": elig}
        res = block_c_random(by_block[bid], meta)
        if res is None:
            continue
        vals = np.array([res[k] for k in SEEDS])
        out[bid] = {"per_seed": res, "mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                    "min": float(vals.min()), "max": float(vals.max()),
                    "frozen_seed0": res[0],
                    "reported": scores.get(bid, {}).get("C_random")}
        print(f"{bid:28} seed0={res[0]:.4f} reported={out[bid]['reported']!s:.6}  "
              f"mean={vals.mean():.4f} sd={vals.std(ddof=1):.4f} "
              f"range=[{vals.min():.4f}, {vals.max():.4f}]", flush=True)

    blocks = sorted(out)
    M = np.array([[out[b]["per_seed"][s] for b in blocks] for s in SEEDS])  # seeds x blocks
    per_seed_mean = M.mean(axis=1)
    agg = {
        "n_blocks": len(blocks), "n_seeds": len(SEEDS),
        "across_block_mean_per_seed": {s: float(per_seed_mean[i]) for i, s in enumerate(SEEDS)},
        "mean_of_means": float(per_seed_mean.mean()),
        "sd_of_means": float(per_seed_mean.std(ddof=1)),
        "min_of_means": float(per_seed_mean.min()),
        "max_of_means": float(per_seed_mean.max()),
        "frozen_seed0_mean": float(per_seed_mean[SEEDS.index(0)]),
        "max_within_block_sd": float(max(out[b]["sd"] for b in blocks)),
        "median_within_block_sd": float(np.median([out[b]["sd"] for b in blocks])),
    }

    # the derived quantity the manuscript interprets
    MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
    best = {}
    for b in blocks:
        v = [scores[b].get(f"{m}:R2") for m in MODELS]
        v = [x for x in v if x is not None and x == x]
        if v:
            best[b] = max(v)
    if best:
        bb = [b for b in blocks if b in best]
        diff = np.array([[best[b] - out[b]["per_seed"][s] for b in bb] for s in SEEDS])
        dm = diff.mean(axis=1)
        agg["best_r2_minus_c_random"] = {
            "mean_over_seeds": float(dm.mean()), "sd_over_seeds": float(dm.std(ddof=1)),
            "min": float(dm.min()), "max": float(dm.max()),
            "at_frozen_seed0": float(dm[SEEDS.index(0)]),
        }

    print("\n" + json.dumps(agg, indent=2))
    dest = os.path.abspath(os.path.join(HERE, "..", "results",
                                        "random_projection_seed_audit.json"))
    json.dump({"per_block": out, "aggregate": agg}, open(dest, "w"), indent=2)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
