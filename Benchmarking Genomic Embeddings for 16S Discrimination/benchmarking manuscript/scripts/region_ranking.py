"""Task T2: region ranking — can a method pick the right 16S region in advance?

Distinct from the classification task (T1) reported everywhere else. T1 asks "how accurate
are you AT a given region". T2 asks "which region should be sequenced at all". A method can
be good at one and poor at the other.

Procedure, per panel, per outer fold (protocol §13 and the pseudocode in
04_Implementation_Schemas_and_Checklists.md §4):

  1. Using TRAINING isolates only, each arm scores every eligible region by inner CV.
  2. Each arm commits to its own top-ranked region. No held-out label is involved.
  3. That arm is evaluated on the held-out isolates at its chosen region.
  4. Separately, every region is evaluated on the held-out isolates to find the best one in
     hindsight (the oracle). Used for metrics ONLY, never for selection.

  regret = oracle best held-out accuracy − held-out accuracy of the region the arm chose

Reported: held-out accuracy of the chosen region, regret, and how often the arm's choice
matches the hindsight-best region. Same isolate-grouped splits, same probe, same seed as
every other analysis.
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
import run_readouts_full as rr

SEED = 0
N_OUTER = 5
N_INNER = 3
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
ARMS = ["C_kmer"] + [f"{m}:R2" for m in MODELS]

SETS = {
    "locked_test": dict(
        blocks=os.path.join(ROOT, "phase4_locked_test", "test_blocks.json"),
        summary=os.path.join(ROOT, "phase4_locked_test", "test_block_summary.json"),
        emb=os.path.join(ROOT, "phase4_locked_test", "emb_test_{m}.npz")),
    "hard_panels": dict(
        blocks=os.path.join(ROOT, "phase7_extensions", "hard_blocks.json"),
        summary=os.path.join(ROOT, "phase7_extensions", "hard_block_summary.json"),
        emb=os.path.join(ROOT, "phase7_extensions", "emb_hard_{m}.npz")),
}


def probe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=2000, class_weight="balanced",
                                            random_state=SEED))


def fit_score(Xtr, ytr, Xte, yte):
    p = probe(); p.fit(Xtr, ytr)
    return rr.balanced_accuracy(yte, p.predict(Xte))


def features(arm, slices_tr, slices_te, lut):
    """Return (Xtr, Xte) for one arm. k-mer vocabulary is fit on training only."""
    if arm == "C_kmer":
        return rr.kmer_fit_transform(slices_tr, slices_te)
    return (np.vstack([lut[s] for s in slices_tr]),
            np.vstack([lut[s] for s in slices_te]))


def eval_panel(bid, per_region, eligible_regions, emb):
    y_by_iso = {}
    for region, d in per_region.items():
        for iso, (_, lab) in d.items():
            y_by_iso[iso] = lab
    isolates = sorted(y_by_iso)
    y_all = np.array([y_by_iso[i] for i in isolates])
    if len(np.unique(y_all)) < 2 or min(np.bincount(y_all)) < N_OUTER:
        return None
    regions = [r for r in eligible_regions if r in per_region]
    if len(regions) < 2:
        return None

    skf = StratifiedKFold(n_splits=N_OUTER, shuffle=True, random_state=SEED)
    acc = defaultdict(list)       # arm -> held-out acc of its chosen region
    regret = defaultdict(list)
    hit = defaultdict(list)       # arm picked the hindsight-best region?
    picks = defaultdict(list)

    for tr_i, te_i in skf.split(np.zeros(len(isolates)), y_all):
        tr_iso = {isolates[i] for i in tr_i}
        te_iso = {isolates[i] for i in te_i}

        # held-out accuracy of every region for every arm (oracle table; metrics only)
        heldout = defaultdict(dict)
        inner = defaultdict(dict)
        for region in regions:
            d = per_region[region]
            tr = [i for i in isolates if i in tr_iso and i in d]
            te = [i for i in isolates if i in te_iso and i in d]
            ytr = np.array([d[i][1] for i in tr])
            yte = np.array([d[i][1] for i in te])
            if len(set(ytr)) < 2 or len(set(yte)) < 2:
                continue
            s_tr = [d[i][0] for i in tr]
            s_te = [d[i][0] for i in te]

            for arm in ARMS:
                lut = emb.get((arm.split(":")[0], region)) if arm != "C_kmer" else None
                if arm != "C_kmer" and lut is None:
                    continue
                try:
                    Xtr, Xte = features(arm, s_tr, s_te, lut)
                except KeyError:
                    continue
                heldout[arm][region] = fit_score(Xtr, ytr, Xte, yte)

                # inner CV on TRAINING isolates only -> this arm's own ranking
                if min(np.bincount(ytr)) < N_INNER:
                    continue
                isk = StratifiedKFold(n_splits=N_INNER, shuffle=True, random_state=SEED)
                sc = []
                for a, b in isk.split(np.zeros(len(ytr)), ytr):
                    ia = [s_tr[i] for i in a]; ib = [s_tr[i] for i in b]
                    try:
                        A, B = features(arm, ia, ib, lut)
                    except KeyError:
                        continue
                    sc.append(fit_score(A, ytr[a], B, ytr[b]))
                if sc:
                    inner[arm][region] = float(np.mean(sc))

        for arm in ARMS:
            if not inner.get(arm) or not heldout.get(arm):
                continue
            chosen = max(inner[arm], key=inner[arm].get)
            if chosen not in heldout[arm]:
                continue
            best_region = max(heldout[arm], key=heldout[arm].get)
            got = heldout[arm][chosen]
            best = heldout[arm][best_region]
            acc[arm].append(got)
            regret[arm].append(best - got)
            hit[arm].append(float(chosen == best_region))
            picks[arm].append(chosen)

    if not acc:
        return None
    return {
        "block": bid,
        "n_regions": len(regions),
        **{f"{a}|acc": float(np.mean(acc[a])) for a in acc},
        **{f"{a}|regret": float(np.mean(regret[a])) for a in regret},
        **{f"{a}|hit": float(np.mean(hit[a])) for a in hit},
        **{f"{a}|picks": Counter(picks[a]).most_common() for a in picks},
    }


def run_set(name, cfg):
    rows = json.load(open(cfg["blocks"]))
    summ = json.load(open(cfg["summary"]))
    emb = {}
    for m in MODELS:
        f = cfg["emb"].format(m=m)
        if not os.path.exists(f):
            continue
        z = np.load(f, allow_pickle=True)
        for key in z.files:
            if key.endswith("__keys"):
                r = key.split("__")[0]
                emb[(m, r)] = {k: v for k, v in zip(z[key], z[f"{r}__vecs"])}

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    jobs = []
    for bid, per_region in by_block.items():
        info = summ.get(bid, {})
        if not info.get("block_eligible"):
            continue
        jobs.append((bid, per_region, info.get("eligible_regions", []), emb))

    print(f"\n=== {name}: {len(jobs)} panels ===", flush=True)
    out = [r for r in Parallel(n_jobs=4, verbose=5)(
        delayed(eval_panel)(*j) for j in jobs) if r]
    json.dump(out, open(os.path.join(HERE, "..", "results",
                                     f"region_ranking_{name}.json"), "w"),
              indent=2, default=str)

    print(f"\n{'arm':14}{'chosen-region acc':>19}{'regret':>10}{'picked best':>13}")
    summary = {}
    for arm in ARMS:
        a = [r[f"{arm}|acc"] for r in out if f"{arm}|acc" in r]
        g = [r[f"{arm}|regret"] for r in out if f"{arm}|regret" in r]
        h = [r[f"{arm}|hit"] for r in out if f"{arm}|hit" in r]
        if not a:
            continue
        summary[arm] = dict(n=len(a), acc=float(np.mean(a)),
                            regret=float(np.mean(g)), hit=float(np.mean(h)))
        print(f"{arm:14}{np.mean(a):>19.4f}{np.mean(g):>10.4f}{np.mean(h):>13.1%}")

    # paired contrast on regret, the metric that defines this task
    if "C_kmer" in summary and "dnabert2:R2" in summary:
        common = [r for r in out if "C_kmer|regret" in r and "dnabert2:R2|regret" in r]
        d = np.array([r["dnabert2:R2|regret"] - r["C_kmer|regret"] for r in common])
        rng = np.random.default_rng(1729)
        dr = d[rng.integers(0, len(d), size=(2000, len(d)))].mean(axis=1)
        lo, hi = np.percentile(dr, [2.5, 97.5])
        print(f"\nregret difference (DNABERT-2 − k-mer), n={len(d)} panels: "
              f"{d.mean():+.4f} [{lo:+.4f}, {hi:+.4f}]")
        print("  negative = the model picks better regions")
        summary["_contrast_regret"] = dict(n=len(d), mean=float(d.mean()), ci=[lo, hi])
    return summary


if __name__ == "__main__":
    allsum = {}
    for name, cfg in SETS.items():
        allsum[name] = run_set(name, cfg)
    json.dump(allsum, open(os.path.join(HERE, "..", "results",
                                        "region_ranking_summary.json"), "w"),
              indent=2, default=float)
    print("\nwrote results/region_ranking_*.json")
