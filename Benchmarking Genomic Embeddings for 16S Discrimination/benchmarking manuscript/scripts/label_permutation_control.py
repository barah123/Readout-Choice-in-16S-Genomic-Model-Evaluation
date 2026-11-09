"""Label-permutation control — the test that decides what C_random actually showed.

WHY THIS EXISTS
---------------
The report previously described C_random as k-mer counts "scrambled through a random
projection, destroying their biological meaning", and concluded from its high score that
most measured accuracy reflects classifier flexibility rather than biological content.

That reasoning is wrong, and a reviewer caught it. A Gaussian random projection approximately
preserves pairwise distances (Johnson-Lindenstrauss), so C_random is a distance-preserving
dimensionality reduction of the k-mer features. It RETAINS k-mer information; it does not
destroy it. Its high score therefore shows only that the specific identity of individual
k-mer features is not needed — a statement about redundancy, not about biology.

The control that can actually test "is this accuracy real signal or classifier flexibility?"
is a label permutation run through the complete supervised pipeline: same features, same
probe, same cross-validation, labels shuffled. If permuted accuracy collapses to chance, the
observed accuracy reflects genuine sequence-label association.

Labels are permuted WITHIN each block, preserving class balance, at the isolate level — the
unit of independence used everywhere else.
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
N_PERM = 20          # permutations per block
ARMS = ["C_kmer", "C_random", "dnabert2:R2"]


def probe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=2000, class_weight="balanced",
                                            random_state=SEED))


def one_block(bid, slices, y, lut, n_perm=N_PERM):
    """Return observed and permuted accuracy for each arm on one block."""
    y = np.asarray(y)
    rng = np.random.default_rng(abs(hash(bid)) % (2 ** 32))
    out = {"block": bid}

    def run(labels):
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
        acc = defaultdict(list)
        for tr, te in skf.split(np.zeros(len(labels)), labels):
            ytr, yte = labels[tr], labels[te]
            if len(set(ytr)) < 2 or len(set(yte)) < 2:
                continue
            s_tr = [slices[i] for i in tr]
            s_te = [slices[i] for i in te]
            Ktr, Kte = rr.kmer_fit_transform(s_tr, s_te)
            c = probe(); c.fit(Ktr, ytr)
            acc["C_kmer"].append(rr.balanced_accuracy(yte, c.predict(Kte)))
            # random projection, fit on TRAINING feature space only
            rp = np.random.default_rng(SEED).standard_normal(
                (Ktr.shape[1], 256)) / np.sqrt(max(Ktr.shape[1], 1))
            cr = probe(); cr.fit(Ktr @ rp, ytr)
            acc["C_random"].append(rr.balanced_accuracy(yte, cr.predict(Kte @ rp)))
            if lut is not None:
                try:
                    Xtr = np.vstack([lut[s] for s in s_tr])
                    Xte = np.vstack([lut[s] for s in s_te])
                    p = probe(); p.fit(Xtr, ytr)
                    acc["dnabert2:R2"].append(
                        rr.balanced_accuracy(yte, p.predict(Xte)))
                except KeyError:
                    pass
        return {k: float(np.mean(v)) for k, v in acc.items() if v}

    obs = run(y)
    for a, v in obs.items():
        out[f"{a}|observed"] = v

    perms = defaultdict(list)
    for _ in range(n_perm):
        yp = rng.permutation(y)          # preserves class balance exactly
        if len(np.unique(yp)) < 2:
            continue
        for a, v in run(yp).items():
            perms[a].append(v)
    for a, v in perms.items():
        out[f"{a}|perm_mean"] = float(np.mean(v))
        out[f"{a}|perm_p95"] = float(np.percentile(v, 95))
    return out


def main():
    rows = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_blocks.json")))
    scores = json.load(open(os.path.join(ROOT, "phase4_locked_test",
                                         "locked_test_scores.json")))
    blocks = {k: v for k, v in scores.items()
              if v.get("_status") == "ok" and v.get("_block_eligible")}

    emb = {}
    z = np.load(os.path.join(ROOT, "phase4_locked_test", "emb_test_dnabert2.npz"),
                allow_pickle=True)
    for key in z.files:
        if key.endswith("__keys"):
            r = key.split("__")[0]
            emb[r] = {k: v for k, v in zip(z[key], z[f"{r}__vecs"])}

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"],
                                                                    r["is_target"])
    jobs = []
    for bid, v in blocks.items():
        region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        d = by_block[bid][region]
        iso = sorted(d)
        y = np.array([d[i][1] for i in iso])
        if len(np.unique(y)) < 2 or min(np.bincount(y)) < 5:
            continue
        jobs.append((bid, [d[i][0] for i in iso], y, emb.get(region)))

    print(f"label-permutation control: {len(jobs)} blocks x {N_PERM} permutations",
          flush=True)
    out = Parallel(n_jobs=4, verbose=5)(delayed(one_block)(*j) for j in jobs)

    print("\n" + "=" * 92)
    print("OBSERVED vs PERMUTED-LABEL accuracy (locked test)")
    print("=" * 92)
    print(f"{'arm':14}{'observed':>11}{'permuted':>11}{'perm 95th':>12}{'gap':>9}")
    summary = {}
    for a in ARMS:
        o = [r[f"{a}|observed"] for r in out if f"{a}|observed" in r]
        p = [r[f"{a}|perm_mean"] for r in out if f"{a}|perm_mean" in r]
        q = [r[f"{a}|perm_p95"] for r in out if f"{a}|perm_p95" in r]
        if not o:
            continue
        summary[a] = dict(observed=float(np.mean(o)), permuted=float(np.mean(p)),
                          perm_p95=float(np.mean(q)),
                          gap=float(np.mean(o) - np.mean(p)), n=len(o))
        print(f"{a:14}{np.mean(o):>11.4f}{np.mean(p):>11.4f}{np.mean(q):>12.4f}"
              f"{np.mean(o)-np.mean(p):>+9.4f}")

    print("\nInterpretation:")
    print("  permuted ~0.50  -> observed accuracy reflects genuine sequence-label signal,")
    print("                     NOT classifier flexibility.")
    print("  permuted >>0.50 -> the pipeline can fit noise; the flexibility concern stands.")

    json.dump({"per_block": out, "summary": summary, "n_perm": N_PERM},
              open(os.path.join(HERE, "..", "results",
                                "label_permutation_control.json"), "w"),
              indent=2, default=float)
    print("\nwrote results/label_permutation_control.json")


if __name__ == "__main__":
    main()
