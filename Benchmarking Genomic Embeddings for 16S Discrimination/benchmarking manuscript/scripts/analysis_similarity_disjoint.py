"""POST HOC: 99%-identity similarity-disjoint evaluation, with a direct cross-fold audit.

NOT part of the frozen confirmatory analysis.

Operational definition, stated in full because the guarantee depends on it:

  identity      the project's canonical measure, 2 * matches / (len_a + len_b), matches
                counted over aligned blocks only, from Biopython global pairwise alignment
                with open_gap_score -5, extend_gap_score -2, match 2, mismatch -1. Reused
                verbatim from run_readouts_full.identity_dist, the same function C_align uses.
  coverage      all sequences compared are slices of the SAME region, so alignment spans the
                full length of both; per-block length spread is reported.
  strand        slices are forward-oriented as extracted. The 22 minus-strand records in the
                corpus (0.09%) are a disclosed, unrepaired pipeline defect and are NOT
                reoriented here, so this analysis inherits that defect.
  ambiguity     N is the only non-ACGT character present. The aligner scores it as a
                mismatch, so an N-bearing slice is LESS similar to its neighbours and is
                therefore, if anything, split off rather than merged. The direction is
                conservative for disjointness.
  clustering    SINGLE LINKAGE: connected components of the graph whose edges are pairs at
                identity >= 0.99. This is the choice that gives the guarantee. Greedy or
                centroid methods (CD-HIT, UCLUST) do NOT ensure every cross-cluster pair
                falls below the threshold; single linkage does, because two sequences in
                different components are by construction not joined by any edge.

Audit, not assumption: after splitting, the maximum identity between any training and any
test sequence is computed per fold and must fall below the threshold. Reported per block.

A threshold caveat that matters here and is reported per block: a fixed identity threshold is
NOT a fixed mismatch budget across regions of 31-107 nt. At a 31 nt V1 slice a single
mismatch gives 96.8% identity, below 99%, so 99%-clustering collapses to exact clustering; at
a 107 nt V4 slice a single mismatch gives 99.1%, above the threshold, so such pairs merge.
The effective mismatch tolerance is printed for every block.

Writes results/similarity_disjoint_results.json.
"""
import json
import os
import sys
import warnings
from collections import Counter, defaultdict

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
import run_readouts_full as rr

warnings.filterwarnings("ignore", category=ConvergenceWarning)

SEED, N_SPLITS, THRESH = 0, 5, 0.99
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "similarity_disjoint_results.json"))


def probe():
    return make_pipeline(StandardScaler(), LogisticRegression(
        max_iter=2000, class_weight="balanced", random_state=SEED))


def components(adj, n):
    """Connected components = single-linkage clusters."""
    comp, cid = [-1] * n, 0
    for s in range(n):
        if comp[s] != -1:
            continue
        stack, comp[s] = [s], cid
        while stack:
            u = stack.pop()
            for v in np.nonzero(adj[u])[0]:
                if comp[v] == -1:
                    comp[v] = cid
                    stack.append(int(v))
        cid += 1
    return comp, cid


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

    results, skipped = {}, {}
    for b in sorted(scores):
        v = scores[b]
        if v.get("_status") != "ok" or not v.get("_block_eligible"):
            continue
        region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        d = per[b][region]
        isolates = sorted(d)
        slices = [d[i][0] for i in isolates]
        y = np.array([d[i][1] for i in isolates])

        uniq = sorted(set(slices))
        D = rr.identity_dist(uniq)                 # distance = 1 - identity
        I = 1.0 - D
        np.fill_diagonal(I, 1.0)
        adj = (I >= THRESH)
        np.fill_diagonal(adj, False)
        comp, n_clust = components(adj, len(uniq))
        cl_of = {s: comp[k] for k, s in enumerate(uniq)}
        groups = np.array([cl_of[s] for s in slices])

        L = [len(s) for s in uniq]
        # effective mismatch tolerance at this length: largest m with 2(L-m)/2L >= THRESH
        Lmed = int(np.median(L))
        tol = int(np.floor(Lmed * (1 - THRESH)))

        cl_cls = defaultdict(lambda: [0, 0])
        for g, lab in zip(groups, y):
            cl_cls[g][0 if lab == 1 else 1] += 1
        n_t = sum(1 for c in cl_cls.values() if c[0] > 0)
        n_n = sum(1 for c in cl_cls.values() if c[1] > 0)
        base = dict(region=region, n_isolates=len(isolates), n_haplotypes=len(uniq),
                    n_clusters_99=n_clust, clusters_target=n_t, clusters_neighbour=n_n,
                    median_slice_len=Lmed, mismatch_tolerance_at_99pct=tol,
                    collapses_to_exact=bool(tol == 0 and n_clust == len(uniq)))
        if n_t < N_SPLITS or n_n < N_SPLITS:
            skipped[b] = {**base, "reason": f"fewer than {N_SPLITS} clusters in a class"}
            print(f"  [skip] {b:26} {region} clusters={n_clust:3} "
                  f"(T {n_t}, N {n_n})  tol={tol}", flush=True)
            continue

        sp = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
        idx_of = {s: k for k, s in enumerate(uniq)}
        arms, max_cross, n_folds = defaultdict(list), 0.0, 0
        for tr, te in sp.split(np.zeros(len(y)), y, groups):
            ytr, yte = y[tr], y[te]
            if len(set(ytr)) < 2 or len(set(yte)) < 2:
                continue
            n_folds += 1
            str_tr = [slices[i] for i in tr]
            str_te = [slices[i] for i in te]
            ti = sorted({idx_of[s] for s in str_tr})
            ei = sorted({idx_of[s] for s in str_te})
            max_cross = max(max_cross, float(I[np.ix_(ti, ei)].max()))

            Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
            c = probe(); c.fit(Ktr, ytr)
            arms["C_kmer"].append(rr.balanced_accuracy(yte, c.predict(Kte)))
            for m in MODELS:
                lut = emb.get((m, region))
                if lut is None or not all(s in lut for s in str_tr + str_te):
                    continue
                p = probe()
                p.fit(np.vstack([lut[s] for s in str_tr]), ytr)
                arms[f"{m}:R2"].append(
                    rr.balanced_accuracy(yte, p.predict(np.vstack([lut[s] for s in str_te]))))

        res = {k: float(np.mean(x)) for k, x in arms.items() if x}
        res["_n_folds"] = n_folds
        res["_max_cross_fold_identity"] = max_cross
        res["_audit_below_threshold"] = bool(max_cross < THRESH)
        results[b] = {**base, "similarity_disjoint": res}
        print(f"  {b:26} {region} clusters={n_clust:3} (T {n_t}, N {n_n}) tol={tol}  "
              f"kmer={res['C_kmer']:.4f} d2={res.get('dnabert2:R2', float('nan')):.4f}  "
              f"maxXid={max_cross:.4f} audit={'PASS' if res['_audit_below_threshold'] else 'FAIL'}",
              flush=True)

    names = sorted(results)
    summary = {"threshold": THRESH, "linkage": "single (connected components)",
               "n_evaluable": len(results), "n_not_evaluable": len(skipped)}
    if names:
        x = np.array([results[b]["similarity_disjoint"]["dnabert2:R2"]
                      - results[b]["similarity_disjoint"]["C_kmer"] for b in names])
        rng = np.random.default_rng(1729)
        dr = x[rng.integers(0, len(x), size=(20000, len(x)))].mean(axis=1)
        lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
        summary["primary_contrast_dnabert2_minus_kmer"] = {
            "n": len(x), "mean": float(x.mean()), "sd": float(x.std(ddof=1)),
            "ci95": [lo, hi],
            "verdict": ("EQUIVALENT" if (lo > -0.010 and hi < 0.010)
                        else "INCONCLUSIVE at this margin")}
        summary["mean_accuracy"] = {
            a: float(np.mean([results[b]["similarity_disjoint"][a] for b in names]))
            for a in ("C_kmer", *[f"{m}:R2" for m in MODELS])}
        summary["audit"] = {
            "max_cross_fold_identity_over_all_blocks": float(max(
                results[b]["similarity_disjoint"]["_max_cross_fold_identity"] for b in names)),
            "all_blocks_below_threshold": all(
                results[b]["similarity_disjoint"]["_audit_below_threshold"] for b in names)}
    json.dump({"summary": summary, "blocks": results, "not_evaluable": skipped},
              open(OUT, "w"), indent=2, default=float)
    print(f"\n  evaluable {len(results)} / not evaluable {len(skipped)}")
    if names:
        c = summary["primary_contrast_dnabert2_minus_kmer"]
        print(f"  DNABERT-2 - k-mer: {100*c['mean']:+.2f} pp "
              f"[{100*c['ci95'][0]:+.2f}, {100*c['ci95'][1]:+.2f}] "
              f"SD {100*c['sd']:.2f}  {c['verdict']}")
        print(f"  AUDIT max cross-fold identity {summary['audit']['max_cross_fold_identity_over_all_blocks']:.4f} "
              f"(threshold {THRESH}) -> all below: {summary['audit']['all_blocks_below_threshold']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
