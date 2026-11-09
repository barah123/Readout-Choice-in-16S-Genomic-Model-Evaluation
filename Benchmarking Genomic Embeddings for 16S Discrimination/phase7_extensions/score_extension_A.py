"""Extension A: embed the hard panels and run the frozen analysis.

Reuses score_test.score_block unchanged - the exact function that produced the Phase 5
confirmatory result - so the only thing that differs is which panels are scored.

Run:  conda run -n glm-benchmark python score_extension_A.py embed dnabert2 dnabert_s hyenadna
      conda run -n glm-ntv2      python score_extension_A.py embed ntv2
      conda run -n glm-benchmark python score_extension_A.py score
"""
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
sys.path.insert(0, os.path.abspath("../phase4_locked_test"))

MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
REGIONS = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"]
MARGIN = 0.010


def do_embed(models):
    from embed_all_models import embed_unique, load
    rows = json.load(open("hard_blocks.json"))
    by_region = defaultdict(set)
    for r in rows:
        by_region[r["region"]].add(r["slice"])
    print({k: len(v) for k, v in sorted(by_region.items())}, flush=True)
    for name in models:
        tok, model, spec = load(name)
        store = {}
        t0 = time.time()
        for region in REGIONS:
            uniq = sorted(by_region.get(region, []))
            if not uniq:
                continue
            E = embed_unique(uniq, tok, model, spec)
            assert np.isfinite(E).all(), f"{name} {region}: non-finite"
            store[f"{region}__keys"] = np.array(uniq, dtype=object)
            store[f"{region}__vecs"] = E
            print(f"{name} {region}: {len(uniq)} unique ({time.time()-t0:.0f}s)", flush=True)
        np.savez_compressed(f"emb_hard_{name}.npz", **store)
        print(f"saved emb_hard_{name}.npz\n", flush=True)


def do_score():
    from score_test import score_block
    from joblib import Parallel, delayed

    rows = json.load(open("hard_blocks.json"))
    summ = json.load(open("hard_block_summary.json"))

    emb = {}
    for m in MODELS:
        f = f"emb_hard_{m}.npz"
        if not os.path.exists(f):
            print(f"[skip] {f} missing", file=sys.stderr)
            continue
        z = np.load(f, allow_pickle=True)
        for r in REGIONS:
            if f"{r}__keys" in z:
                emb[(m, r)] = {k: v for k, v in zip(z[f"{r}__keys"], z[f"{r}__vecs"])}
    if not emb:
        sys.exit("no embeddings")

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    jobs = []
    for bid, per_region in by_block.items():
        info = summ.get(bid, {})
        meta = {"_family": info.get("family"), "_rank": info.get("rank"),
                "_expected_class": "HARD",
                "_eligible_regions": info.get("eligible_regions", []),
                "_block_eligible": info.get("block_eligible", False),
                "_screen_kmer": info.get("screen_kmer"),
                "_ceilings": {r: v.get("ceiling_bal_isolate")
                              for r, v in info.get("regions", {}).items()}}
        jobs.append((bid, per_region, emb, meta))

    print(f"scoring {len(jobs)} hard panels under the frozen plan", flush=True)
    out = Parallel(n_jobs=4, verbose=5)(delayed(score_block)(*j) for j in jobs)
    res = {r["_block"]: r for r in out}
    json.dump(res, open("hard_scores.json", "w"), indent=2, default=float)

    blocks = {k: v for k, v in res.items()
              if v.get("_status") == "ok" and v.get("_block_eligible")}
    names = sorted(blocks)
    n = len(names)

    print("\n" + "=" * 100)
    print("EXTENSION A - PRE-REGISTERED TEST ON SCREENED-HARD PANELS")
    print("=" * 100)
    print(f"eligible hard panels: n = {n}  (target >=12; <8 => underpowered)")
    if n < 8:
        print("*** UNDERPOWERED - reported as inconclusive per PREREGISTRATION.md ***")
    print()
    print(f"{'panel':34}{'ceil':>7}{'C_kmer':>9}{'DNABERT-2':>11}{'delta':>9}"
          f"{'best-of-4':>11}")
    for b in names:
        v = blocks[b]
        ch = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        ceils = [c for c in (v.get("_ceilings") or {}).values() if c == c]
        ceil = (v.get("_ceilings") or {}).get(ch, max(ceils) if ceils else float("nan"))
        best = max(v[f"{m}:R2"] for m in MODELS if f"{m}:R2" in v)
        print(f"{b:34}{ceil:7.3f}{v['C_kmer']:9.3f}{v['dnabert2:R2']:11.3f}"
              f"{v['dnabert2:R2']-v['C_kmer']:+9.3f}{best:11.3f}")

    def boot(d, seed=1729):
        rng = np.random.default_rng(seed)
        dr = d[rng.integers(0, len(d), size=(2000, len(d)))].mean(axis=1)
        return float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))

    prim = np.array([blocks[b]["dnabert2:R2"] - blocks[b]["C_kmer"] for b in names])
    lo, hi = boot(prim)
    print()
    print("PRIMARY (directional H1: models BETTER on hard panels)")
    print(f"  mean DNABERT-2:R2 - C_kmer = {prim.mean():+.5f}")
    print(f"  95% bootstrap CI           = [{lo:+.5f}, {hi:+.5f}]")
    print(f"  panels favouring model / classical / tie(+-0.01): "
          f"{(prim>0.01).sum()} / {(prim<-0.01).sum()} / "
          f"{((prim>=-0.01)&(prim<=0.01)).sum()}")
    if lo > 0:
        verdict = "H1 SUPPORTED - models better on hard panels"
    elif hi < 0:
        verdict = "models WORSE on hard panels"
    else:
        verdict = "H1 NOT SUPPORTED - CI contains zero"
    print(f"  VERDICT: {verdict}")

    bestd = np.array([max(blocks[b][f"{m}:R2"] for m in MODELS) - blocks[b]["C_kmer"]
                      for b in names])
    blo, bhi = boot(bestd, 1730)
    print(f"\n  secondary, best-of-4 minus C_kmer = {bestd.mean():+.5f} "
          f"[{blo:+.5f}, {bhi:+.5f}]")

    # leave-one-out stability, the check Extension B showed is essential
    loo = [np.delete(prim, i).mean() for i in range(len(prim))]
    print(f"\n  leave-one-out mean range: [{min(loo):+.5f}, {max(loo):+.5f}]  "
          f"sign stable: {not (min(loo) < 0 < max(loo))}")

    json.dump({"n": n, "mean": float(prim.mean()), "ci": [lo, hi], "verdict": verdict,
               "best_of_4": {"mean": float(bestd.mean()), "ci": [blo, bhi]},
               "loo_range": [float(min(loo)), float(max(loo))],
               "per_panel": {b: {"C_kmer": blocks[b]["C_kmer"],
                                 "dnabert2": blocks[b]["dnabert2:R2"],
                                 "delta": blocks[b]["dnabert2:R2"] - blocks[b]["C_kmer"]}
                             for b in names}},
              open("extension_A_results.json", "w"), indent=2, default=float)
    print("\nwrote extension_A_results.json")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "score"
    if mode == "embed":
        do_embed(sys.argv[2:] or ["dnabert2"])
    else:
        do_score()
