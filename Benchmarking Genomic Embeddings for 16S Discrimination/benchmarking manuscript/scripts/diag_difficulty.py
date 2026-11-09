"""Phase 6 diagnostic A (EXPLORATORY): does the equivalence hold on hard blocks,
or is it driven by easy blocks where every method saturates near the ceiling?

This is the main threat to interpreting Phase 5. If nearly all blocks sit at a ceiling of
~0.99 then no method can differ, and "equivalent" would say more about the panels than
about the representations.
"""
import json

import numpy as np

SCORES = "../phase4_locked_test/locked_test_scores.json"
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]

scores = json.load(open(SCORES))
blocks = {k: v for k, v in scores.items()
          if v.get("_status") == "ok" and v.get("_block_eligible")}

rows = []
for b, v in blocks.items():
    ceils = [c for c in (v.get("_ceilings") or {}).values() if c == c]
    from collections import Counter
    chosen = Counter(v["_chosen_regions"]).most_common(1)[0][0]
    ceil = (v.get("_ceilings") or {}).get(chosen, max(ceils) if ceils else float("nan"))
    ck = v["C_kmer"]
    rows.append(dict(block=b, ceiling=ceil, kmer=ck,
                     headroom=ceil - ck,
                     delta=v["dnabert2:R2"] - ck,
                     best_delta=max(v[f"{m}:R2"] for m in MODELS) - ck,
                     err_kmer=1 - ck))

rows.sort(key=lambda r: -r["err_kmer"])

print("=" * 96)
print("EXPLORATORY DIAGNOSTIC A - is equivalence an artifact of easy blocks?")
print("=" * 96)
print(f"{'block':26}{'ceiling':>9}{'C_kmer':>9}{'kmer error':>12}{'headroom':>10}"
      f"{'delta(D2)':>11}{'delta(best)':>12}")
for r in rows:
    print(f"{r['block']:26}{r['ceiling']:9.3f}{r['kmer']:9.3f}{r['err_kmer']:12.3f}"
          f"{r['headroom']:10.3f}{r['delta']:+11.3f}{r['best_delta']:+12.3f}")

err = np.array([r["err_kmer"] for r in rows])
d = np.array([r["delta"] for r in rows])
bd = np.array([r["best_delta"] for r in rows])

print()
print(f"blocks where classical error > 5% (real headroom): "
      f"{(err > 0.05).sum()} / {len(err)}")
print(f"blocks where classical error < 1% (saturated)   : "
      f"{(err < 0.01).sum()} / {len(err)}")

hard = err > 0.05
easy = err < 0.01
mid = ~hard & ~easy


def boot(x, seed=99):
    if len(x) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    dr = x[rng.integers(0, len(x), size=(2000, len(x)))].mean(axis=1)
    return float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))


print()
print("mean DNABERT-2 minus C_kmer, stratified by how much room there is to improve:")
for name, m in [("HARD   (kmer error >5%)", hard),
                ("MIDDLE (1-5%)", mid),
                ("EASY   (kmer error <1%)", easy)]:
    if m.sum():
        lo, hi = boot(d[m])
        print(f"  {name:26} n={m.sum():2}  mean={d[m].mean():+.4f}  "
              f"CI=[{lo:+.4f}, {hi:+.4f}]")

print()
print("Spearman(classical error, delta) - does the gap grow where the task is harder?")
from scipy.stats import spearmanr
rho, p = spearmanr(err, d)
print(f"  DNABERT-2 : rho={rho:+.3f}  p={p:.3f}")
rho2, p2 = spearmanr(err, bd)
print(f"  best-of-4 : rho={rho2:+.3f}  p={p2:.3f}")

print()
print("Interpretation guide (stated before reading the numbers):")
print("  - If the HARD stratum also shows ~0 difference, equivalence is a real finding.")
print("  - If HARD shows a large |difference| hidden by saturated EASY blocks, the pooled")
print("    equivalence verdict is misleading and must be qualified.")
