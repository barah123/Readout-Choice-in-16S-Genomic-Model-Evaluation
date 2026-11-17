"""Does incompleteness-aware prediction of resolvability beat the naive indicator?

This is the kill-criterion test for regioneer's proposed primary contribution, and it is
runnable without new curation.

The problem. The ambiguity ceiling is computable when observed haplotype diversity is
near-complete. The claimed contribution is a CALIBRATED ceiling when it is not -- a target
with 3 sequenced genomes against a neighbour with 500. An earlier argument held that the
contribution was already demonstrated, because a properly weighted ceiling explains R^2 = 0.96
of block accuracy where a raw collision indicator explains 0.23. That compares two quantities
both measured on sequences in hand; it is a statement about weighting, not extrapolation.

Design. Treat each block's full reference as ground truth. Subsample it to d genomes per
species, recompute the ceiling from the reduced reference (the NAIVE estimate, which assumes
unobserved haplotypes do not collide), and compare against the full-reference ceiling. Then
ask whether sample coverage -- Good-Turing, C = 1 - f1/n on haplotype frequencies -- PREDICTS
the resulting bias. If it does, a coverage-aware correction is learnable and the contribution
is real. If it does not, the honest deliverable is a well-weighted calculator with a coverage
warning, which is a resource rather than a statistical contribution.

No estimator is invented and tuned on this data; the test is whether the bias is predictable.

Writes results/step2_reference_incompleteness.json.
"""
import json, os, sys
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
CACHE = os.path.join(RES, "step1c_slices.json")
OUT = os.path.join(RES, "step2_reference_incompleteness.json")

SPANS = ["V4", "A_V4_515F_806R", "A_V3V4_341F_806R", "A_V1V3_27F_534R"]
DEPTHS = [2, 3, 5, 10, 20]
N_REP = 20
RESOLVABLE = 0.99


def ceiling_bal(seqs, y):
    y = np.asarray(y)
    n_t, n_n = int((y == 1).sum()), int((y == 0).sum())
    if not n_t or not n_n:
        return float("nan")
    c = defaultdict(lambda: [0, 0])
    for s, lab in zip(seqs, y):
        c[s][0 if lab == 1 else 1] += 1
    return sum(max(a / (2 * n_t), b / (2 * n_n)) for a, b in c.values())


def coverage(seqs):
    """Good-Turing sample coverage: 1 - (singleton haplotypes / n)."""
    n = len(seqs)
    if n == 0:
        return float("nan")
    f1 = sum(1 for v in Counter(seqs).values() if v == 1)
    return 1.0 - f1 / n


def main():
    data = json.load(open(CACHE))
    # isolate_group -> species, from the raw records (no re-slicing needed)
    raw = os.path.join(ROOT, "phase4_locked_test", "raw_species_test")
    sp_of = {}
    for f in os.listdir(raw):
        if f.endswith(".json"):
            for r in json.load(open(os.path.join(raw, f))):
                sp_of[r["isolate_group"]] = r.get("species_query", "?")

    rows = []
    for b, per in data.items():
        for span in SPANS:
            v = per.get(span)
            if not v:
                continue
            iso, seqs, y = v["iso"], v["seq"], np.array(v["y"])
            if len(set(y)) < 2:
                continue
            full = ceiling_bal(seqs, y)
            by_sp = defaultdict(list)
            for i, g in enumerate(iso):
                by_sp[sp_of.get(g, "?")].append(i)
            for d in DEPTHS:
                for rep in range(N_REP):
                    rng = np.random.default_rng(hash((b, span, d, rep)) % (2**32))
                    keep = []
                    for s_, idxs in by_sp.items():
                        keep.extend(idxs if len(idxs) <= d
                                    else list(rng.choice(idxs, d, replace=False)))
                    ys = y[keep]
                    if len(set(ys)) < 2:
                        continue
                    ss = [seqs[i] for i in keep]
                    naive = ceiling_bal(ss, ys)
                    rows.append(dict(block=b, span=span, depth=d, rep=rep,
                                     n_sub=len(keep), n_full=len(iso),
                                     ceiling_full=full, ceiling_naive=naive,
                                     bias=naive - full,
                                     coverage=coverage(ss),
                                     cov_target=coverage([ss[i] for i in range(len(ss)) if ys[i] == 1]),
                                     resolvable_full=int(full >= RESOLVABLE),
                                     resolvable_naive=int(naive >= RESOLVABLE)))
    print(f"{len(rows)} subsample evaluations\n")

    bias = np.array([r["bias"] for r in rows])
    inc = np.array([1 - r["coverage"] for r in rows])          # incompleteness
    ok = ~(np.isnan(bias) | np.isnan(inc))
    bias, inc = bias[ok], inc[ok]

    print("1. IS THE NAIVE ESTIMATE BIASED?  (ceiling_naive - ceiling_full)")
    print(f"{'depth':>7}{'n':>7}{'mean bias pp':>14}{'SD':>8}{'% optimistic':>14}")
    for d in DEPTHS:
        bb = np.array([r["bias"] for r in rows if r["depth"] == d and r["bias"] == r["bias"]])
        print(f"{d:>7}{len(bb):>7}{100*bb.mean():>14.2f}{100*bb.std():>8.2f}"
              f"{100*(bb>0).mean():>14.1f}")

    print("\n2. DOES COVERAGE PREDICT THE BIAS?  (the kill criterion)")
    sl, ic, r, p, se = stats.linregress(inc, bias)
    print(f"   bias ~ (1 - coverage):  R^2 = {r**2:.3f}   p = {p:.2e}   slope = {sl:.3f}")
    print("   per depth:")
    for d in DEPTHS:
        sub = [(1 - x["coverage"], x["bias"]) for x in rows
               if x["depth"] == d and x["bias"] == x["bias"] and x["coverage"] == x["coverage"]]
        if len(sub) > 5 and len({a for a, _ in sub}) > 2:
            a = np.array([x[0] for x in sub]); bq = np.array([x[1] for x in sub])
            rr_ = stats.linregress(a, bq)
            print(f"     d={d:<3} R^2 = {rr_.rvalue**2:.3f}  p = {rr_.pvalue:.1e}")

    print("\n3. OPERATIONAL: does the naive rule misclassify resolvability?")
    print(f"{'depth':>7}{'naive says resolvable':>24}{'truly resolvable':>19}{'FALSE POSITIVE rate':>21}")
    for d in DEPTHS:
        sub = [x for x in rows if x["depth"] == d]
        fp = [x for x in sub if x["resolvable_naive"] == 1 and x["resolvable_full"] == 0]
        pos = [x for x in sub if x["resolvable_naive"] == 1]
        print(f"{d:>7}{np.mean([x['resolvable_naive'] for x in sub]):>24.3f}"
              f"{np.mean([x['resolvable_full'] for x in sub]):>19.3f}"
              f"{(len(fp)/max(len(pos),1)):>21.3f}")

    summary = dict(
        n_rows=len(rows),
        bias_by_depth={str(d): float(np.nanmean([r["bias"] for r in rows if r["depth"] == d]))
                       for d in DEPTHS},
        coverage_predicts_bias=dict(r2=float(r**2), p=float(p), slope=float(sl)),
        false_positive_by_depth={
            str(d): float(len([x for x in rows if x["depth"] == d and x["resolvable_naive"] == 1
                               and x["resolvable_full"] == 0])
                          / max(len([x for x in rows if x["depth"] == d
                                     and x["resolvable_naive"] == 1]), 1))
            for d in DEPTHS})
    json.dump({"summary": summary, "rows": rows}, open(OUT, "w"), indent=2, default=float)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
