"""Phase 5: the PRESPECIFIED analysis of the locked test set. Run once.

Implements FREEZE.md §3, §4, §7 exactly and adds nothing:
  PRIMARY   (1 test, uncorrected): DNABERT-2:R2 - C_kmer, paired over eligible blocks,
            bootstrap 2000 draws seed 1729, 95% percentile CI, equivalence at +-0.010
            declared only if the ENTIRE CI lies inside the margin.
  SECONDARY (4 tests, Holm-corrected): DNABERT-S:R2, HyenaDNA:R2, NTv2:R2, and
            best-of-4-R2, each vs C_kmer.
  SENSITIVITY: the primary CI also read against the protocol's illustrative +-0.03
            (deviation D-001).
  EXPLORATORY: everything else, labelled as such and never promoted.
"""
import json
from collections import Counter, defaultdict

import numpy as np

SEED = 1729          # frozen
N_BOOT = 2000        # frozen
MARGIN = 0.010       # frozen
MARGIN_SENS = 0.030  # protocol's illustrative margin, reported as sensitivity
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]

SCORES = "../phase4_locked_test/locked_test_scores.json"
SUMMARY = "../phase4_locked_test/test_block_summary.json"


def boot_ci(d, seed, n_boot=N_BOOT):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    draws = d[idx].mean(axis=1)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    p = 2 * min((draws <= 0).mean(), (draws >= 0).mean())
    return float(lo), float(hi), float(min(p, 1.0)), draws


def verdict(lo, hi, margin):
    if lo > margin:
        return "SUPERIOR (model)"
    if hi < -margin:
        return "SUPERIOR (classical)"
    if lo > -margin and hi < margin:
        return "EQUIVALENT (practically negligible difference)"
    return "INCONCLUSIVE at this margin"


def holm(pvals, labels):
    order = np.argsort(pvals)
    m = len(pvals)
    out, running = {}, 0.0
    for rank, i in enumerate(order):
        adj = min(1.0, (m - rank) * pvals[i])
        running = max(running, adj)      # enforce monotonicity
        out[labels[i]] = running
    return out


def main():
    scores = json.load(open(SCORES))
    summ = json.load(open(SUMMARY))

    blocks = {k: v for k, v in scores.items()
              if v.get("_status") == "ok" and v.get("_block_eligible")}
    names = sorted(blocks)
    n = len(names)

    print("=" * 108)
    print("PHASE 5 - PRESPECIFIED ANALYSIS OF THE LOCKED TEST SET")
    print("=" * 108)
    print(f"eligible blocks entering confirmatory analysis: n = {n}")
    print("freeze precondition: >=14 for planned precision, <10 => report underpowered")
    if n < 10:
        print("*** UNDERPOWERED - primary reported as inconclusive per FREEZE.md section 5 ***")
    print()

    arms = [f"{m}:R2" for m in MODELS] + ["C_kmer", "C_exact", "C_align", "C_random"]
    missing = {a: [b for b in names if a not in blocks[b] or blocks[b][a] != blocks[b][a]]
               for a in arms}
    tot_missing = sum(len(v) for v in missing.values())
    print(f"MISSING CELLS: {tot_missing} across {len(arms)} arms x {n} blocks "
          f"({'none' if tot_missing == 0 else missing})")
    print()

    print("-" * 108)
    print("PER-BLOCK RESULTS (balanced accuracy; region chosen support-only, "
          "same region for both arms)")
    print("-" * 108)
    print(f"{'block':26}{'ceil':>7}{'C_kmer':>9}{'DNABERT-2':>11}{'delta':>9}"
          f"{'DNA-S':>8}{'Hyena':>8}{'NTv2':>8}{'C_rand':>8}")
    rows = []
    for b in names:
        v = blocks[b]
        ceils = [c for c in (v.get("_ceilings") or {}).values() if c == c]
        chosen = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        ceil = (v.get("_ceilings") or {}).get(chosen, max(ceils) if ceils else float("nan"))
        d = v["dnabert2:R2"] - v["C_kmer"]
        rows.append((b, ceil, v["C_kmer"], v["dnabert2:R2"], d))
        print(f"{b:26}{ceil:7.3f}{v['C_kmer']:9.3f}{v['dnabert2:R2']:11.3f}{d:+9.3f}"
              f"{v['dnabert_s:R2']:8.3f}{v['hyenadna:R2']:8.3f}{v['ntv2:R2']:8.3f}"
              f"{v['C_random']:8.3f}")

    prim = np.array([blocks[b]["dnabert2:R2"] - blocks[b]["C_kmer"] for b in names])
    lo, hi, p, _ = boot_ci(prim, SEED)
    print()
    print("=" * 108)
    print("PRIMARY ENDPOINT (prespecified, uncorrected): DNABERT-2:R2  minus  C_kmer")
    print("=" * 108)
    print(f"  n blocks                       : {n}")
    print(f"  mean paired difference         : {prim.mean():+.5f}")
    print(f"  95% bootstrap CI (2000, s=1729): [{lo:+.5f}, {hi:+.5f}]   "
          f"half-width {(hi-lo)/2:.5f}")
    print(f"  two-sided bootstrap p (H0: 0)  : {p:.4f}")
    print(f"  blocks model>classical / < / tie(+-0.01): "
          f"{(prim>0.01).sum()} / {(prim<-0.01).sum()} / "
          f"{((prim>=-0.01)&(prim<=0.01)).sum()}")
    print()
    print(f"  VERDICT at frozen margin +-{MARGIN:.3f} : {verdict(lo, hi, MARGIN)}")
    print(f"  sensitivity at protocol's +-{MARGIN_SENS:.3f}  : {verdict(lo, hi, MARGIN_SENS)}")

    print()
    print("=" * 108)
    print("SECONDARY ENDPOINTS (4 tests, Holm-corrected within this family)")
    print("=" * 108)
    sec_labels, sec_p, sec_stats = [], [], {}
    for m in ["dnabert_s", "hyenadna", "ntv2"]:
        d = np.array([blocks[b][f"{m}:R2"] - blocks[b]["C_kmer"] for b in names])
        l, h, pp, _ = boot_ci(d, SEED + 7)
        sec_labels.append(f"{m}:R2 - C_kmer"); sec_p.append(pp)
        sec_stats[f"{m}:R2 - C_kmer"] = (d.mean(), l, h, pp)
    best = np.array([max(blocks[b][f"{m}:R2"] for m in MODELS) - blocks[b]["C_kmer"]
                     for b in names])
    l, h, pp, _ = boot_ci(best, SEED + 7)
    sec_labels.append("best-of-4 R2 - C_kmer"); sec_p.append(pp)
    sec_stats["best-of-4 R2 - C_kmer"] = (best.mean(), l, h, pp)

    adj = holm(np.array(sec_p), sec_labels)
    print(f"{'contrast':26}{'mean':>10}{'95% CI':>24}{'p raw':>9}{'p Holm':>9}"
          f"  verdict @ +-0.010")
    for lab in sec_labels:
        mu, l, h, pp = sec_stats[lab]
        ci = f"[{l:+.4f},{h:+.4f}]"
        print(f"{lab:26}{mu:+10.4f}{ci:>24}{pp:9.4f}{adj[lab]:9.4f}  "
              f"{verdict(l, h, MARGIN)}")

    print()
    print("=" * 108)
    print("EXPLORATORY (not confirmatory; never promoted to a confirmatory claim)")
    print("=" * 108)
    for a in arms:
        print(f"  {a:16} mean balanced accuracy over {n} blocks = "
              f"{np.mean([blocks[b][a] for b in names]):.4f}")
    rnd = np.array([max(blocks[b][f"{m}:R2"] for m in MODELS) - blocks[b]["C_random"]
                    for b in names])
    l, h, _, _ = boot_ci(rnd, SEED + 11)
    print(f"\n  flexibility reference: best-of-4 R2 minus C_random = {rnd.mean():+.4f} "
          f"[{l:+.4f}, {h:+.4f}]")

    byc = defaultdict(list)
    for b in names:
        v = blocks[b]
        ceils = [c for c in (v.get("_ceilings") or {}).values() if c == c]
        if ceils:
            byc[v.get("_expected_class")].append(max(ceils))
    print("\n  a-priori expected class vs measured ceiling (calibration):")
    for c in "ABC":
        if byc.get(c):
            print(f"    expected {c}: n={len(byc[c]):2}  mean max-ceiling="
                  f"{np.mean(byc[c]):.3f}  range=[{min(byc[c]):.3f}, {max(byc[c]):.3f}]")

    allregions = Counter(r for b in names for r in blocks[b]["_chosen_regions"])
    print(f"\n  regions chosen by support-only selection across {5*n} folds: "
          f"{dict(allregions.most_common())}")

    print()
    print("=" * 108)
    print("BLOCK ACCOUNTING")
    print("=" * 108)
    inelig = [b for b, v in summ.items() if not v.get("block_eligible")]
    print(f"  declared locked-test blocks        : 22")
    print(f"  eligible and scored (confirmatory) : {n}")
    print(f"  descriptive-only (failed gate)     : {len(inelig)}  {inelig}")

    json.dump({
        "n_blocks": n,
        "primary": {"contrast": "dnabert2:R2 - C_kmer", "mean": float(prim.mean()),
                    "ci95": [lo, hi], "p_two_sided": p, "margin": MARGIN,
                    "verdict": verdict(lo, hi, MARGIN),
                    "verdict_sensitivity_0.03": verdict(lo, hi, MARGIN_SENS)},
        "secondary": {k: {"mean": v[0], "ci95": [v[1], v[2]], "p_raw": v[3],
                          "p_holm": adj[k], "verdict": verdict(v[1], v[2], MARGIN)}
                      for k, v in sec_stats.items()},
        "per_block": {b: {"ceiling": r[1], "C_kmer": r[2], "dnabert2_R2": r[3],
                          "delta": r[4]} for b, r in zip(names, rows)},
        "missing_cells": tot_missing,
        "ineligible_blocks": inelig,
    }, open("locked_test_results.json", "w"), indent=2, default=float)
    print("\nwrote locked_test_results.json")


if __name__ == "__main__":
    main()
