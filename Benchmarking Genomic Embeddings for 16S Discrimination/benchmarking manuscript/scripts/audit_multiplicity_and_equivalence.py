"""Audit of the confirmatory multiplicity structure and of the equivalence family.

Addresses four pre-submission review items without re-running any model:

  1. MULTIPLICITY. Re-derives the two families actually declared in the frozen plan
     (`S2 Text` / FREEZE.md section 7) and shows the adjusted results for each.
  2. SECONDARY EQUIVALENCE. Reports the two component one-sided tests per secondary
     contrast, so that no interval is presented as a formal equivalence test without
     its components.
  3. BOOTSTRAP DIFFERENCE TEST. Reproduces the two-sided p-values exactly as
     `analyze_locked_test.boot_ci` computes them, and records the attainable
     resolution at B = 2000.
  4. BEST-OF-FOUR. Separates the per-block oracle maximum from the operational
     pre-named rule, and reports the selection gap.

Reads only frozen artefacts. Writes results/multiplicity_equivalence_audit.json.
"""
import json
import os

import numpy as np
from scipy import stats

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
SCORES = os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results",
                   "multiplicity_equivalence_audit.json")

SEED_PRIMARY = 1729      # frozen; analyze_locked_test.SEED
SEED_SECONDARY = 1736    # frozen; SEED + 7
N_BOOT = 2000            # frozen
N_BOOT_TOST = 20000      # post hoc corroboration only
MARGIN = 0.010           # frozen equivalence margin, balanced-accuracy units
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]


def boot_ci(d, seed, n_boot=N_BOOT):
    """Verbatim reimplementation of analyze_locked_test.boot_ci.

    The p-value is the two-sided percentile-bootstrap achieved significance level:
    the resampling distribution is NOT shifted to the null, so this p inverts the
    percentile interval rather than testing a shifted-null hypothesis.
    """
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    draws = d[idx].mean(axis=1)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    p = 2 * min((draws <= 0).mean(), (draws >= 0).mean())
    lo90, hi90 = np.percentile(draws, [5, 95])
    return dict(mean=float(d.mean()), ci95=[float(lo), float(hi)],
                ci90=[float(lo90), float(hi90)], p_two_sided=float(min(p, 1.0)))


def tost_paired_t(d, margin=MARGIN):
    """Two one-sided paired-t tests against +-margin, with both components."""
    m, sd, n = float(d.mean()), float(d.std(ddof=1)), len(d)
    se = sd / np.sqrt(n)
    t_lower = (m + margin) / se           # H0: delta <= -margin
    t_upper = (m - margin) / se           # H0: delta >= +margin
    p_lower = float(stats.t.sf(t_lower, n - 1))
    p_upper = float(stats.t.cdf(t_upper, n - 1))
    half = float(stats.t.ppf(0.95, n - 1) * se)
    return dict(n=n, mean=m, sd=sd, se=se,
                t_lower=float(t_lower), p_lower=p_lower,
                t_upper=float(t_upper), p_upper=p_upper,
                p_tost=max(p_lower, p_upper),
                ci90=[m - half, m + half])


def tost_bootstrap(d, seed=SEED_PRIMARY, margin=MARGIN, n_boot=N_BOOT_TOST):
    """Distribution-free corroboration. Tail counts have a floor of 1/n_boot."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    draws = d[idx].mean(axis=1)
    p_lower = float((draws <= -margin).mean())
    p_upper = float((draws >= margin).mean())
    lo, hi = np.percentile(draws, [5, 95])
    return dict(n_boot=n_boot, p_lower=p_lower, p_upper=p_upper,
                p_tost=max(p_lower, p_upper), resolution_floor=1.0 / n_boot,
                ci90=[float(lo), float(hi)])


def holm(pvals, labels):
    order = np.argsort(pvals)
    m = len(pvals)
    out, running = {}, 0.0
    for rank, i in enumerate(order):
        adj = min(1.0, (m - rank) * pvals[i])
        running = max(running, adj)
        out[labels[i]] = float(running)
    return out


def verdict_ci(ci, margin=MARGIN):
    lo, hi = ci
    if lo > margin:
        return "SUPERIOR (model)"
    if hi < -margin:
        return "SUPERIOR (classical)"
    if lo > -margin and hi < margin:
        return "EQUIVALENT (whole CI inside margin)"
    return "INCONCLUSIVE at this margin"


def main():
    scores = json.load(open(SCORES))
    blocks = {k: v for k, v in scores.items()
              if v.get("_status") == "ok" and v.get("_block_eligible")}
    names = sorted(blocks)
    n = len(names)

    primary_lab = "dnabert2:R2 - C_kmer"
    deltas = {primary_lab: np.array([blocks[b]["dnabert2:R2"] - blocks[b]["C_kmer"]
                                     for b in names])}
    for m in ["dnabert_s", "hyenadna", "ntv2"]:
        deltas[f"{m}:R2 - C_kmer"] = np.array(
            [blocks[b][f"{m}:R2"] - blocks[b]["C_kmer"] for b in names])
    deltas["best-of-4 R2 - C_kmer"] = np.array(
        [max(blocks[b][f"{mm}:R2"] for mm in MODELS) - blocks[b]["C_kmer"]
         for b in names])
    secondary = [k for k in deltas if k != primary_lab]

    # ---- family 1: difference tests
    diff = {}
    for lab, d in deltas.items():
        seed = SEED_PRIMARY if lab == primary_lab else SEED_SECONDARY
        B = boot_ci(d, seed)
        B["seed"] = seed
        B["p_paired_t"] = float(stats.ttest_1samp(d, 0.0).pvalue)
        diff[lab] = B
    holm_boot = holm(np.array([diff[k]["p_two_sided"] for k in secondary]), secondary)
    holm_t = holm(np.array([diff[k]["p_paired_t"] for k in secondary]), secondary)
    for k in secondary:
        diff[k]["p_holm_bootstrap"] = holm_boot[k]
        diff[k]["p_holm_paired_t"] = holm_t[k]
    diff[primary_lab]["correction"] = "none (single prespecified primary)"

    # ---- family 2: equivalence tests
    eq = {}
    for lab, d in deltas.items():
        T = tost_paired_t(d)
        T["bootstrap_corroboration"] = tost_bootstrap(d)
        T["verdict_frozen_rule_95pct_CI_containment"] = verdict_ci(diff[lab]["ci95"])
        T["verdict_tost_alpha_0.05"] = ("EQUIVALENT" if T["p_tost"] < 0.05
                                        else "NOT DECLARED")
        eq[lab] = T
    holm_eq = holm(np.array([eq[k]["p_tost"] for k in secondary]), secondary)
    for k in secondary:
        eq[k]["p_tost_holm"] = holm_eq[k]
        eq[k]["verdict_tost_holm"] = ("EQUIVALENT" if holm_eq[k] < 0.05
                                      else "NOT DECLARED")
    eq[primary_lab]["correction"] = "none (single prespecified primary)"

    # ---- best-of-four accounting
    per_block_argmax = {}
    for b in names:
        vals = {m: blocks[b][f"{m}:R2"] for m in MODELS}
        per_block_argmax[b] = max(vals, key=vals.get)
    counts = {m: sum(1 for v in per_block_argmax.values() if v == m) for m in MODELS}
    mean_individual = float(np.mean([deltas[f"{m}:R2 - C_kmer"].mean()
                                     for m in MODELS]))
    best_of_four = {
        "selection_input": "per-block held-out cross-validated balanced accuracy "
                           "(the same quantity the contrast is computed from)",
        "selection_uses_heldout_outcomes": True,
        "character": "oracle (post hoc per-block maximum), not an operational rule",
        "per_block_argmax": per_block_argmax,
        "argmax_counts": counts,
        "mean_of_four_individual_contrasts": mean_individual,
        "best_of_four_contrast": float(deltas["best-of-4 R2 - C_kmer"].mean()),
        "selection_gap": float(deltas["best-of-4 R2 - C_kmer"].mean() - mean_individual),
        "operational_prenamed_rule": {
            "rule": "always DNABERT-2:R2, named from development before unsealing",
            "contrast": float(deltas[primary_lab].mean()),
        },
    }

    out = {
        "source": "phase4_locked_test/locked_test_scores.json",
        "n_eligible_blocks": n,
        "blocks": names,
        "units": "balanced-accuracy units (multiply by 100 for percentage points)",
        "frozen_margin": MARGIN,
        "frozen_plan_reference": "S2 Text / FREEZE.md section 7",
        "family_structure": {
            "difference_family": {
                "members": [primary_lab] + secondary,
                "structure": "hierarchical: primary tested unadjusted; Holm applied "
                             "across the four secondary contrasts only",
            },
            "equivalence_family": {
                "members": [primary_lab] + secondary,
                "prespecified_rule": "whole 95% bootstrap CI inside +-0.010; no "
                                     "multiplicity adjustment was prespecified",
                "added_here": "component one-sided TOST p-values, plus Holm across "
                              "the four secondary equivalence tests, reported as a "
                              "post hoc sensitivity",
            },
        },
        "bootstrap_difference_test": {
            "estimator": "mean paired per-block difference",
            "resampling_unit": "block",
            "n_boot": N_BOOT,
            "seed_primary": SEED_PRIMARY,
            "seed_secondary": SEED_SECONDARY,
            "interval": "percentile, 2.5th and 97.5th",
            "p_definition": "2 * min(Pr[mean* <= 0], Pr[mean* >= 0]) over the "
                            "UNSHIFTED resampling distribution; an inversion of the "
                            "percentile interval, not a shifted-null bootstrap test",
            "smallest_nonzero_two_sided_p": 2.0 / N_BOOT,
            "note": "a zero tail count returns exactly 0.0 and must be reported as "
                    "< %.4f, never as a measured value" % (2.0 / N_BOOT),
        },
        "difference_family_results": diff,
        "equivalence_family_results": eq,
        "best_of_four": best_of_four,
    }
    json.dump(out, open(OUT, "w"), indent=2, default=float)

    pp = 100.0
    print(f"n eligible blocks = {n}\n")
    print("DIFFERENCE FAMILY (percentile-bootstrap, as implemented)")
    print(f"{'contrast':26}{'mean pp':>9}{'95% CI pp':>20}{'p':>8}{'p Holm':>9}")
    for lab in [primary_lab] + secondary:
        r = diff[lab]
        ci = f"[{pp*r['ci95'][0]:+.3f},{pp*r['ci95'][1]:+.3f}]"
        h = r.get("p_holm_bootstrap")
        print(f"{lab:26}{pp*r['mean']:+9.3f}{ci:>20}{r['p_two_sided']:8.4f}"
              f"{('-' if h is None else f'{h:.4f}'):>9}")
    print("\nEQUIVALENCE FAMILY (paired-t TOST components, margin +-1.00 pp)")
    print(f"{'contrast':26}{'t_lower':>9}{'p_lower':>11}{'t_upper':>9}{'p_upper':>11}"
          f"{'p TOST':>11}{'p Holm':>11}")
    for lab in [primary_lab] + secondary:
        r = eq[lab]
        h = r.get("p_tost_holm")
        print(f"{lab:26}{r['t_lower']:9.3f}{r['p_lower']:11.3e}{r['t_upper']:9.3f}"
              f"{r['p_upper']:11.3e}{r['p_tost']:11.3e}"
              f"{('-' if h is None else f'{h:.3e}'):>11}")
    print("\nBEST-OF-FOUR")
    print(f"  per-block argmax counts        : {counts}")
    print(f"  mean of four individual        : {pp*mean_individual:+.3f} pp")
    print(f"  per-block oracle maximum       : "
          f"{pp*best_of_four['best_of_four_contrast']:+.3f} pp")
    print(f"  selection gap                  : "
          f"{pp*best_of_four['selection_gap']:+.3f} pp")
    print(f"  operational pre-named rule     : "
          f"{pp*best_of_four['operational_prenamed_rule']['contrast']:+.3f} pp")
    print(f"\nwrote {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()
