"""Block-level analysis of the full benchmark (protocol §13, §16).

Primary estimand (§13): mean PAIRED difference between a prespecified frozen-model strategy
and a prespecified strong classical strategy, across INDEPENDENT BLOCKS.

Prespecified from the 2-panel pilot (phase2_pilot_report.md "Implications for Phase 3"):
  model strategy     = best R2 linear probe across the 4 core models, selected per block
                       (a development-selected rule, declared before this full run)
  classical strategy = best of C_exact / C_align / C_kmer, selected per block
  C_random is reported alongside as a flexibility reference, NOT as the comparator.

Inference unit is the block, not the cell (§16: "Do not t-test correlated folds, count
pairwise distances as independent, or use the number of model-region cells as biological n").
Within a block, the region is chosen by the SAME rule for both arms, then one value per
block enters the bootstrap. Bootstrap resamples blocks (2000 draws, §16).

Strata kept separate: species-rank vs subspecies-rank (§3); confirmatory-eligible vs
descriptive-only (the pre-declared gate); Enterobacteriaceae collapsed in a stringent
sensitivity analysis (§6.3).
"""
import json
from collections import defaultdict

import numpy as np

MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
CLASSICAL = ["C_exact", "C_align", "C_kmer"]
SEED = 1729
N_BOOT = 2000


def ok(v):
    return v is not None and isinstance(v, float) and v == v


def main():
    res = json.load(open("readout_results_full.json"))
    summ = json.load(open("block_summary_full.json"))

    cells = [r for r in res.values() if r.get("_status") == "ok"]
    print(f"{len(cells)} scored cells across {len({c['_block'] for c in cells})} blocks\n")

    # ---------------------------------------------------------------- per-block reduction
    per_block = {}
    for c in cells:
        b = c["_block"]
        if not c.get("_eligible"):
            continue
        cl = {k: c.get(k) for k in CLASSICAL if ok(c.get(k))}
        md = {f"{m}:R2_linear": c.get(f"{m}:R2_linear")
              for m in MODELS if ok(c.get(f"{m}:R2_linear"))}
        if not cl or not md:
            continue
        # region selection: the region where the CLASSICAL arm is strongest, applied
        # identically to both arms so the choice cannot favour either.
        score = max(cl.values())
        prev = per_block.get(b)
        if prev is None or score > prev["_sel_score"]:
            per_block[b] = {
                "_sel_score": score, "region": c["_region"],
                "classical": max(cl.values()), "classical_name": max(cl, key=cl.get),
                "model": max(md.values()), "model_name": max(md, key=md.get),
                "C_random": c.get("C_random"),
                "ceiling": c.get("_ceiling_bal"),
                "rank": c.get("_rank"), "family": c.get("_family"),
                "expected_class": c.get("_expected_class"),
                "target_isolates": c.get("_target_isolates"),
                **{f"{m}:R0": c.get(f"{m}:R0_geometry_auc") for m in MODELS},
                **{f"{m}:R1": c.get(f"{m}:R1_nn") for m in MODELS},
                **{f"{m}:R2": c.get(f"{m}:R2_linear") for m in MODELS},
            }

    species_blocks = {k: v for k, v in per_block.items() if v["rank"] == "species"}
    subsp_blocks = {k: v for k, v in per_block.items() if v["rank"] == "subspecies"}

    print("=" * 112)
    print("CONFIRMATORY-ELIGIBLE BLOCKS (species rank) - best classical vs best R2 probe")
    print("=" * 112)
    print(f"{'block':30}{'reg':5}{'ceil':>7}{'classical':>22}{'best R2 model':>26}"
          f"{'delta':>8}{'rand':>7}")
    for b, v in sorted(species_blocks.items(),
                       key=lambda kv: -(kv[1]["model"] - kv[1]["classical"])):
        d = v["model"] - v["classical"]
        cl_label = "{} {:.3f}".format(v["classical_name"], v["classical"])
        md_label = "{} {:.3f}".format(v["model_name"].replace(":R2_linear", ""), v["model"])
        rnd_val = v["C_random"] if ok(v["C_random"]) else float("nan")
        print(f"{b:30}{v['region']:5}{v['ceiling']:7.3f}{cl_label:>22}{md_label:>26}"
              f"{d:>+8.3f}{rnd_val:>7.3f}")

    # ---------------------------------------------------------------- primary contrast
    deltas = np.array([v["model"] - v["classical"] for v in species_blocks.values()])
    if len(deltas) >= 2:
        rng = np.random.default_rng(SEED)
        boot = np.array([rng.choice(deltas, len(deltas), replace=True).mean()
                         for _ in range(N_BOOT)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        print(f"\nPRIMARY CONTRAST (species-rank, eligible blocks, n={len(deltas)} blocks)")
        print(f"  mean delta (best-R2 minus best-classical) = {deltas.mean():+.4f}")
        print(f"  95% bootstrap CI over blocks              = [{lo:+.4f}, {hi:+.4f}]")
        print(f"  blocks favouring model / classical / tie(+-0.01): "
              f"{(deltas > 0.01).sum()} / {(deltas < -0.01).sum()} / "
              f"{((deltas >= -0.01) & (deltas <= 0.01)).sum()}")

    # flexibility reference
    rnd = np.array([v["model"] - v["C_random"] for v in species_blocks.values()
                    if ok(v["C_random"])])
    if len(rnd) >= 2:
        rng = np.random.default_rng(SEED + 1)
        b2 = np.array([rng.choice(rnd, len(rnd), replace=True).mean() for _ in range(N_BOOT)])
        print(f"\n  reference: best-R2 minus C_random         = {rnd.mean():+.4f} "
              f"[{np.percentile(b2,2.5):+.4f}, {np.percentile(b2,97.5):+.4f}]")

    # ---------------------------------------------------------------- stringent family grouping
    fam = defaultdict(list)
    for v in species_blocks.values():
        fam[v["family"]].append(v["model"] - v["classical"])
    fd = np.array([np.mean(v) for v in fam.values()])
    if len(fd) >= 2:
        rng = np.random.default_rng(SEED + 2)
        b3 = np.array([rng.choice(fd, len(fd), replace=True).mean() for _ in range(N_BOOT)])
        print(f"\nSTRINGENT (family-grouped, n={len(fd)} families): mean delta = {fd.mean():+.4f} "
              f"[{np.percentile(b3,2.5):+.4f}, {np.percentile(b3,97.5):+.4f}]")

    # ---------------------------------------------------------------- readout hierarchy
    print("\n" + "=" * 112)
    print("READOUT HIERARCHY, averaged over eligible species-rank blocks")
    print("=" * 112)
    print(f"{'model':12}{'R0 geometry':>14}{'R1 nearest':>14}{'R2 probe':>12}")
    for m in MODELS:
        r0 = [v[f"{m}:R0"] for v in species_blocks.values() if ok(v.get(f"{m}:R0"))]
        r1 = [v[f"{m}:R1"] for v in species_blocks.values() if ok(v.get(f"{m}:R1"))]
        r2 = [v[f"{m}:R2"] for v in species_blocks.values() if ok(v.get(f"{m}:R2"))]
        print(f"{m:12}{np.mean(r0):>14.3f}{np.mean(r1):>14.3f}{np.mean(r2):>12.3f}")
    cl = [v["classical"] for v in species_blocks.values()]
    print(f"{'(classical)':12}{'-':>14}{'-':>14}{np.mean(cl):>12.3f}")

    # ---------------------------------------------------------------- expectation calibration
    print("\n" + "=" * 112)
    print("A-PRIORI EXPECTED CLASS vs MEASURED CEILING (calibration of our own predictions)")
    print("=" * 112)
    byc = defaultdict(list)
    for v in per_block.values():
        byc[v["expected_class"]].append(v["ceiling"])
    for c in "ABC":
        if byc[c]:
            print(f"  expected {c}: n={len(byc[c]):2}  measured ceiling "
                  f"mean={np.mean(byc[c]):.3f}  range=[{min(byc[c]):.3f}, {max(byc[c]):.3f}]")

    # ---------------------------------------------------------------- coverage accounting
    print("\n" + "=" * 112)
    print("BLOCK ACCOUNTING (protocol §19: distinct codes, never all rendered as zero)")
    print("=" * 112)
    tot = len(summ)
    elig_blocks = {v["block_id"] for v in summ.values() if v.get("block_eligible")}
    print(f"  blocks declared                          : 31")
    print(f"  blocks with fetched records              : {tot}")
    print(f"  blocks with >=1 eligible region          : {len(elig_blocks)}")
    print(f"  blocks entering primary contrast (species): {len(species_blocks)}")
    print(f"  subspecies-stratum blocks (reported separately): {len(subsp_blocks)}")
    ineligible = [v["block_id"] for v in summ.values() if not v.get("block_eligible")]
    if ineligible:
        print(f"  descriptive-only (failed eligibility gate): {len(ineligible)}")
        for b in ineligible:
            r = summ[b]["regions"]
            mx = max((x["target_isolates"] for x in r.values()), default=0)
            print(f"      {b:30} max target isolates in any region = {mx}")

    for b, v in subsp_blocks.items():
        print(f"\nSUBSPECIES STRATUM {b}: ceiling={v['ceiling']:.3f} "
              f"classical={v['classical']:.3f} best_R2={v['model']:.3f} "
              f"(not pooled with species-rank results)")

    json.dump(per_block, open("block_level_results.json", "w"), indent=2, default=float)
    print("\nwrote block_level_results.json")


if __name__ == "__main__":
    main()
