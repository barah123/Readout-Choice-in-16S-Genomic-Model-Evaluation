"""Extension A: full fetch, build, embed and score for the 15 screened HARD panels.

Runs the frozen Phase 3 analysis without modification (same checkpoints, same pooling, same
R2 linear probe, isolate-grouped CV, support-only region selection applied identically to
both arms, primary estimand DNABERT-2:R2 - C_kmer, bootstrap 2000 draws seed 1729).

The only difference from Phase 4/5 is which panels enter: these were selected for difficulty
using classical methods only, before any model was run on them.
"""
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
sys.path.insert(0, os.path.abspath("../phase4_locked_test"))
import build_blocks_full as bb
import fetch_blocks

fetch_blocks.OUTDIR = "raw_species_hard"
fetch_blocks.RETMAX = 300


def selected_panels():
    out = []
    for f in ["screen_results.json", "screen_results_round2.json"]:
        for bid, v in json.load(open(f)).items():
            if v.get("hard"):
                out.append(dict(block_id=bid, family=v["family"], rank="mixed",
                                target=v["target"], neighbours=v["neighbours"],
                                screen_ceiling=v["ceiling"], screen_kmer=v["kmer"]))
    return out


def step_fetch(panels):
    os.makedirs(fetch_blocks.OUTDIR, exist_ok=True)
    for p in panels:
        for sp, is_t in ([(p["target"], True)]
                         + [(n, False) for n in p["neighbours"]]):
            try:
                r, cached = fetch_blocks.fetch_species(p["block_id"], sp, is_t)
                print(f"  {'cached' if cached else 'fetched':7} {p['block_id']:32} "
                      f"{sp:48} {len(r):4} records", flush=True)
            except Exception as e:
                print(f"  FAILED {sp}: {e}", flush=True)


def step_build(panels):
    all_rows, summary = [], {}
    for p in panels:
        recs = []
        for sp in [p["target"]] + p["neighbours"]:
            safe = sp.replace(" ", "_").replace("/", "_").replace(".", "")
            path = os.path.join(fetch_blocks.OUTDIR, f"{p['block_id']}__{safe}.json")
            if os.path.exists(path):
                recs.extend(json.load(open(path)))
        if not recs:
            continue
        sliced = Parallel(n_jobs=8, batch_size=16)(
            delayed(bb._slice_one)(r["sequence"]) for r in recs)
        for r, sl in zip(recs, sliced):
            r["slices"] = sl

        binfo = {"block_id": p["block_id"], "family": p["family"], "rank": p["rank"],
                 "target": p["target"], "expected_class": "HARD",
                 "screen_kmer": p["screen_kmer"], "n_records": len(recs), "regions": {}}
        for region in bb.REGIONS:
            present = [r for r in recs if r["slices"][region]]
            if not present:
                continue
            by_iso = defaultdict(list)
            for r in present:
                by_iso[r["isolate_group"]].append(r)
            iso_rows = []
            for iso, rs in by_iso.items():
                modal = Counter(x["slices"][region] for x in rs).most_common(1)[0][0]
                iso_rows.append((iso, rs[0]["species_query"], rs[0]["is_target"], modal))
            n_t = sum(1 for x in iso_rows if x[2])
            n_n = len(iso_rows) - n_t
            cnt = defaultdict(lambda: [0, 0])
            for _, _, is_t, sl in iso_rows:
                cnt[sl][0 if is_t else 1] += 1
            acc, bal = bb.ceilings(cnt, n_t, n_n)
            per_nb = Counter(sp for _, sp, is_t, _ in iso_rows if not is_t)
            nb_ok = sum(1 for v in per_nb.values() if v >= bb.MIN_NEIGHBOUR_ISOLATES)
            status = ("no_coverage_one_class" if (n_t == 0 or n_n == 0)
                      else "insufficient_target_isolates" if n_t < bb.MIN_TARGET_ISOLATES
                      else "insufficient_neighbour_species" if nb_ok < bb.MIN_NEIGHBOUR_SPECIES
                      else "eligible")
            binfo["regions"][region] = dict(
                n_isolates=len(iso_rows), target_isolates=n_t, neighbour_isolates=n_n,
                ceiling_acc_isolate=acc, ceiling_bal_isolate=bal,
                neighbour_species_ok=nb_ok, status=status,
                eligible=bool(status == "eligible"))
            for iso, sp, is_t, sl in iso_rows:
                all_rows.append(dict(block_id=p["block_id"], region=region,
                                     isolate_group=iso, species=sp,
                                     is_target=int(is_t), slice=sl))
        elig = [r for r, v in binfo["regions"].items() if v["eligible"]]
        binfo["eligible_regions"] = elig
        binfo["block_eligible"] = bool(elig)
        summary[p["block_id"]] = binfo
        mx = max((v["ceiling_bal_isolate"] for v in binfo["regions"].values()
                  if v["ceiling_bal_isolate"] == v["ceiling_bal_isolate"]), default=float("nan"))
        print(f"{p['block_id']:32} recs={len(recs):5} eligible={len(elig)}/9 "
              f"max_ceiling={mx:.3f}", flush=True)

    json.dump(all_rows, open("hard_blocks.json", "w"))
    json.dump(summary, open("hard_block_summary.json", "w"), indent=2)
    n = sum(1 for v in summary.values() if v["block_eligible"])
    print(f"\nwrote hard_blocks.json ({len(all_rows)} rows); eligible blocks: "
          f"{n}/{len(summary)}")


if __name__ == "__main__":
    panels = selected_panels()
    print(f"Extension A: {len(panels)} screened-hard panels\n")
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "fetch"):
        print("--- fetch ---", flush=True)
        step_fetch(panels)
    if what in ("all", "build"):
        print("\n--- build ---", flush=True)
        step_build(panels)
