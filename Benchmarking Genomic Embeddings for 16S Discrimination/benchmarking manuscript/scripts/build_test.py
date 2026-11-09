"""Phase 4 step 2: slice the locked-test blocks and apply the frozen eligibility gate.

Reuses build_blocks_full's slicing and ceiling code unchanged (same code path as
development), pointed at the test block list and test raw data.
"""
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
sys.path.insert(0, os.path.abspath("../phase3_freeze"))

from joblib import Parallel, delayed

import build_blocks_full as bb
from test_blocks_definition import TEST_BLOCKS

RAW = "raw_species_test"


def load_block(b):
    recs = []
    for sp in [b["target"]] + b["neighbours"]:
        safe = sp.replace(" ", "_").replace("/", "_").replace(".", "")
        path = os.path.join(RAW, f"{b['block_id']}__{safe}.json")
        if os.path.exists(path):
            recs.extend(json.load(open(path)))
    return recs


def main():
    all_rows, summary = [], {}
    for b in TEST_BLOCKS:
        recs = load_block(b)
        if not recs:
            print(f"[skip] {b['block_id']}: no records", flush=True)
            continue
        sliced = Parallel(n_jobs=8, batch_size=16)(
            delayed(bb._slice_one)(r["sequence"]) for r in recs)
        for r, sl in zip(recs, sliced):
            r["slices"] = sl

        binfo = {"block_id": b["block_id"], "family": b["family"], "rank": b["rank"],
                 "target": b["target"], "expected_class": b["expected_class"],
                 "n_records": len(recs), "regions": {}}

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

            if n_t == 0 or n_n == 0:
                status = "no_coverage_one_class"
            elif n_t < bb.MIN_TARGET_ISOLATES:
                status = "insufficient_target_isolates"
            elif nb_ok < bb.MIN_NEIGHBOUR_SPECIES:
                status = "insufficient_neighbour_species"
            else:
                status = "eligible"

            binfo["regions"][region] = dict(
                n_isolates=len(iso_rows), target_isolates=n_t, neighbour_isolates=n_n,
                hap_target=len({sl for _, _, t, sl in iso_rows if t}),
                hap_neighbour=len({sl for _, _, t, sl in iso_rows if not t}),
                ambiguous_haplotypes=sum(1 for v in cnt.values() if v[0] and v[1]),
                ceiling_acc_isolate=acc, ceiling_bal_isolate=bal,
                neighbour_species_ok=nb_ok, status=status,
                eligible=bool(status == "eligible"))

            for iso, sp, is_t, sl in iso_rows:
                all_rows.append(dict(block_id=b["block_id"], region=region,
                                     isolate_group=iso, species=sp,
                                     is_target=int(is_t), slice=sl))

        elig = [r for r, v in binfo["regions"].items() if v["eligible"]]
        binfo["eligible_regions"] = elig
        binfo["block_eligible"] = bool(elig)
        summary[b["block_id"]] = binfo
        mx = max((v["ceiling_bal_isolate"] for v in binfo["regions"].values()
                  if v["ceiling_bal_isolate"] == v["ceiling_bal_isolate"]),
                 default=float("nan"))
        print(f"{b['block_id']:28} recs={len(recs):5} eligible_regions={len(elig)}/9 "
              f"max_ceiling={mx:.3f}", flush=True)

    json.dump(all_rows, open("test_blocks.json", "w"))
    json.dump(summary, open("test_block_summary.json", "w"), indent=2)
    nelig = sum(1 for v in summary.values() if v["block_eligible"])
    print(f"\nwrote test_blocks.json ({len(all_rows)} isolate x region rows)")
    print(f"ELIGIBLE BLOCKS: {nelig}/{len(summary)}  "
          f"(freeze requires >=14 for planned precision; <10 => underpowered)")


if __name__ == "__main__":
    main()
