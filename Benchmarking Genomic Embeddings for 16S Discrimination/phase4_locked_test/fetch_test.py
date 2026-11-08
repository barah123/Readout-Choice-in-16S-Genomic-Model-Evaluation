"""Phase 4 step 1: fetch the LOCKED TEST blocks.

Reuses fetch_blocks.fetch_species unchanged (same code path as development) with a separate
output directory, so the test data is produced by exactly the validated pipeline.
"""
import sys, os
sys.path.insert(0, os.path.abspath("../phase2_pilot"))
sys.path.insert(0, os.path.abspath("../phase3_freeze"))

import fetch_blocks
from test_blocks_definition import TEST_BLOCKS, all_species

fetch_blocks.OUTDIR = "raw_species_test"


def main():
    os.makedirs(fetch_blocks.OUTDIR, exist_ok=True)
    targets = all_species()
    print(f"{len(targets)} species queries across {len(TEST_BLOCKS)} locked-test blocks",
          flush=True)
    tally = {}
    for n, (block_id, species, is_target) in enumerate(targets, 1):
        try:
            recs, cached = fetch_blocks.fetch_species(block_id, species, is_target)
        except Exception as e:
            print(f"[{n}/{len(targets)}] FAILED {species}: {e}", flush=True)
            continue
        n_iso = len({r["isolate_group"] for r in recs})
        tally[(block_id, species)] = (len(recs), n_iso)
        print(f"[{n}/{len(targets)}] {'cached' if cached else 'fetched':7} "
              f"{block_id:28} {species:44} records={len(recs):4} isolates={n_iso:4}",
              flush=True)

    print("\n--- per-block target support (pre-slicing) ---", flush=True)
    for b in TEST_BLOCKS:
        t = tally.get((b["block_id"], b["target"]))
        ns = [tally.get((b["block_id"], x)) for x in b["neighbours"]]
        ns = [x for x in ns if x]
        print(f"{b['block_id']:28} target_isolates={(t[1] if t else 0):4}  "
              f"neighbours_ge5={sum(1 for x in ns if x[1] >= 5)}/{len(b['neighbours'])}",
              flush=True)


if __name__ == "__main__":
    main()
