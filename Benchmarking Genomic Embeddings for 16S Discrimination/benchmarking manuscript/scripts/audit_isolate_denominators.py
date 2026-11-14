"""Block-by-block, fold-by-fold reconciliation of the locked-test isolate denominators.

Three counts circulate for the locked test and are NOT interchangeable:

  6,760  isolate groups across all 22 DECLARED blocks (manifest, S7/S8)
  5,977  isolate groups present at the MODAL chosen region, 20 eligible blocks
  5,902  held-out isolates summed over folds, each fold at ITS OWN chosen region

A reviewer argued the last two must be equal, since 5-fold CV holds every group out exactly
once. That is true of the SPLIT population but not of the EVALUATED population. The outer
split is defined over the union of isolates that have a slice in at least one eligible
region; a fold then evaluates only those held-out isolates that have a slice at the region
that fold selected. Region availability varies (V4 100%, V9 48.7% across the corpus), and
folds may select different regions within one block, so the evaluated count is at or below
the split count. This script prints the full chain and writes it as a table.

Reads only frozen artefacts. Writes results/isolate_denominator_audit.csv.
"""
import csv
import json
import os
from collections import Counter, defaultdict

import numpy as np
from sklearn.model_selection import StratifiedKFold

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "isolate_denominator_audit.csv"))


def main():
    rows = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_blocks.json")))
    scores = json.load(open(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json")))
    summ = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_block_summary.json")))

    per = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        per[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    # every isolate group named anywhere in the locked-test manifest, all 22 declared blocks
    declared_groups = {b: {i for d in per[b].values() for i in d} for b in per}
    n_declared = sum(len(v) for v in declared_groups.values())
    inelig = sorted(b for b, v in summ.items() if not v.get("block_eligible"))

    out, tot = [], Counter()
    for b in sorted(scores):
        v = scores[b]
        if v.get("_status") != "ok" or not v.get("_block_eligible"):
            continue
        pr = per[b]
        labels = {i: lab for region, d in pr.items() for i, (_, lab) in d.items()}
        isolates = sorted(labels)                 # the SPLIT population for this block
        y = np.array([labels[i] for i in isolates])
        chosen = v["_chosen_regions"]
        modal = Counter(chosen).most_common(1)[0][0]

        n_split = len(isolates)
        n_modal = len(pr[modal])
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        n_eval, per_fold = 0, []
        for k, (_, te_idx) in enumerate(skf.split(np.zeros(n_split), y)):
            if k >= len(chosen):
                break
            d = pr[chosen[k]]
            held = len(te_idx)
            present = sum(1 for i in te_idx if isolates[i] in d)
            per_fold.append((chosen[k], held, present))
            n_eval += present

        out.append(dict(
            block=b, split_population=n_split, modal_region=modal,
            groups_at_modal_region=n_modal,
            distinct_regions_chosen=len(set(chosen)),
            regions_by_fold="|".join(chosen),
            heldout_total=sum(f[1] for f in per_fold),
            evaluated_total=n_eval,
            dropped_absent_at_fold_region=sum(f[1] for f in per_fold) - n_eval))
        tot["split"] += n_split
        tot["modal"] += n_modal
        tot["held"] += sum(f[1] for f in per_fold)
        tot["eval"] += n_eval

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    print("LOCKED-TEST ISOLATE DENOMINATOR CHAIN")
    print("=" * 78)
    print(f"  (1) isolate groups, all 22 DECLARED blocks        : {n_declared:,}")
    print(f"      less the 2 blocks failing the eligibility gate: "
          f"-{sum(len(declared_groups[b]) for b in inelig):,}  {inelig}")
    print(f"  (2) isolate groups, 20 ELIGIBLE blocks (split pop): {tot['split']:,}")
    print(f"  (3) groups present at each block's MODAL region   : {tot['modal']:,}"
          f"   (-{tot['split']-tot['modal']:,} absent at that region)")
    print(f"  (4) held out across 5 folds (= split population)  : {tot['held']:,}")
    print(f"  (5) EVALUATED: held out AND present at that fold's")
    print(f"      own chosen region                             : {tot['eval']:,}"
          f"   (-{tot['held']-tot['eval']:,} absent at that fold's region)")
    print()
    print(f"  (2) == (4): {tot['split'] == tot['held']}  "
          "— every group is held out exactly once, as expected")
    print(f"  (3) != (5) by {tot['modal']-tot['eval']:,}: the modal-region count and the "
          "per-fold evaluated")
    print("      count are different populations whenever folds pick different regions.")
    multi = [r for r in out if r["distinct_regions_chosen"] > 1]
    print(f"\n  blocks where folds chose >1 distinct region: {len(multi)} of {len(out)}")
    for r in multi:
        print(f"    {r['block']:26} {r['regions_by_fold']:20} "
              f"modal={r['groups_at_modal_region']:4}  evaluated={r['evaluated_total']:4}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
