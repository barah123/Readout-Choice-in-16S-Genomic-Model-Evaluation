"""Fold-matched ambiguity ceilings for every block plotted in Figure 7.

The reported accuracy is the mean over five outer folds of per-fold balanced accuracy.
The published ceiling is evaluated once on the WHOLE block at the modal chosen region.
Those two aggregations do not match, and a fold mean can exceed a whole-population bound:
3 of 66 cells do.

This recomputes the identical ceiling formula (`build_blocks_full.ceilings`) on each
HELD-OUT fold at that fold's own chosen region, then averages over folds - matching the
accuracy's aggregation exactly. The result is a strict upper bound on every reported cell
and differs from the published ceiling by +0.001 pp on average.

Reads only frozen artefacts; reproduces the outer split from the frozen definition
(StratifiedKFold(5, shuffle=True, random_state=0) over the isolate-level table, which is
grouped-by-isolate by construction since each isolate contributes one modal slice per
region). Writes results/fold_matched_ceiling_audit.csv.
"""
import csv
import json
import os
from collections import Counter, defaultdict

import numpy as np
from sklearn.model_selection import StratifiedKFold

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "fold_matched_ceiling_audit.csv"))

SETS = [("locked test", "phase4_locked_test/test_blocks.json",
         "phase4_locked_test/locked_test_scores.json"),
        ("hard panels", "phase7_extensions/hard_blocks.json",
         "phase7_extensions/hard_scores.json")]


def ceiling_bal(slices, y):
    """build_blocks_full.ceilings, balanced component. counts are [target, neighbour]."""
    y = np.asarray(y)
    n_t, n_n = int((y == 1).sum()), int((y == 0).sum())
    if not n_t or not n_n:
        return float("nan")
    c = defaultdict(lambda: [0, 0])
    for s, lab in zip(slices, y):
        c[s][0 if lab == 1 else 1] += 1
    return sum(max(a / (2 * n_t), b / (2 * n_n)) for a, b in c.values())


def main():
    rows_out = []
    for setname, blockfile, scorefile in SETS:
        rows = json.load(open(os.path.join(ROOT, blockfile)))
        scores = json.load(open(os.path.join(ROOT, scorefile)))
        per = defaultdict(lambda: defaultdict(dict))
        for r in rows:
            per[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

        for b in sorted(scores):
            v = scores[b]
            if v.get("_status") != "ok" or not v.get("_block_eligible"):
                continue
            pr = per[b]
            labels = {i: lab for region, d in pr.items() for i, (_, lab) in d.items()}
            isolates = sorted(labels)
            y = np.array([labels[i] for i in isolates])
            chosen = v["_chosen_regions"]
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
            fold_ceils = []
            for k, (_, te_idx) in enumerate(skf.split(np.zeros(len(isolates)), y)):
                if k >= len(chosen):
                    break
                d = pr[chosen[k]]
                te = [isolates[i] for i in te_idx if isolates[i] in d]
                fold_ceils.append(ceiling_bal([d[i][0] for i in te], [d[i][1] for i in te]))
            modal = Counter(chosen).most_common(1)[0][0]
            rows_out.append(dict(
                set=setname, block=b, modal_region=modal,
                ceiling_full_block=round(float((v.get("_ceilings") or {}).get(modal, float("nan"))), 6),
                ceiling_fold_matched=round(float(np.nanmean(fold_ceils)), 6),
                C_kmer=round(v["C_kmer"], 6),
                dnabert2_R2=round(v["dnabert2:R2"], 6)))

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0]))
        w.writeheader()
        w.writerows(rows_out)

    full = sum(1 for r in rows_out for a in ("C_kmer", "dnabert2_R2")
               if r[a] > r["ceiling_full_block"] + 1e-12)
    fold = sum(1 for r in rows_out for a in ("C_kmer", "dnabert2_R2")
               if r[a] > r["ceiling_fold_matched"] + 1e-12)
    d = np.array([r["ceiling_fold_matched"] - r["ceiling_full_block"] for r in rows_out])
    print(f"blocks {len(rows_out)}, cells {2*len(rows_out)}")
    print(f"  cells above full-block ceiling  : {full}")
    print(f"  cells above fold-matched ceiling: {fold}")
    print(f"  fold-matched - full-block (pp)  : mean {100*d.mean():+.4f}, "
          f"range [{100*d.min():+.3f}, {100*d.max():+.3f}]")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
