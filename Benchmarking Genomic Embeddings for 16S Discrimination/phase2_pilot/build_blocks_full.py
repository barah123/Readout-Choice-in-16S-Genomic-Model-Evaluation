"""Slice V1-V9 for every block, compute isolate-weighted ambiguity ceilings, apply the
pre-declared eligibility gate, and emit the analysis-ready table.

Same design decisions as build_blocks.py (isolate is the unit of independence; records are
not string-deduplicated; ceilings are isolate-weighted to match the metric), extended to the
full 31-block set and parallelised across records because alignment dominates at this scale.
"""
import json
import os
from collections import Counter, defaultdict

import numpy as np
from Bio.Align import PairwiseAligner
from joblib import Parallel, delayed

from blocks_definition import BLOCKS

REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
           "V8": (1243, 1294), "V9": (1435, 1465)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}

# pre-declared eligibility gate (blocks_definition.py rule 6)
MIN_TARGET_ISOLATES = 10
MIN_NEIGHBOUR_ISOLATES = 5
MIN_NEIGHBOUR_SPECIES = 2

_EC = json.load(open("../../glm_prototype_bundle (1)/ecoli_16s.txt"))["seq"]


def _slice_one(seq):
    al = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                         match_score=2, mismatch_score=-1)
    aln = al.align(_EC, seq)[0]
    ref_to_query = {}
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1 = int(t_blk[0]), int(t_blk[1]); q0 = int(q_blk[0])
        for i in range(t0, t1):
            ref_to_query[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in REGIONS.items():
        qs = [ref_to_query[p] for p in range(s - 1, e) if p in ref_to_query]
        sl = seq[min(qs):max(qs) + 1] if qs else ""
        out[name] = sl if sl and 0.75 * EXPECTED[name] <= len(sl) <= 1.30 * EXPECTED[name] else None
    return out


def ceilings(counts, n_pos, n_neg):
    total = n_pos + n_neg
    if not total or not n_pos or not n_neg:
        return float("nan"), float("nan")
    err = sum(min(a, b) for a, b in counts.values())
    acc = 1 - err / total
    bal = sum(max(a / (2 * n_pos), b / (2 * n_neg)) for a, b in counts.values())
    return acc, bal


def load_block(b):
    recs = []
    for sp in [b["target"]] + b["neighbours"]:
        safe = sp.replace(" ", "_").replace("/", "_").replace(".", "")
        path = os.path.join("raw_species", f"{b['block_id']}__{safe}.json")
        if os.path.exists(path):
            recs.extend(json.load(open(path)))
    return recs


def main():
    all_rows, summary = [], {}

    for b in BLOCKS:
        recs = load_block(b)
        if not recs:
            print(f"[skip] {b['block_id']}: no records", flush=True)
            continue

        seqs = [r["sequence"] for r in recs]
        sliced = Parallel(n_jobs=8, batch_size=16)(delayed(_slice_one)(s) for s in seqs)
        for r, sl in zip(recs, sliced):
            r["slices"] = sl

        binfo = {"block_id": b["block_id"], "family": b["family"], "rank": b["rank"],
                 "target": b["target"], "expected_class": b["expected_class"],
                 "n_records": len(recs), "regions": {}}

        for region in REGIONS:
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
            acc, bal = ceilings(cnt, n_t, n_n)

            per_nb = Counter(sp for _, sp, is_t, _ in iso_rows if not is_t)
            nb_ok = sum(1 for v in per_nb.values() if v >= MIN_NEIGHBOUR_ISOLATES)

            binfo["regions"][region] = dict(
                n_isolates=len(iso_rows), target_isolates=n_t, neighbour_isolates=n_n,
                hap_target=len({sl for _, _, t, sl in iso_rows if t}),
                hap_neighbour=len({sl for _, _, t, sl in iso_rows if not t}),
                ambiguous_haplotypes=sum(1 for v in cnt.values() if v[0] and v[1]),
                ceiling_acc_isolate=acc, ceiling_bal_isolate=bal,
                neighbour_species_ok=nb_ok,
                eligible=bool(n_t >= MIN_TARGET_ISOLATES and nb_ok >= MIN_NEIGHBOUR_SPECIES),
            )
            for iso, sp, is_t, sl in iso_rows:
                all_rows.append(dict(block_id=b["block_id"], region=region,
                                     isolate_group=iso, species=sp,
                                     is_target=int(is_t), slice=sl))

        elig = [r for r, v in binfo["regions"].items() if v["eligible"]]
        binfo["eligible_regions"] = elig
        binfo["block_eligible"] = bool(elig)
        summary[b["block_id"]] = binfo
        mx = max((v["ceiling_bal_isolate"] for v in binfo["regions"].values()
                  if v["ceiling_bal_isolate"] == v["ceiling_bal_isolate"]), default=float("nan"))
        print(f"{b['block_id']:32} recs={len(recs):5} eligible_regions={len(elig)}/9 "
              f"max_ceiling={mx:.3f}", flush=True)

    json.dump(all_rows, open("blocks_full.json", "w"))
    json.dump(summary, open("block_summary_full.json", "w"), indent=2)
    nelig = sum(1 for v in summary.values() if v["block_eligible"])
    print(f"\nwrote blocks_full.json ({len(all_rows)} isolate x region rows)")
    print(f"blocks with >=1 confirmatory-eligible region: {nelig}/{len(summary)}")


if __name__ == "__main__":
    main()
