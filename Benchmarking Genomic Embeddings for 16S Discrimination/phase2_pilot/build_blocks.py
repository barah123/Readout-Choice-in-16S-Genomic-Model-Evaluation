"""Phase 2 step 3: slice V1-V9, compute matched ambiguity ceilings, and emit the
analysis-ready block table. Shared by every downstream readout/baseline script.

Design decisions declared here, before any model result is seen:

- Unit of independence = isolate_group (protocol §12.1). All records from one isolate stay
  together in any split. This is strictly stronger than the Phase 0-3a string-dedup.
- Records are NOT deduplicated by sequence string. Identical strings from *different* isolates
  are biologically real repeated observations and are kept, with the dependence handled by
  grouping, per protocol §8.3 ("duplicating a shared string across biological genomes does not
  create new sequence information; the dependence must be handled").
- Ceilings are computed isolate-weighted (each isolate contributes once per region, via its
  modal slice) AND record-weighted, reported separately, because protocol §8.1 requires the
  ceiling to match the metric's weighting rather than assuming one universal ceiling.
- Slicing logic copied verbatim from phase3a_staph_baseline.py per the CLAUDE.md convention.
"""
import json
from collections import Counter, defaultdict

import numpy as np
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                          match_score=2, mismatch_score=-1)
REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
           "V8": (1243, 1294), "V9": (1435, 1465)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}

PANELS = {
    "staph_CoNS": dict(path="staph_panel_meta.json", target="S_aureus"),
    "cereus_group": dict(path="cereus_panel_meta.json", target="B_anthracis"),
}


def slice_regions(seq, ec_ref):
    aln = aligner.align(ec_ref, seq)[0]
    ref_to_query = {}
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1 = int(t_blk[0]), int(t_blk[1]); q0, q1 = int(q_blk[0]), int(q_blk[1])
        for i in range(t0, t1):
            ref_to_query[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in REGIONS.items():
        qs = [ref_to_query[p] for p in range(s - 1, e) if p in ref_to_query]
        sl = seq[min(qs):max(qs) + 1] if qs else ""
        out[name] = sl if sl and 0.75 * EXPECTED[name] <= len(sl) <= 1.30 * EXPECTED[name] else None
    return out


def ceilings(labels_counts, n_pos, n_neg):
    """labels_counts: {sequence -> [n_positive, n_negative]} under one weighting scheme."""
    total = n_pos + n_neg
    if not total or not n_pos or not n_neg:
        return float("nan"), float("nan")
    err = sum(min(a, b) for a, b in labels_counts.values())
    acc = 1 - err / total
    bal = sum(max(a / (2 * n_pos), b / (2 * n_neg)) for a, b in labels_counts.values())
    return acc, bal


def main():
    ec_ref = json.load(open("../../glm_prototype_bundle (1)/ecoli_16s.txt"))["seq"]
    out_rows = []
    summary = {}

    for panel, cfg in PANELS.items():
        recs = json.load(open(cfg["path"]))
        target = cfg["target"]
        print("=" * 100)
        print(f"PANEL {panel}   target={target}   records={len(recs)}")
        print("=" * 100)

        # slice every record once
        for r in recs:
            r["slices"] = slice_regions(r["sequence"], ec_ref)

        summary[panel] = {}
        hdr = (f"{'reg':5}{'N_rec':>7}{'N_iso':>7}{'tgt_iso':>8}{'neg_iso':>8}"
               f"{'hap_tgt':>9}{'hap_neg':>9}{'ambig':>7}"
               f"{'ceil_acc_iso':>14}{'ceil_bal_iso':>14}{'ceil_bal_rec':>14}")
        print(hdr)

        for region in REGIONS:
            present = [r for r in recs if r["slices"][region]]
            if not present:
                continue

            # --- isolate-weighted view: one modal slice per isolate ---
            by_iso = defaultdict(list)
            for r in present:
                by_iso[r["isolate_group"]].append(r)
            iso_rows = []
            for iso, rs in by_iso.items():
                modal = Counter(x["slices"][region] for x in rs).most_common(1)[0][0]
                iso_rows.append((iso, rs[0]["species_group"], modal))

            tgt_iso = [x for x in iso_rows if x[1] == target]
            neg_iso = [x for x in iso_rows if x[1] != target]

            cnt_iso = defaultdict(lambda: [0, 0])
            for _, sp, sl in iso_rows:
                cnt_iso[sl][0 if sp == target else 1] += 1
            ceil_acc_iso, ceil_bal_iso = ceilings(cnt_iso, len(tgt_iso), len(neg_iso))

            # --- record-weighted view (what Phase 0-3a implicitly used) ---
            cnt_rec = defaultdict(lambda: [0, 0])
            n_t_rec = n_n_rec = 0
            for r in present:
                is_t = r["species_group"] == target
                cnt_rec[r["slices"][region]][0 if is_t else 1] += 1
                n_t_rec += is_t
                n_n_rec += not is_t
            _, ceil_bal_rec = ceilings(cnt_rec, n_t_rec, n_n_rec)

            hap_t = len({sl for _, sp, sl in iso_rows if sp == target})
            hap_n = len({sl for _, sp, sl in iso_rows if sp != target})
            ambig = sum(1 for v in cnt_iso.values() if v[0] and v[1])

            summary[panel][region] = dict(
                n_records=len(present), n_isolates=len(iso_rows),
                target_isolates=len(tgt_iso), neighbor_isolates=len(neg_iso),
                hap_target=hap_t, hap_neighbor=hap_n, ambiguous_haplotypes=ambig,
                ceiling_acc_isolate=ceil_acc_iso, ceiling_bal_isolate=ceil_bal_iso,
                ceiling_bal_record=ceil_bal_rec,
            )
            print(f"{region:5}{len(present):7}{len(iso_rows):7}{len(tgt_iso):8}{len(neg_iso):8}"
                  f"{hap_t:9}{hap_n:9}{ambig:7}"
                  f"{ceil_acc_iso:14.3f}{ceil_bal_iso:14.3f}{ceil_bal_rec:14.3f}")

            for iso, sp, sl in iso_rows:
                out_rows.append(dict(panel=panel, region=region, isolate_group=iso,
                                     species_group=sp, is_target=int(sp == target),
                                     slice=sl))
        print()

    with open("blocks.json", "w") as f:
        json.dump(out_rows, f)
    with open("block_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"wrote blocks.json ({len(out_rows)} isolate x region rows) and block_summary.json")


if __name__ == "__main__":
    main()
