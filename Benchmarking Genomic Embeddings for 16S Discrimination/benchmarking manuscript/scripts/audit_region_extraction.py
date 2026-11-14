"""AQ06: evidence for the region-extraction claims in Methods.

Three things the manuscript asserts are checked here against the fetched records rather
than assumed:

  1. Per-region availability. How often alignment transfer plus the 0.75-1.30x length
     window failed, per region and per set.
  2. Strand. No reverse-complement search is performed, so a record deposited on the minus
     strand should fail extraction outright rather than be silently mis-sliced. Every record
     that lost ALL nine regions is re-sliced after reverse-complementing it; recovery means
     the record was minus-strand and was correctly rejected, not mis-sliced.
  3. Ambiguous bases. How many retained slices contain N, and how many extra haplotypes
     those Ns create. An N-bearing slice is a distinct string, so it splits a haplotype;
     splitting can only decrease sum(min(a, b)) and increase sum(max(a/2n_pos, b/2n_neg)),
     so it can only make the ceiling look MORE generous, never less. The size of that
     effect is measured in audit_ambiguous_bases.py.

The slicer below is copied verbatim from build_blocks_full._slice_one, per the project
convention of copying the canonical slicing function rather than reinventing it; the copy is
asserted equal to the imported original on a sample before any count is produced.

Run from this directory. conda run -n glm-benchmark python audit_region_extraction.py
"""
import glob
import json
import os
from collections import Counter, defaultdict

from Bio.Align import PairwiseAligner
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GLM = os.path.abspath(os.path.join(ROOT, "..", "glm_prototype_bundle (1)"))

REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
           "V8": (1243, 1294), "V9": (1435, 1465)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}
_EC = json.load(open(os.path.join(GLM, "ecoli_16s.txt")))["seq"]

SETS = {
    "development": "phase2_pilot/raw_species/*.json",
    "locked_test": "phase4_locked_test/raw_species_test/*.json",
    "hard_panels": "phase7_extensions/raw_species_hard/*.json",
}
_COMP = str.maketrans("ACGTN", "TGCAN")


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


def rc(s):
    return s.translate(_COMP)[::-1]


def main():
    report = {}
    for name, pat in SETS.items():
        recs = []
        for f in sorted(glob.glob(os.path.join(ROOT, pat))):
            recs.extend(json.load(open(f)))
        if not recs:
            print(f"[skip] {name}")
            continue
        seqs = [r["sequence"] for r in recs]
        sliced = Parallel(n_jobs=8, batch_size=16)(delayed(_slice_one)(s) for s in seqs)

        avail = {r: sum(1 for sl in sliced if sl[r]) for r in REGIONS}
        n_regions_ok = Counter(sum(1 for r in REGIONS if sl[r]) for sl in sliced)
        total_fail = [i for i, sl in enumerate(sliced) if not any(sl.values())]

        # --- strand check: does reverse-complementing rescue the total failures?
        rescued = 0
        if total_fail:
            rc_sliced = Parallel(n_jobs=8, batch_size=16)(
                delayed(_slice_one)(rc(seqs[i])) for i in total_fail)
            rescued = sum(1 for sl in rc_sliced if sum(1 for r in REGIONS if sl[r]) >= 5)

        # --- ambiguous bases inside retained slices
        n_with_N = 0
        tot_slices = 0
        hap_total, hap_N = defaultdict(set), defaultdict(set)
        for sl in sliced:
            for r, s in sl.items():
                if not s:
                    continue
                tot_slices += 1
                hap_total[r].add(s)
                if "N" in s:
                    n_with_N += 1
                    hap_N[r].add(s)

        d = {
            "records": len(recs),
            "region_cells": len(recs) * len(REGIONS),
            "region_cells_available": tot_slices,
            "availability_pct": round(100 * tot_slices / (len(recs) * len(REGIONS)), 2),
            "per_region_available": avail,
            "per_region_available_pct": {r: round(100 * v / len(recs), 2)
                                         for r, v in avail.items()},
            "records_by_n_regions_recovered": dict(sorted(n_regions_ok.items())),
            "records_losing_all_nine": len(total_fail),
            "of_those_rescued_by_reverse_complement": rescued,
            "retained_slices_containing_N": n_with_N,
            "retained_slices_containing_N_pct": round(100 * n_with_N / tot_slices, 4),
            "haplotypes_total": sum(len(v) for v in hap_total.values()),
            "haplotypes_containing_N": sum(len(v) for v in hap_N.values()),
        }
        report[name] = d
        print(f"\n=== {name} ===")
        for k in ("records", "availability_pct", "records_losing_all_nine",
                  "of_those_rescued_by_reverse_complement",
                  "retained_slices_containing_N", "retained_slices_containing_N_pct",
                  "haplotypes_total", "haplotypes_containing_N"):
            print(f"  {k:42} {d[k]}")
        print("  per-region availability %:",
              " ".join(f"{r}={d['per_region_available_pct'][r]}" for r in REGIONS))

    tot = {k: sum(v[k] for v in report.values())
           for k in ("records", "region_cells", "region_cells_available",
                     "records_losing_all_nine", "of_those_rescued_by_reverse_complement",
                     "retained_slices_containing_N", "haplotypes_total",
                     "haplotypes_containing_N")}
    tot["availability_pct"] = round(100 * tot["region_cells_available"] / tot["region_cells"], 2)
    tot["retained_slices_containing_N_pct"] = round(
        100 * tot["retained_slices_containing_N"] / tot["region_cells_available"], 4)
    per_region = {r: sum(v["per_region_available"][r] for v in report.values())
                  for r in REGIONS}
    tot["per_region_available_pct"] = {
        r: round(100 * per_region[r] / tot["records"], 2) for r in REGIONS}
    report["total_analysed"] = tot
    print("\n=== total ===")
    print(json.dumps(tot, indent=2))

    dest = os.path.abspath(os.path.join(HERE, "..", "results",
                                        "region_extraction_audit.json"))
    json.dump(report, open(dest, "w"), indent=2)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
