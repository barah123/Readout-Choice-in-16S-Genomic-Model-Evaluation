"""AQ06 strand audit: find minus-strand deposits, and test whether the locked-test primary
result depends on them.

No reverse-complement search is performed anywhere in the pipeline. The assumption recorded
during the study was that a minus-strand deposit would fail the per-region length window and
be dropped. That assumption is WRONG and this script is the evidence: global alignment of a
reverse-complemented 16S against the forward reference still yields a contiguous query span
for most regions, and those spans usually fall inside the 0.75-1.30x window, so the record
enters the analysis carrying mis-sliced sequence.

Detection: a record is minus-strand if its reverse complement scores higher against the
E. coli J01695 reference than the deposited orientation, under the same aligner and scoring
used for extraction.

Sensitivity: the frozen locked-test scoring is re-run for every affected block with the
minus-strand records dropped, reusing the existing embedding cache (no re-embedding, and no
re-run of the frozen analysis itself). The primary estimand is then recomputed with those
block values substituted, so the question "does the confirmatory conclusion depend on these
records?" is answered numerically rather than argued.

conda run -n glm-benchmark python audit_strand_sensitivity.py
"""
import glob
import json
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
from Bio.Align import PairwiseAligner
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
sys.path.insert(0, os.path.join(ROOT, "phase4_locked_test"))
sys.path.insert(0, HERE)
import audit_region_extraction as are    # noqa: E402  (slicer + reference)
import score_test as st                  # noqa: E402

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

TEST = os.path.join(ROOT, "phase4_locked_test")
SETS = {
    "development": "phase2_pilot/raw_species/*.json",
    "locked_test": "phase4_locked_test/raw_species_test/*.json",
    "hard_panels": "phase7_extensions/raw_species_hard/*.json",
}
_COMP = str.maketrans("ACGTN", "TGCAN")


def scores(seq):
    al = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                         match_score=2, mismatch_score=-1)
    return float(al.score(are._EC, seq)), float(al.score(are._EC, seq.translate(_COMP)[::-1]))


def find_minus_strand():
    found = {}
    for name, pat in SETS.items():
        recs = []
        for f in sorted(glob.glob(os.path.join(ROOT, pat))):
            recs.extend(json.load(open(f)))
        sc = Parallel(n_jobs=8, batch_size=32)(delayed(scores)(r["sequence"]) for r in recs)
        bad = []
        for r, (fwd, rev) in zip(recs, sc):
            if rev > fwd:
                sl = are._slice_one(r["sequence"])
                bad.append({"accession": r["accession"], "block_id": r["block_id"],
                            "species_query": r["species_query"],
                            "isolate_group": r["isolate_group"], "slen": r["slen"],
                            "score_forward": fwd, "score_reverse_complement": rev,
                            "regions_retained_despite_wrong_strand":
                                [k for k, v in sl.items() if v]})
        found[name] = bad
        print(f"{name:12} records={len(recs):6} minus_strand={len(bad)} "
              f"blocks={sorted({b['block_id'] for b in bad})}")
    return found


def rerun_locked_test(drop_accessions):
    """Re-run score_test.score_block for affected blocks, with and without the drops."""
    rows = json.load(open(os.path.join(TEST, "test_blocks.json")))
    summ = json.load(open(os.path.join(TEST, "test_block_summary.json")))
    raw = []
    for f in sorted(glob.glob(os.path.join(ROOT, SETS["locked_test"]))):
        raw.extend(json.load(open(f)))

    # isolate groups that contain a dropped accession and no other record
    by_group = defaultdict(list)
    for r in raw:
        by_group[r["isolate_group"]].append(r["accession"])
    drop_groups = {g for g, accs in by_group.items()
                   if set(accs) & drop_accessions and not (set(accs) - drop_accessions)}
    mixed = {g for g, accs in by_group.items()
             if set(accs) & drop_accessions and (set(accs) - drop_accessions)}
    print(f"\nisolate groups fully removed by the drop: {len(drop_groups)}; "
          f"groups that also hold a plus-strand record: {len(mixed)}")

    emb = {}
    for m in st.MODELS:
        p = os.path.join(TEST, f"emb_test_{m}.npz")
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        for region in st.REGIONS:
            ks, vs = f"{region}__keys", f"{region}__vecs"
            if ks in z:
                emb[(m, region)] = dict(zip(z[ks].tolist(), z[vs]))

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    affected = sorted({r["block_id"] for r in raw if r["accession"] in drop_accessions})
    out = {}
    for bid in affected:
        info = summ.get(bid, {})
        meta = {"_family": info.get("family"), "_rank": info.get("rank"),
                "_expected_class": info.get("expected_class"),
                "_eligible_regions": info.get("eligible_regions", []),
                "_block_eligible": info.get("block_eligible", False),
                "_ceilings": {r: v.get("ceiling_bal_isolate")
                              for r, v in info.get("regions", {}).items()}}
        base = st.score_block(bid, by_block[bid], emb, dict(meta))
        kept = {r: {k: v for k, v in d.items() if k not in drop_groups}
                for r, d in by_block[bid].items()}
        drop = st.score_block(bid, kept, emb, dict(meta))
        out[bid] = {"as_run": base, "minus_strand_dropped": drop}
        print(f"  {bid:26} C_kmer {base.get('C_kmer', float('nan')):.4f} -> "
              f"{drop.get('C_kmer', float('nan')):.4f}   "
              f"dnabert2:R2 {base.get('dnabert2:R2', float('nan')):.4f} -> "
              f"{drop.get('dnabert2:R2', float('nan')):.4f}")
    return out, drop_groups


def main():
    found = find_minus_strand()
    drop_acc = {b["accession"] for v in found.values() for b in v}

    print("\n--- locked-test sensitivity: frozen scoring re-run without minus-strand records ---")
    rerun, drop_groups = rerun_locked_test(drop_acc)

    # exactly phase5's confirmatory set and bootstrap settings
    scores_all = json.load(open(os.path.join(TEST, "locked_test_scores.json")))
    blocks = sorted(b for b, v in scores_all.items()
                    if v.get("_status") == "ok" and v.get("_block_eligible"))
    def diff(d):
        return d["dnabert2:R2"] - d["C_kmer"]
    as_run = {b: diff(scores_all[b]) for b in blocks}
    adj = dict(as_run)
    for b, v in rerun.items():
        if v["minus_strand_dropped"].get("_status") == "ok":
            adj[b] = diff(v["minus_strand_dropped"])

    a = np.array([as_run[b] for b in blocks]); d = np.array([adj[b] for b in blocks])
    def boot(x):
        rng = np.random.default_rng(1729)          # phase5 frozen seed
        draws = x[rng.integers(0, len(x), size=(2000, len(x)))].mean(axis=1)
        lo, hi = np.percentile(draws, [2.5, 97.5])
        p = 2 * min((draws <= 0).mean(), (draws >= 0).mean())
        return float(lo), float(hi), float(min(p, 1.0))

    def verdict(lo, hi, margin=0.01):
        if lo > margin: return "SUPERIOR (model)"
        if hi < -margin: return "SUPERIOR (classical)"
        if lo > -margin and hi < margin: return "EQUIVALENT (practically negligible difference)"
        return "INCONCLUSIVE at this margin"
    la, ha, pa = boot(a)
    ld, hd, pd_ = boot(d)
    res = {
        "n_blocks": len(blocks),
        "blocks_changed": sorted(b for b in blocks if adj[b] != as_run[b]),
        "primary_as_run": {"mean": float(a.mean()), "ci95": [la, ha], "p": pa,
                           "verdict": verdict(la, ha)},
        "primary_minus_strand_dropped": {"mean": float(d.mean()), "ci95": [ld, hd], "p": pd_,
                                         "verdict": verdict(ld, hd)},
        "largest_single_block_shift_pp": float(np.abs(d - a).max() * 100),
        "shift_in_primary_pp": float((d.mean() - a.mean()) * 100),
        "margin": 0.01,
    }
    print("\n" + json.dumps(res, indent=2))

    dest = os.path.abspath(os.path.join(HERE, "..", "results", "strand_audit.json"))
    json.dump({"minus_strand_records": found,
               "isolate_groups_dropped": sorted(drop_groups),
               "locked_test_block_rerun": rerun,
               "primary_sensitivity": res}, open(dest, "w"), indent=2, default=str)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
