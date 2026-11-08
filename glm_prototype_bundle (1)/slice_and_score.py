"""Phase 0-1 prototype: slice V1-V9 by E. coli J01695 coordinates, compute per-region
identity margins for B. anthracis (target) vs other B. cereus-group members."""
import json
from Bio import Align
from Bio.Align import PairwiseAligner

records = json.load(open("cereus_group_16s.json"))
groups = [r for r in records if r[0] != "E_coli_ref"]
ec_ref = json.load(open("ecoli_16s.txt"))["seq"]
assert len(ec_ref) > 1500, f"E. coli ref too short: {len(ec_ref)}"

# Chakravorty et al. 2007 V-region coordinates (E. coli J01695 numbering, 1-based inclusive)
REGIONS = {
    "V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
    "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
    "V8": (1243, 1294), "V9": (1435, 1465),
}

aligner = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                          match_score=2, mismatch_score=-1)

def slice_regions(seq):
    """Align seq to E. coli ref; return {region: subsequence}."""
    aln = aligner.align(ec_ref, seq)[0]
    # aln.aligned shape = (2, n_blocks, 2): [target_blocks, query_blocks]
    ref_to_query = {}
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1 = int(t_blk[0]), int(t_blk[1])
        q0, q1 = int(q_blk[0]), int(q_blk[1])
        for i in range(t0, t1):
            ref_to_query[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in REGIONS.items():
        qs = [ref_to_query[p] for p in range(s - 1, e) if p in ref_to_query]
        out[name] = seq[min(qs):max(qs) + 1] if qs else ""
    return out

# dedupe identical full sequences; QA: require near-full-length genes
seen = set()
seqs = []  # (group, region_dict)
for grp, hdr, seq in groups:
    if seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    seqs.append((grp, slice_regions(seq)))
print(f"{len(seqs)} unique near-full-length sequences after QA", flush=True)

# QA region slices: drop slices whose length deviates >25% from expected
EXPECTED_LEN = {k: e - s + 1 for k, (s, e) in REGIONS.items()}
def ok_slice(sl, name):
    exp = EXPECTED_LEN[name]
    return sl and 0.75 * exp <= len(sl) <= 1.30 * exp

def region_identity(a, b):
    if not a or not b:
        return None
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return None
    aln = aligner.align(a, b)[0]
    matches = 0
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        seg_t = a[int(t_blk[0]):int(t_blk[1])]
        seg_q = b[int(q_blk[0]):int(q_blk[1])]
        matches += sum(1 for x, y in zip(seg_t, seg_q) if x == y)
    return 2 * matches / (la + lb)

import numpy as np
results = {}
for name in REGIONS:
    an = [r[name] for grp, r in seqs if grp == "B_anthracis" and ok_slice(r[name], name)]
    other = [r[name] for grp, r in seqs if grp in ("B_cereus", "B_thuringiensis", "B_weihenstephanensis") and ok_slice(r[name], name)]
    sub = [r[name] for grp, r in seqs if grp == "B_subtilis" and ok_slice(r[name], name)]
    # within-anthracis
    within = [region_identity(x, y) for i, x in enumerate(an) for y in an[i+1:]]
    within = [v for v in within if v is not None]
    # anthracis vs other cereus-group
    between = [region_identity(x, y) for x in an for y in other]
    between = [v for v in between if v is not None]
    # anthracis vs outgroup (B. subtilis)
    vs_out = [region_identity(x, y) for x in an for y in sub]
    vs_out = [v for v in vs_out if v is not None]
    if within and between:
        results[name] = {
            "n_anthracis": len(an), "n_other_cereus": len(other), "n_subtilis": len(sub),
            "within_min": min(within), "within_mean": float(np.mean(within)),
            "between_max": max(between), "between_mean": float(np.mean(between)),
            "margin": max(between) - min(within),  # <= 0 means cross-group identities fall below within-group -> separable
            "vs_subtilis_mean": float(np.mean(vs_out)) if vs_out else None,
        }

header = f"{'region':6}{'nA':>4}{'nO':>4}{'within_min':>11}{'within_mean':>12}{'btwn_max':>10}{'btwn_mean':>11}{'sep_margin':>11}{'vs_subtilis':>12}"
print(header)
for name, r in results.items():
    print(f"{name:6}{r['n_anthracis']:>4}{r['n_other_cereus']:>4}"
          f"{r['within_min']:>11.4f}{r['within_mean']:>12.4f}{r['between_max']:>10.4f}"
          f"{r['between_mean']:>11.4f}{r['margin']:>11.4f}"
          f"{(r['vs_subtilis_mean'] if r['vs_subtilis_mean'] is not None else float('nan')):>12.4f}")

json.dump(results, open("region_margins.json", "w"), indent=2)
