"""Phase 3 pre-check for the S. aureus / CoNS target (MASTER_PLAN.md §4 item 2).

Before fetching more data or spending any embedding compute on this target, apply the
haplotype gate Phase 2c established: slice V1-V9 with the canonical slice_regions() function,
dedup to unique near-full-length sequences, and count unique slice HAPLOTYPES per class per
region (not raw sequence count - see CLAUDE.md "Before scoring any classifier..."). A region
is scorable only if both classes have >=2 unique haplotypes (and ideally >=5 per the
phase2cd_report.md §5.3 suggested minimum).
"""
import json
from collections import defaultdict
from Bio.Align import PairwiseAligner

records = json.load(open("staph_group_16s.json"))
ec_ref = json.load(open("ecoli_16s.txt"))["seq"]
aligner = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                          match_score=2, mismatch_score=-1)
REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
           "V8": (1243, 1294), "V9": (1435, 1465)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}


def slice_regions(seq):
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


seqs, seen = [], set()
for grp, hdr, seq in records:
    if seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    seqs.append((grp, hdr.split()[0], slice_regions(seq)))

counts = defaultdict(int)
for grp, _, _ in seqs:
    counts[grp] += 1
print(f"{len(seqs)} unique near-full-length sequences after QA (of {len(records)} raw)")
for g, n in sorted(counts.items(), key=lambda kv: -kv[1]):
    print(f"  {g:18} n={n}")

target = "S_aureus"
an = [i for i, (g, _, _) in enumerate(seqs) if g == target]
cg = [i for i, (g, _, _) in enumerate(seqs) if g != target]
print(f"\ntarget = {target} (n={len(an)}) vs CoNS pooled (n={len(cg)})\n")

print(f"{'reg':4} {'N':>4} {'n_tgt':>6} {'n_cons':>7} {'hap_tgt':>8} {'hap_cons':>9} "
      f"{'hap_ambig':>10} {'scorable':>9} {'ample':>6} {'oracle_acc':>10} {'oracle_bal':>10}")
for name in REGIONS:
    an_p = [seqs[i][2][name] for i in an if seqs[i][2][name] is not None]
    cg_p = [seqs[i][2][name] for i in cg if seqs[i][2][name] is not None]
    hap_an, hap_cg = set(an_p), set(cg_p)
    ambig = len(hap_an & hap_cg)
    scorable = len(hap_an) >= 2 and len(hap_cg) >= 2
    ample = len(hap_an) >= 5 and len(hap_cg) >= 5

    # oracle ceiling: per Bayes-optimal rule given haplotype->label counts (phase2c convention)
    hap = defaultdict(lambda: [0, 0])  # haplotype -> [n_target, n_cons]
    for s in an_p:
        hap[s][0] += 1
    for s in cg_p:
        hap[s][1] += 1
    N = len(an_p) + len(cg_p)
    err = sum(min(a, c) for a, c in hap.values())
    oracle_acc = 1 - err / N if N else float("nan")
    oracle_bal = sum(max(a / (2 * len(an_p)), c / (2 * len(cg_p))) for a, c in hap.values()) if an_p and cg_p else float("nan")

    print(f"{name:4} {len(an_p)+len(cg_p):4} {len(an_p):6} {len(cg_p):7} {len(hap_an):8} "
          f"{len(hap_cg):9} {ambig:10} {str(scorable):>9} {str(ample):>6} "
          f"{oracle_acc:10.3f} {oracle_bal:10.3f}")
