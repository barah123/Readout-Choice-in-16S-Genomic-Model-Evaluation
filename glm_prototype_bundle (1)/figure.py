"""Dump pair-level region identities and build the prototype summary figure."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner

records = json.load(open("cereus_group_16s.json"))
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

def identity(a, b):
    if not a or not b:
        return None
    aln = aligner.align(a, b)[0]
    matches = 0
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        seg_t = a[int(t_blk[0]):int(t_blk[1])]
        seg_q = b[int(q_blk[0]):int(q_blk[1])]
        matches += sum(1 for x, y in zip(seg_t, seg_q) if x == y)
    return 2 * matches / (len(a) + len(b))

seqs = []
seen = set()
for grp, hdr, seq in records:
    if grp == "E_coli_ref" or seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    seqs.append((grp, hdr.split()[1] if "strain" not in hdr else hdr.split()[3], slice_regions(seq)))

rows = []  # region, pair_type, identity
for name in REGIONS:
    an = [(h, r[name]) for g, h, r in seqs if g == "B_anthracis" and r[name]]
    ot = [(h, r[name]) for g, h, r in seqs if g in ("B_cereus", "B_thuringiensis", "B_weihenstephanensis") and r[name]]
    for i in range(len(an)):
        for j in range(i + 1, len(an)):
            v = identity(an[i][1], an[j][1])
            if v is not None:
                rows.append((name, "within B. anthracis", v))
    for (h1, s1) in an:
        for (h2, s2) in ot:
            v = identity(s1, s2)
            if v is not None:
                rows.append((name, "B. anthracis vs cereus group", v))

import csv
with open("pair_identities.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["region", "comparison", "identity"])
    w.writerows(rows)
print(f"wrote {len(rows)} pair-level rows")

# ---- figure ----
order = list(REGIONS)
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [2.1, 1]})
ax = axes[0]
regions_x = np.arange(len(order))
for k, (label, color) in enumerate([("within B. anthracis", "#1f77b4"),
                                     ("B. anthracis vs cereus group", "#d62728")]):
    pts = {r_: [] for r_ in order}
    for reg, comp, v in rows:
        if comp == label:
            pts[reg].append(v)
    for i, reg in enumerate(order):
        vs = pts[reg]
        if vs:
            xs = np.full(len(vs), i + (k - 0.5) * 0.35)
            ax.scatter(xs, vs, s=14, alpha=0.45, color=color, label=label if i == 0 else None)
ax.axhline(1.0, color="grey", lw=0.8, ls="--", alpha=0.6)
ax.set_xticks(regions_x, order)
ax.set_xlabel("16S hypervariable region")
ax.set_ylabel("Pairwise sequence identity")
ax.set_title("A. Identity by region: B. anthracis vs B. cereus group (NCBI sampling)")
ax.legend(frameon=False, loc="lower left", fontsize=9)
ax.set_ylim(0.95, 1.012)
ax.grid(axis="y", alpha=0.25)

# Panel B: V6 haplotype sharing
ax2 = axes[1]
anth_v6 = Counter = {}
cereus_v6 = {}
for g, h, r in seqs:
    if not r["V6"]:
        continue
    if g == "B_anthracis":
        anth_v6[r["V6"]] = anth_v6.get(r["V6"], 0) + 1
    elif g in ("B_cereus", "B_thuringiensis", "B_weihenstephanensis"):
        cereus_v6[r["V6"]] = cereus_v6.get(r["V6"], 0) + 1
exclusive = sum(v for h, v in anth_v6.items() if h not in cereus_v6)
shared = sum(v for h, v in anth_v6.items() if h in cereus_v6)
ax2.bar([0, 1], [exclusive, shared], color=["#2ca02c", "#d62728"], width=0.55)
ax2.set_xticks([0, 1], ["V6 haplotype\nunique to B. anthracis", "V6 haplotype\nshared with cereus group"])
ax2.set_ylabel("No. of B. anthracis sequences")
ax2.set_title("B. V6 does not separate B. anthracis\nfrom all cereus-group strains")
for i, v in enumerate([exclusive, shared]):
    ax2.text(i, v + 0.2, str(v), ha="center", fontsize=11)
ax2.set_ylim(0, max(exclusive, shared) + 2)

plt.tight_layout()
plt.savefig("prototype_summary.png", dpi=180, bbox_inches="tight")
print("saved prototype_summary.png")
