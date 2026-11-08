"""Phase 3a summary figure. Reads phase3a_staph_summary.json; writes phase3a_summary.png."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = json.load(open("phase3a_staph_summary.json"))
REG = [f"V{i}" for i in range(1, 10)]
x = np.arange(len(REG))
nan = float("nan")


def g(r, k):
    v = S.get(r, {}).get(k)
    return v if v is not None else nan


fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# ---- A: three arms vs oracle ceiling ----
ax = axes[0, 0]
b1 = [g(r, "b1_bal") for r in REG]
b2 = [g(r, "rf_kmer_bal") for r in REG]
ge = [g(r, "rf_embed_bal") for r in REG]
orc = [g(r, "oracle_bal_acc") for r in REG]
w = 0.25
ax.bar(x - w, b1, w, label="B1: identity 1-NN", color="#7f8c8d")
ax.bar(x, b2, w, label="B2: k-mer + RF", color="#2980b9")
ax.bar(x + w, ge, w, label="gLM: DNABERT-2 + RF", color="#27ae60")
ax.plot(x, orc, "k_", ms=22, mew=1.6, label="oracle ceiling")
ax.axhline(0.5, ls=":", c="k", lw=1.2)
ax.annotate("chance", (8.45, 0.515), ha="right", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(REG); ax.set_ylim(0, 1.05)
ax.set_ylabel("balanced accuracy (LOHO-CV)")
ax.set_title("A  All three arms vs. the Bayes ceiling\n(BH-significant cells only - see panel D)", fontsize=10)
ax.legend(fontsize=8, loc="lower right"); ax.grid(axis="y", alpha=.3)

# ---- B: V1 close-up - the headline ----
ax = axes[0, 1]
labs = ["B1\nidentity", "B2\nk-mer+RF", "gLM\nembed+RF"]
vals = [g("V1", "b1_bal"), g("V1", "rf_kmer_bal"), g("V1", "rf_embed_bal")]
cols = ["#7f8c8d", "#2980b9", "#27ae60"]
bars = ax.bar(labs, vals, color=cols)
ax.axhline(0.5, ls=":", c="k", lw=1.2)
ax.axhline(g("V1", "oracle_bal_acc"), ls="--", c="k", lw=1.4)
ax.annotate("oracle ceiling = 1.000", (1.5, 1.015), ha="center", fontsize=8)
ax.annotate("chance", (2.35, 0.47), fontsize=7)
for b, v in zip(bars, vals):
    ax.annotate(f"{v:.3f}", (b.get_x() + b.get_width() / 2, v + 0.02), ha="center", fontsize=9, fontweight="bold")
ax.set_ylim(0, 1.1)
ax.set_ylabel("balanced accuracy (LOHO-CV)")
ax.set_title("B  V1 (published 'best' region, oracle = 1.000):\nk-mer nearly solves it, gLM is at exact chance", fontsize=10)
ax.grid(axis="y", alpha=.3)

# ---- C: haplotype structure (why V1/V9 are fair tests) ----
ax = axes[1, 0]
ht = [g(r, "hap_tgt") for r in REG]
hc = [g(r, "hap_cns") for r in REG]
amb = [g(r, "hap_ambig") for r in REG]
ax.bar(x - 0.2, ht, 0.4, label="S. aureus (target)", color="#c0392b")
ax.bar(x + 0.2, hc, 0.4, label="CoNS pooled", color="#2980b9")
for i, a in enumerate(amb):
    if a == a and a > 0:
        ax.annotate(f"{int(a)} ambig", (i, max(ht[i], 2) + 1), ha="center", fontsize=7, color="#c0392b")
ax.set_xticks(x); ax.set_xticklabels(REG)
ax.set_ylabel("unique slice haplotypes")
ax.set_title("C  Haplotype structure: V1 and V9 have ZERO\ncross-class ambiguity (oracle = 1.000)", fontsize=10)
ax.legend(fontsize=8); ax.grid(axis="y", alpha=.3)

# ---- D: BH q-values ----
ax = axes[1, 1]
qv = S["_bh"]["qvalues"]
cells = sorted(qv.items(), key=lambda kv: kv[1])
lbl = [c for c, _ in cells]
q = [v for _, v in cells]
cols2 = {"B1_identity": "#7f8c8d", "B2_kmer": "#2980b9", "gLM_embed": "#27ae60"}
bar_cols = [cols2[c.split(":")[1]] for c in lbl]
y = np.arange(len(lbl))
ax.barh(y, q, color=bar_cols)
ax.axvline(0.05, ls="--", c="#c0392b", lw=1.5)
ax.annotate("alpha = 0.05", (0.06, len(lbl) - 1.2), fontsize=8, color="#c0392b")
ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("BH q-value"); ax.set_xlim(0, max(q) * 1.1)
ax.set_title(f"D  {S['_bh']['n_tests']} (region x arm) tests, BH-corrected:\n"
             f"{len(S['_bh']['survivors'])} survive (gLM: V2, V9 only)", fontsize=10)
ax.grid(axis="x", alpha=.3)

fig.suptitle("Phase 3a - S. aureus vs. coagulase-negative staphylococci: the first real (non-null) head-to-head\n"
             "At the published-best region (V1), k-mer+RF nearly solves it (0.997); DNABERT-2 embedding is at exact chance (0.500)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("phase3a_summary.png", dpi=160)
print("wrote phase3a_summary.png")
