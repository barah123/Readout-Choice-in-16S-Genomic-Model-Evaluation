"""Phase 2c+2d summary figure. Reads phase2c_b2_rf_summary.json; writes phase2cd_summary.png."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = json.load(open("phase2c_b2_rf_summary.json"))
REG = [f"V{i}" for i in range(1, 10)]
x = np.arange(len(REG))
nan = float("nan")

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# ---- A: reported LOO-NN vs the Bayes ceiling it exceeds -------------------
ax = axes[0, 0]
rep = [S[r]["reported_embed_loo"] for r in REG]
orc = [S[r]["oracle_acc"] for r in REG]
ax.bar(x - 0.2, rep, 0.4, label="Phase 2 reported LOO-NN", color="#c0392b")
ax.bar(x + 0.2, orc, 0.4, label="Bayes-optimal ceiling", color="#7f8c8d")
for i, r in enumerate(REG):
    if S[r]["reported_exceeds_oracle"]:
        ax.annotate("impossible", (i, rep[i]), textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=7, color="#c0392b", fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(REG); ax.set_ylim(0, 1)
ax.set_ylabel("accuracy")
ax.set_title("A  Erratum A: reported LOO-NN exceeds what any\nclassifier could achieve (V1, V4, V7)", fontsize=10)
ax.legend(fontsize=8); ax.grid(axis="y", alpha=.3)

# ---- B: haplotype collapse ------------------------------------------------
ax = axes[0, 1]
han = [S[r]["hap_an"] for r in REG]
hcg = [S[r]["hap_cg"] for r in REG]
ax.bar(x - 0.2, han, 0.4, label="B. anthracis (target)", color="#c0392b")
ax.bar(x + 0.2, hcg, 0.4, label="B. cereus group", color="#2980b9")
ax.axhline(2, ls="--", c="k", lw=1)
ax.annotate("below this line: not estimable", (8.4, 2.15), ha="right", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(REG)
ax.set_ylabel("unique slice haplotypes")
ax.set_title("B  The real constraint: 1-3 effective observations\nin the target class, not 9 sequences", fontsize=10)
ax.legend(fontsize=8); ax.grid(axis="y", alpha=.3)

# ---- C: corrected region-level result -------------------------------------
ax = axes[1, 0]
km = [S[r]["rf_kmer_bal"] if S[r]["rf_kmer_bal"] is not None else nan for r in REG]
em = [S[r]["rf_embed_bal"] if S[r]["rf_embed_bal"] is not None else nan for r in REG]
cl = [S[r]["oracle_bal_acc"] for r in REG]
ax.bar(x - 0.2, km, 0.4, label="B2: k-mer + RF", color="#2980b9")
ax.bar(x + 0.2, em, 0.4, label="gLM: DNABERT-2 + RF", color="#27ae60")
ax.plot(x, cl, "k_", ms=22, mew=1.6, label="ceiling")
ax.axhline(0.5, ls=":", c="k", lw=1.2)
ax.annotate("chance", (8.45, 0.515), ha="right", fontsize=7)
for i, r in enumerate(REG):
    if not S[r]["estimable"]:
        ax.annotate("n/e", (i, 0.03), ha="center", fontsize=8, color="#555")
ax.set_xticks(x); ax.set_xticklabels(REG); ax.set_ylim(0, 1)
ax.set_ylabel("balanced accuracy (leave-one-haplotype-out CV)")
ax.set_title("C  Corrected: gLM arm is flat at chance;\nnothing survives BH correction", fontsize=10)
ax.legend(fontsize=8, loc="upper right"); ax.grid(axis="y", alpha=.3)

# ---- D: BH q-values --------------------------------------------------------
ax = axes[1, 1]
qv = S["_bh"]["qvalues"]
cells = sorted(qv.items(), key=lambda kv: kv[1])
lbl = [c for c, _ in cells]
q = [v for _, v in cells]
cols = ["#27ae60" if "embed" in c else "#2980b9" for c in lbl]
y = np.arange(len(lbl))
ax.barh(y, q, color=cols)
ax.axvline(0.05, ls="--", c="#c0392b", lw=1.5)
ax.annotate("alpha = 0.05", (0.057, len(lbl) - 1.2), fontsize=8, color="#c0392b")
ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=8); ax.invert_yaxis()
ax.set_xlabel("BH q-value"); ax.set_xlim(0, 1)
ax.set_title(f"D  {S['_bh']['n_tests']} (region x arm) tests, BH-corrected:\n"
             f"{len(S['_bh']['survivors'])} survivors (min q = {min(q):.2f})", fontsize=10)
ax.grid(axis="x", alpha=.3)

fig.suptitle("Phase 2c+2d - real B2 baseline, and the two errata that invalidate the Phase 2/2b region-level results\n"
             "B. anthracis vs B. cereus group (Panel A expected-null target)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("phase2cd_summary.png", dpi=160)
print("wrote phase2cd_summary.png")
