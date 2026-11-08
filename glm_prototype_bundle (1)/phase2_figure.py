"""Phase 2 figure: embedding distance vs identity; per-region separability."""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("pair_identity_vs_embedding.csv")))
summary = json.load(open("phase2_summary.json"))

order = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# Panel A: identity vs embedding distance
ax = axes[0]
for ptype, color, label in [("within", "#1f77b4", "within B. anthracis"),
                            ("between", "#d62728", "B. anthracis vs cereus group")]:
    x = [float(r["sequence_identity"]) for r in rows if r["pair_type"] == ptype]
    y = [float(r["embedding_cosine_distance"]) for r in rows if r["pair_type"] == ptype]
    ax.scatter(x, y, s=8, alpha=0.35, color=color, label=label)
from numpy import corrcoef
sp = 0
# quick spearman
def rank(v):
    v = np.asarray(v); o = np.argsort(v); r = np.empty(len(v)); r[o] = np.arange(len(v)); return r
xi = [float(r["sequence_identity"]) for r in rows]
yi = [float(r["embedding_cosine_distance"]) for r in rows]
sp = float(np.corrcoef(rank(xi), rank(yi))[0, 1])
ax.set_xlabel("Pairwise sequence identity (region slice)")
ax.set_ylabel("DNABERT-2 embedding cosine distance")
ax.set_title(f"A. Embedding distance tracks identity\nSpearman = {sp:.2f} (n = {len(rows)} pairs)")
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.grid(alpha=0.25)

# Panel B: per-region LOO nearest-neighbor accuracy
ax = axes[1]
accs = [summary[r]["loo_acc"] for r in order]
colors = ["#2ca02c" if a > 0.5 else "#d62728" for a in accs]
ax.bar(np.arange(len(order)), accs, color=colors, width=0.6)
ax.axhline(0.5, color="grey", ls="--", lw=1, label="chance (0.5)")
ax.set_xticks(np.arange(len(order)), order)
ax.set_ylabel("LOO nearest-neighbor accuracy\n(B. anthracis vs cereus group)")
ax.set_title("B. Embedding separability by region")
ax.legend(frameon=False, fontsize=9)
ax.set_ylim(0, 1)
for i, a in enumerate(accs):
    ax.text(i, a + 0.02, f"{a:.2f}", ha="center", fontsize=8)

# Panel C: region ranking comparison — identity vs embedding
ax = axes[2]
btwn_ident = {r: None for r in order}
for r in rows:
    if r["pair_type"] == "between":
        btwn_ident[r["region"]] = float(btwn_ident[r["region"]] or 0) + float(r["sequence_identity"])
counts = {r: sum(1 for x in rows if x["region"] == r and x["pair_type"] == "between") for r in order}
btwn_ident = {r: btwn_ident[r] / counts[r] for r in order}
btwn_embed = {r: summary[r]["between_embed"] for r in order}
x = [1 - btwn_ident[r] for r in order]  # divergence
y = [btwn_embed[r] for r in order]
ax.scatter(x, y, color="#1f77b4")
for r in order:
    ax.annotate(r, (1 - btwn_ident[r], btwn_embed[r]), textcoords="offset points", xytext=(4, 3), fontsize=9)
rr = float(np.corrcoef(rank(x), rank(y))[0, 1])
ax.set_xlim(left=min(x) * 0.9, right=max(x) + (max(x) - min(x)) * 0.25)
ax.set_xlabel("Mean between-group divergence (1 - identity)")
ax.set_ylabel("Mean between-group\nembedding distance")
ax.set_title(f"C. Regions rank consistently\nSpearman = {rr:.2f}")
ax.grid(alpha=0.25)

plt.tight_layout()
plt.savefig("phase2_summary.png", dpi=180, bbox_inches="tight")
print("saved phase2_summary.png")
