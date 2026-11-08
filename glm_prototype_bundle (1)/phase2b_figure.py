"""Phase 2b figure: is the embedding signal redundant with k-mer/GC/length composition?"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pairs = pd.read_csv("pair_identity_vs_embedding_vs_kmer.csv")
summ = json.load(open("phase2b_kmer_diagnostic_summary.json"))
phase2 = json.load(open("phase2_summary.json"))
kmer_region = summ["region_level_kmer_summary"]
order = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"]

def rank(v):
    v = np.asarray(v); o = np.argsort(v); r = np.empty(len(v)); r[o] = np.arange(len(v)); return r

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# Panel A: kmer_dist vs embedding_dist, colored by pair_type
ax = axes[0]
for ptype, color, label in [("within", "#1f77b4", "within B. anthracis"),
                            ("between", "#d62728", "B. anthracis vs cereus group")]:
    sub = pairs[pairs["pair_type"] == ptype]
    ax.scatter(sub["kmer_dist"], sub["embedding_cosine_distance"], s=8, alpha=0.3, color=color, label=label)
sp = summ["sp_kmer_embed_pooled"]
ax.set_xlabel("k-mer (4-6mer) frequency distance")
ax.set_ylabel("DNABERT-2 embedding cosine distance")
ax.set_title(f"A. Embedding distance tracks k-mer distance\nSpearman = {sp:.2f} (n = {len(pairs)} pairs)")
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.grid(alpha=0.25)

# Panel B: head-to-head R^2 predicting true identity
ax = axes[1]
labels = ["k-mer+GC+len\nalone", "embedding\nalone", "both\ncombined"]
vals = [summ["r2_ident_from_kmer_gc_len"], summ["r2_ident_from_embed"], summ["r2_ident_from_both"]]
colors = ["#7f7f7f", "#2ca02c", "#1f77b4"]
ax.bar(np.arange(3), vals, color=colors, width=0.55)
for i, v in enumerate(vals):
    ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
ax.set_xticks(np.arange(3), labels)
ax.set_ylabel("R$^2$ predicting true pairwise identity")
ax.set_ylim(0, 1)
ax.set_title(f"B. Embedding ≈ composition at recovering identity\n(+{summ['incremental_r2_embed_over_kmer']:.3f} R$^2$ from embedding on top of composition)")

# Panel C: region-level LOO-NN accuracy, k-mer vs embedding
ax = axes[2]
kmer_acc = [kmer_region[r]["loo_acc"] for r in order]
embed_acc = [phase2[r]["loo_acc"] for r in order]
x = np.arange(len(order))
w = 0.35
ax.bar(x - w/2, kmer_acc, width=w, color="#7f7f7f", label="k-mer NN (proxy B2)")
ax.bar(x + w/2, embed_acc, width=w, color="#2ca02c", label="DNABERT-2 embedding")
ax.axhline(0.5, color="grey", ls="--", lw=1)
ax.set_xticks(x, order)
ax.set_ylim(0, 1)
ax.set_ylabel("LOO nearest-neighbor accuracy")
ax.set_title("C. Region-level: embedding beats raw k-mer NN,\nbut ranking doesn't correlate (Spearman = %.2f)" % summ["sp_region_loo_kmer_vs_embed"])
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig("phase2b_summary.png", dpi=180, bbox_inches="tight")
print("saved phase2b_summary.png")
