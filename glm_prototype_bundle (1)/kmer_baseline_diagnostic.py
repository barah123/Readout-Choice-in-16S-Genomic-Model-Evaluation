"""Phase 2b diagnostic: is the DNABERT-2 embedding separability signal just k-mer
composition in disguise?

Reuses the identical slicing/QA logic from slice_and_score.py / embed_and_score.py so
that every pair here lines up by (region, id1, id2) with the existing
pair_identity_vs_embedding.csv rows. Adds:
  - B2: k-mer (4,5,6-mer, pooled) frequency-vector distance per pair
  - B3: per-slice Shannon entropy (order-0) difference per pair
  - GC-content and length difference per pair (covariates flagged in phase2_report.md)

Then:
  1. Pair-level (N=2328): Spearman(kmer_dist, embedding_dist), Spearman(kmer_dist, identity),
     OLS embedding_dist ~ kmer_dist + GC_diff + len_diff -> R^2 and residual,
     partial correlation of that residual against true identity.
  2. Region-level (N=9): k-mer-based LOO-NN accuracy & silhouette, computed the same way
     Phase 2 computed the embedding-based ones, for a head-to-head region ranking check.
"""
import json
import csv
import itertools
import numpy as np
import pandas as pd
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

# ---- same QA as Phase 0-2: dedupe, near-full-length only ----
seqs = []  # (group, id, region_dict)
seen = set()
for grp, hdr, seq in records:
    if grp == "E_coli_ref" or seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    acc = hdr.split()[0]
    seqs.append((grp, acc, slice_regions(seq)))
print(f"{len(seqs)} sequences (matches Phase 2 QA count)", flush=True)
by_id = {acc: (grp, r) for grp, acc, r in seqs}

# ---- B2: pooled 4/5/6-mer frequency vectors per slice ----
def kmer_vector(seq, ks=(4, 5, 6)):
    vec = {}
    for k in ks:
        n = len(seq) - k + 1
        if n <= 0:
            continue
        counts = {}
        for i in range(n):
            kmer = seq[i:i + k]
            if "N" in kmer:
                continue
            counts[(k, kmer)] = counts.get((k, kmer), 0) + 1
        tot = sum(counts.values())
        if tot:
            for key, c in counts.items():
                vec[key] = c / tot  # normalize within each k separately, then concat
    return vec

def kmer_dist(va, vb):
    """Euclidean distance over the union of keys (sparse dict vectors)."""
    keys = set(va) | set(vb)
    return float(np.sqrt(sum((va.get(k, 0.0) - vb.get(k, 0.0)) ** 2 for k in keys)))

def shannon_entropy(seq):
    if not seq:
        return np.nan
    counts = {}
    for b in seq:
        counts[b] = counts.get(b, 0) + 1
    n = len(seq)
    ps = np.array([c / n for c in counts.values()])
    return float(-(ps * np.log2(ps)).sum())

def gc_content(seq):
    if not seq:
        return np.nan
    return sum(1 for b in seq if b in "GC") / len(seq)

# ---- load existing Phase 2 pair table and join k-mer/entropy/GC/len features ----
pairs = pd.read_csv("pair_identity_vs_embedding.csv")

kmer_cache = {}   # (region, id) -> kmer vector
scalar_cache = {} # (region, id) -> (entropy, gc, length)
def get_slice(region, acc):
    grp, r = by_id[acc]
    return r[region]

def get_kmer(region, acc):
    key = (region, acc)
    if key not in kmer_cache:
        s = get_slice(region, acc)
        kmer_cache[key] = kmer_vector(s) if s else {}
    return kmer_cache[key]

def get_scalars(region, acc):
    key = (region, acc)
    if key not in scalar_cache:
        s = get_slice(region, acc)
        scalar_cache[key] = (shannon_entropy(s) if s else np.nan,
                              gc_content(s) if s else np.nan,
                              len(s) if s else np.nan)
    return scalar_cache[key]

kmer_d, ent_diff, gc_diff, len_diff = [], [], [], []
for _, row in pairs.iterrows():
    region, id1, id2 = row["region"], row["id1"], row["id2"]
    v1, v2 = get_kmer(region, id1), get_kmer(region, id2)
    kmer_d.append(kmer_dist(v1, v2) if v1 and v2 else np.nan)
    e1, g1, l1 = get_scalars(region, id1)
    e2, g2, l2 = get_scalars(region, id2)
    ent_diff.append(abs(e1 - e2) if not (np.isnan(e1) or np.isnan(e2)) else np.nan)
    gc_diff.append(abs(g1 - g2) if not (np.isnan(g1) or np.isnan(g2)) else np.nan)
    len_diff.append(abs(l1 - l2) if not (np.isnan(l1) or np.isnan(l2)) else np.nan)

pairs["kmer_dist"] = kmer_d
pairs["entropy_diff"] = ent_diff
pairs["gc_diff"] = gc_diff
pairs["len_diff"] = len_diff
pairs = pairs.dropna(subset=["kmer_dist", "embedding_cosine_distance", "sequence_identity"]).copy()
pairs.to_csv("pair_identity_vs_embedding_vs_kmer.csv", index=False)
print(f"wrote pair_identity_vs_embedding_vs_kmer.csv ({len(pairs)} rows)", flush=True)

def spearman(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    def rank(a):
        order = np.argsort(a)
        r = np.empty(len(a))
        r[order] = np.arange(1, len(a) + 1)
        return r
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return float("nan")
    return float(np.corrcoef(rank(x), rank(y))[0, 1])

def ols(X, y):
    """X: (n,p) without intercept; adds intercept. Returns coef, R2, residuals, predicted."""
    Xd = np.column_stack([np.ones(len(X)), X])
    coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    pred = Xd @ coef
    resid = y - pred
    ss_res = (resid ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return coef, r2, resid, pred

# ======================================================================
# 1. PAIR-LEVEL DIAGNOSTIC (N=2328, robust)
# ======================================================================
print("\n=== Pair-level (N=%d) ===" % len(pairs))

sp_kmer_embed = spearman(pairs["kmer_dist"], pairs["embedding_cosine_distance"])
sp_kmer_ident = spearman(pairs["kmer_dist"], pairs["sequence_identity"])
sp_embed_ident = spearman(pairs["embedding_cosine_distance"], pairs["sequence_identity"])
print(f"Spearman(kmer_dist, embedding_dist)     pooled = {sp_kmer_embed:.3f}")
print(f"Spearman(kmer_dist, true identity)      pooled = {sp_kmer_ident:.3f}")
print(f"Spearman(embedding_dist, true identity) pooled = {sp_embed_ident:.3f}  (Phase 2 reported -0.80)")

per_region_pair = {}
for region, grp in pairs.groupby("region"):
    per_region_pair[region] = dict(
        n=len(grp),
        sp_kmer_embed=spearman(grp["kmer_dist"], grp["embedding_cosine_distance"]),
        sp_kmer_ident=spearman(grp["kmer_dist"], grp["sequence_identity"]),
        sp_embed_ident=spearman(grp["embedding_cosine_distance"], grp["sequence_identity"]),
    )
print("\nregion   n     Sp(kmer,embed)  Sp(kmer,ident)  Sp(embed,ident)")
for region in REGIONS:
    if region in per_region_pair:
        d = per_region_pair[region]
        print(f"{region:6} {d['n']:5d} {d['sp_kmer_embed']:14.3f} {d['sp_kmer_ident']:15.3f} {d['sp_embed_ident']:16.3f}")

# OLS: embedding_dist ~ kmer_dist + gc_diff + len_diff  -> how much of embedding distance
# is linearly explained by composition-based features alone?
X = pairs[["kmer_dist", "gc_diff", "len_diff"]].to_numpy(dtype=float)
y = pairs["embedding_cosine_distance"].to_numpy(dtype=float)
coef, r2, resid, pred = ols(X, y)
print(f"\nOLS embedding_dist ~ kmer_dist + gc_diff + len_diff: R^2 = {r2:.3f}")
print(f"  coefficients (intercept, kmer_dist, gc_diff, len_diff) = {np.round(coef, 4)}")

# Partial check: does the *residual* embedding distance (unexplained by composition)
# still track true identity? If yes, gLM carries identity-relevant signal beyond
# k-mer/GC/length. If the partial correlation is ~0, the embedding's identity-tracking
# is fully redundant with composition.
sp_resid_ident = spearman(resid, pairs["sequence_identity"])
print(f"Spearman(residual embedding_dist, true identity) = {sp_resid_ident:.3f}  "
      f"(vs raw Spearman(embed_dist, identity) = {sp_embed_ident:.3f})")

# Reverse check: how much does kmer_dist alone (+gc/len) explain of true identity,
# vs how much embedding_dist alone explains -- direct head-to-head on the ground truth.
Xk = pairs[["kmer_dist", "gc_diff", "len_diff"]].to_numpy(dtype=float)
yident = pairs["sequence_identity"].to_numpy(dtype=float)
_, r2_kmer_ident, _, _ = ols(Xk, yident)
Xe = pairs[["embedding_cosine_distance"]].to_numpy(dtype=float)
_, r2_embed_ident, _, _ = ols(Xe, yident)
Xboth = pairs[["kmer_dist", "gc_diff", "len_diff", "embedding_cosine_distance"]].to_numpy(dtype=float)
_, r2_both_ident, _, _ = ols(Xboth, yident)
print(f"\nR^2 predicting true identity from kmer+gc+len alone   = {r2_kmer_ident:.3f}")
print(f"R^2 predicting true identity from embedding_dist alone = {r2_embed_ident:.3f}")
print(f"R^2 predicting true identity from BOTH combined        = {r2_both_ident:.3f}  "
      f"(incremental R^2 from adding embedding on top of kmer+gc+len = {r2_both_ident - r2_kmer_ident:.3f})")

# ======================================================================
# 2. REGION-LEVEL DIAGNOSTIC (N=9, exploratory) -- k-mer-based LOO-NN/silhouette,
#    computed exactly like embed_and_score.py did for the embedding arm.
# ======================================================================
print("\n=== Region-level (N=9 regions) ===")
groups_arr = np.array([s[0] for s in seqs])
ids_arr = np.array([s[1] for s in seqs])
region_dicts = [r for _, _, r in seqs]

kmer_region_summary = {}
for name in REGIONS:
    sl = [r[name] for r in region_dicts]
    an = [i for i in range(len(seqs)) if groups_arr[i] == "B_anthracis" and sl[i]]
    ot = [i for i in range(len(seqs)) if groups_arr[i] in ("B_cereus", "B_thuringiensis", "B_weihenstephanensis") and sl[i]]
    if not an or not ot:
        continue
    idx = an + ot
    labels = [0] * len(an) + [1] * len(ot)
    vecs = [kmer_vector(sl[i]) for i in idx]
    n = len(idx)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = kmer_dist(vecs[i], vecs[j])
            D[i, j] = D[j, i] = d
    accs = []
    for i in range(n):
        d = D[i].copy()
        d[i] = np.inf
        nn = int(np.argmin(d))
        accs.append(1 if labels[nn] == labels[i] else 0)
    loo_acc = float(np.mean(accs))
    own, other = [], []
    for i in range(n):
        same = [D[i, j] for j in range(n) if labels[j] == labels[i] and j != i]
        diff = [D[i, j] for j in range(n) if labels[j] != labels[i]]
        own.append(np.mean(same)); other.append(np.mean(diff))
    own, other = np.array(own), np.array(other)
    silh = float(np.mean((other - own) / np.maximum(other, own)))
    kmer_region_summary[name] = dict(loo_acc=loo_acc, silhouette=silh)

json.dump(kmer_region_summary, open("phase2b_kmer_region_summary.json", "w"), indent=2)

phase2_summary = json.load(open("phase2_summary.json"))
print("region   kmer_LOO  embed_LOO  kmer_silh  embed_silh")
rows_for_sp = []
for name in REGIONS:
    if name in kmer_region_summary and name in phase2_summary:
        k = kmer_region_summary[name]; e = phase2_summary[name]
        print(f"{name:6} {k['loo_acc']:9.3f} {e['loo_acc']:10.3f} {k['silhouette']:10.3f} {e['silhouette']:11.3f}")
        rows_for_sp.append((k["loo_acc"], e["loo_acc"], k["silhouette"], e["silhouette"]))

rows_for_sp = np.array(rows_for_sp)
sp_loo = spearman(rows_for_sp[:, 0], rows_for_sp[:, 1])
sp_silh = spearman(rows_for_sp[:, 2], rows_for_sp[:, 3])
print(f"\nSpearman(kmer LOO-acc, embedding LOO-acc) across 9 regions  = {sp_loo:.3f}")
print(f"Spearman(kmer silhouette, embedding silhouette) across 9 regions = {sp_silh:.3f}")

with open("phase2b_kmer_diagnostic_summary.json", "w") as f:
    json.dump(dict(
        n_pairs=len(pairs),
        sp_kmer_embed_pooled=sp_kmer_embed,
        sp_kmer_ident_pooled=sp_kmer_ident,
        sp_embed_ident_pooled=sp_embed_ident,
        per_region_pair=per_region_pair,
        ols_embed_on_kmer_gc_len_r2=r2,
        ols_coef=coef.tolist(),
        sp_residual_embed_vs_ident=sp_resid_ident,
        r2_ident_from_kmer_gc_len=r2_kmer_ident,
        r2_ident_from_embed=r2_embed_ident,
        r2_ident_from_both=r2_both_ident,
        incremental_r2_embed_over_kmer=r2_both_ident - r2_kmer_ident,
        region_level_kmer_summary=kmer_region_summary,
        sp_region_loo_kmer_vs_embed=sp_loo,
        sp_region_silhouette_kmer_vs_embed=sp_silh,
    ), f, indent=2)
print("\nsaved phase2b_kmer_diagnostic_summary.json")
