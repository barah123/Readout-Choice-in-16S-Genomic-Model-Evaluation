"""Phase 2: DNABERT-2 embeddings on the same V1-V9 slices; compare with identity baseline.

CPU inference with mean-pooled last-layer token embeddings (special tokens excluded).
"""
import json
import numpy as np
import torch
from Bio.Align import PairwiseAligner
from transformers import AutoModel, AutoTokenizer

torch.manual_seed(0)
np.random.seed(0)

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

# ---- collect sequences with same QA as Phase 0-1 ----
seqs = []  # (group, id, region_dict)
seen = set()
for grp, hdr, seq in records:
    if grp == "E_coli_ref" or seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    acc = hdr.split()[0]
    seqs.append((grp, acc, slice_regions(seq)))
print(f"{len(seqs)} sequences", flush=True)

# ---- load model ----
tok = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
model = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
model.eval()
for name, module in model.named_modules():
    if hasattr(module, "rebuild_alibi_tensor") and getattr(module, "alibi", "x") is None:
        module.rebuild_alibi_tensor(size=module._current_alibi_size)
special = tok("ACGT")["input_ids"]
print("tokenizer special-token check:", special[:3], "...", flush=True)

def embed_batch(seqlist, batch_size=16):
    """Mean-pooled last-layer embeddings, excluding first/last (special) tokens."""
    outs = []
    for i in range(0, len(seqlist), batch_size):
        batch = seqlist[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=256)
        with torch.no_grad():
            out = model(**enc)[0]  # (B, L, 768)
        mask = torch.ones_like(enc["input_ids"], dtype=torch.bool)
        mask[:, 0] = False
        mask[:, -1] = False
        mask[enc["input_ids"] == tok.pad_token_id] = False
        m = mask.unsqueeze(-1).float()
        emb = (out * m).sum(1) / m.sum(1).clamp(min=1)
        outs.append(emb.numpy())
    return np.vstack(outs)

# ---- embed each region ----
embeddings = {}   # region -> (N, 768) aligned with `seqs` order (NaN rows where slice missing)
for name in REGIONS:
    present = [r[name] if r[name] is not None else "" for _, _, r in seqs]
    ok = [i for i, s in enumerate(present) if s]
    if ok:
        E = np.full((len(seqs), 768), np.nan)
        E[ok] = embed_batch([present[i] for i in ok])
    else:
        E = np.full((len(seqs), 768), np.nan)
    embeddings[name] = E
    print(f"embedded {name}: {len(ok)}/{len(seqs)} slices", flush=True)

np.savez_compressed("dnabert2_embeddings.npz",
                    **{f"{k}": v for k, v in embeddings.items()},
                    ids=np.array([s[1] for s in seqs]),
                    groups=np.array([s[0] for s in seqs]))

# ---- pair-level join: embedding cosine distance vs sequence identity ----
groups = np.array([s[0] for s in seqs])
rows = []
for name in REGIONS:
    E = embeddings[name]
    sl = [r[name] for _, _, r in seqs]
    an_idx = [i for i in range(len(seqs)) if groups[i] == "B_anthracis" and sl[i]]
    ot_idx = [i for i in range(len(seqs)) if groups[i] in ("B_cereus", "B_thuringiensis", "B_weihenstephanensis") and sl[i]]
    En = E / np.linalg.norm(E, axis=1, keepdims=True)
    for pair_type, idx_a, idx_b in [("within", an_idx, an_idx), ("between", an_idx, ot_idx)]:
        for a in idx_a:
            for b in idx_b:
                if pair_type == "within" and b <= a:
                    continue
                cos_sim = float(En[a] @ En[b])
                ident = identity(sl[a], sl[b])
                rows.append((name, seqs[a][1], seqs[b][1], groups[b] if pair_type == "between" else "B_anthracis",
                             pair_type, ident, 1 - cos_sim))

import csv
with open("pair_identity_vs_embedding.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["region", "id1", "id2", "group2", "pair_type", "sequence_identity", "embedding_cosine_distance"])
    w.writerows(rows)
print(f"wrote {len(rows)} pair rows", flush=True)

# ---- per-region separability metrics ----
def spearman(x, y):
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    if n < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return float("nan")
    rx = stats_rank(x); ry = stats_rank(y)
    return float(np.corrcoef(rx, ry)[0, 1])

def stats_rank(x):
    x = np.asarray(x)
    order = np.argsort(x)
    ranks = np.empty(len(x))
    ranks[order] = np.arange(1, len(x) + 1)
    return ranks

print("\nregion  LOO_NN_acc  silh   within_embed  between_embed  btwn_min_ident  spearman(ident,embeddist)")
summary = {}
for name in REGIONS:
    E = embeddings[name]
    sl = [r[name] for _, _, r in seqs]
    an = [i for i in range(len(seqs)) if groups[i] == "B_anthracis" and sl[i]]
    ot = [i for i in range(len(seqs)) if groups[i] in ("B_cereus", "B_thuringiensis", "B_weihenstephanensis") and sl[i]]
    if not an or not ot:
        continue
    idx = an + ot
    labels = [0] * len(an) + [1] * len(ot)
    En = E[idx] / np.linalg.norm(E[idx], axis=1, keepdims=True)
    S = En @ En.T
    # LOO nearest-neighbor: label of nearest neighbor (exclude self)
    accs = []
    for i in range(len(idx)):
        sims = S[i].copy()
        sims[i] = -2
        nn = int(np.argmax(sims))
        pred = labels[nn]
        # nearest neighbor of the *other* class:
        accs.append(1 if pred == labels[i] else 0)
    loo_acc = float(np.mean(accs))
    # silhouette
    dist = 1 - S
    own = []
    other = []
    for i in range(len(idx)):
        same = [dist[i, j] for j in range(len(idx)) if labels[j] == labels[i] and j != i]
        diff = [dist[i, j] for j in range(len(idx)) if labels[j] != labels[i]]
        own.append(np.mean(same)); other.append(np.mean(diff))
    own = np.array(own); other = np.array(other)
    silh = float(np.mean((other - own) / np.maximum(other, own)))
    # embedding distances (global indices)
    En_full = E / np.linalg.norm(E, axis=1, keepdims=True)
    within_d = [1 - float(En_full[a] @ En_full[b]) for ai, a in enumerate(an) for b in an[ai+1:]]
    between_d = [1 - float(En_full[a] @ En_full[b]) for a in an for b in ot]
    # identity from pair rows
    pr = [r for r in rows if r[0] == name]
    btwn_min_ident = min(r[5] for r in pr if r[4] == "between")
    sp = spearman([r[5] for r in pr], [r[6] for r in pr])
    summary[name] = dict(loo_acc=loo_acc, silhouette=silh,
                         within_embed=float(np.mean(within_d)),
                         between_embed=float(np.mean(between_d)),
                         btwn_min_ident=btwn_min_ident, spearman_ident_embed=sp)
    print(f"{name:6} {loo_acc:8.3f} {silh:7.3f} {np.mean(within_d):11.4f} {np.mean(between_d):13.4f} {btwn_min_ident:13.4f} {sp:9.3f}")

# overall pooled correlation
sp_all = spearman([r[5] for r in rows], [r[6] for r in rows])
print(f"\npooled Spearman(identity, embedding-distance) over all {len(rows)} pairs: {sp_all:.3f}")

json.dump(summary, open("phase2_summary.json", "w"), indent=2)
print("saved phase2_summary.json")
