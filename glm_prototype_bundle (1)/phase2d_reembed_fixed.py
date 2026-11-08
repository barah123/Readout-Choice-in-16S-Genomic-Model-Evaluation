"""Phase 2d: re-embed the V1-V9 slices with a BATCH-INVARIANT pooling function.

Why
---
The Phase 2 cache (`dnabert2_embeddings.npz`) is contaminated: string-identical slices
received different embeddings in 6 of 9 regions. Proof (see phase2d_report.md §2): in V4 a
107-nt slice shared by 28 sequences has |emb| = 2.492214 in batches 0 and 2 but 2.531707 in
batch 1 - the value is constant within a batch and changes exactly at batch boundaries.
Because `embed_and_score.py` batched in `seqs` order, which is grouped by species, batch
membership correlates with class label (batch 1 was entirely B. thuringiensis). That is
label leakage: a classifier can separate groups using a padding artifact rather than biology.

The bug is in the mean-pooling mask:

    mask[:, 0]  = False   # drop [CLS]  - correct for every row
    mask[:, -1] = False   # drop [SEP]  - only correct for the LONGEST row in the batch

[SEP] sits at each row's own final real-token position, not at the last column of a padded
batch. For every shorter row the real [SEP] is wrongly kept in the mean while a pad column is
dropped instead, so the pooled vector depends on the batch's max token count - i.e. on batch
composition.

Fix: build the mask from the tokenizer's own attention_mask and clear each row's first and
own-last real token individually. The result is provably independent of batching, which this
script asserts rather than assumes (§ verification below).

Writes `dnabert2_embeddings_fixed.npz`. The original cache is left untouched.
"""
import json
import numpy as np
import torch
from Bio.Align import PairwiseAligner
from transformers import AutoModel, AutoTokenizer

torch.manual_seed(0)
np.random.seed(0)

MODEL_DIR = "dnabert2_local"  # patched local copy - NOT the bare HF repo id (see CLAUDE.md)

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


seqs, seen = [], set()
for grp, hdr, seq in records:
    if grp == "E_coli_ref" or seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    seqs.append((grp, hdr.split()[0], slice_regions(seq)))
print(f"{len(seqs)} sequences after QA", flush=True)

tok = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModel.from_pretrained(MODEL_DIR, trust_remote_code=True)
model.eval()
for _n, module in model.named_modules():
    if hasattr(module, "rebuild_alibi_tensor") and getattr(module, "alibi", "x") is None:
        module.rebuild_alibi_tensor(size=module._current_alibi_size)


def embed_batch(seqlist, batch_size=16):
    """Mean-pool last-layer token embeddings, excluding [CLS], each row's own [SEP], and pads.

    Unlike the Phase 2 version this uses the tokenizer's attention_mask and each row's own
    final real-token index, so the pooled vector does not depend on what else is in the batch.
    """
    outs = []
    for i in range(0, len(seqlist), batch_size):
        batch = seqlist[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=256)
        with torch.no_grad():
            out = model(**enc)[0]
        am = enc["attention_mask"].bool()
        mask = am.clone()
        lengths = am.sum(1)
        mask[:, 0] = False                                    # [CLS]
        rows = torch.arange(mask.size(0))
        mask[rows, lengths - 1] = False                       # each row's own [SEP]
        m = mask.unsqueeze(-1).float()
        outs.append(((out * m).sum(1) / m.sum(1).clamp(min=1)).numpy())
    return np.vstack(outs)


# Embed each UNIQUE slice once, then broadcast to every sequence carrying it. This makes
# the embedding a pure function of the slice *by construction* - with per-row embedding the
# residual float32 GEMM non-determinism across differing batch shapes still leaves identical
# strings differing by ~3e-6, which is harmless but needlessly non-reproducible.
embeddings = {}
for name in REGIONS:
    present = [r[name] if r[name] is not None else "" for _, _, r in seqs]
    ok = [i for i, s in enumerate(present) if s]
    E = np.full((len(seqs), 768), np.nan)
    if ok:
        uniq = sorted({present[i] for i in ok})
        U = embed_batch(uniq)
        lookup = {s: U[j] for j, s in enumerate(uniq)}
        for i in ok:
            E[i] = lookup[present[i]]
    embeddings[name] = E
    print(f"embedded {name}: {len(ok)}/{len(seqs)} slices "
          f"({len({present[i] for i in ok})} unique)", flush=True)

# ---------------------------------------------------------------
# VERIFICATION 1 - string-identical slices must now be bit-identical
# ---------------------------------------------------------------
print("\nVerification 1: identical slices -> identical embeddings")
worst = 0.0
for name in REGIONS:
    E = embeddings[name]
    byslice = {}
    for i, (_, _, r) in enumerate(seqs):
        if r[name]:
            byslice.setdefault(r[name], []).append(i)
    for s, idxs in byslice.items():
        if len(idxs) > 1:
            sub = E[idxs]
            worst = max(worst, float(np.abs(sub - sub[0]).max()))
    print(f"  {name}: max within-haplotype deviation = "
          f"{max([float(np.abs(E[v] - E[v][0]).max()) for v in byslice.values() if len(v) > 1], default=0.0):.3e}")
assert worst == 0.0, f"still batch-dependent: max deviation {worst:.3e}"
print(f"PASS - worst deviation across all regions = {worst:.3e} (exactly 0 by construction)")

# ---------------------------------------------------------------
# VERIFICATION 2 - explicit batch-invariance: re-embed one region with a
# different batch size and a shuffled order; results must match.
# ---------------------------------------------------------------
print("\nVerification 2: batch-size / ordering invariance (region V4)")
present = [r["V4"] for _, _, r in seqs if r["V4"]]
a = embed_batch(present, batch_size=16)
b = embed_batch(present, batch_size=7)
perm = np.random.default_rng(0).permutation(len(present))
c = embed_batch([present[i] for i in perm], batch_size=5)
c_unperm = np.empty_like(c)
c_unperm[perm] = c
d1 = float(np.abs(a - b).max())
d2 = float(np.abs(a - c_unperm).max())
print(f"  batch 16 vs batch 7           : max abs diff = {d1:.3e}")
print(f"  batch 16 vs shuffled batch 5  : max abs diff = {d2:.3e}")
# Tolerance 1e-4, not 0: this checks the raw per-row path, where differing batch shapes change
# float32 GEMM reduction order. That residual is ~1e-6 on vectors of norm ~2.5, i.e. five
# orders of magnitude below the ~0.7 L2 systematic offset the original masking bug produced.
# The saved cache is exactly deterministic anyway because unique slices are embedded once.
assert d1 < 1e-4 and d2 < 1e-4, "embedding still systematically depends on batching"
print("PASS - no systematic batch dependence (residual is float32 reduction-order noise)")

np.savez_compressed("dnabert2_embeddings_fixed.npz",
                    **{k: v for k, v in embeddings.items()},
                    ids=np.array([s[1] for s in seqs]),
                    groups=np.array([s[0] for s in seqs]))
print("\nsaved dnabert2_embeddings_fixed.npz (original cache left untouched)")

# ---------------------------------------------------------------
# How much did the bug move things?
# ---------------------------------------------------------------
old = np.load("dnabert2_embeddings.npz", allow_pickle=True)
print("\nOld vs new cache, per region (cosine distance between old and new vector):")
print(f"{'reg':5} {'mean':>10} {'max':>10}")
for name in REGIONS:
    O, N = old[name], embeddings[name]
    ok = ~np.isnan(O).any(1) & ~np.isnan(N).any(1)
    On = O[ok] / np.linalg.norm(O[ok], axis=1, keepdims=True)
    Nn = N[ok] / np.linalg.norm(N[ok], axis=1, keepdims=True)
    cd = 1 - (On * Nn).sum(1)
    print(f"{name:5} {cd.mean():10.3e} {cd.max():10.3e}")
