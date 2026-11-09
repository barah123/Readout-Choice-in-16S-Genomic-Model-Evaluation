"""Phase 6 diagnostic B (EXPLORATORY): on the HARD blocks only, does a different
layer or pooling rule do better than the frozen last-layer mean pooling?

Phase 5's stated limitation was that only one layer and one pooling rule were tested, so
"R2 is the best accessible readout" was never established - only "R2 beats R0/R1". Diagnostic
A showed the hard blocks are where any real difference would live. This tests both at once,
on the 5 blocks where classical error exceeds 5%.

EXPLORATORY. Selected after seeing Phase 5. Cannot be promoted to a confirmatory claim.
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from transformers import AutoModel, AutoTokenizer

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
import run_readouts_full as rr

HARD = ["brucella_ochrobactrum", "veillonella", "francisella",
        "methylobacterium", "acetobacteraceae"]
MODEL_DIR = "../../glm_prototype_bundle (1)/dnabert2_local"
SEED = 0

rows = json.load(open("../phase4_locked_test/test_blocks.json"))
scores = json.load(open("../phase4_locked_test/locked_test_scores.json"))

tok = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModel.from_pretrained(MODEL_DIR, trust_remote_code=True)
model.eval()
for _n, mod in model.named_modules():
    if hasattr(mod, "rebuild_alibi_tensor") and getattr(mod, "alibi", "x") is None:
        mod.rebuild_alibi_tensor(size=mod._current_alibi_size)


def embed_variants(seqs, batch_size=16):
    """Return {variant_name: matrix} for several layer/pooling choices."""
    acc = defaultdict(list)
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=256)
        with torch.no_grad():
            out = model(**enc)
        # NOTE on scope: intermediate layers are NOT tested here. MosaicBERT runs its encoder
        # on UNPADDED sequences internally (output_all_encoded_layers returns tensors shaped
        # (total_real_tokens, hidden), not (batch, seq, hidden)), so recovering per-layer
        # padded representations needs its unpad/repad index bookkeeping. Reimplementing that
        # risks a silently wrong representation, which is worse than not running the check.
        # Final-layer pooling variants are tested instead; the layer question is recorded as
        # an explicit open limitation.
        h = out[0]
        am = enc["attention_mask"].bool()
        mask = am.clone()
        lengths = am.sum(1)
        mask[:, 0] = False
        r = torch.arange(mask.size(0))
        mask[r, (lengths - 1).clamp(min=0)] = False
        empty = mask.sum(1) == 0
        if empty.any():
            mask[empty] = am[empty]
        m = mask.unsqueeze(-1).float()

        acc["mean_frozen"].append(
            ((h * m).sum(1) / m.sum(1).clamp(min=1)).float().numpy())
        acc["cls_token"].append(h[:, 0].float().numpy())
        masked = h.masked_fill(~mask.unsqueeze(-1), -1e9)
        acc["max_pool"].append(masked.max(1).values.float().numpy())
        mean_v = (h * m).sum(1) / m.sum(1).clamp(min=1)
        acc["mean_plus_max"].append(
            torch.cat([mean_v, masked.max(1).values], dim=-1).float().numpy())
    return {k: np.vstack(v) for k, v in acc.items()}


def probe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=2000, class_weight="balanced",
                                            random_state=SEED))


print("=" * 100)
print("EXPLORATORY DIAGNOSTIC B - layer / pooling sensitivity on the 5 HARD blocks")
print("=" * 100)
print("Same isolate split and same region as the frozen analysis, so only the")
print("representation changes.\n")

by_block = defaultdict(lambda: defaultdict(dict))
for r in rows:
    by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

results = {}
for b in HARD:
    v = scores[b]
    region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
    d = by_block[b][region]
    isolates = sorted(d)
    y = np.array([d[i][1] for i in isolates])
    slices = [d[i][0] for i in isolates]

    uniq = sorted(set(slices))
    variants = embed_variants(uniq)
    lut = {k: {s: M[j] for j, s in enumerate(uniq)} for k, M in variants.items()}

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    out = defaultdict(list)
    for tr, te in skf.split(np.zeros(len(y)), y):
        ytr, yte = y[tr], y[te]
        if len(set(ytr)) < 2 or len(set(yte)) < 2:
            continue
        str_tr = [slices[i] for i in tr]
        str_te = [slices[i] for i in te]
        Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
        c = probe(); c.fit(Ktr, ytr)
        out["C_kmer"].append(rr.balanced_accuracy(yte, c.predict(Kte)))
        for name, L in lut.items():
            Xtr = np.vstack([L[s] for s in str_tr])
            Xte = np.vstack([L[s] for s in str_te])
            p = probe(); p.fit(Xtr, ytr)
            out[name].append(rr.balanced_accuracy(yte, p.predict(Xte)))
    results[b] = {k: float(np.mean(x)) for k, x in out.items()}
    results[b]["_region"] = region
    results[b]["_frozen_dnabert2"] = v["dnabert2:R2"]

variant_names = [k for k in results[HARD[0]] if not k.startswith("_") and k != "C_kmer"]
print(f"{'block':24}{'reg':5}{'C_kmer':>9}{'frozen':>9}" + "".join(f"{n:>14}" for n in variant_names))
for b in HARD:
    r = results[b]
    print(f"{b:24}{r['_region']:5}{r['C_kmer']:9.3f}{r['_frozen_dnabert2']:9.3f}"
          + "".join(f"{r[n]:>14.3f}" for n in variant_names))

print()
print("mean over the 5 hard blocks, as a difference from C_kmer:")
for n in ["_frozen_dnabert2"] + variant_names:
    d = np.mean([results[b][n] - results[b]["C_kmer"] for b in HARD])
    label = "frozen last-layer mean" if n == "_frozen_dnabert2" else n
    print(f"  {label:24} {d:+.4f}")

json.dump(results, open("diag_layers_results.json", "w"), indent=2, default=float)
print("\nwrote diag_layers_results.json")
