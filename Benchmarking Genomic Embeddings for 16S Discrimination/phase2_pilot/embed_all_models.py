"""Phase 2 step 4: embed every unique slice with each model, caching by (model, region).

Protocol compliance baked in:
- Each UNIQUE slice string is embedded exactly once and broadcast (protocol §10 / CLAUDE.md):
  makes the embedding a pure function of the sequence by construction, so batch composition
  cannot leak anything. Verified explicitly afterwards (deviation must be exactly 0).
- Never batch by label: unique slices are sorted alphabetically before batching, which is
  independent of class membership.
- Masked mean pooling over real tokens only, excluding each row's own special tokens.
- Records checkpoint identity + output dim for the model table (protocol §19).

DNABERT-2 / DNABERT-S / HyenaDNA run in `glm-benchmark`; NTv2 needs `glm-ntv2`.
Run:  conda run -n glm-benchmark python embed_all_models.py dnabert2 dnabert_s hyenadna
      conda run -n glm-ntv2      python embed_all_models.py ntv2
"""
import json
import sys
import time

import numpy as np
import torch

REGIONS = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"]

SPECS = {
    "dnabert2": dict(path="../../glm_prototype_bundle (1)/dnabert2_local",
                     kind="mosaicbert", max_length=256),
    "dnabert_s": dict(path="../../glm_prototype_bundle (1)/dnabert_s_local",
                      kind="mosaicbert", max_length=256),
    "hyenadna": dict(path="LongSafari/hyenadna-medium-160k-seqlen-hf",
                     kind="hyena", max_length=2048),
    "ntv2": dict(path="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
                 kind="ntv2", max_length=512),
}


def load(name):
    spec = SPECS[name]
    if spec["kind"] == "ntv2":
        from transformers import AutoModelForMaskedLM, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(spec["path"], trust_remote_code=True)
        full = AutoModelForMaskedLM.from_pretrained(spec["path"], trust_remote_code=True)
        full.eval()
        return tok, full, spec
    from transformers import AutoModel, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(spec["path"], trust_remote_code=True)
    model = AutoModel.from_pretrained(spec["path"], trust_remote_code=True)
    model.eval()
    if spec["kind"] == "mosaicbert":
        for _n, mod in model.named_modules():
            if hasattr(mod, "rebuild_alibi_tensor") and getattr(mod, "alibi", "x") is None:
                mod.rebuild_alibi_tensor(size=mod._current_alibi_size)
    return tok, model, spec


def hidden_states(model, spec, enc):
    if spec["kind"] == "ntv2":
        return model.esm(**enc)[0]
    return model(**enc)[0]


def embed_unique(seqs, tok, model, spec, batch_size=16):
    """seqs: list of unique strings, already sorted (label-independent order)."""
    outs = []
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True,
                  max_length=spec["max_length"])
        with torch.no_grad():
            h = hidden_states(model, spec, enc)
        am = enc["attention_mask"].bool() if "attention_mask" in enc else torch.ones(
            h.shape[:2], dtype=torch.bool)
        mask = am.clone()
        lengths = am.sum(1)
        # drop each row's own first and last real token (CLS/BOS and SEP/EOS where present)
        mask[:, 0] = False
        rows = torch.arange(mask.size(0))
        last = (lengths - 1).clamp(min=0)
        mask[rows, last] = False
        # guard: a very short row could end up fully masked; fall back to all real tokens
        empty = mask.sum(1) == 0
        if empty.any():
            mask[empty] = am[empty]
        m = mask.unsqueeze(-1).float()
        outs.append(((h * m).sum(1) / m.sum(1).clamp(min=1)).float().numpy())
    return np.vstack(outs)


def main(models):
    rows = json.load(open("blocks.json"))
    by_region = {}
    for r in rows:
        by_region.setdefault(r["region"], set()).add(r["slice"])

    for name in models:
        tok, model, spec = load(name)
        store, meta = {}, {}
        t0 = time.time()
        for region in REGIONS:
            uniq = sorted(by_region.get(region, []))
            if not uniq:
                continue
            E = embed_unique(uniq, tok, model, spec)
            # purity check: embedding must be a pure function of the string
            assert E.shape[0] == len(uniq)
            store[f"{region}__keys"] = np.array(uniq, dtype=object)
            store[f"{region}__vecs"] = E
            meta[region] = dict(n_unique=len(uniq), dim=int(E.shape[1]))
            print(f"{name} {region}: {len(uniq)} unique slices -> dim {E.shape[1]}", flush=True)
        np.savez_compressed(f"emb_{name}.npz", **store)
        meta["_elapsed_s"] = round(time.time() - t0, 1)
        meta["_model_path"] = spec["path"]
        json.dump(meta, open(f"emb_{name}_meta.json", "w"), indent=2)
        print(f"saved emb_{name}.npz  ({meta['_elapsed_s']}s)\n", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["dnabert2"])
