"""Phase 1 env audit: real inference smoke test for DNABERT-S, patched via the same 4 fixes
as DNABERT-2 (see glm_prototype_bundle (1)/dnabert_s_local/, patches transferred from
dnabert2_patched/README_PATCHES.md - same MosaicBERT codebase, confirmed identical bert_padding.py).
"""
import time
import torch
from transformers import AutoTokenizer, AutoModel

MODEL_DIR = "../../glm_prototype_bundle (1)/dnabert_s_local"

t0 = time.time()
tok = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModel.from_pretrained(MODEL_DIR, trust_remote_code=True)
model.eval()
for _n, module in model.named_modules():
    if hasattr(module, "rebuild_alibi_tensor") and getattr(module, "alibi", "x") is None:
        module.rebuild_alibi_tensor(size=module._current_alibi_size)
load_s = time.time() - t0
print(f"load time: {load_s:.1f}s")

test_seqs = {
    "V1_like_29nt": "ACGTACGTACGTACGTACGTACGTACGTA"[:29],
    "V4_like_107nt": "".join(["ACGT"[i % 4] for i in range(107)]),
    "full_length_1500nt": "".join(["ACGT"[i % 4] for i in range(1500)]),
}

results = {}
for name, seq in test_seqs.items():
    enc = tok(seq, return_tensors="pt", truncation=True, max_length=256)
    with torch.no_grad():
        out1 = model(**enc)[0]
    with torch.no_grad():
        out2 = model(**enc)[0]
    finite = torch.isfinite(out1).all().item()
    deterministic = torch.allclose(out1, out2, atol=0)
    results[name] = dict(input_len=len(seq), token_len=enc["input_ids"].shape[1],
                          output_shape=tuple(out1.shape), finite=finite,
                          deterministic_repeat=deterministic)
    print(f"{name}: input_len={len(seq)} token_len={enc['input_ids'].shape[1]} "
          f"output_shape={tuple(out1.shape)} finite={finite} deterministic={deterministic}")

assert all(r["finite"] for r in results.values()), "non-finite values in DNABERT-S output"
assert all(r["deterministic_repeat"] for r in results.values()), "DNABERT-S not deterministic on repeated identical input"
print("\nPASS: DNABERT-S real inference smoke test (tokenize + forward pass + determinism + finite check)")
