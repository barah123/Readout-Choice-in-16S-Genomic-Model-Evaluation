"""Phase 1 env audit: real inference smoke test for Nucleotide Transformer v2 (not an import
test). Per CLAUDE.md: must use AutoModelForMaskedLM + trust_remote_code=True (not plain
AutoModel - that silently falls through to transformers' built-in ESM impl and fails on a
real weight-shape mismatch), hidden states via model.esm(...).
"""
import time
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

MODEL = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"

t0 = time.time()
tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True)
model.eval()
load_s = time.time() - t0
print(f"load time: {load_s:.1f}s")

test_seqs = {
    "V1_like_29nt": "ACGTACGTACGTACGTACGTACGTACGTA"[:29],
    "V4_like_107nt": "".join(["ACGT"[i % 4] for i in range(107)]),
    "full_length_1500nt": "".join(["ACGT"[i % 4] for i in range(1500)]),
}

results = {}
for name, seq in test_seqs.items():
    enc = tok(seq, return_tensors="pt")
    with torch.no_grad():
        out1 = model.esm(**enc)[0]
    with torch.no_grad():
        out2 = model.esm(**enc)[0]
    finite = torch.isfinite(out1).all().item()
    deterministic = torch.allclose(out1, out2, atol=0)
    results[name] = dict(input_len=len(seq), token_len=enc["input_ids"].shape[1],
                          output_shape=tuple(out1.shape), finite=finite,
                          deterministic_repeat=deterministic)
    print(f"{name}: input_len={len(seq)} token_len={enc['input_ids'].shape[1]} "
          f"output_shape={tuple(out1.shape)} finite={finite} deterministic={deterministic}")

assert all(r["finite"] for r in results.values()), "non-finite values in NTv2 output"
assert all(r["deterministic_repeat"] for r in results.values()), "NTv2 not deterministic on repeated identical input"
print("\nPASS: NTv2 real inference smoke test (tokenize + forward pass + determinism + finite check)")
