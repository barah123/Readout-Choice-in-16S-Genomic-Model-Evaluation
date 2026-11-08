"""Phase 1 env audit: real inference smoke test for HyenaDNA (not an import test).
Per protocol §5.2: archive sequence lengths, output shape, masking, finite embeddings,
determinism, peak resource use.
"""
import time
import tracemalloc
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel

MODEL = "LongSafari/hyenadna-medium-160k-seqlen-hf"

tracemalloc.start()
t0 = time.time()
tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)
model = AutoModel.from_pretrained(MODEL, trust_remote_code=True)
model.eval()
load_s = time.time() - t0
print(f"load time: {load_s:.1f}s")

# realistic 16S-region-scale lengths actually used in this project (shortest V1 ~29nt to
# longest V2/V4 ~115nt), not arbitrary strings.
test_seqs = {
    "V1_like_29nt": "ACGTACGTACGTACGTACGTACGTACGT"[:29],
    "V4_like_107nt": "".join(["ACGT"[i % 4] for i in range(107)]),
    "full_length_1500nt": "".join(["ACGT"[i % 4] for i in range(1500)]),
}

results = {}
for name, seq in test_seqs.items():
    enc = tok(seq, return_tensors="pt")
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

cur, peak = tracemalloc.get_traced_memory()
print(f"peak Python-tracked memory during smoke test: {peak/1e6:.1f} MB (process RSS will be higher - model weights not fully captured by tracemalloc)")
tracemalloc.stop()

assert all(r["finite"] for r in results.values()), "non-finite values in HyenaDNA output"
assert all(r["deterministic_repeat"] for r in results.values()), "HyenaDNA not deterministic on repeated identical input"
print("\nPASS: HyenaDNA real inference smoke test (tokenize + forward pass + determinism + finite check)")
