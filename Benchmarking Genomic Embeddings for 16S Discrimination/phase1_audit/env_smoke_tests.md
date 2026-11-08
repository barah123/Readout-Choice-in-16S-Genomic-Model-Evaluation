# Environment smoke tests — re-verified this session, protocol §5.2's last checklist item

**Why re-run these:** `CLAUDE.md`/`MASTER_PLAN.md` already described HyenaDNA as "confirmed
working" and NTv2 as passing "a determinism check (cosine similarity = 1.0)" from an earlier
session. Per the audit principle, a prior session's self-description is a claim to verify, not
a fact to inherit — especially the distinction protocol §5.2 draws explicitly: *"Passing an
import test is not model readiness."* All three models below were re-run fresh this session
with real tokenize→forward-pass→determinism→finite-value checks, not import checks.

## Scope of what "smoke test" means here

Per §5.2: sequence lengths spanning what this project actually uses (shortest region ~29nt,
a mid-length region ~107nt, a full-length 16S gene ~1500nt), output shape, finite values,
determinism (identical input twice → bit-identical output), and load time. Not yet covered:
peak RSS/GPU memory under a full embedding sweep (deferred to Phase 2's pilot, where it's
actually exercised at scale) and masking-behavior parity across models (deferred to Phase 2
once a common embedding step exists for all models).

## Results

| model | env | load time | lengths tested | output shape (768/1024/256-dim) | finite | deterministic |
|---|---|---|---|---|---|---|
| DNABERT-2 | `glm-benchmark` | (prior sessions — real embedding sweeps already run at scale on both panels) | 29–1500nt | (1, L, 768) | yes | yes (purity-checked, §6 of audit_report.md) |
| HyenaDNA | `glm-benchmark` | 5.4s | 28/107/1500nt | (1, L+1, 256) | yes | yes |
| NTv2 | `glm-ntv2` | 5.0s | 29/107/1500nt | (1, ~L/6, 1024) | yes | yes |
| DNABERT-S | `glm-benchmark` (new: `dnabert_s_local/`) | 1.3s | 29/107/1500nt (truncated to 256 tokens) | (1, L, 768) | yes | yes |

Scripts: `hyenadna_smoke_test.py`, `ntv2_smoke_test.py`, `dnabert_s_smoke_test.py` (this
directory). DNABERT-2 was not re-run here since its determinism/purity was already verified at
full scale across both real panels (Phase 2d, Phase 3a) — a toy smoke test would be strictly
weaker evidence than what's already archived.

**Note on HyenaDNA's token count:** `token_len = input_len + 1` for every test case — consistent
with a single BOS/CLS-style token prepended, not an error. **Note on NTv2's token count:**
roughly `input_len / 6` — consistent with its published 6-mer tokenizer, not truncation or
information loss within the tested lengths.

## New this session: DNABERT-S is runnable, via the same 4 patches as DNABERT-2

The protocol's critical review (§2.2, §9) flags DNABERT-S as "essential prior art" and a
required Phase 2 model — it is the one model in the core four actually trained with a
species-aware objective, unlike DNABERT-2/HyenaDNA/NTv2's generic pretraining. It was not
previously in either conda environment. Checked and set up this session:

- `zhihan1996/DNABERT-S`'s `trust_remote_code` files are the same MosaicBERT codebase as
  DNABERT-2 (`bert_padding.py` is byte-identical; `bert_layers.py`/`configuration_bert.py`/
  `flash_attn_triton.py` differ only in exactly the places the existing 4 documented patches
  touch — confirmed via direct diff before writing anything).
- Applied the same 4 patches surgically to DNABERT-S's own files (not copied wholesale — its
  `configuration_bert.py` inherits from the real `BertConfig` rather than `PretrainedConfig`,
  and it carries a `BertForNextSentencePrediction` stub DNABERT-2's patched copy doesn't; left
  both alone rather than risk silently dropping something DNABERT-S-specific).
- Loaded cleanly on the first attempt with all 4 patches applied — no new errors encountered,
  consistent with CLAUDE.md's "try it in glm-benchmark first, patch forward if shallow" guidance.
- `config.json`'s `_name_or_path: "zhihan1996/DNABERT-2-117M"` confirms DNABERT-S is literally
  initialized from the same DNABERT-2 checkpoint and contrastively fine-tuned from there, which
  is *why* the same patches transfer — not a coincidence to rely on for any future model.
- License: Apache-2.0 (from the repo's README frontmatter; no separate `LICENSE` file bundled,
  unlike DNABERT-2-117M's repo which does include one — worth a second check directly on the
  HF repo page before any redistribution, not just this snapshot's metadata).
- New local copy at `glm_prototype_bundle (1)/dnabert_s_local/` (git-ignored: `pytorch_model.bin`
  only, same as `dnabert2_local/` — patched `.py` files are tracked).

**Verdict: DNABERT-S is ready for Phase 2.** No new environment, no new patch category needed.
