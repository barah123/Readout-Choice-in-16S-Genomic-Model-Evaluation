# DNABERT-2 custom code, patched for transformers 5.x / torch 2.15 / CPU

These are the trust_remote_code files from zhihan1996/DNABERT-2-117M (revision 7bce263b)
with four fixes applied so the model runs on CPU with current libraries:

1. bert_layers.py: `padding_idx=getattr(config, "pad_token_id", 0)` — transformers 5.x
   BertConfig no longer defines pad_token_id.
2. bert_layers.py: ALiBi tensor construction wrapped in try/except — transformers 5.x
   initializes weights under a meta-device context, which breaks the mixed CPU/meta
   tensor math; alibi is rebuilt after load (see embed_and_score.py,
   the `rebuild_alibi_tensor` loop).
3. bert_layers.py: attention fallback gated on `not qkv.is_cuda` — the custom forward
   otherwise calls the GPU-only triton flash-attention kernel unconditionally.
4. **flash_attn_triton.py: replaced with a stub that raises `ImportError` immediately,
   with no `import triton` anywhere in the file** (macOS / Apple Silicon only — added
   2026-11-06, verified on an M4 MacBook Air). `triton` has no published wheel for
   macOS/arm64 at all (`pip install triton` fails outright there), unlike on the
   original sandbox where patches 1–3 were made (triton was importable there, just
   unable to run its CUDA kernels — patch 3 already handled that). The problem isn't
   runtime — bert_layers.py already wraps its import of this module in
   `try: ... except ImportError: flash_attn_qkvpacked_func = None` and every call site
   is additionally gated on `qkv.is_cuda`, so the real kernel was never going to run on
   CPU regardless of platform. The problem is `transformers`' *static*
   `check_imports` pre-flight check (`dynamic_module_utils.py`), which recursively
   inspects this file's own top-level imports before any code actually executes, and
   refuses to even cache the module if `triton` isn't importable — independent of
   whether that code path would ever run. It does tolerate imports wrapped in
   `try/except` in the file doing the importing (bert_layers.py), but not an
   unconditional `import triton` inside the dependency file itself. Replacing the
   kernel with a no-op `ImportError` reproduces the exact "flash attention
   unavailable" state bert_layers.py already knows how to handle, without requiring
   triton to be importable at all. No numerical effect — this file's real content was
   never reachable on CPU in the first place.

No numerical changes to the model. **How to use (verified approach):** don't rely on
overwriting transformers' HF-hub module cache directly — `get_cached_module_file` can
re-fetch from the Hub and clobber local edits. Instead: `snapshot_download` the model
repo to a local directory (e.g. `dnabert2_local/`, real files not symlinks:
`local_dir_use_symlinks=False`), copy these four patched files over the downloaded
ones, then `AutoModel.from_pretrained("./dnabert2_local", trust_remote_code=True)` —
loading from an explicit local path means every "fetch" is just reading your already-
patched source. See `phase2_report.md` §1 for the original three-patch context and
`MASTER_PLAN.md` for when patch 4 was added.
