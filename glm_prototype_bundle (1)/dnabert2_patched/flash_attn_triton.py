# Patched stub (4th CPU-compatibility patch, macOS/Apple Silicon specific).
#
# Original content: MosaicML's Triton implementation of Flash Attention
# (Copyright 2022 MosaicML Examples authors / Tri Dao, Apache-2.0).
#
# `triton` has no published wheel for macOS / Apple Silicon (pip install triton
# fails outright), so the real kernel in this file can never be imported, let
# alone run, on this platform. bert_layers.py already wraps its import of this
# module in `try: ... except ImportError: flash_attn_qkvpacked_func = None`,
# and every call site is additionally gated on `qkv.is_cuda` (see the 3rd
# documented patch in README_PATCHES.md) -- so this path is never taken on
# CPU regardless. The only problem was transformers' static `check_imports`
# pre-flight check, which recursively inspects this file's own top-level
# imports (ignoring try/except in the *caller*, bert_layers.py, but not able
# to tolerate an unconditional `import triton` *here*) and refuses to even
# cache the module if `triton` isn't importable -- independent of whether the
# code path that needs it would ever execute.
#
# Replacing the real kernel with a no-op ImportError reproduces exactly the
# "flash attention unavailable" state bert_layers.py already knows how to
# handle, without requiring triton to be importable at all. The 3 patches in
# README_PATCHES.md were made where `triton` was pip-installable (evidently a
# Linux sandbox) even though its CUDA kernels couldn't run there either, so
# the hard-dependency-at-check_imports-time problem never surfaced; it does
# on macOS/arm64.
raise ImportError("triton has no macOS/Apple Silicon wheel; flash-attention kernel unavailable on this platform")
