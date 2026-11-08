# gLM — 16S Region Benchmark

## What this project is

A benchmark study testing whether a pretrained genomic language model (gLM — DNABERT-2,
and later HyenaDNA / Nucleotide Transformer v2) can predict, from sequence embeddings alone,
which 16S rRNA hypervariable region (V1–V9) best discriminates a target bacterial organism
from its near neighbors — before anyone sequences anything. It is positioned as a follow-on
to the user's systematic review on breastfeeding vs. formula feeding and infant gut
microbiome method-dependence (16S-inferred vs. shotgun-measured function): this study tests
a different instance of the same "method choice shapes conclusions" thesis, applied to
region choice rather than profiling technology. Working title for the eventual R package:
`regioneer`.

Full protocol: **`gLM 16S Region Benchmark — Analysis Plan v0.1.md`** (project root). Read
that before making any methodological change — it carries the hypotheses (H1: gLM vs.
classical baselines; H2: embedding scores predict empirical taxonomic resolution), the
panel definitions, the statistical plan, and pre-registered kill-switch criteria (§9).

**Current status and next steps live in `MASTER_PLAN.md`, not here — read that first for
"what's done / what's next."** This file is oriented to structure and conventions.

## Directory layout

```
gLM/
├── CLAUDE.md                                    # this file
├── MASTER_PLAN.md                               # status, phase log, recommendations — read first
├── gLM 16S Region Benchmark — Analysis Plan v0.1.md   # the protocol (source of truth for design)
├── literature.pdf                               # Perplexity export of the planning conversation
│                                                 # that produced the analysis plan; background only
└── glm_prototype_bundle (1)/                    # all code, data, and reports produced so far
    ├── fetch_seqs.py              # NCBI retrieval (B. cereus group 16S records)
    ├── slice_and_score.py         # Phase 0-1: align to E. coli ref, slice V1-V9, B1 identity baseline
    ├── figure.py                  # Phase 0-1 figure
    ├── embed_and_score.py         # Phase 2: DNABERT-2 embeddings on the same slices + scoring
    ├── phase2_figure.py           # Phase 2 figure
    ├── kmer_baseline_diagnostic.py # Phase 2b: k-mer/GC/length composition diagnostic (reviewer-requested)
    ├── phase2b_figure.py          # Phase 2b figure
    ├── phase2c_b2_rf_baseline.py  # Phase 2c: real B2 (k-mer+RF, LOHO-CV, permutation null, BH correction) + B3/B4
    ├── phase2d_reembed_fixed.py   # Phase 2d: re-embeds V1-V9 with a batch-invariant pooling fix (see below)
    ├── phase2cd_figure.py         # Phase 2c+2d figure
    ├── fetch_staph.py             # Phase 3a: NCBI retrieval, S. aureus + 6 CoNS species (retmax=300 — see below)
    ├── haplotype_gate_staph.py    # Phase 3a: haplotype pre-check (MASTER_PLAN §4 item 2 gate), run BEFORE
    │                               # committing embed compute to any new target
    ├── phase3a_staph_baseline.py  # Phase 3a: full B1-B4 + gLM baseline on S. aureus vs CoNS, built with every
    │                               # Phase 2c/2d lesson applied from the start (LOHO-CV, per-arm permutation
    │                               # null, BH correction, embedding purity check)
    ├── phase3a_figure.py          # Phase 3a figure
    ├── dnabert2_local/            # patched local copy of zhihan1996/DNABERT-2-117M (snapshot_download +
    │                               # the 4 patches applied) — load DNABERT-2 from here, not the bare HF repo id
    ├── cereus_group_16s.json      # 83 fetched sequences, B. cereus group (cached — avoids re-hitting NCBI)
    ├── staph_group_16s.json       # 787 fetched sequences, S. aureus + CoNS (cached — avoids re-hitting NCBI)
    ├── ecoli_16s.txt              # E. coli J01695 coordinate reference
    ├── pair_identities.csv        # 2,328 pair-level identities (corrected, post-erratum)
    ├── pair_identity_vs_embedding.csv        # + embedding cosine distance per pair
    ├── pair_identity_vs_embedding_vs_kmer.csv # + k-mer/GC/length features per pair
    ├── region_margins.json, phase2_summary.json, phase2b_kmer_region_summary.json,
    │   phase2b_kmer_diagnostic_summary.json, phase2c_b2_rf_summary.json, phase3a_staph_summary.json
    │                                                  # per-region / per-pair metric dumps
    ├── dnabert2_embeddings.npz    # ORIGINAL Phase 2 cache (B. cereus group) — leaks the class label via a
    │                               # batch-composition bug (phase2cd_report.md §3). Do not use for scoring;
    │                               # kept as historical record.
    ├── dnabert2_embeddings_fixed.npz  # CORRECTED cache (Phase 2d), B. cereus group — use for embedding scoring
    ├── staph_dnabert2_embeddings.npz  # Phase 3a embedding cache, S. aureus/CoNS — built pure-function-of-
    │                               # sequence from the start (dedupe-then-broadcast), no erratum needed
    ├── prototype_report.md, phase2_report.md, phase2b_kmer_diagnostic_report.md,
    │   phase2cd_report.md, phase3a_report.md           # the five reports written so far, in order —
    │                                                    # phase2cd_report.md supersedes the region-level
    │                                                    # numbers in the two reports before it; phase3a_report.md
    │                                                    # is the first report with a real (non-null) result
    ├── prototype_summary.png, phase2_summary.png, phase2b_summary.png, phase2cd_summary.png,
    │   phase3a_summary.png                             # matching figures
    └── dnabert2_patched/           # DNABERT-2 trust_remote_code files with CPU/transformers-5.x
                                     # patches applied (no numerical changes) + README_PATCHES.md
```

Reports are numbered in the order they were written (`prototype_report.md` = Phase 0–1,
`phase2_report.md` = Phase 2, `phase2b_...` = Phase 2b) and each one's §4/§5 ("Interpretation
for the study design") is the authoritative record of what was decided and why — prefer
reading those over re-deriving conclusions from the raw CSV/JSON.

## Environment

**There is no project-local Python environment — don't assume `python3` on PATH has anything
installed.** The system Python (3.14, Homebrew) has only `pip`/`wheel`. Three conda envs are
relevant here, each with a specific, verified role — set up and smoke-tested 2026-11-06 on
this machine (MacBook Air, Apple M4, 10 cores, 16GB RAM, arm64):

- **`glm-benchmark`** (Python 3.11, Biopython 1.88, numpy/pandas/scipy/matplotlib/scikit-learn,
  torch 2.15.0, transformers 5.22.1, einops, accelerate) — the main project env. Use for
  everything: slicing, identity scoring, the k-mer/RF baseline, figures, **and** DNABERT-2 /
  HyenaDNA embedding. `conda run -n glm-benchmark python <script>.py`
  - **DNABERT-2** needs the patched local copy at `dnabert2_local/` (not the bare HF repo id —
    see below). **HyenaDNA** (`LongSafari/hyenadna-medium-160k-seqlen-hf`) loads directly from
    the hub with no patches, confirmed working on this env/transformers version.
- **`glm-ntv2`** (Python 3.11, transformers **4.34.1** pinned, torch 2.15.0, Biopython, einops,
  accelerate) — a **separate** env, used only for Nucleotide Transformer v2
  (`InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`). Load via
  `AutoModelForMaskedLM.from_pretrained(..., trust_remote_code=True)` (not plain `AutoModel` —
  this repo's `auto_map` has no custom-code entry for it, so plain `AutoModel` silently falls
  through to transformers' own *built-in*, now-architecturally-diverged ESM implementation and
  fails on a real weight-shape mismatch, not just a config error); get hidden states via
  `model.esm(...)`. Do not try to run this model in `glm-benchmark` — see "Why two environments
  for embedding models" below.
- **`qiime2-amplicon-2024.5`** (pre-existing on this machine, shared with other projects) —
  no longer needed for this project now that `glm-benchmark` exists, but was used for Phase 0,
  1, and 2b before `glm-benchmark` was created. Fine to keep using it for one-off CSV/JSON-only
  scripts if ever convenient; don't add project-specific packages to it (it's shared).

### Why two environments for embedding models

DNABERT-2 (revision `7bce263b`) and Nucleotide Transformer v2 are both older (~2023,
pre-dating transformers 5.x) `trust_remote_code` repos, and each broke against current
`transformers` (5.22.1) in a different way:

- DNABERT-2's breakage was shallow (config fallback, a meta-device init ordering issue, a
  CUDA-only kernel branch, and — macOS-specific — `triton` has no wheel here at all) — four
  small, targeted patches fixed it cleanly. See `dnabert2_patched/README_PATCHES.md` for all
  four, and the **"How to use"** section there for why patches must go into a local
  `snapshot_download`-ed copy (`dnabert2_local/`), not the HF module cache directly (the cache
  gets silently re-fetched/clobbered otherwise).
- NTv2's breakage went progressively deeper — missing `rope_theta`, then `is_decoder`, then
  `add_cross_attention` on its config object, then a `transformers.pytorch_utils` function
  removed upstream, and finally a failure *inside* `transformers`' own internal weight-tying
  finalization code (`all_tied_weights_keys`) that has nothing to do with this repo's code at
  all. That last one is the tell: patching forward means reimplementing part of current
  `transformers`' internals, which risks silently-wrong tied weights — a correctness risk, not
  just an inconvenience. The actual fix was to stop patching and pin transformers back to
  **4.34.1** (close to the `4.32.0.dev0` this repo's `config.json` was saved under) in its own
  env — the unmodified repo code then loads and runs with zero patches, determinism check
  cosine similarity = 1.0 exactly.

**Lesson for any future model (DNABERT-S, Evo, etc. if ever added):** try it in `glm-benchmark`
first. If the failures stay shallow (config field, import error) patch it like DNABERT-2. If a
traceback ever originates *inside* `transformers/modeling_utils.py` or similar framework
internals rather than the model repo's own files, stop patching — pin an older transformers in
a new dedicated env instead, the way `glm-ntv2` does.

## Conventions

- The analysis plan's own workflow spec (§8) calls for **Python for model inference only,
  R for all statistics/visualization/packaging**. The prototype phases (0, 1, 2, 2b) did not
  follow this — everything so far is Python (numpy/scipy/pandas/matplotlib), for speed during
  prototyping. This is a known, acknowledged deviation, not an oversight — flag it in
  `MASTER_PLAN.md` if continuing to ignore it, don't silently "fix" it mid-phase.
- Every phase so far reuses the **identical slicing/QA logic** (same V1–V9 coordinates, same
  alignment parameters, same 1,450–1,600 nt full-length filter, same per-region length QA
  band) copy-pasted across scripts rather than factored into a shared module. When adding a
  new script that needs the slices, copy the `slice_regions()` function from
  `embed_and_score.py` (the most recent, canonical version) rather than reinventing it —
  this is what keeps every pair-level CSV joinable by `(region, id1, id2)`.
- Identity = `2 × matches / (len_a + len_b)` computed from Biopython's global pairwise
  alignment, matches counted over aligned blocks only. **A previous version of this function
  counted aligned positions instead of matching positions and silently inflated all
  equal-length gapless identities to 1.0** — see the Phase 2 erratum in `prototype_report.md`
  §2. If you ever touch this function again, the identical-sequence sanity check (does a
  sequence compared to itself, or to a known-identical strain, return exactly 1.0 only when
  truly identical?) is mandatory before trusting any output.
- Dataset used in Phase 0–2d is *B. anthracis* vs. the *B. cereus* group — this is Panel A's
  **expected-null case** (anthracis is a monophyletic subclade of cereus; no 16S region
  separates them). Don't mistake a null result on this target for a pipeline bug; it's the
  documented, correct outcome. The actual test of whether the gLM arm adds value needs a
  target with real partial separation — see `MASTER_PLAN.md` for which panel is next.
- **Embeddings must be a pure function of the input sequence — assert this before trusting any
  embedding-based score.** The Phase 2 DNABERT-2 cache (`dnabert2_embeddings.npz`) mean-pooled
  with `mask[:, -1] = False` to drop `[SEP]`, which is only correct for the longest row in a
  padded batch; shorter rows kept their real `[SEP]` while dropping a pad position instead, so
  the pooled vector depended on **batch composition**. Because Phase 2 batched sequences in
  species-grouped order, this leaked the class label into embeddings that were supposed to be
  content-only — see the Phase 2d erratum in `phase2cd_report.md` §3. The fix
  (`phase2d_reembed_fixed.py`): build the pooling mask from the tokenizer's own
  `attention_mask` plus each row's own sequence length, and — simplest and most robust — embed
  each **unique** sequence/slice once and broadcast to every row sharing it, so the embedding
  is a pure function of the string by construction. Before scoring any new embedding cache
  (HyenaDNA, NTv2, or a re-run of DNABERT-2), assert that string-identical inputs produce
  (near-)identical vectors — this is the embedding-arm analogue of the identity-function check
  above, and would have caught this bug immediately.
- **Before scoring any classifier on a per-region panel, count unique slice haplotypes per
  class per region — not raw sequence count.** On the anthracis/cereus panel, 9 sequences per
  class collapsed to 1–3 unique haplotypes per region (V1 and V8: exactly 1, identical across
  the class boundary). This caps achievable accuracy well below 1.0 (compute and report the
  Bayes-optimal ceiling — see `phase2c_b2_rf_baseline.py`'s oracle calculation) and makes
  nearest-neighbor accuracy meaningless under heavy ties (LOO-NN in `embed_and_score.py` /
  `kmer_baseline_diagnostic.py` reports tie-break order, not separability, once >~70% of
  samples are tied at distance zero — three regions reported accuracy *above* their own Bayes
  ceiling, a literal impossibility; see `phase2cd_report.md` §2). Use grouped CV (leave-one-
  haplotype-out) and a label-permutation null instead of a plain LOO nearest-neighbor score.
  **Apply this gate *before* fetching the full panel or spending embed compute**, not after —
  `haplotype_gate_staph.py` is the template: slice, count haplotypes per class per region, and
  only proceed if the target clears it (Phase 3a's *S. aureus*/CoNS panel did; it's why that
  phase could show a real result instead of another null).
- **NCBI `esearch` `retmax` matters a lot for full-length-sequence yield — don't reuse
  `fetch_seqs.py`'s `retmax=40` as a template without reconsidering it.** A first attempt at
  fetching *S. aureus*/CoNS with `retmax=40` per species yielded only 5 unique *S. aureus*
  sequences after the standard 1450–1600nt QA filter — most NCBI records titled "16S ribosomal
  RNA" for these species are short clinical-typing amplicons (a few hundred nt), not
  full-length genes, so a small `retmax` draw is mostly filtered out. This looked like genuine
  scarcity until checked against NCBI's true total hit count (`esearch` with `retmax=0`, read
  `.esearchresult.count`): 5,534 for *S. aureus* alone. `retmax=300` fixed it (195 unique
  sequences survived QA). Check the true NCBI count for a candidate species **before**
  concluding a panel is too thin to pursue on retmax=40-scale sampling.
