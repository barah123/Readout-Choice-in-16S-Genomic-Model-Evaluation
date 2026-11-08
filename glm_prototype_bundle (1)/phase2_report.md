# Phase 2 Report: DNABERT-2 Embeddings on the Same 16S Region Slices

**Date:** November 2026
**Scope:** Phase 2 of the analysis plan — first gLM arm. DNABERT-2-117M embeddings of the same V1–V9 slices used in Phase 0–1, scored for separability of *B. anthracis* vs the *B. cereus* group, and compared head-to-head with the identity baseline.

---

## 1. What was done

1. **Model.** DNABERT-2-117M (BPE-tokenized, 12 layers, 768-dim), run on CPU. Mean-pooled last-layer token embeddings, with [CLS]/[SEP]/[PAD] tokens excluded from pooling. Determinism check passed: identical sequences produce identical embeddings (cosine ≈ 1.0) regardless of batch padding.
2. **Data.** The identical 43-sequence QA'd set and V1–V9 slices from Phase 0–1 (V9 partially represented: 25/43 slices pass QA). 387 slices embedded.
3. **Scoring.** Per region: leave-one-out nearest-neighbor (LOO-NN) classification accuracy for *B. anthracis* vs cereus group; silhouette width; mean within- and between-group embedding cosine distance; and pair-level Spearman correlation between embedding distance and true sequence identity.
4. **Compatibility fixes (documented for reproducibility).** DNABERT-2's custom code required three patches to run under transformers 5.x / torch 2.15 / CPU: (i) `config.pad_token_id` fallback, (ii) deferred ALiBi tensor construction past the meta-device weight-init context, (iii) CPU fallback in the attention forward (the PyTorch fallback path already existed; the GPU triton branch was gated on `qkv.is_cuda`). No numerical changes to the model.

## 2. Methodological finding: Phase 2 caught a baseline bug

The gLM arm's identical-sequence sanity check exposed a bug in the Phase 0–1 identity function: it counted *aligned* positions rather than *matching* positions, so any equal-length gapless pair scored 1.0 regardless of SNPs. The function was fixed and all Phase 0–1 numbers recomputed (`prototype_report.md` carries an erratum). The qualitative conclusions survived: btwn_max = 1.0 in every region is now a genuine string-identity claim, and V6 — while still non-separating — turns out to carry the largest mean between-group divergence of any region, directionally consistent with Chakravorty's heterogeneity claim.

**Protocol recommendation:** make the identical-sequence embedding check a standing QA step in the two-arm design. The gLM arm validated the baseline arm — an argument for keeping both arms in the benchmark even if the gLM never beats the baseline.

## 3. Results

| Region | LOO-NN accuracy | Silhouette | Within-embed dist | Between-embed dist | Spearman(identity, embed dist) |
|---|---|---|---|---|---|
| V1 | 0.800 | −0.161 | 0.0000 | 0.1170 | −0.749 |
| V2 | 0.750 | −0.217 | 0.0043 | 0.0350 | −0.892 |
| V3 | 0.389 | −0.290 | 0.0132 | 0.0364 | −0.788 |
| V4 | 0.806 | −0.264 | 0.0012 | 0.0149 | −0.649 |
| V5 | 0.778 | −0.220 | 0.0105 | 0.0602 | −0.789 |
| V6 | 0.528 | −0.092 | 0.0505 | 0.0873 | −0.715 |
| V7 | 0.778 | −0.211 | 0.0097 | 0.0635 | −0.783 |
| V8 | 0.528 | −0.301 | 0.0000 | 0.0517 | −0.708 |
| V9 | 0.667 | −0.090 | 0.0566 | 0.1332 | −0.869 |

Figure: `phase2_summary.png`.

## 4. Findings

1. **The gLM captures substitution-level signal in short slices.** Embedding cosine distance tracks true sequence identity strongly: pooled Spearman −0.80 over all 2,328 pairs (per-region −0.65 to −0.89). This is the necessary validity check for the gLM arm — DNABERT-2 "sees" single SNPs in 30–110 nt slices (e.g., the two anthracis V6 haplotypes, one SNP apart, are clearly separated in embedding space).
2. **Region rankings agree across arms.** Regions ranked by mean between-group divergence (identity) vs mean between-group embedding distance correlate at Spearman 0.73 (panel C). The gLM arm does not contradict the baseline's region hierarchy.
3. **V6 does not separate — the gLM agrees with the baseline.** LOO-NN accuracy in V6 is 0.528, indistinguishable from chance, and every region's silhouette is negative: *B. anthracis* sits inside the cereus-group cloud in embedding space for all nine regions. The Phase 0–1 conclusion (expected-null case) is confirmed by the gLM arm.
4. **No region is a clean separator, but V1/V4 carry the most partial signal** (LOO-NN 0.80–0.81). Note this is a *statistical* separation only — exact string identity between at least one anthracis–cereus pair in every region (btwn_max = 1.0) means no region can support a perfect assay, which is the operationally relevant claim.

## 5. Interpretation for the study design

- The two-arm design is working as intended: the gLM arm both validated (identity tracking) and stress-tested (QA bug) the baseline. Keep both arms.
- On this *expected-null* target, the gLM and the alignment baseline agree — as they should, since some pairs are string-identical and no representation can separate identical inputs. The discriminating test of the gLM's value must come from a target where partial separation exists: **Phase 3 should prioritize *S. aureus* vs CoNS (published V1 claim) and *B. longum* subsp. *infantis* vs *longum* (the infant-gut target)**, where embedding geometry could plausibly beat single-metric identity.
- Two additions before Phase 3: the k-mer RF baseline (B2, currently missing) and a hold-out check that DNABERT-2's embeddings are not simply recovering GC content or slice length (regress embedding distance on those covariates).

## 6. Deliverables in this directory

- `embed_and_score.py` — embedding + scoring (includes the CPU compatibility fixes)
- `phase2_figure.py` — figure generation
- `dnabert2_embeddings.npz` — 387 slice embeddings (43 sequences × 9 regions, NaN where QA-dropped)
- `pair_identity_vs_embedding.csv` — 2,328 pair-level rows (corrected identity vs embedding distance)
- `phase2_summary.json` — per-region metrics
- `phase2_summary.png` — Figure 2

## 7. Sources

- DNABERT-2 model card and weights: https://huggingface.co/zhihan1996/DNABERT-2-117M
- Zhou et al. 2023, DNABERT-2 (arXiv:2306.15006): https://arxiv.org/abs/2306.15006
- Chakravorty et al. 2007: https://pmc.ncbi.nlm.nih.gov/articles/PMC2562909/
