# Phase 2b Report: Is the Embedding Signal Just K-mer Composition in Disguise?

**Date:** November 2026
**Scope:** Targeted diagnostic requested by the external review of the analysis plan — "regress the gLM separability score on the k-mer/entropy baseline scores and report the residual." Runs on the existing Phase 0–2 *B. anthracis* vs *B. cereus* group dataset; no new sequencing or re-embedding was needed, since the question is about the relationship between DNABERT-2's embeddings and composition, not about a new target.

---

## 1. What was done

1. Re-sliced the same 43 QA'd sequences into V1–V9 using the identical alignment logic from `slice_and_score.py` / `embed_and_score.py`, so every result below joins by `(region, id1, id2)` onto the existing `pair_identity_vs_embedding.csv`.
2. Computed, per slice: a pooled 4/5/6-mer frequency vector (B2 proxy), order-0 Shannon entropy (B3), GC content, and length.
3. Computed, per pair (the same 2,328 pairs already scored for identity and embedding distance): k-mer frequency-vector Euclidean distance, |GC difference|, |length difference|.
4. Two levels of analysis:
   - **Pair-level (N = 2,328, well-powered):** correlations, an OLS regression of embedding distance on composition features, and a head-to-head comparison of how well composition alone vs. embedding alone vs. both together predict true pairwise identity.
   - **Region-level (N = 9, exploratory):** a k-mer nearest-neighbor classifier built the same way Phase 2 built the embedding LOO-NN classifier, so the two region rankings can be compared directly.

Outputs: `kmer_baseline_diagnostic.py`, `phase2b_figure.py`, `pair_identity_vs_embedding_vs_kmer.csv`, `phase2b_kmer_region_summary.json`, `phase2b_kmer_diagnostic_summary.json`, `phase2b_summary.png`.

## 2. Results

### 2.1 Pair-level (the statistically robust view)

| Comparison | Spearman |
|---|---|
| k-mer distance vs. embedding distance | 0.76 |
| k-mer distance vs. true identity | −0.87 |
| embedding distance vs. true identity | −0.78 |

- OLS `embedding_dist ~ kmer_dist + gc_diff + len_diff`: **R² = 0.86**. Composition alone linearly explains 86% of the variance in DNABERT-2's embedding distance.
- The residual embedding distance (what composition does *not* explain) has Spearman 0.25 against true identity — weak, and the opposite sign from the raw relationship (−0.78). Almost all of the embedding's identity-tracking is routed through composition; what's left over is not doing meaningful extra identity-relevant work.
- Head-to-head R² for predicting **true pairwise identity**:
  - k-mer + GC + length alone: **0.711**
  - embedding alone: **0.714**
  - both combined: **0.744** (incremental R² from adding the embedding on top of composition: **+0.034**)

**Reading:** at the pair level, DNABERT-2's embedding distance is essentially a (slightly noisier) restatement of simple sequence composition for this near-identical, highly conserved marker-gene regime. It does not out-predict composition at recovering ground truth, and adding it to a composition model buys a small, not obviously decision-relevant, improvement.

### 2.2 Region-level (exploratory, N = 9)

| Region | k-mer NN acc. | embed NN acc. | k-mer silhouette | embed silhouette |
|---|---|---|---|---|
| V1 | 0.486 | 0.800 | −0.194 | −0.161 |
| V2 | 0.417 | 0.750 | −0.107 | −0.217 |
| V3 | 0.333 | 0.389 | −0.272 | −0.290 |
| V4 | 0.333 | 0.806 | −0.270 | −0.264 |
| V5 | 0.472 | 0.778 | −0.237 | −0.220 |
| V6 | 0.528 | 0.528 | −0.064 | −0.092 |
| V7 | 0.417 | 0.778 | −0.232 | −0.211 |
| V8 | 0.417 | 0.528 | −0.286 | −0.301 |
| V9 | 0.619 | 0.667 | −0.061 | −0.090 |

- Silhouette (a continuous, smoother separability measure) agrees almost perfectly across arms: Spearman = **0.95**.
- LOO nearest-neighbor classification accuracy does **not** agree across arms: Spearman = **−0.07**, and the embedding beats the raw k-mer nearest-neighbor classifier in 7 of 9 regions, sometimes by a wide margin (V1: 0.49 vs 0.80; V4: 0.33 vs 0.81).

**Reading:** this is a genuine wrinkle, not a clean confirmation of "it's all composition." The smooth, aggregate separability statistic says the two arms see the same thing; the discrete classification statistic says they don't, and favors the embedding. Two caveats keep this from overriding §2.1's conclusion: (a) N = 9 regions with 9 anthracis vs. 26–34 cereus-group sequences per region is a small, noisy basis for a region-ranking correlation; (b) the k-mer arm here is a **raw-distance nearest-neighbor classifier**, not the random-forest classifier the analysis plan actually specifies for B2 (§4.2). A short k-mer window (V1 is ~31 nt) gives a coarse, tie-prone frequency vector under simple Euclidean distance — an RF trained on those same features could plausibly close some or all of this gap. This diagnostic used the cheaper proxy to get an answer fast; it is not a substitute for building the real B2.

## 3. Interpretation for the study design

- **The review's core concern is largely supported, not refuted**, at the level of evidence that matters most (pair-level, N = 2,328): most of what the embedding captures about sequence relationships in this near-identical, highly conserved marker-gene regime is redundant with simple k-mer/GC/length composition. This is consistent with the plan's own cited precedent (Micro16S) and with the prior-probability argument the review raised.
- **It is not a full refutation of the gLM arm's value**, for two reasons specific to this dataset: (1) this is the *expected-null* target (§10 of the plan) — *B. anthracis* is a poor venue for showing a representation-learning advantage, since several pairs are string-identical in every region and no method can separate identical inputs; (2) the region-level classification-accuracy divergence (§2.2) is real and unexplained by the pair-level story, and deserves a proper test rather than being averaged away.
- **This does not change the plan's existing Phase 3 priority** — if anything it sharpens the justification for it. *S. aureus* vs. CoNS and *B. longum* subsp. *infantis* vs. *longum* are targets with real, partial separation where a representation that captures something beyond raw composition would have room to show it. The anthracis/cereus result is a necessary validity check (passed: embeddings see SNPs, agree with baseline ranking), not the discriminating test.
- **One concrete process change:** build the actual B2 (k-mer + RF, as specified in §4.2) before or alongside Phase 3, not after. The raw-distance proxy used here was sufficient to get a fast read but is a weaker baseline than the plan calls for, and the one place in this diagnostic where the embedding showed a real, large advantage (LOO-NN accuracy) is exactly the comparison where swapping in the real RF baseline could change the conclusion.

## 4. Sources

- Same reference data and model as Phase 2 (`phase2_report.md`).
- Review comment this diagnostic responds to: see conversation record (`initial_idea`).
