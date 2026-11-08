# Phase 2c + 2d Report — the real B2 baseline, and two errata that invalidate the Phase 2/2b region-level results

**Date:** 2026-11-06
**Scripts:** `phase2c_b2_rf_baseline.py`, `phase2d_reembed_fixed.py`
**Outputs:** `phase2c_b2_rf_summary.json`, `dnabert2_embeddings_fixed.npz`
**Supersedes:** the region-level numbers in `phase2_report.md` and `phase2b_kmer_diagnostic_report.md`

---

## 1. What this phase set out to do, and what it found instead

`MASTER_PLAN.md` §4 item 1 called for replacing the cheap raw-distance k-mer
nearest-neighbour proxy used in Phase 2b with the random-forest classifier the analysis plan
§4.2 actually specifies (B2), and for promoting B3/B4 from diagnostic covariates to formal
baselines. §5 flagged the proxy-vs-real-classifier gap as "the single most load-bearing
unresolved question right now."

Building that baseline required looking at the per-region class structure, and that exposed
**two independent defects that invalidate every region-level number in Phases 2 and 2b** —
including the one the open question was about. The open question is now answered, but the
answer is that it was malformed: Phase 2b's region-level comparison was between two
artifacts, not two feature spaces.

Both defects are documented below with the evidence that establishes them, followed by the
corrected results.

---

## 2. Erratum A — the region-level LOO-NN metric was measuring tie-break order

**The structural fact.** After the standard QA (43 unique near-full-length sequences), the
target class does not have 9 independent members. It has **1–3 unique slice haplotypes per
region**:

| region | n_anthracis | unique anthracis haplotypes | n_cereus-group | unique cg haplotypes | haplotypes shared across the class boundary |
|---|---|---|---|---|---|
| V1 | 9 | **1** | 26 | 9 | 1 |
| V2 | 9 | 3 | 27 | 16 | 2 |
| V3 | 9 | 2 | 27 | 8 | 1 |
| V4 | 9 | 2 | 27 | 8 | 1 |
| V5 | 9 | 2 | 27 | 11 | 1 |
| V6 | 9 | 2 | 27 | 13 | 2 |
| V7 | 9 | 2 | 27 | 10 | 2 |
| V8 | 9 | **1** | 27 | 9 | 1 |
| V9 | 6 | 2 | 15 | 8 | 1 |

In V1 and V8 all nine *B. anthracis* sequences are **one identical string**. And every region
contains at least one haplotype that occurs in *both* classes, so the label is not a function
of the slice.

**The consequence.** 72–90% of samples therefore have **tied** nearest neighbours at distance
zero, and `np.argmin`/`np.argmax` resolves those ties by index order. Because the sequence
list is ordered by species, that index order correlates with the class label.

**The proof that this is fatal, not cosmetic.** Three regions report an accuracy *above the
Bayes-optimal accuracy* — the best any classifier could achieve given that identical strings
carry both labels:

| region | Phase 2 reported LOO-NN | Bayes-optimal accuracy | verdict |
|---|---|---|---|
| V1 | 0.800 | 0.743 | **impossible** |
| V4 | 0.806 | 0.778 | **impossible** |
| V7 | 0.778 | 0.750 | **impossible** |

The metric is also wildly unstable to tie policy. Switching only the tie rule — first-index
vs majority-vote, same data, same distances — moves V4 from 0.361 to 0.750, V7 from 0.361 to
0.722, and V3 from 0.278 to 0.694. The *region ranking*, which is the study's entire output,
reorders depending on an arbitrary implementation detail.

Recomputed tie-aware (expected balanced accuracy averaged over uniform random tie-breaking),
every region falls to 0.45–0.62 in both arms.

## 3. Erratum B — the Phase 2 embedding cache leaks the class label through batch composition

**The observation.** String-identical slices received **different embeddings** in 6 of 9
regions (max within-haplotype L2 deviation up to 0.71).

**The mechanism.** In `embed_and_score.py` the mean-pooling mask does:

```python
mask[:, 0]  = False   # drop [CLS] — correct for every row
mask[:, -1] = False   # drop [SEP] — only correct for the LONGEST row in the batch
```

`[SEP]` sits at each row's own final real-token position, not at the last column of a padded
batch. For every shorter row the real `[SEP]` is kept in the mean while a pad column is
dropped instead, so the pooled vector depends on the batch's max token count — i.e. on
**batch composition**.

**The proof.** In V4, a 107-nt slice shared by 28 sequences has `|emb|` **exactly constant
within each batch and changing exactly at batch boundaries**:

| batch (size 16, in input order) | members | `|emb|` |
|---|---|---|
| 0 | 8 anthracis + 7 cereus | 2.492214 |
| 1 | **all 10 *B. thuringiensis*** | 2.531707 |
| 2 | 3 weihenstephanensis | 2.492214 |

`embed_and_score.py` batched in `seqs` order, which is grouped by species — so batch
membership correlates with class label. A classifier can separate *B. thuringiensis* from the
rest using a **padding artifact**, not biology. This is the direct cause of the
above-the-ceiling scores in Erratum A's table.

**The fix** (`phase2d_reembed_fixed.py`): build the mask from the tokenizer's own
`attention_mask` and clear each row's first and own-last real token individually; then embed
each *unique* slice once and broadcast to the sequences carrying it, making the embedding a
pure function of the slice **by construction**. Verified: max within-haplotype deviation is
now exactly `0.000e+00` in all 9 regions, and re-embedding at batch sizes 16 / 7 / shuffled-5
is bit-identical. Written to `dnabert2_embeddings_fixed.npz`; the original cache is left in
place but **must not be used for scoring**.

---

## 4. Corrected results

### 4.1 Pair-level (N = 2,328) — Phase 2b's conclusions survive intact

This arm was never affected by the tie-breaking defect, and is only mildly affected by the
embedding defect. Re-running it on the corrected cache reproduces Phase 2b's published values
on the old cache (0.862 / 0.711 / 0.714 / +0.28 — confirming the comparison is like-for-like):

| quantity | old cache | **corrected cache** | Phase 2b reported |
|---|---|---|---|
| Spearman(embedding distance, true identity) | −0.771 | **−0.877** | −0.80 |
| R² embedding distance ~ k-mer + GC + length | 0.862 | **0.862** | 0.86 |
| R² identity ~ composition alone | 0.711 | **0.711** | 0.711 |
| R² identity ~ embedding alone | 0.714 | **0.707** | 0.714 |
| incremental R² from adding embedding | +0.034 | **+0.031** | +0.034 |
| Spearman(residual embedding distance, identity) | +0.283 | **+0.207** | +0.25 |

The embedding's raw identity-tracking actually *improves* once the batch artifact is removed
(−0.771 → −0.877): the artifact was adding noise. But **Phase 2b's verdict is unchanged** —
composition explains 86% of embedding distance, and the embedding adds only +0.031 R² beyond
composition in predicting true identity.

### 4.2 Region-level, leakage-controlled — the real B2, and a reversal

Random forest (500 trees, balanced class weights) on k-mer (4/5/6-mer) features vs on
DNABERT-2 embeddings, under **identical leave-one-haplotype-out CV** so that string-identical
slices never straddle a train/test split, scored by balanced accuracy against a 200-draw
label-permutation null:

| region | anthracis haplotypes | **RF k-mer** | **RF embedding** | p(k-mer) | p(embed) | balanced-accuracy ceiling |
|---|---|---|---|---|---|---|
| V1 | 1 | *not estimable* | *not estimable* | — | — | 0.712 |
| V2 | 3 | **0.722** | 0.481 | 0.015 | 0.592 | 0.815 |
| V3 | 2 | 0.537 | 0.481 | 0.232 | 0.584 | 0.630 |
| V4 | 2 | 0.519 | 0.500 | 0.455 | 0.294 | 0.630 |
| V5 | 2 | **0.685** | 0.500 | 0.020 | 0.439 | 0.685 |
| V6 | 2 | 0.444 | 0.500 | 0.687 | 0.473 | 0.778 |
| V7 | 2 | 0.463 | 0.500 | 0.616 | 0.470 | 0.667 |
| V8 | 1 | *not estimable* | *not estimable* | — | — | 0.648 |
| V9 | 2 | 0.433 | 0.500 | 0.657 | 0.414 | 0.800 |

**After Benjamini-Hochberg across all 14 (region × arm) cells (plan §6): nothing survives**
(smallest q = 0.143).

Three things changed relative to Phase 2b, all in the same direction:

1. **The embedding arm shows no separability whatever** — 0.481–0.500 in every estimable
   region, every p > 0.29. On the *contaminated* cache the same analysis gave V2 = 0.778
   (p = 0.010) and V4 = 0.759 (p = 0.005); both vanish entirely once the batch leak is closed.
   Those two "significant" cells were pure artifact.
2. **The direction of the gLM-vs-baseline comparison reverses.** Phase 2b reported the
   embedding beating raw k-mer nearest-neighbour in 7 of 9 regions (V1 0.49 vs 0.80; V4 0.33
   vs 0.81). Leakage-controlled, the k-mer RF is the *better* arm wherever there is any
   nominal signal at all (V2, V5), and the gLM arm is flat.
3. **V1 and V8 are not estimable** — with one anthracis haplotype, holding it out removes the
   class from training. Reporting a number there would be fabrication; Phase 2/2b reported
   0.800 and 0.528.

### 4.3 B3 and B4 as formal baselines

| region | B3 per-site Shannon entropy (all) | B3 (anthracis only) | B4 mean length | B4 variable-site density |
|---|---|---|---|---|
| V1 | 0.603 | 0.000 | 28.9 | 1.000 |
| V2 | 0.183 | 0.009 | 115.0 | 0.557 |
| V3 | 0.209 | 0.008 | 66.0 | 0.708 |
| V4 | 0.148 | 0.005 | 107.0 | 0.679 |
| V5 | 0.589 | 0.009 | 59.4 | 0.863 |
| V6 | 0.295 | 0.016 | 58.4 | 0.655 |
| V7 | 0.309 | 0.009 | 56.1 | 0.893 |
| V8 | 0.281 | 0.000 | 52.0 | 0.647 |
| V9 | 0.313 | 0.084 | 31.5 | 0.750 |

Within-anthracis entropy is ~0 in every region except V9 — the quantitative statement of the
clonality that drives the haplotype collapse in §2.

---

## 5. Interpretation for the study design

**5.1 The headline is now a clean, defensible negative on this panel.** With leakage
controlled and multiple testing corrected, no region and no arm separates *B. anthracis* from
the *B. cereus* group. That is the *correct* result for this target — anthracis is a
monophyletic subclade of cereus *sensu lato* — and Phase 1 already reclassified it as the
expected-null case. What Phase 2c adds is that the null is now **quantified and bounded**: the
balanced-accuracy ceiling is 0.63–0.82, not 1.0, because identical strings carry both labels.
The expected-null call is no longer an assertion, it is a measured ceiling.

**5.2 The gLM-vs-baseline question remains open, but the Phase 2b evidence for the gLM side
is withdrawn.** MASTER_PLAN §5 warned not to cite the Phase 2b region-level numbers without
rerunning them against the real B2. Rerun, they reverse: the gLM arm is flat and the k-mer RF
is nominally better. Combined with the pair-level result (§4.1: +0.031 incremental R²), both
arms of the Phase 2b diagnostic now point the *same* way — toward the reviewer's concern and
the Micro16S precedent that the analysis plan §1.2 already pre-registers as a reportable
outcome. This is not yet a kill-switch trigger (§9 concerns Panel A groups that B1 separates;
B1 separates none here), but it removes the one piece of evidence that was pointing the other
way.

**5.3 The binding constraint on this design is unique haplotypes per region per class, not
sequence count.** 43 QA'd sequences sounds adequate; 1–3 effective observations in the target
class is not. MASTER_PLAN §4.2 already flagged checking per-group sequence counts before
committing compute — this sharpens it into a hard gate, and the unit must be **unique slice
haplotypes within each class, per region**, counted *after* slicing. A region where either
class has <2 haplotypes cannot support any cross-validated classifier and should be reported
as not estimable rather than scored. Suggested gate: ≥5 unique haplotypes per class per
region for a scorable cell, and report the cross-class-ambiguous fraction alongside every
score.

**5.4 Two protocol changes adopted, both generalising existing conventions.**
- *Embeddings must be a pure function of the sequence.* Any future embedding step (HyenaDNA,
  NTv2, DNABERT-S) must assert that string-identical inputs produce bit-identical vectors
  before its output is scored. This is the embedding-arm analogue of the existing
  identical-sequence QA check that caught the Phase 0–1 identity bug, and it would have caught
  this one immediately. The dedupe-then-broadcast pattern in `phase2d_reembed_fixed.py`
  satisfies it by construction and is also cheaper.
- *Nearest-neighbour metrics are banned on panels with duplicate slices.* Where ties are this
  prevalent, LOO-NN reports tie-break order. Use a grouped-CV classifier with an explicit
  permutation null, and always report the Bayes ceiling next to the score — the ceiling is
  what caught both bugs.

**5.5 What this does not change.** The pair-level arm (N = 2,328) is unaffected and its Phase
2b conclusions stand (§4.1). The Phase 0–1 identity baseline is unaffected. The environment
work is unaffected. The recommendation to move to a partial-separation target is unchanged and
now better motivated: on an expected-null panel with a 0.63–0.82 ceiling and 1–3 haplotypes
per class, *no* method can demonstrate value, so this panel cannot adjudicate H1 either way.

---

## 6. Recommended next step

Unchanged in priority from MASTER_PLAN §4 item 2, but now with a prerequisite:

1. **Apply the §5.3 haplotype gate to candidate Panel A/B targets *before* fetching or
   embedding anything.** For *S. aureus* vs CoNS and *B. longum* subspp., slice first, then
   count unique haplotypes per class per region. This is cheap and will determine whether
   either target can support the analysis at all.
2. Then proceed to *S. aureus* vs coagulase-negative staphylococci (published V1 claim, real
   partial separation expected), re-using the corrected embedding path.
3. Re-run the full B1–B4 + gLM comparison there, with the permutation null and ceiling
   reporting established here.
