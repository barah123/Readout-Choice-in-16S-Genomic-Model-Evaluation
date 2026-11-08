# Phase 3 — FROZEN ANALYSIS PLAN

**Protocol:** `../05_Third_Idea_Critical_Review_and_Benchmark_Protocol.md` §18 Phase 3
**Frozen on:** 2026-11-08
**Status:** binding. Any change after this file is committed requires a dated entry in
`DEVIATIONS.md` using the §12 protocol-deviation form, stating whether it was made before or
after locked-test inspection.

---

## 0. The problem this freeze must solve, stated plainly

**All 31 blocks built in Phase 2 are now DEVELOPMENT data and are burned for confirmatory
use.** Every decision below — that R2 is the primary readout, that `C_random` belongs in the
comparator set, that DNABERT-2 and `C_kmer` are the single pre-named arms, the equivalence
margin, the precision calculation — was derived by looking at results on those 31 blocks.
Re-using them as a test set would be exactly the failure protocol §18's stopping rules name:
*"Repeated test-guided model selection: acquire a new test or label analysis exploratory."*

Therefore Phase 3 declares:

| set | blocks | role |
|---|---|---|
| **Development** | the 31 blocks in `blocks_definition.py` | burned. All tuning/selection done. Results are exploratory/descriptive forever. |
| **Locked test** | the new blocks in `test_blocks_definition.py` | never fetched, sliced, embedded or scored before this file was committed. Single use. |

Nothing in Phase 2's reports changes status: they remain valid **development** findings and
are reported as such. The Phase 2-full headline (frozen ≈ classical, CI ±0.005) is a
*hypothesis* that Phase 5 will test on locked data, not a confirmed result.

---

## 1. Models and checkpoints (§9, §19)

Family-representative rule, fixed: one checkpoint per model family, the one already used in
development, no size sweep.

| model | source | revision / checkpoint identity | dim | env |
|---|---|---|---|---|
| DNABERT-2 | `zhihan1996/DNABERT-2-117M` → local patched `dnabert2_local/` | weights sha256 `7ff39ec77a484dd01070a41bfd6e95cd…`, 468,354,983 B | 768 | `glm-benchmark` |
| DNABERT-S | `zhihan1996/DNABERT-S` → local patched `dnabert_s_local/` | weights sha256 `f3cfc3d0541859df64759e758cc4bc40…`, 468,320,053 B | 768 | `glm-benchmark` |
| HyenaDNA | `LongSafari/hyenadna-medium-160k-seqlen-hf` | hub revision `7ebf71773d22c0ede2cc55cb2be15ee8c289e1ce` | 256 | `glm-benchmark` |
| NTv2 | `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` | hub revision `06615c1660c892fc199840c18123f8385b3542a8` | 1024 | `glm-ntv2` |

All four carry Apache-2.0 or equivalent permissive terms. Patches to the two MosaicBERT repos
are CPU/transformers-compatibility only, no numerical change
(`dnabert2_patched/README_PATCHES.md`). **No model may be added, swapped or re-versioned
after this point**; if a checkpoint becomes unavailable, record it as a technical exclusion
rather than substituting another.

Pooling is frozen as **masked mean over real tokens, excluding each row's own first and last
special token, last hidden layer**, embedding each unique slice once and broadcasting
(guarantees purity by construction). No layer or pooling sweep in the confirmatory analysis;
any such analysis is exploratory and labelled so (§17.2).

## 2. Readout (§10)

**Primary readout: R2, a regularized linear probe** — `StandardScaler` → `LogisticRegression`
(`max_iter=2000`, `class_weight="balanced"`, L2, default C=1.0, `random_state=0`), fit on
training folds only.

Justification from development: spread across readouts within a model is 0.21–0.27 balanced
accuracy, versus 0.007 across models within R2. R0 (unsupervised geometry) demonstrably
misranks — it scored 0.79–0.83 on cereus regions whose ceiling is 0.55, i.e. separating
something that cannot be the class label.

R0 and R1 are still computed and reported, as **diagnostics only**. R3 (bounded nonlinear
probe) and R4 (layer/pooling variants) are **not** in the confirmatory family; if run, they
are exploratory.

## 3. Primary endpoint (§13)

**Metric:** isolate-weighted, species-balanced accuracy. Each isolate contributes its modal
slice once per region; target and pooled-neighbour halves weighted equally. Matched
isolate-weighted ceiling reported beside every cell.

**Cross-validation:** 5-fold `StratifiedGroupKFold`, grouped on **isolate**, `random_state=0`.
No isolate's records cross a fold.

**Region selection — support-only, no query labels.** For each locked-test block, the region
is chosen by inner grouped cross-validation *within the training folds using the classical arm
only*, then applied **identically to both arms**. This isolates classification ability at a
common region and keeps query labels out of the selection, per §13 ("No query-derived
performance may enter the ranking procedure").

**Primary estimand — a single pre-named pair, no selection inside the statistic:**

> mean paired difference, across locked-test blocks, of
> **DNABERT-2 : R2_linear** minus **C_kmer**

Both arms are pre-named from development (DNABERT-2 best at R2 in 20/25 blocks; `C_kmer`
best classical in 17/25). Naming them removes the model-selection and control-selection steps
from the statistic, which §14.2 would otherwise require any permutation reference to
replicate.

**Inference:** block is the unit. Bootstrap over blocks, 2000 draws, seed 1729, 95% percentile
CI. No t-test over folds; folds, pairwise distances and seeds never count as n.

## 4. Equivalence margin (§15)

**Margin: ±0.010 absolute balanced-accuracy units.**

Justified on application grounds, not on observed precision: a difference below 0.01 balanced
accuracy does not change which region or amplicon a recommendation tool would return, and is
smaller than the fold-to-fold variation seen within single cells in development. A method
needing 0.01 to justify ~1000× the compute of k-mer counting is not practically preferable.

**Procedure:** two one-sided tests via the bootstrap CI — equivalence is declared only if the
**entire 95% CI lies inside (−0.010, +0.010)**. A non-significant superiority test alone will
not be reported as equivalence.

Possible outcomes, all pre-committed as reportable:
- CI entirely inside ±0.010 → *practically negligible difference*
- CI entirely above +0.010 → frozen representations superior
- CI entirely below −0.010 → classical superior
- CI straddling a boundary → *inconclusive at this margin and n*

## 5. Precision plan and required n (§16)

Development block-level SD of the paired delta = **0.01138** (n=25).

| target 95% CI half-width | blocks needed |
|---|---|
| 0.0150 | 2.2 |
| 0.0100 | 5.0 |
| 0.0075 | 8.8 |
| 0.0050 | 19.9 |

**Target: ≥14 eligible locked-test blocks**, giving an expected half-width ≈ 0.006 —
comfortably inside the ±0.010 margin, so an equivalence verdict is attainable rather than
foredoomed to "inconclusive". In development 5/31 declared blocks (16%) failed eligibility;
the test set deliberately includes rarer taxa, so assume a higher ~35% failure rate and
**declare 22 test blocks** to land ≈14 eligible. If fewer than 10 prove eligible, the primary
analysis is reported as **underpowered and inconclusive** rather than spun.

## 6. Eligibility gate (unchanged from development, restated as binding)

A region-cell is confirmatory-eligible iff the target has **≥10 independent isolates** and
**≥2 neighbour species each have ≥5**. A block is eligible if ≥1 region qualifies. Ineligible
blocks are **retained and reported as descriptive-only**, never deleted. Per-region status
codes stay distinct: `eligible` / `insufficient_target_isolates` /
`insufficient_neighbour_species` / `no_coverage_one_class` (§19).

## 7. Multiplicity (§16)

Declared confirmatory family, **5 tests total**:

| # | test | correction |
|---|---|---|
| 1 | **PRIMARY**: DNABERT-2:R2 − C_kmer, equivalence at ±0.010 | none (single primary) |
| 2–5 | **SECONDARY**: each of DNABERT-S, HyenaDNA, NTv2, and best-of-4-R2 vs C_kmer | Holm across these 4 |

Everything else — R0/R1 results, per-region tables, the T2 region-ranking endpoint, ceiling
calibration, `C_random`, subspecies stratum, family-grouped sensitivity, failure-mode
diagnostics — is **exploratory** and will be labelled exploratory in every table and in the
manuscript. No exploratory result may be promoted to a confirmatory claim.

## 8. Secondary endpoint T2 (region ranking), exploratory

Each arm selects its own region support-only; report held-out accuracy of that choice and
regret against the evaluation-only best region. Reported as exploratory because the primary
region rule (§3) deliberately fixes a common region for both arms.

## 9. Comparator set (frozen)

`C_exact` (exact slice lookup with ambiguity) · `C_align` (tie-aware 1-NN on pairwise
identity) · **`C_kmer`** (canonical 4/5/6-mer, vocabulary fit per training fold, same probe
as R2) · `C_random` (matched-dim random projection of the k-mer space).

`C_random` is retained on the strength of development, where it reached 0.980 against
`C_kmer`'s 0.981 — small deltas between high-dimensional feature sets must be read against
it. It is a flexibility reference, **not** the primary comparator.

## 10. Compute budget

Development cost, measured: fetch ≈ 65 min; slicing 12,356 records ≈ 25 min (8 cores);
embedding 7,878 unique slices = 45 s (DNABERT-2), 62 s (DNABERT-S), 474 s (HyenaDNA), 187 s
(NTv2); readouts 277 cells ≈ 1 min (6 workers). All CPU, MacBook Air M4, 16 GB.

Locked-test budget: **≤3 hours wall-clock total**, same hardware, no GPU. If the test set
exceeds this, reduce the number of declared blocks *before* fetching — never after seeing
results. Cold-start and cached inference recorded separately.

## 11. Unsealing conditions (§9 of the schemas doc)

The locked test is scored **once**. Before unsealing:
- [ ] this FREEZE.md committed and unmodified
- [ ] `test_blocks_definition.py` committed before any test fetch
- [ ] split manifest hash and accession list archived
- [ ] all analysis code able to run end-to-end on development data without modification
- [ ] a negative or equivalent result is pre-accepted as the reportable outcome

If the primary result prompts a change of method, that change is **development work requiring
a new locked test set** — it may not be reported as a confirmatory finding on this one.
