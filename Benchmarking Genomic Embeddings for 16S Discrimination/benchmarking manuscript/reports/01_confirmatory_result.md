# Phase 5 — confirmatory results on the locked test set

**Analysis plan:** `../phase3_freeze/FREEZE.md`, committed 2026-11-08 07:20:42, unmodified
since (verified by `git diff HEAD`).
**Locked test data:** first record fetched 2026-11-08 08:54:19 — **94 minutes after** the
plan was frozen. No test taxon was fetched, sliced, embedded or scored before the freeze.
**Manifest:** 6,933 records / 6,926 accessions / 6,760 isolate groups,
`accession_manifest.tsv` sha256 `a4f1aea327a5d053…`, `test_blocks.json` sha256
`29cebfcf8ab16f1d…`.
**Run once.** No tuning, no re-runs, no post-hoc changes.

---

## 1. Primary endpoint — prespecified, uncorrected

**DNABERT-2:R2 − C_kmer**, paired over 20 independent eligible blocks, bootstrap 2000 draws,
seed 1729.

| quantity | value |
|---|---|
| n blocks | **20** (freeze required ≥14) |
| mean paired difference | **−0.00211** |
| 95% bootstrap CI | **[−0.00595, +0.00163]** |
| CI half-width | 0.00379 |
| two-sided bootstrap p | 0.253 |
| blocks model > / < / tie (±0.01) | 1 / 2 / **17** |
| **verdict at the frozen ±0.010 margin** | **EQUIVALENT — practically negligible difference** |
| sensitivity at the protocol's illustrative ±0.030 | EQUIVALENT |

The entire confidence interval lies inside ±0.006. This is a **positive equivalence finding**,
not a failure to detect a difference: the data exclude any advantage for the frozen
representation larger than about 0.002, and any deficit larger than about 0.006, on the
metric and population specified in advance.

Observed precision matched the plan almost exactly — predicted half-width ≈0.005 from the
development SD, achieved 0.0038.

## 2. Secondary endpoints — 4 tests, Holm-corrected

| contrast | mean | 95% CI | p raw | p Holm | verdict @ ±0.010 |
|---|---|---|---|---|---|
| DNABERT-S:R2 − C_kmer | −0.0011 | [−0.0047, +0.0021] | 0.580 | 0.580 | EQUIVALENT |
| **HyenaDNA:R2 − C_kmer** | **−0.0088** | **[−0.0148, −0.0038]** | **0.001** | **0.004** | INCONCLUSIVE at this margin |
| NTv2:R2 − C_kmer | −0.0028 | [−0.0068, +0.0011] | 0.168 | 0.336 | EQUIVALENT |
| best-of-4 R2 − C_kmer | +0.0029 | [+0.0003, +0.0057] | 0.025 | 0.075 | EQUIVALENT |

**HyenaDNA is the one arm that is reliably, if slightly, worse than k-mer composition**
(Holm-adjusted p = 0.004; CI entirely below zero). The effect is −0.009 — statistically clear
but below the ±0.010 practical margin, so it is formally "inconclusive" under the frozen
equivalence rule while being a confirmed small deficit under the superiority test. Both facts
are reported; neither is suppressed.

The best-of-4 arm is nominally positive (+0.0029) but does not survive Holm correction
(p = 0.075), and its advantage is within the margin. Picking the best of four models per
block does not buy a practically meaningful gain.

## 3. What this confirms, and what it does not

**Confirmed (confirmatory, locked test, 20 independent blocks from 20 families disjoint from
all development data):** on fine-scale 16S target-vs-neighbour discrimination, a frozen
pretrained genomic language model with a linear probe performs **equivalently to 4/5/6-mer
composition with the same probe**, within ±0.010 balanced accuracy.

**This replicates the development-phase finding** (+0.0014 [−0.0032, +0.0057] over 25
development blocks) on entirely independent taxa, with a point estimate of the opposite sign
and a comparably narrow interval. Two independent block sets, 45 blocks total, same
conclusion.

**Not supported by this result** (protocol §22's unsupported extensions):
- *"Genomic language models contain no 16S discriminatory information."* False — every model
  reaches 0.95–0.96 mean balanced accuracy, against a mean ceiling near 0.98.
- *"Genomic language models are useless for 16S."* Not tested. This is one task (closed-panel
  target-vs-neighbour), one readout family, four checkpoints, clean reference sequences.
- *"Classical methods are superior."* The primary contrast is equivalence, not superiority.
  Only HyenaDNA shows a reliable deficit, and it is below the practical margin.

## 4. Exploratory (never promoted to confirmatory)

Mean balanced accuracy over the 20 blocks:

| arm | mean |
|---|---|
| C_random | 0.9616 |
| C_kmer | 0.9609 |
| DNABERT-S:R2 | 0.9598 |
| DNABERT-2:R2 | 0.9588 |
| NTv2:R2 | 0.9581 |
| HyenaDNA:R2 | 0.9521 |
| C_align | 0.9429 |
| C_exact | 0.9028 |

**The most uncomfortable number in this table is `C_random` at the top.** A random projection
of the k-mer space scores 0.9616 — nominally above every learned representation and above
tuned k-mer itself. Best-of-4 minus C_random is +0.0021 [−0.0005, +0.0052], i.e.
indistinguishable. This is the §11 flexibility reference doing exactly the job it was kept
for: at this difficulty level, most of the measured performance reflects a flexible
classifier over a high-dimensional feature space, not the specific content of any
representation. Any future claim that a representation is "better" on this task must clear
this bar first.

**Calibration of the a-priori classes was poor on the test set** (expected A 0.998, B 0.999,
C 0.914): the A/B distinction collapsed entirely. Only the "highly ambiguous" class C carried
information, driven by brucella (0.786) and veillonella (0.849). Predicting *degree* of
16S separability in advance is evidently hard; predicting the extremes is feasible.

**Region selection** (support-only, 100 folds): V2 chosen 34×, then V6 14×, V1 12×, V8 11×,
V4 9×, V5 7×, V3 6×, V7 6×, V9 1×. V2's dominance is consistent with it being the longest
window in the panel (~115 nt).

## 5. Accounting

| | |
|---|---|
| declared locked-test blocks | 22 |
| eligible and scored (confirmatory) | 20 |
| descriptive-only (failed pre-declared gate) | 2 — `actinomyces` (7 target isolates), `bartonella` (1/3 neighbour species ≥5) |
| **missing or NaN cells** | **0** across 8 arms × 20 blocks |

## 6. Files

`analyze_locked_test.py` → `locked_test_results.json`, `phase5_output.txt`.
Inputs: `../phase4_locked_test/locked_test_scores.json`, `test_block_summary.json`,
`accession_manifest.tsv`, `manifest_hashes.json`.
