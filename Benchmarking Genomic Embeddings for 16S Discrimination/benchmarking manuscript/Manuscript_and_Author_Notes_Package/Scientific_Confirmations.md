# Answers to the ten scientific confirmations

Every answer below was extracted from the code or result files, not from recollection. Two
items required correcting the manuscript; they are marked **CORRECTION**.

---

## 1. Reconcile 42 versus 53 screened candidates, and the record counts

**53 counts screening attempts; 42 is the number of distinct candidates.** Round 1 screened
30 candidates. Round 2 contained 23 entries, of which **11 were re-screens** of round-1
candidates that had failed on insufficient data at the cheap sampling depth, and 12 were
new. So 30 + 23 = 53 attempts, but 30 + 12 = 42 unique candidates. **Use 42 for candidates
and state the 11 re-screens separately.** The manuscript currently says 42; any text saying
53 refers to attempts and should say so.

**Record counts** (all verifiable from the manifests):

| set | records | unique accessions | isolate groups |
|---|---|---|---|
| Development | 12,356 | 12,243 | 12,188 |
| Locked test | 6,933 | 6,926 | 6,760 |
| Hard panels | 4,720 | 4,720 | 4,657 |
| Screening only | 6,788 | 6,304 | 6,150 |
| **Analysed total** | **24,009** | 23,889 | 23,605 |

Records exceed accessions because the same accession can be returned for more than one
species query where NCBI organism annotation overlaps. Accessions exceed isolate groups
because several deposits share a strain identifier.

**The development subspecies stratum** is `bifidobacterium_longum_subsp`. It is included in
the 31 declared and 26 eligible development blocks and in the 12,356 records, but **excluded
from the 25 analysed blocks** and from every pooled estimate, because species-rank and
subspecies-rank panels are separate strata. It is reported on its own in the development
report.

## 2. Region-ranking oracle, tie tolerance, repeated selections

**The oracle is method-specific.** For each outer fold and each arm, held-out accuracy is
computed for every eligible region *for that arm*, and the oracle is that arm's own maximum
(`best_region = max(heldout[arm], key=heldout[arm].get)`, `region_ranking.py` line 144).
Regret is therefore each method's shortfall against its own best attainable region, not
against a shared oracle. This is the conservative choice: it does not penalise a method for
being weaker overall, only for choosing badly among its own options.

**Tie tolerance is zero.** The oracle uses a strict `max`; where several regions tie at the
maximum, Python returns the first encountered. This is why a method can record **zero regret
with a hit rate below 100 percent**: selecting a tied-best region costs nothing but does not
match the single argmax used for the hit flag. Four locked-test blocks show this
(`rothia_micrococcaceae` 0.000 regret / 0.00 hit; `brucella_ochrobactrum` 0.000 / 0.40;
`porphyromonas` 0.000 / 0.20; `fusobacterium` 0.000 / 0.80).

**Five selections per block**, one per outer fold of the five-fold split; reported values are
means over those five. Selection used inner three-fold cross-validation on training isolates
only.

**Recommended manuscript addition:** state that the hit rate is a strict-argmax criterion and
is therefore pessimistic under ties, and that regret is the better-behaved metric.

## 3. Strain-identifier normalization and neighbor-selection sources

**Normalization rule, in full** (`fetch_blocks.clean_strain`): take the `strain` field from
esummary, or failing that the `strain` key parsed from the `subname`/`subtype` pair; strip
whitespace; treat as missing the values `-`, `na`, `n/a`, `none`, `unknown`, `missing`
(case-insensitive); treat as missing any value that exactly matches a bare culture-collection
prefix from the set {ATCC, DSM, NCTC, CECT, JCM, KCTC, CIP, LMG, NBRC, CGMCC, BCRC, NRRL,
CCUG, VTT, CCM, IFO, KACC, BAA}. The isolate group is then `species|strain`, or
`species|ACC:accession` when the strain is missing. No further normalization was applied: no
case folding, no punctuation stripping, no cross-referencing to culture-collection catalogues.

**This is deposited metadata, not verified provenance.** Two deposits of the same physical
isolate under differently formatted strain strings would be treated as independent, and
standalone 16S records carry no link to a parent genome assembly. The manuscript states this
as a limitation.

**Neighbor selection** was made on a documented rule (most frequently confusable congeners,
or closest relatives by published phylogeny) with a one-line rationale recorded per panel in
the committed definition files. **Per-panel primary citations were not recorded.** This is a
genuine gap: the rationales are author judgement informed by the literature, not traceable
to a specific source for each panel. Either add citations per panel or describe the rule as
expert selection rather than literature-derived.

## 4. Region coordinates, lengths, strand, N handling, tie rules, processing order

**Coordinates** are *Escherichia coli* J01695 positions, inclusive, transferred through a
global pairwise alignment:

| region | start–end | expected length (nt) |
|---|---|---|
| V1 | 69–99 | 31 |
| V2 | 137–242 | 106 |
| V3 | 433–497 | 65 |
| V4 | 576–682 | 107 |
| V5 | 822–879 | 58 |
| V6 | 986–1043 | 58 |
| V7 | 1117–1173 | 57 |
| V8 | 1243–1294 | 52 |
| V9 | 1435–1465 | 31 |

**Retention window:** an extracted region is kept only if its length lies within 0.75 to 1.30
times the expected length; otherwise it is recorded as unavailable for that sequence and the
isolate contributes nothing to that region. Unavailable regions are counted and reported, not
imputed.

**Processing order at fetch** (`fetch_blocks.py` lines 96–97): whitespace removal, uppercase,
**then** U to T conversion, **then** the 1,350–1,600 nt length filter and the alphabet check.
So uracil conversion precedes both filters, and a sequence written in RNA alphabet is not
rejected by the alphabet check.

**N handling:** N is permitted in the retained alphabet {A, C, G, T, N} at the sequence level.
In k-mer featurization any k-mer containing N is skipped, so Ns contribute no features. They
are not excluded from length counts and are not imputed.

**Strand:** no strand normalization or reverse-complement search was performed. Records are
used in the orientation deposited; the global alignment to a forward-strand reference means a
reverse-complemented deposit would fail the length window and be dropped rather than silently
mis-sliced. This was not separately audited.

**Modal-sequence tie rule:** where an isolate contributes several records to one region, the
isolate is represented by `Counter(...).most_common(1)[0][0]`, which under a tie returns the
**first-inserted** sequence, i.e. the first encountered in file order. This is an arbitrary
but deterministic rule. It affects only isolates with multiple records and a tie, which is a
small minority (6 records across 3 isolates in the staphylococcal panel).

## 5. Projection seed and scaling; alignment and exact-match ties

**Scaling: 1/p is the VARIANCE.** The implementation is
`rng.standard_normal((p, 256)) / np.sqrt(p)`, so each entry has standard deviation
1/√p and variance 1/p. Verified empirically: for p = 1,299 the realized variance is
7.718 × 10⁻⁴ against 1/p = 7.698 × 10⁻⁴. The manuscript's "N(0, 1/p)" is correct with 1/p
read as variance; say so explicitly.

**CORRECTION — the seed differs between two scripts, and the Methods text is accurate for
only one of them.**

- `run_readouts_full.py` (main analysis) draws the projection from the **shared, advancing**
  generator, so a **different** projection is used for every fold and every cell.
- `label_permutation_control.py` creates a **fresh generator from a fixed seed**, so the
  **same** projection is used in every fold.

Both are legitimate random projections and neither invalidates a result, but the Methods
sentence stating a single fixed seed applies only to the permutation control. The main-run
behaviour should be described as an independent draw per fold, which is arguably the stronger
design because it averages over projections rather than depending on one.

**Alignment and exact-match ties and abstentions.** `C_align` is tie-aware: all training
neighbours within 1 × 10⁻⁹ of the minimum distance are collected and the score is the
fraction carrying the true label, so a tie yields fractional rather than arbitrary credit.
`C_exact` assigns the proportion of identical training sequences carrying the true label and
**abstains at 0.5 when the test sequence is absent from training**; abstentions are scored,
not dropped, so they lower the reported accuracy rather than shrinking the denominator.

## 6. Bootstrap TOST — **CORRECTION**

**The previously reported p = 0.0001 was under-resolved and was not a boundary-null
construction. The correct value is p = 4.1 × 10⁻⁴ from a paired-t TOST.**

What was actually computed before: the proportion of 20,000 percentile bootstrap means
falling beyond each margin. Two draws fell below −1.0 point (2/20,000 = 1.0 × 10⁻⁴) and zero
above +1.0 point. A count of zero gives a floor of 1/B = 5 × 10⁻⁵, not a measured value, and
percentile tail counts are not a shifted-null bootstrap test. The reviewer is right on both
points.

**Correct analysis**, n = 20 blocks, mean −0.2108 points, SD 0.8892 points, SE 0.1988:

| quantity | value |
|---|---|
| t for H₀: δ ≤ −1.0 | 3.969, one-sided p = 4.109 × 10⁻⁴ |
| t for H₀: δ ≥ +1.0 | −6.089, one-sided p = 3.715 × 10⁻⁶ |
| **TOST p (the maximum)** | **4.1 × 10⁻⁴** |
| Parametric 90% CI | [−0.555, +0.133] points |
| Bootstrap 90% CI (B = 20,000) | [−0.529, +0.105] points |

Both 90% intervals lie entirely inside ±1.0 point, which is the TOST decision rule at
α = 0.05, so equivalence is declared by either route. **Report the paired-t TOST as the
primary test and the bootstrap interval as distribution-free corroboration**, rather than
quoting a bootstrap tail probability.

## 7. Multiplicity against the frozen plan — confirmed, no change

The frozen plan (`S2 Text`, section 7) declares a confirmatory family of exactly five tests:
one primary contrast, uncorrected, and four secondary contrasts Holm-corrected within their
family. The manuscript and Table 2 follow this exactly. No silent change was made. All other
analyses, including the hard panels, degradation, layers and region selection, are labelled
exploratory in both the plan and the manuscript.

## 8. Simulated variation: rate, distribution, seed, timing

| | substitution-error conditions | `copy_var` |
|---|---|---|
| Model | Binomial(L, r) substitutions, r = 0.001 or 0.01 | a fixed count drawn Uniform{1, 2, 3}, independent of length |
| Positions | sampled without replacement, uniform over the sequence | same |
| Replacement base | uniform over the three bases other than the original | same |
| Length dependence | **yes**, expected count scales with L | **no**, 1–3 regardless of length |

**Seeding:** one generator per (block, condition), seeded as
`abs(hash((bid, cond))) % 2**32`. The same generator supplies all perturbations within that
cell. Note that Python's `hash` on strings is salted per process unless `PYTHONHASHSEED` is
fixed, so exact reproduction requires setting it; the qualitative result does not depend on
the draw.

**Timing:** perturbation is applied **after** region extraction and **before** the
cross-validation split, to the isolate-level modal sequence, and both arms receive the
identical perturbed string.

**Distinction to state in the manuscript:** `copy_var` is not a length-scaled error model and
is not derived from observed multiple 16S copies; it introduces 1–3 substitutions into one
deposited sequence. As noted previously, the name overstates what it represents and a rename
to *within-isolate sequence variation* is recommended.

## 9. Software versions and checkpoints — verified against live environments

Queried directly from the two conda environments:

- **glm-benchmark:** Python 3.11.18, NumPy 2.4.6, SciPy 1.17.1, scikit-learn 1.9.1,
  PyTorch 2.15.0, Transformers 5.22.1, Biopython 1.88, joblib 1.6.0
- **glm-ntv2:** Python 3.11.18, PyTorch 2.15.0, Transformers 4.34.1

Checkpoints: DNABERT-2 weights SHA-256 begins `7ff39ec77a48`, 468,354,983 bytes; DNABERT-S
SHA-256 begins `f3cfc3d05418`, 468,320,053 bytes; HyenaDNA hub revision
`7ebf71773d22c0ede2cc55cb2be15ee8c289e1ce`; NTv2 hub revision
`06615c1660c892fc199840c18123f8385b3542a8`.

**Compatibility patches are not independently validated.** Four patches were applied to the
two MosaicBERT repositories for CPU and library compatibility. They are documented and are
asserted to make no numerical change, but no formal numerical-equivalence test against an
unpatched reference was run, because the unpatched code does not execute in this environment.
This should be stated rather than implied.

## 10. Plots and numerical precision

**Perturbation intervals come from `robustness_results.json`**, produced by the single
`robustness_perturbation.py` run. No new model run was performed when the manuscript was
assembled, and no interval was recomputed by hand.

**Figures** were rebuilt after the R0 correction. Figure 2 now places R0 on a separate axis
labelled as a pair AUC and places only balanced-accuracy quantities on the shared-axis panel.
No figure presents an AUC-to-accuracy change as a single effect. Figure 3 shows both the
permuted-label floor and the empirical ceiling.

**Precision:** all differences are reported in percentage points of balanced accuracy to two
decimals, matching the precision in the source JSON. The one place to watch is the primary
estimate, which is −0.2108 points and appears as −0.21 in tables and −0.2 in prose; these
should be harmonized to two decimals throughout.

---

## Actions arising

1. **Correct the TOST value** from 0.0001 to 4.1 × 10⁻⁴ and describe it as a paired-t TOST
   with the bootstrap interval as corroboration. *(item 6)*
2. **Correct the random-projection seeding description** to distinguish the main run from the
   permutation control. *(item 5)*
3. Add per-panel citations for neighbour selection, or describe it as expert selection.
   *(item 3)*
4. State that the strict-argmax hit rate is pessimistic under ties. *(item 2)*
5. State that compatibility patches were not numerically validated against unpatched code.
   *(item 9)*
6. Fix `PYTHONHASHSEED` or replace the hash-based perturbation seed for exact
   reproducibility. *(item 8)*
7. Harmonize the primary estimate to two decimals throughout. *(item 10)*
