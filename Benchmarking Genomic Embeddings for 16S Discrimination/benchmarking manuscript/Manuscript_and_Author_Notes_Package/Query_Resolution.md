# Resolution of author queries AQ01–AQ24 and the nine specific checks

Each query is marked **RESOLVED** (the information is now in the PLOS manuscript or a named
file), **PARTIAL** (resolved as far as the data allow, with the residual stated), or
**AUTHOR** (requires information only the authors hold). Nothing is marked resolved unless
the supporting file exists.

---

## Numbered queries

| # | Status | Resolution |
|---|---|---|
| AQ01 | **RESOLVED** | Table 1 plus its note. Records / unique accessions / isolate groups: development 12,356 / 12,243 / 12,188; locked test 6,933 / 6,926 / 6,760; hard panels 4,720 / 4,720 / 4,657; screening-only 6,788 / 6,304 / 6,150. Development reports 25 analysed blocks because one is a subspecies stratum reported separately. |
| AQ02 | **AUTHOR** | Repository URL, release DOI and licences. Internal provenance exists: 25 commits, tagged per phase, freeze commit timestamped 94 minutes before the first locked-test record. `S1–S4 Text`. |
| AQ03 | **RESOLVED** | Block, family and rank recorded for every panel in `S1 Data` (`S5`/`S6` panel definitions). Development and locked-test families are disjoint by construction; the locked-test definition file asserts zero family overlap. |
| AQ04 | **RESOLVED** | Manifest with per-record species, strain, isolate group and sequence checksum is `S2 Data`. Methods now give the grouping rule in full (esummary `strain`, else the `strain` key parsed from `subtype`/`subname`; whitespace stripped; the placeholder set `-`, `na`, `n/a`, `none`, `unknown`, `missing`; the 18 bare culture-collection prefixes; key = `species|strain` or `species|ACC:accession`), state explicitly that **no** normalisation was applied beyond them (no case folding, no punctuation stripping, no catalogue cross-reference), and quantify the consequence: strain coverage 23,459/24,009 (97.7%); aggressive renormalisation would merge only 45 of 23,605 groups (0.19%); 313 groups hold >1 record (717 records, max 6); 120 accessions appear under two species queries, 112 within one block, and these lower the ceiling rather than inflate accuracy. Reproduced by `audit_isolate_grouping.py` into `results/isolate_grouping_audit.json`. The residual is now carried as a stated limitation in both Methods and Discussion: the key is submitter-supplied metadata, not verified provenance; standalone 16S deposits have no parent-genome link; the error direction is under-merging, hence optimistic. |
| AQ05 | **RESOLVED** | 42 candidates screened; 33 returned a score and 9 had insufficient data; 15 selected. Thresholds were ceiling ≤ 0.97 and k-mer balanced accuracy ≤ 0.95, fixed before screening. Full list and per-candidate values in `screen_results.json` and `screen_results_round2.json`. |
| AQ06 | **RESOLVED** | Methods give the aligner, its scoring parameters, the coordinate reference (*E. coli* J01695), the 0.75–1.30× length window and the handling of unavailable regions, and now add: (a) why coordinate transfer was chosen over an HMM extractor such as V-Xtractor [9] — every arm must see the same span, so a boundary set by a learned criterion would be confounded with the representation under test; (b) measured per-region availability (187,889/216,081 cells, 87.0%; V4 100%, V5–V8 ≥ 99.97%, V2 98.5%, V3 77.2%, V1 58.3%, V9 48.7%; no record failed all nine); (c) strand — no orientation check was performed, and an audit found 22/24,009 records (0.09%) deposited reverse-complemented and **retained** rather than rejected, contradicting the earlier assumption; a leave-out sensitivity analysis moves the primary estimand by 0.003 pp and changes no verdict, and the defect is now stated in Methods and Limitations; (d) ambiguous bases — N is the only non-ACGT character present (0.55% of slices, 4.7% of haplotypes), treated as a mismatch by the aligner, skipped by the k-mer featurizer, passed literally to tokenizers, with a measured ceiling effect of ≤0.03 pp. |
| AQ07 | **RESOLVED** | Perturbation algorithms, rates, trim length, seeding per (block, condition) and the fact that perturbation follows extraction and precedes splitting are in Methods and `robustness_perturbation.py`. Described as stress tests, not simulations. |
| AQ08 | **RESOLVED** | Checkpoint provenance in Methods: SHA-256 prefixes and byte sizes for the two local models, hub revision hashes for the other two, plus the statement that patches are compatibility-only. |
| AQ09 | **RESOLVED** | Pair construction and subsampling for R0, tie-aware fractional credit for R1, and the exact probe specification (StandardScaler, LogisticRegression, L2, C = 1.0, balanced class weights, max_iter 2000, seed 0) are in Methods. No hyperparameter search was performed; the probe was fixed in the frozen plan. |
| AQ10 | **RESOLVED** | Five depths per model at 0.25, 0.50, 0.75, 0.90 and 1.00 of encoder depth. The manuscript states that depth was chosen on the same blocks used for the comparison, making the best-depth figure optimistic rather than a bound. |
| AQ11 | **PARTIAL** | Permutation unit (isolate), stratification (class balance preserved within block), count (20 per block), full refit per permutation and seeding are in Methods; per-block values in `label_permutation_control.json`. **Residual**: 20 permutations give coarse tail resolution. The manuscript describes it as an initial null check and notes the 46-point separation is not sensitive to the count. |
| AQ12 | **RESOLVED** | Margin rationale in Methods: 1 point is below the level that would change a region recommendation and below observed within-block fold-to-fold variation. Fixed before locked-test collection; recorded in `S2 Text` and deviation D-001 in `S3 Text`. |
| AQ13 | **RESOLVED** | Ceiling equation is given explicitly in Methods, with unit, weighting and region stated as matched to the metric. The residual mismatch is quantified: 2 of 40 locked-test cells exceed the modal-region ceiling, maximum 0.30 points, both in Veillonella. Cause identified below. |
| AQ14 | **RESOLVED** | Region-selection algorithm (inner 3-fold on training isolates, C_kmer only, applied identically to all arms) and the evaluation-only oracle are in Methods and `region_ranking.py`. |
| AQ15 | **PARTIAL** | TOST is a bootstrap implementation on the same paired block-level estimand as the interval, 20,000 draws; overall p = 0.0001 as the maximum of the two one-sided values. **Residual**: the script and raw output should be deposited, and the small bootstrap-versus-t difference explained in the supplement. |
| AQ16 | **RESOLVED** | All scripts in `S1 Code`; unit-level predictions in `S3 Data`; bootstrap code inline in the analysis scripts; multiplicity family declared in `S2 Text`; missing cells reported as zero in the locked test. Family-clustered sensitivity was run on the development set (22 families, mean −0.0001). |
| AQ17 | **RESOLVED** | Both primary arms were named in advance from development results (DNABERT-2 best at R2 in 20 of 25 blocks; C_kmer best classical in 17 of 25), which removes selection from the statistic. The best-of-four rule is reported only as a Holm-corrected secondary and is described as upward-biased. |
| AQ18 | **RESOLVED** | The 20.5–26.9 point and 30-fold statements are withdrawn. R0 is a pair-level AUC and is never subtracted from R1 or R2. The same-endpoint statement is R1 → R2 of +2.8 to +3.0 points against a 0.73-point model spread, roughly fourfold. Table 4 labels each metric. |
| AQ19 | **RESOLVED** | Oracle definition, tie handling, fold structure and the use of a common probe across arms are in Methods. The zero-regret-with-low-hit-rate pattern is explained below. |
| AQ20 | **AUTHOR** | Repository URL, archived DOI, access date, licences. |
| AQ21 | **AUTHOR** | Institutional determination. All data are public records; no human or animal subjects. |
| AQ22 | **AUTHOR** | Authors, affiliations, corresponding author, CRediT roles, funding, competing interests. |
| AQ23 | **AUTHOR** | Verify all 13 references against publisher records before submission. |
| AQ24 | **RESOLVED** | Johnson–Lindenstrauss is cited as the primary source for distance preservation. The manuscript now states correctly that the **locked test** and the permutation control both use one fixed projection (seed 0) — only the development readouts draw a new projection per fold — and makes no invariance claim on that basis. A 20-seed recomputation is reported as a post hoc check: seed 0 reproduces all 20 block values exactly; across-block mean 96.21% (SD 0.09 pp, range 96.02–96.36); best-of-four R2 − C_random is +0.21 pp at seed 0 and +0.02 to +0.36 pp across seeds, every seed inside the ±1 pp margin. |

---

## The nine specific checks

**1. Veillonella ceiling mismatch — explained, and it is not a labelling problem.**
All five folds selected V8, so region switching does not account for it. The true cause is a
unit mismatch: the ceiling is a full-population bound, whereas the reported figure is the
mean of per-fold balanced accuracies, each computed on a fold with its own class
composition. A fold can exceed the population bound, so their mean can too. Across the
locked test this occurred in 2 of 40 cells, maximum excess 0.30 points. Now stated in
Methods rather than relabelled.

**2. TOST specification.** Bootstrap TOST on the same paired block-level estimand, 20,000
resamples, seed 1729. Component one-sided p-values 0.0001 and <0.00005. The 90% interval is
[−0.53, +0.11] and the 95% interval [−0.60, +0.16]; both lie inside ±1.0, so the conclusion
does not depend on which is used. Script deposit flagged as a residual.

**3. Random-projection seeds — now measured.** The confirmatory `C_random` arm uses one
fixed projection (seed 0): `score_test.py` reseeds inside the fold loop, as does the
permutation control. Only the development readouts draw a new projection per fold, and an
earlier note had this backwards; the Methods text has been corrected. The manuscript states
the single seed and claims no invariance on its basis. A 20-seed recomputation holding folds,
chosen regions and k-mer features fixed gives: seed 0 reproduces all 20 reported block values
exactly; across-block mean 96.215% (SD 0.087 pp, range 96.020–96.358); median within-block SD
0.155 pp, max 1.174 pp (francisella, 55 units); best-of-four R2 − C_random +0.214 pp at seed 0
and +0.020 to +0.358 pp across seeds, every seed inside the ±1 pp margin with no sign reversal.
Reported in Methods as post hoc, outside the frozen analysis. Label permutation remains
described as an initial null check.

**4. Zero regret with a low best-region hit rate — consistent, not contradictory.**
Several regions frequently tie at the maximum held-out accuracy. A method that selects a
tied region incurs zero regret while not matching the single argmax used for the hit flag.
Four locked-test blocks show this pattern (`rothia_micrococcaceae` regret 0.000, hit 0.00;
`brucella_ochrobactrum` 0.000 / 0.40; `porphyromonas` 0.000 / 0.20; `fusobacterium`
0.000 / 0.80). Now explained in Results and in the Table 7 note.

**5. "Copy variation" is misnamed — accepted.** The condition introduces one to three
substitutions into a single deposited sequence. It does not use multiple verified 16S copies
from one genome assembly, because standalone 16S deposits carry no parent-genome link. The
manuscript flags this and proposes renaming it within-isolate sequence variation.

**6. Isolate identity — documented and quantified.** The group key is the deposited strain
string: useful metadata, not verified provenance. Methods now state the full derivation and
the normalisation that was *not* performed (no case folding, no punctuation stripping, no
culture-collection catalogue cross-reference), so deposits of one isolate held under an ATCC
number and its DSM equivalent remain separate groups. The formatting component is bounded at
45 of 23,605 groups (0.19%) by re-deriving groups after upper-casing and stripping all
non-alphanumeric characters; the catalogue-equivalence component cannot be bounded from the
records. Coverage is 23,459/24,009 records (97.7%), and grouping binds on 313 groups holding
717 records, so grouped CV departs from record-level CV only for those. The limitation and
its direction (under-merging, therefore optimistic) are stated in Methods and in Discussion
Limitations. Numbers reproduced by `audit_isolate_grouping.py`.

**7. Figure axes.** Figure 2 has been rebuilt: R0 now sits on its own axis explicitly
labelled as a pair AUC and not comparable, and only balanced-accuracy quantities share the
common-axis panel. No figure now presents an AUC-to-accuracy change as a single effect.

**8. Hard-panel wording.** The manuscript states that the interval includes zero and also
spans the margin, so it neither demonstrates a benefit nor establishes equivalence. It is
not described as proving no benefit.

**9. Units.** All differences are reported in percentage points of balanced accuracy, with
the metric named at first use in each section.

---

## Residual blockers before submission

1. Author, affiliation, funding, competing-interest and ethics statements (AQ02, AQ20–AQ23).
2. Public repository and archived DOI.
3. Deposit the TOST and permutation scripts and raw outputs (AQ15, AQ11).
4. Decide whether to rename the copy-variation condition (check 5).
5. Verify all references against publisher records (AQ23).

Closed since the last revision: strand handling and ambiguous-base treatment (AQ06) and the projection-seed question (AQ24, check 3) are now measured and written into Methods. AQ06 surfaced one unrepaired pipeline defect (minus-strand records are mis-sliced, not rejected); it is disclosed with its sensitivity analysis rather than silently fixed, because the locked test is scored once under a frozen plan.
