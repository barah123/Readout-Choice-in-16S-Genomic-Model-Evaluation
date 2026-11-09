# Author Completion and Evidence Notes

Prepared November 8, 2026. This file accompanies the manuscript draft. The draft is complete in structure but not ready for submission until these details are resolved. The latest pasted report supersedes the earlier narrative for its stated corrections, but it does not replace missing source outputs.

## Evidence rules used

- Uploaded numeric result files support the primary contrast and existing extension estimates.
- The new TOST p = 0.0001, 90% CI, 256-dimensional Gaussian control, and label-permutation means are accepted as author-supplied findings. Their test output/scripts were not among the reviewed files.
- The primary equivalence result is retained. It is not generalized to HyenaDNA, every hard panel, or all perturbation conditions.
- R0 pair AUC and R2 balanced accuracy are kept distinct. The 30-fold and 50-fold ratios are not presented as effects on a common accuracy endpoint.
- Random projection is not described as removing biological information. Similar performance does not by itself establish feature-level redundancy.
- End trimming is not treated as an empirical primer-binding experiment.
- Best layer/model selected on the evaluation data is an optimistic estimate, not a formal upper bound.
- Development, locked-test, and hard-panel comparisons use different selection contexts and are not described as three independent formal equivalence replications.

## Numbered author queries

### AQ01

reconcile record totals across unique accessions, records, isolate groups, and repeated phase use; explain the development stratum reported separately.

### AQ02

cite and provide the protocol, commit hash, repository, timestamps, panel definitions, and accession manifests.

### AQ03

provide the family/block assignments and overlap matrix across all phases.

### AQ04

supply the complete accession manifest and describe how strain identifiers, isolate groups, taxonomic labels, and duplicate records were defined for each phase.

### AQ05

provide the full screened-candidate list, selection thresholds, taxonomy source/version, and exclusions.

### AQ06

specify the extraction software and version, coordinate/region definition method, strand normalization, sequence-length criteria, handling of ambiguous bases, and the reason each region was unavailable.

### AQ07

document the exact perturbation algorithms, random seeds, number of replicates, whether perturbations were applied before or after splitting, and how read or copy weights were handled. The report correctly notes that uniform substitutions and fixed end trimming are stress tests, not full simulations of sequencing and PCR.

### AQ08

provide the model configuration files and checkpoint provenance.

### AQ09

define the pair construction, AUC calculation, tie treatment, probe implementation, preprocessing, regularization, hyperparameter selection, and training folds.

### AQ10

specify which layer depths were selected and whether any layer or pooling decisions were made using the same panels later used for the layer comparison.

### AQ11

provide the permutation-level results, exchangeability unit, within-group stratification, number of permutations, full-refitting procedure, seeds, and uncertainty. Twenty permutations provide coarse tail resolution and should be described as an initial null check unless a larger analysis is supplied.

### AQ12

explain the application-based rationale for this smallest effect of interest and provide the original margin-selection documentation.

### AQ13

provide the exact ceiling equation, the per-block values, how shared sequences and multiple copies were weighted, and verify that each reported accuracy and ceiling uses identical labels, weights, and query units.

### AQ14

provide the region-selection algorithm and evaluation-only oracle calculation.

### AQ15

provide the TOST code/output, exact test statistic and component one-sided p-values, and clarify whether the reported TOST uses the same paired block-level estimand and assumptions as the bootstrap interval.

### AQ16

supply analysis scripts, unit-level predictions, bootstrap resampling code, multiplicity family, missing-data handling, and any family-clustered sensitivity analysis.

### AQ17

specify the exact selection rule for “best model” and “best classical” and ensure no test information entered this comparison.

### AQ18

revise the “20.5–26.9 point within-model spread” and “30-fold” statements unless all compared quantities are placed on a common, interpretable scale.

### AQ19

specify the oracle region definition, tie-aware hit criterion, fold/repeat structure, and whether all methods were evaluated with a common downstream classifier.

### AQ20

insert public repository URL, archived release/DOI, access date, licenses, and confirm that the described supplementary materials are included in the submission package.

### AQ21

provide institutional determination or state whether formal ethics review was not applicable, consistent with institutional policy.

### AQ22

provide funding sources, competing-interest statement, author names and affiliations, and contribution roles.

### AQ23

verify complete author list, article number, and citation formatting from the publisher record before submission.

### AQ24

verify the citation and preferred primary source for the random-projection statement.

## Specific checks that remain important

1. Explain the Veillonella ceiling mismatch in the existing aggregate file, or provide corrected per-split ceilings. Relabeling alone does not establish matched measurement.
2. Confirm the TOST nulls, numerical method, degrees of freedom or bootstrap rule, component p-values, unit of inference, and why the reported 90% interval differs slightly from a simple t interval over the uploaded differences. Differences can be legitimate but must be explained.
3. Verify random projection across seeds if interpreting invariance; treat label permutation as diagnostic rather than a strong tail test with only 20 permutations.
4. Reconcile zero regret with zero best-region hit rates and distinguish common versus method-specific oracle regions.
5. Verify whether the perturbation named copy variation represents true genomic copies or alternative database records linked by a strain field. Do not infer parent-genome copy completeness from standalone 16S deposits.
6. Confirm the unit of isolate identity. A reported strain string is useful metadata, but independent isolate provenance and normalized collection identifiers should be documented.
7. If the plotted figures still call R0 AUC 'accuracy' or visualize AUC-to-accuracy changes as a common quantitative effect, revise their axes/legends before submission.
8. Do not claim that a non-significant hard-panel result proves no benefit. Describe the estimated contrast and its uncertainty.
9. Check exact source manuscript/table precision and harmonize percentage points versus balanced-accuracy units.

## Submission metadata

Author names/order; affiliations; corresponding author; CRediT roles; funding; competing interests; ethics/institutional determination; public repository/archive identifiers; data/model licenses; target journal. None should be supplied by inference from personal background.

## Source traceability

Primary endpoint and locked test: [file:432][file:469]. Development and pilot readouts: [file:466][file:470][file:444]. Hard panels: [file:464][file:433]. Robustness: [file:468][file:435]. Layers: [file:467][file:430][file:453]. Region selection: [file:434][file:436][file:431]. Study counts and method table: [file:441][file:451].

Literature records used: Johnson et al. [web:527]; DNABERT-2 [web:506]; DNABERT-S [web:331]; HyenaDNA [web:497]; Nucleotide Transformer [web:517][web:518]; DART-Eval [web:516]; V-Xtractor [web:548][web:549]; mopo16S [web:332]; Qscore [web:540]; equivalence testing [web:487]; random projection [web:489].

The BibTeX library contains verified identifiers and titles. Entries using 'and others' deliberately avoid inventing unchecked coauthor lists; import the verified DOI in a reference manager for the complete publisher metadata. The manuscript is not final typeset copy until its author-date citations are linked to that library and the journal style is applied.
