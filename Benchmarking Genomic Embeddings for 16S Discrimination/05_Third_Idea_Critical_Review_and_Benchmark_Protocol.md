# Protocol 3: Benchmarking pretrained genomic representations for fine-scale 16S discrimination

Version 1.0 | Prepared 7 November 2026 | Prospective implementation guide

## Executive assessment

This is a feasible and potentially valuable comparative evaluation, not a conventional literature review. Adding models creates primary computational results. A systematic literature search can accompany the work, but does not convert the experimental benchmark into a review article.

Do not begin from the conclusion that existing pretrained models do not capture the signal. The defensible research question is whether, under specified representations and readouts, pretrained models preserve information useful for fine-scale discrimination and ranking of 16S regions. Failure of mean-pooled cosine distances is not proof that the encoder contains no relevant information. DART-Eval and recent representation analyses illustrate the importance of readout choice and the distinction between frozen probes and adaptation. [web:393][web:390]

This third idea substantially overlaps the scientific core of Protocol 1. It can be the first, narrower manuscript: evaluate representation utility before developing a recommender. Protocol 1 becomes downstream assay/software translation; Protocol 2 remains adaptation. Do not split the same core model-by-panel experiments into multiple nominally separate papers without distinct questions, data, analyses, and transparent cross-referencing.

The quoted claims about an existing pipeline, above-ceiling results, fixed leakage, and verified model environments are user-provided descriptions, not independently audited findings. Audit scripts, manifests, logs, and immutable outputs before relying on them in a manuscript.

## 1. Working title and central questions

Provisional neutral title:

Do pretrained genomic representations preserve fine-scale 16S discriminatory information? A sequence-ambiguity-aware benchmark across organism panels and readout strategies.

Alternate title:

Biological information limits and representation accessibility in 16S-based taxonomic discrimination.

Primary question: On biologically independent query genomes, how does the discrimination accessible from frozen genomic representations compare with strong sequence-based baselines?

Secondary question: Do those representations correctly rank candidate regions by held-out discriminatory performance?

Mechanistic question: If performance is poor, is the main limitation observed-sequence ambiguity, raw embedding geometry, pooling/readout, support-reference density, or task/pretraining mismatch?

This protocol does not require fine-tuning or a deployed recommender. A small readout diagnostic is in scope; model adaptation remains a separate hypothesis-driven extension.

## 2. Literature context and contribution

### 2.1 Biological limits

Johnson et al. compare full-length and subregion 16S resolution and emphasize intragenomic copy variation. This supports explicit biological information controls and prevents equating ASV/haplotype count with species count. [web:297]

### 2.2 Relevant prior representation work

DNABERT-S directly targets species-aware DNA embeddings and already compares general models with tetra-nucleotide features. It is essential prior art and an important positive comparator. Its published training uses 10,000-bp sequence pairs; transfer to short conserved 16S observations is not guaranteed. [web:331]

Micro16S learns 16S phylogenetic embeddings. Its six downstream benchmark tasks involve microbiome sample classification, not six demonstrations of failure on target-specific region discrimination. It separately evaluates taxonomy, where lower-rank performance is limited. Do not conflate these tasks. The reviewed March 2026 version is a preprint. [web:318]

DART-Eval compares zero-shot, probed, and fine-tuned settings on regulatory DNA, finding task-dependent limitations and strong non-foundation-model baselines. It also shows that pooled embedding distance can miss information accessible through other model uses. This is methodological precedent, not direct evidence about 16S. [web:393]

A 2026 representation-analysis preprint reports that pooling/readout can obscure local biological information. It supports appropriately narrow interpretation of a negative embedding benchmark, while not proving any particular model works on 16S. [web:390]

### 2.3 Distinct contribution to pursue

Combine fine-scale target-neighbor panels, copy-aware observed-sequence ambiguity, genome/taxon-held-out evaluation, multiple frozen model families, readout diagnostics, and honest null/equivalence interpretation.

Neither comparing genomic models nor reporting a classical-method win is automatically novel. Demonstrate why the fine-scale marker task and biological controls expose a meaningful, reproducible gap not answered by broad species classification.

## 3. Task definitions and scope

Task T1, primary: closed-panel target-versus-neighbor classification of one observed 16S sequence. A labeled support library is available; hidden queries come from independent genome/isolate groups.

Task T2, secondary: rank candidate regions or amplicons according to the classification performance they support on independent query genomes. No query-derived performance may enter the ranking procedure.

Task T3, diagnostic: evaluate whether raw embedding geometry reflects within-versus-between-species separation, conditional on the task and available reference variation.

Optional T4: open-panel rejection with genuinely unlisted taxa. Only include if a separate unknown-taxon set and calibration policy are available.

Not primary: community abundance prediction, pathogenicity, clinical identification, strain-level attribution, or de novo foundation-model training.

Species and subspecies experiments must be separate strata. A Bifidobacterium longum subspecies panel answers a subspecies question and may face different annotation reliability; do not count it as an independent species-level replication without explanation.

A Bacillus anthracis-containing panel, if used, is a public-sequence biological stress test, not an assertion that all its 16S sequences/regions are identical or that the panel is categorically null. Measure its empirical ambiguity under the selected taxonomy and reference snapshot. No wet-laboratory handling of high-risk organisms is required or proposed.

## 4. Hypotheses and claim hierarchy

H1: frozen representation performance is model-, region-, panel-, and readout-dependent.
H2: some model/readout combinations underperform alignment or k-mer baselines even where the observed sequence has demonstrable held-out discriminating value.
H3: geometry-only failure is not equivalent to failure under a trained probe.
H4: some hard panels are limited chiefly by shared observed sequences, not by representation quality.
H5: performance on species classification does not necessarily imply accurate ranking of regions.

Claim ladder:
1. No detectable separation under the tested geometry/statistic.
2. No practically useful discrimination under the tested frozen readouts, with adequately narrow intervals.
3. Underperformance relative to the specified classical reference.
4. No added value beyond classical features under a defined comparison.
5. Absence of signal anywhere in the model.

The first four may be testable within this protocol. The fifth is not supported by a finite collection of pooled embeddings and probes. Use model/readout/task-specific language.

## 5. Audit existing work before expansion

### 5.1 Preserve history

Archive the current code commit, environment, initial outputs, and patch log. Recompute all reported results after fixes. Keep pre-fix results as engineering diagnostics, not empirical evidence for model performance.

### 5.2 Audit checklist

- How many independent genomes, isolates, full-length genes, copies, unique region sequences, and labels exist?
- Are raw sequences or haplotypes collapsed? Which weights survive collapse?
- Can the same observed sequence map to multiple species/subspecies?
- Does nearest-neighbor evaluation include the query itself, another copy from its parent genome, or an identical technical duplicate?
- What is the tie policy when multiple species are equally near?
- Do batches, normalization, pooling, file order, or model mode change embeddings?
- Does any feature derive from labels, target-specific batching, filenames, or test-set summary statistics?
- Do all models observe precisely the same nucleotide bases?
- Are the purported environment checks only imports, or successful tokenization and inference on realistic inputs?

Passing an import test is not model readiness. Archive an inference smoke test with sequence lengths, output shape, masking behavior, finite embeddings, deterministic settings, and peak resource use.

## 6. Reference data and panel design

### 6.1 Data hierarchy

Use genome-backed full-length 16S with stable assembly/sequence accessions and isolate/copy relationships. Auxiliary curated 16S references can expand diversity but must be labeled when genome linkage is unavailable. Freeze the taxonomy and source snapshot.

Required fields: assembly accession/version, isolate group, source/release, species and subspecies labels, taxonomy system/version, copy ID, coordinates/strand, sequence checksum, QC, inclusion/exclusion reason, and source study.

Retain conflicting labels in an audit. Exclude uncertain labels from confirmatory supervised estimates according to a rule fixed without inspecting model performance.

### 6.2 Panel classes

Assemble panels prospectively across an observed-sequence separability spectrum:
A. informative: sequence baselines and ambiguity audit suggest useful discriminatory information;
B. partial: some strain/species pairs share observations while others differ;
C. highly ambiguous: substantial cross-label identical-sequence overlap;
D. low-reference: deliberately evaluate limited support, separately from well-supported panels.

Select panels using source diversity, taxonomy, and sequence criteria, not whichever panels make embeddings look weak. A positive result must be as publishable as a negative result.

Pilot: 3–5 blocks can validate engineering and detect major errors. Main planning target: approximately 15–30 independent taxonomic blocks if feasible; more may be needed for heterogeneous effects. These are pragmatic scope ranges, not a completed sample-size calculation. A smaller 3–5-block study may be a useful focused case series but cannot support universal model claims.

Aim for enough independent genomes per label to form meaningful support/query splits, e.g. at least 5 and preferably 10+, while acknowledging that even 10 may not represent population diversity. Genome counts, not amplicon expansion or repeated reads, define biological support.

### 6.3 Independence and overlap

If several target panels share species or genomes, connect them into one taxonomic block. Broader related genera may be grouped for stringent inference. Report the number of blocks after grouping, not merely the number of target-versus-neighbor panels.

All copies, regions, read realizations, and repeated assemblies of a genome/isolate remain linked. These are repeated observations, not independent replicates.

## 7. Region and observed-input construction

Perform a two-level evaluation:
- Intrinsic sequence benchmark: validated variable-region windows and full-length reference sequence.
- Practical candidate benchmark: primer-defined inserts under explicit read coverage, trimming, and compatibility rules.

If resources are limited, complete the intrinsic benchmark first and state clearly that it does not establish laboratory assay utility.

Use validated alignment/HMM extraction rather than identical raw nucleotide coordinates across taxa. V-Xtractor is an established starting point for region extraction. [web:282]

Candidate set may include V1–V9 where boundaries are valid, plus common multi-region windows and full-length sequences. Do not imply that every theoretical variable region has an equivalent feasible primer assay.

Primer sequences must be verified from source publications; clip oligo-derived bases before scoring. Do not supply an unsequenced middle segment to a model processing paired reads. Track extraction failure and unsupported inputs explicitly.

Fit neither token truncation nor length-based exclusions to favor one model. Report performance on the common valid input set alongside coverage/failure results in the intended universe.

## 8. Observed-sequence ambiguity and ceilings

### 8.1 Definition

A sequence type/haplotype is one unique observed insert sequence under a declared strand policy. It is not equivalent to one genome, strain, species, or independently sampled unit.

For each sequence s and class y aggregate declared biological weights w(s,y). For the fixed empirical population, the maximum accuracy of a deterministic classifier whose only input is s is:

A_ceiling = sum_s max_y w(s,y) / sum_s sum_y w(s,y).

For species-balanced or target-balanced accuracy, encode the corresponding weights before applying the formula. Do not compare a raw-frequency ceiling with a balanced metric. The formula is for accuracy, not a universal ceiling for AUROC, macro-F1, or distance-separation statistics.

### 8.2 Interpretation

It is a finite-reference, weight-dependent, observed-sequence ceiling. It changes when sequence coverage, class priors, labels, or available contextual inputs change. A full-length-template ceiling is not automatically the correct ceiling for a shorter observed insert, or vice versa.

For a deterministic sequence-only predictor evaluated on exactly the same weighted population, surpassing the matched ceiling is a red flag. Possible causes include metric mismatch, weighting errors, query-specific ties/self-matches, stochastic record-specific predictions, leakage, wrong observed input, or an incorrect ceiling implementation. Exceeding it does not identify one particular bug by itself.

### 8.3 Collapsed and uncollapsed analyses

Primary: genome-aware evaluation with copy information and weights retained.
Sensitivity: unique-haplotype analysis with explicit weights and multi-label membership.
Sensitivity: same-input ambiguity subset, reported independently.

Do not collapse a haplotype to a single arbitrary species or majority label and then describe the result as unbiased species classification. Conversely, duplicating a shared string across biological genomes does not create new sequence information; the dependence must be handled.

### 8.4 Positive information controls

Include full-length observations, direct sequence matching/alignment, and a simulated fixture with known diagnostic bases. A real informative panel should demonstrate held-out baseline performance, not merely have unique training sequences.

Synthetic fixtures test implementation and sensitivity; they do not establish biological generalization or replace natural panels.

## 9. Model selection and inclusion rules

Primary feasible core:
1. DNABERT-2: general multi-species backbone.
2. DNABERT-S: species-aware comparator.
3. HyenaDNA: representative sequence architecture.
4. Nucleotide Transformer v2: representative transformer/tokenization family.

Optional, contingent on resources and checkpoint suitability: Caduceus; one Evo-family checkpoint; a compatible Micro16S checkpoint as a specialized non-generic representation control. Models must have usable code, licensed weights, supported input context, and a successful smoke test. Availability and environment readiness must be verified rather than assumed from the quoted proposal.

Use a family-representative selection rule fixed before main testing. An example is the smallest documented checkpoint capable of the input lengths, plus one size-sensitivity comparison. Do not select each checkpoint by test accuracy. Different parameter counts, training corpora, and architectures mean the study is a practical comparison, not a clean causal architecture experiment.

Record checkpoint revision/hash, tokenizer revision, license, pretraining scope, model size, context length, precision, pooling, orientation, and dependency versions. Treat undocumented pretraining overlap as unknown.

## 10. Readout hierarchy: essential for interpreting negative findings

R0, intrinsic geometry: masked mean-pooled embeddings and a prespecified distance (e.g. cosine), no supervised encoder updates.
R1, reference lookup: nearest-neighbor classification using the same frozen embeddings, explicit support references, and a fixed tie policy.
R2, linear probe: regularized logistic regression or equivalent trained on support-training data only.
R3, bounded nonlinear probe: a small MLP with fixed development-selected capacity and regularization.
R4, diagnostic readout sensitivity: a limited alternative layer/pooling or token-level local readout, where architecture permits, selected only on development data.

Primary readout can be R2, with R0/R1 testing natural geometry. R3/R4 are diagnostic checks before concluding broad representation failure. A model with poor geometry but good linear-probe performance has accessible discriminating information under supervision; a model requiring a nonlinear head supports a different conclusion.

A frozen encoder plus trained probe is not zero-shot learning. Nearest-reference classification uses labeled support; label the supervision level precisely. Likelihood/pseudo-likelihood may be an exploratory diagnostic, but tokenization/objective differences prevent naive direct comparison of raw likelihoods across models.

For every model: test padding invariance, batch order, batch size within numerical tolerance, deterministic inference, strand handling, and same-input consistency. Never batch by the true query label. Fit any scaler/PCA/normalization on support-training only.

## 11. Strong classical controls

Mandatory controls:
- prevalence/majority and balanced-chance references, interpreted under task weights;
- exact sequence lookup with shared-label ambiguity;
- global/local alignment or edit distance appropriate to the extracted inputs;
- canonical k-mer composition with a small development-selected k grid;
- k-mer + the same regularized linear probe used for embeddings;
- a simple diagnostic-position/alignment-feature classifier where valid.

Length-only and GC-only classifiers are shortcut diagnostics, not the main classical competitor. A random embedding with matched dimension can reveal classifier flexibility and chance structure; it is not a biological baseline.

Use comparable support, splits, and reasonable development budgets. No embedding/classical comparison may use a different query population without showing the effect of that difference.

## 12. Experimental splitting

### 12.1 Primary biological holdout

Partition independent blocks into development and locked evaluation groups. Within every panel partition isolate/genome groups into support and hidden query. Copies/regions from the same genome cannot cross roles.

Global hyperparameters/readout recipes are selected on development blocks. On a locked panel, local support labels may train the fixed probe recipe. Query labels remain inaccessible until predictions are frozen.

### 12.2 Different novelty regimes

Report separately:
A. new genome, possibly known observed haplotype: realistic support-library usage;
B. new observed haplotype: stricter sequence-novelty challenge;
C. globally new taxonomic block with local support: representation-transfer test.

Exact overlap between distinct genomes is biologically real and should not be erased from all analyses. Report it; use sequence-novelty analyses as an additional selected subset with its changed population stated.

### 12.3 Pretraining contamination

Downstream holdout does not prove absence from pretraining. Audit documented sources and sequence overlap where feasible, but record uncertainty when the corpus is not available. A date cutoff alone does not guarantee a novel sequence.

## 13. Metrics and primary contrasts

Primary metric: genome-weighted, species-balanced target-versus-neighbor accuracy on hidden queries. For each genome average its copy/read contributions under a fixed policy; average genomes within labels; weight target and pooled-neighbor halves equally, with equal neighbor-species weights inside the latter.

Primary estimand: mean paired difference between a prespecified frozen-model strategy and a prespecified strong classical strategy across independent blocks. Define whether the strategy is a named model or a development-selected best-model rule. A test-selected best model is descriptive only.

Secondary: macro-F1, per-neighbor sensitivity/specificity, AUROC/AUPRC with prevalence, calibration if probabilities exist, ceiling gap, rank stability, computational cost, and open-set rejection if separately designed.

Ceiling gap: A_ceiling minus observed accuracy, using matched weighting and input. It describes unrealized performance under that empirical population; it is not an estimate of hidden information anywhere in the model. Avoid dividing by tiny ceiling-minus-chance denominators near unresolvable cases.

Geometry metric: within-versus-between-label distance contrast or pairwise separation AUC with declared pair weighting. Pairs are dependent. Infer at genome/block level, not with a naive test treating n(n-1)/2 distances as independent observations.

Region-ranking metric: held-out accuracy of support-selected region; regret relative to the evaluation-only highest-performing feasible region; rank correlation with uncertainty/ties. Geometry-derived region rankers must be fixed before query inspection.

## 14. Permutation design: do not assume sequence-level shuffling is universally correct

The valid randomization unit depends on exchangeability and the null hypothesis. A label-shuffling method can be wrong at either the haplotype or sequence-row level. Repeated copies/regions of a genome, replicate reads, and related isolates are not automatically exchangeable. Permutation validity requires preserving the dependence relevant to the null. [web:392][web:394]

### 14.1 Define the null before choosing the shuffle

Possible diagnostic null: among the allowed biological units within specified strata, features have no association with the assigned target labels. This is a conditional computational reference, not a literal biological model in which species labels are randomized in nature.

If the unit is a genome/isolate, permute one label per allowed unit and propagate it consistently to all copies and regions. Preserve class counts within the declared exchangeability strata. If labels cannot be exchanged within the plausible strata, a meaningful unrestricted permutation test may not be available; do not force one.

If unique haplotypes with class-frequency weights are the actual analysis unit, derive the null for that weighted dataset rather than substituting an unrelated row shuffle. Neither level is intrinsically always correct.

### 14.2 Fit/select inside each permutation when the statistic requires it

If the statistic includes fitting a probe or selecting the best region/model, the permutation reference must repeat the corresponding training/selection steps under the null. Otherwise it does not represent the same statistic and can be anti-conservative.

A simple geometry statistic with a completely frozen encoder can often reuse label-independent distance calculations, but must recompute label-dependent aggregation/weights consistently. Keep one genome-level label assignment consistent across regions in a permutation when testing a shared region family.

A practical Monte Carlo p-value is (1 + number of null statistics at least as extreme as observed)/(B + 1). Choose B to resolve the relevant corrected significance level; 999 permutations may be too coarse for a large family. Start with development profiling, then declare the final budget.

### 14.3 Recommendation

Use effect sizes and block-aware confidence intervals as primary evidence. Keep permutation analyses as well-defined diagnostics, not as the sole basis for 'signal exists/does not exist.' Consult a statistician before making taxonomy-related exchangeability claims across strongly related genomes.

## 15. 'No signal' versus inconclusive versus equivalence

A non-significant p-value means insufficient evidence under that test, not proof of absent signal. Classify results using performance estimates and uncertainty:

- Supported useful discrimination: performance exceeds a declared practical criterion with an adequately narrow interval.
- Evidence of practically negligible improvement: the uncertainty interval excludes gains larger than a prespecified meaningful margin under the stated inferential design.
- Inconclusive: interval spans both meaningful benefit and little/no benefit.
- Information-limited: sequence ambiguity materially limits all sequence-only methods under the evaluated population.
- Implementation failure: missing/invalid outputs; never encode as a biological no-signal cell.

For formal equivalence, prespecify margins and a valid equivalence procedure. Example planning margin: plus/minus 0.03 absolute balanced-accuracy units, chosen with biological/application justification before final testing. This is not a universal accepted threshold. Equivalence requires an appropriate interval/test, not both one-sided superiority tests being non-significant.

## 16. Statistical analysis and multiplicity

Block-level paired contrasts are the primary unit. Bootstrap independent blocks, averaging nested panels and technical seeds first; optionally resample genomes inside blocks to reflect within-block sampling. Provide 95% confidence intervals and individual-block plots.

Do not t-test correlated folds, count pairwise distances as independent, or use the number of model-region cells as biological n.

Multiplicity: declare the entire confirmatory comparison family. Correcting across regions alone is insufficient when also testing many models, panels, metrics, layers, pooling rules, and readouts. Options include one primary pooled contrast and Holm correction for a limited secondary family, or a prespecified hierarchical testing scheme. Mark exploratory grids as exploratory.

Use shared bootstrap/permutation draws where appropriate to preserve paired structure. If reporting simultaneous cell-level claims, use a declared familywise strategy or clearly separate FDR discovery claims from confirmatory evidence.

Plan sample size from development-only block-level paired variation and scientifically meaningful effects/interval widths. Simulate prospective uncertainty across feasible block counts rather than claim that adding more models increases power. More models increase breadth and multiplicity, not independent organism replication.

## 17. Failure-mode diagnostics

### 17.1 Local diagnostic variation

Use training/support-derived diagnostic positions to examine natural allele pairs and carefully designed counterfactual substitutions. Record whether a change is naturally observed or synthetic. Quantify embedding/probe sensitivity, controlling for length and position where possible.

Do not choose sites using test labels and then present the probe as independent discovery. A perturbation response is mechanistic evidence only; it does not replace performance on real held-out genomes.

### 17.2 Pooling and tokenization

Compare one primary pooling rule with a limited development-selected alternative; map token representation back to bases where possible. Single-base changes can alter byte-pair tokenization beyond one token, so token-index attribution is not a direct nucleotide-causal explanation. DART-Eval discusses this complication. [web:393]

### 17.3 Support density and imbalance

Evaluate controlled support budgets and equal-class weighting. Learning curves can show whether poor performance reflects limited references rather than universally uninformative embeddings. Preserve biological grouping during subsampling.

### 17.4 Batch and metadata controls

The same sequence should not gain taxonomic information from its batch companions, filenames, order, or query-label-dependent preprocessing. Inference-mode normalization and deterministic pooling should be tested with randomized batching. Numerical variation is not automatically label leakage; document actual mechanism and magnitude.

### 17.5 Errors and strand

Add a bounded secondary error/strand robustness experiment after the clean-input benchmark is stable. Uniform synthetic errors are stress tests, not complete PCR/sequencing realism. Full assay validation belongs mainly to the regioneer translation project.

## 18. Execution plan and stopping rules

Phase 1, audit: archive code and fix ceiling/tie/weight issues; generate synthetic fixtures and current-panel summaries.
Phase 2, pilot: DNABERT-2, HyenaDNA, NTv2, DNABERT-S if runnable; classical controls; 3–5 development blocks; primary and linear-probe readouts.
Phase 3, freeze: define final model/checkpoint list, panel eligibility, primary metric, readout selection rule, multiplicity plan, practical margins, and compute budget.
Phase 4, scale: run the locked model-by-panel-by-region matrix on audited data; cache label-independent embeddings by exact settings.
Phase 5, evaluate: freeze predictions, unseal labels, execute prespecified summaries, and document missing cells.
Phase 6, diagnose/report: perform clearly labeled exploratory layer/readout analyses, prepare reproducible release, and write conclusions conditional on results.

Stopping rules must distinguish:
- Invalid data or leakage: stop that analysis until repaired.
- Biological ambiguity: report it; do not chase impossible separation.
- Resource incompatibility: mark technical exclusion with reasons.
- Null result: continue the prespecified benchmark, not indefinite model hunting.
- Repeated test-guided model selection: acquire a new test or label analysis exploratory.

A kill-switch used to stop investment in one deployment method does not have to stop a planned comparative benchmark. However, there should still be a fixed resource/model budget and data-quality gates. Keep null cells; do not expand until a desired narrative emerges.

## 19. Data/output schemas

Reference table: genome/isolate ID; accession; taxonomy label/rank; copy ID; region/assay; observed string/hash; source; QC.
Model table: family; checkpoint hash; tokenizer; pooling/layer; strand; model mode; input token coverage; license.
Split table: block; genome group; global development/test; local support/query; fold/seed.
Prediction table: model/readout; panel/block; region; query genome/copy; truth; prediction; probability if valid; failure status.
Metric table: weighting; metric; estimate; interval; eligible n_genomes/n_blocks; ceiling where defined; multiplicity family.
Audit table: bug/limitation; mechanism; affected outputs; fix commit; rerun verification.

Keep biological no-signal, insufficient data, technical failure, excluded-by-protocol, and not-run-for-resources as distinct codes. A heatmap must not render all of them as zero.

## 20. Suggested figures and manuscript structure

Figure 1: biological hierarchy, input observations, and independent splits.
Figure 2: shared-sequence ambiguity and matched accuracy ceilings by panel/region.
Figure 3: model-by-panel performance matrix, explicitly coding missing/technical cells.
Figure 4: paired block-level comparisons against alignment/k-mer controls.
Figure 5: geometry versus linear/nonlinear readout accessibility.
Figure 6: region-ranking accuracy/regret and model/classical agreement.
Figure 7: support-density, local variant, and resource trade-off diagnostics.

Tables: panel provenance and rank; models/checkpoints; readout/hyperparameter budget; primary contrasts and intervals; failure mechanisms and corrected audit tests; limitations.

Manuscript sections: biological motivation; related work; data/splits; ambiguity framework; model/readout comparison; region ranking; diagnostics; computational cost; limitations; reproducibility. A systematic literature-search appendix can be included if actually performed and documented.

## 21. Ethical and computational considerations

Use public/authorized data, respect licenses, avoid unnecessary host sequence, and do not upload controlled data to hosted APIs. This is not a clinical diagnostic validation. Public in-silico analysis of a challenging pathogenic taxon does not require laboratory culture; any future wet-lab work requires institutional biosafety review and is outside this protocol.

Measure actual runtime/memory before selecting expensive optional model families. Prefer an executable four-family benchmark with strong controls and biological diversity to an incomplete ten-family sweep. GPUs, seeds, sequence hashes, and environment locks must be recorded.

## 22. Decision table for the final interpretation

| Result | Supported interpretation | Unsupported extension |
|---|---|---|
| Raw geometry weak; probe strong | Useful information is accessible with the tested supervised readout | The pretrained model lacks the signal |
| Embeddings weak; alignment/k-mers strong on independent queries | Tested frozen representations/readouts fail to use available discrimination effectively | Every genomic model fails universally |
| All methods poor; shared sequences frequent | Observed marker ambiguity limits the task in this reference population | Fine-tuning will recover absent sequence information |
| All methods poor; classical controls not validated | Study may be underpowered, mis-specified, or information-limited; investigate | Embeddings uniquely fail |
| Embeddings and classical methods both strong | Representations preserve useful information; compare cost and ranking | Expensive models are necessary |
| Model improves classification but not region ranking | Classification benefit does not translate to selection benefit here | A superior region recommender was established |
| Intervals wide | Inconclusive | Evidence of no signal |

## 23. Annotated source reading list

1. Johnson et al. Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome analysis. Nature Communications (2019). DOI 10.1038/s41467-019-13036-1. Biological resolution/copy variation anchor. [web:297]
2. Zhou et al. DNABERT-S: pioneering species differentiation with species-aware DNA embeddings. Bioinformatics 41(Suppl 1), i255–i264 (2025). DOI 10.1093/bioinformatics/btaf188. Species-aware prior art, model comparator, long-fragment-to-short-amplicon transfer question. [web:331]
3. Bishop et al. Micro16S: Universal Phylogenetic 16S rRNA Gene Representations for Deep Learning of the Microbiome. bioRxiv (reviewed v1, 2026). DOI 10.64898/2026.03.21.713432. Learned 16S phylogenetic representations; distinguish taxonomic classification from six sample-classification tasks. [web:318]
4. Patel, Singhal, Wang et al. DART-Eval: A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA. Reviewed arXiv version 2412.05430v2; verify final venue/version before manuscript submission. Read for strong baselines, task alignment, geometry/probing distinctions, and variant sensitivity. [web:393]
5. A Representation Analysis of Genomic Language Models. Reviewed arXiv record 2608.05329 (2026); preprint. Read for layer/pooling/accessibility limits, not direct 16S evidence. [web:390]
6. The Exchangeability Assumption for Permutation Tests of Multiple…. Reviewed arXiv record 2406.07756v2. Read for dependence and valid randomization; export complete bibliographic title before submission. [web:392]
7. Permutation tests for experimental data. PMC record PMC10066020 (2023). Read for stratification and matching the shuffle to the actual null. [web:394]
8. Benchmarking taxonomic assignments based on 16S rRNA gene profiling of the microbiota from commonly sampled environments. PubMed record 29762668 (2018). Related taxonomic benchmark and ecological reference dependence. [web:388]
9. Comprehensive Assessment of 16S rRNA Gene Amplicon Sequencing for Microbiome Profiling across Multiple Habitats. Microbiology Spectrum (2023). DOI 10.1128/spectrum.00563-23. Qscore prior art for any downstream region-selection tool claim. [web:300]
10. V-Xtractor: An open-source, high-throughput software tool to identify and extract hypervariable regions of small subunit (16S/18S) ribosomal RNA gene sequences (2010). Established extraction precedent. [web:282]

This reading list is targeted, not exhaustive. The earlier shared annotated review supplies additional database, primer, simulation, and adaptation readings. Verify current publication status, exact checkpoint releases, and bibliographic exports before preregistration/submission.

## 24. Final recommendation

Proceed with a preregistered frozen-representation benchmark, not a paper premised on proving that models fail. Make biological ambiguity, meaningful positive controls, readout hierarchy, independent genomes, and careful uncertainty central. Treat regioneer as the subsequent translation and custom adaptation as a contingent extension. The strongest outcome is an accurate map of where representations do and do not help—not a uniformly negative heatmap.
