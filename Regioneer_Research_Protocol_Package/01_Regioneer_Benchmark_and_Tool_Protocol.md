# Protocol 1: Regioneer benchmark and target-conditioned 16S amplicon recommender

Version 1.0 | Prospective protocol | 7 November 2026

## 1. Title and study design

Working title: Target-conditioned selection of 16S rRNA amplicons for species discrimination: a leakage-resistant benchmark of classical sequence methods and frozen genomic representations.

Design: computational comparative-methods study, followed by independent read-level validation and optional prospective laboratory validation. Software development is a translational deliverable, not a substitute for validation.

Unit of scientific inference: independent organism groups or taxonomic blocks. Genomes/isolates are biological sampling units within groups. Copies, regions, reads, pairwise distances, random seeds, and cross-validation folds are not independent biological replicates.

## 2. Rationale and literature context

16S contains region-dependent discriminatory information; full-length sequences can improve resolution, but intragenomic copies complicate the mapping of sequence variants to organisms. Johnson et al. directly examined full-length versus subregion resolution and copy variation. Do not generalize their specific clustering thresholds as universal species boundaries. [web:297]

Primer sequence and amplification compatibility affect which organisms are observed. Klindworth et al. evaluated 175 primers and 512 pairs, while mopo16S explicitly optimized efficiency, coverage, and bias. This argues for assay-level recommendations rather than rankings of abstract region labels alone. [web:360][web:332]

Qscore already integrates amplification, annotation, and cost and offers taxa-/habitat-informed configuration recommendations. Regioneer must be compared with this conceptual and, if executable, software baseline. Its contribution cannot be the mere existence of an online region selector. [web:300]

Genomic representations are plausible candidates, not guaranteed improvements. DNABERT-S is species-aware, whereas Micro16S optimizes phylogenetic representation and partial region invariance. These objectives differ from choosing an amplicon with high local discriminatory information. [web:331][web:318]

## 3. Research questions and hypotheses

Primary question: Across independent organism blocks, does a prespecified recommendation strategy improve the downstream discrimination achieved by its selected feasible amplicon relative to a strong classical strategy?

Secondary questions:
1. How often is target-versus-neighbor ambiguity unavoidable under the observed sequence and reference population?
2. Are rankings stable across isolates, reference databases, primer mismatch rules, sequencing errors, and candidate-neighbor definitions?
3. When does full-length 16S materially outperform feasible short amplicons?
4. Does a hybrid representation improve performance enough to justify runtime and operational complexity?
5. Can uncertainty and an explicit abstention policy prevent unsupported recommendations?

Hypotheses:
- H1: the preferred assay varies across target panels rather than following one universal region ranking.
- H2: performance estimates deteriorate under genome-held-out or taxonomic-block-held-out evaluation compared with naive sequence-row splits.
- H3: frozen embeddings may improve some panels, but superiority over strong classical methods is not assumed.
- H4: amplification constraints can change the preferred assay relative to a sequence-only ranking.

H1/H2 can be supported even if H3 fails. All hypotheses and primary comparisons must be locked before independent evaluation.

## 4. Intended use and explicit non-goals

Primary intended use: planning a research assay to distinguish a target species from a declared neighbor panel using a supported set of primer-defined amplicons.

Distinct use cases to label:
- Closed-panel sequence assignment: the query is known to come from one listed species.
- Open-panel sequence assignment: queries may come from unlisted taxa; rejection is evaluated.
- Community presence/absence: the target may be absent in a mixed sample; evaluate sample-level detection separately.

Non-goals for first release: clinical diagnosis; pathogenicity inference; absolute abundance; universal primer design from scratch; strain identification; genome-derived functional prediction; prediction of species from identical observed sequences.

## 5. Operational definitions

Target: a species identifier under a locked taxonomy, not merely a free-text name.
Near neighbor: a species selected by a preregistered combination of genome-based relatedness, 16S similarity, and application relevance. Maintain both biological and sequence-confusability definitions.
Panel: one target plus its declared neighbor set and associated reference sequences.
Taxonomic block: a connected component of overlapping panels, or a higher taxonomic grouping, used to prevent shared organisms leaking across panel splits.
Region: a positional/hypervariable designation used for mechanistic analysis.
Assay: primer pair + strand/orientation policy + trimming + read layout + usable read length + downstream processing.
Observed sequence: bases actually available after primer removal and processing, not the full template used to simulate them.
Recommendation: an assay ranked under constraints, with evidence and uncertainty, not an unconditional biological optimum.
Unresolvable: no candidate meets a prespecified criterion under the evaluated reference and measurement model. This is not proof that every possible marker or future assay fails.

## 6. Data sources and governance

### 6.1 Reference tiers

Tier A: complete or sufficiently well-validated isolate genomes with stable assembly accessions, linked taxonomy, and extracted 16S copies. Prefer complete genomes for copy-complete analyses; do not equate absence of assembled rRNA in a draft/MAG with biological absence.
Tier B: curated full-length reference 16S sequences with reliable species attribution; useful for diversity coverage but frequently lack strain/copy provenance.
Tier C: independent empirical isolate/mock-community read datasets with independently established composition or isolate identity.

Candidate sources: NCBI assembly/reference records, SILVA, a chosen GTDB release, and rrnDB for copy-number contextual information. These are not interchangeable taxonomic authorities; preserve each source label. SILVA supports quality-checked rRNA reference use; rrnDB documents copy-number variability. [web:368][web:359]

### 6.2 Required manifest

For each genome: accession/version; BioSample/isolate link when available; assembly level; source; retrieval date; declared species/taxid; reconciled taxon; reference taxonomy version; study/source; assembly quality flags; checksum; inclusion reason.
For each 16S: parent genome; copy ID; coordinates/strand; extractor/version; gene coverage; ambiguity count; canonical sequence hash; taxonomy conflict flag.
For each external dataset: accession; composition ground truth; primer sequences; chemistry; read length; extraction/library protocol if available; negatives; independent biological replicates; license/access conditions.

### 6.3 Taxonomy reconciliation

Maintain a versioned crosswalk between NCBI, GTDB, SILVA, and names used in manuscripts. Preserve synonyms and conflicting assignments. Never resolve disagreements by whichever label improves model accuracy. Use genome-derived species labels for the primary genome-backed set, selected before analysis; run a sensitivity analysis using an alternative defensible taxonomy.

Ambiguous labels are excluded from species-supervised primary analyses but retained in an explicit ambiguity audit. Sequence labels cannot override taxonomic uncertainty. Manual curation decisions require an audit record and ideally a second reviewer.

### 6.4 Sequence QC

Proposed starting screen for a full-length primary set: approximately 1,200–1,800 nt, no unresolved bases in the evaluated observed sequence, consistent orientation, credible SSU annotation, and no unresolved contamination/chimera flag. This length screen is pragmatic, not a universal definition of bacterial 16S; audit taxa excluded by it.

Preserve genuine intragenomic variants. Remove redundant downloaded records without deleting their genome/copy/taxonomy memberships. Separate sequence deduplication from biological provenance deduplication.

Use one locked extraction method for primary processing and a second method/subset to audit boundary accuracy. V-Xtractor supplies a documented region-extraction precedent. Alignments and extraction must not silently substitute E. coli integer coordinates for all taxa. [web:282]

## 7. Organism-panel construction

### 7.1 Panel selection

A first pilot may include 2–4 organism groups relevant to microbiome research. These are engineering panels, not evidence of broad generalization. Candidate groups could include Bifidobacterium, Bacteroides-related taxa, Lactobacillaceae, and one deliberately challenging species complex; choose actual targets only after accession and taxonomy auditing.

For the main study, aim initially to assemble approximately 30–60 non-overlapping taxonomic blocks if feasible, spanning variable sequence divergence, reference density, and expected discriminability. This is a planning target, not a power claim. Fewer defensible blocks are preferable to many pseudo-independent panels. If only 1–2 blocks exist, frame the result as a case study.

For each target select a prespecified number of nearest species, e.g. 5–10, and retain all biologically justified ties. Use a genome-based taxonomy/relatedness criterion where available, then a 16S-confusability criterion as a separate sensitivity definition. Avoid defining neighbors from embeddings under evaluation.

Include easy, intermediate, and intrinsically ambiguous panels. Do not exclude hard panels because methods perform poorly. Cap overrepresented species at a declared genome count, e.g. 50 in pilot analyses, sampling by isolate/source diversity rather than first accession order. Retain sampling weights and repeat cap sensitivities.

### 7.2 Eligibility for inferential evaluation

A pragmatic target is at least 5 distinct genome/isolate groups per species for within-panel assessment, preferably 10 or more. Counts are not guarantees of representativeness. With too few groups to form train/development/test partitions, mark a panel as descriptive only or use it solely in a separately labeled low-reference experiment.

All copies and candidate regions from one genome remain linked. Identical observed sequences with different species labels remain visible as contradictory/ambiguous evidence; never relabel them to make training easier.

## 8. Candidate regions and assays

### 8.1 Two-track analysis

Track R: accurately mapped V1–V9 and selected multi-region windows; diagnoses where information resides.
Track A: real published primer pairs plus supported read configurations; supplies actionable recommendations.

Examples to consider, subject to exact primer validation: V1–V2, V1–V3, V3–V4, V4, V4–V5, V6–V8, and full-length 16S. There is no assumption that all nine single-variable regions have practical validated assays.

Build `primers.tsv` with oligos in 5'-3' orientation, IUPAC degeneracies, original citation, intended domain, expected span, and supported platform. Verify oligos against primary records before use. Primer names alone are insufficient because variants share names.

### 8.2 In-silico PCR

Search both strands, expand/match degeneracies correctly, reverse-complement the reverse primer for template matching, and record every accepted product. Define a primary strict compatibility rule and sensitivity rules allowing e.g. 0, 1, and 2 mismatches per primer with an explicit terminal 3'-end policy. These are matching scenarios, not validated PCR efficiency models.

For each product record binding coordinates, primer mismatches by position, multiple-product count, amplicon length, and whether the observed insert is usable. Remove primer-derived bases before information scoring: the oligo sequence is not independently observed organism evidence.

Handle no product, multiple products, and ambiguous length as explicit outcomes. A missing amplicon is not an empty sequence, and silently dropping it will inflate performance.

### 8.3 Platform/read constraints

For paired-end assays derive the actual overlap after primer removal, quality trimming, and minimum merging requirements. If the insert is not covered, do not supply the unsequenced middle to the model. Either reject the configuration or define a separately validated paired-read representation without fabricating contiguous sequence.

Full-length sequence is an information reference under error-free evaluation, not automatically the practical performance ceiling under different sequencing errors and costs. Report both sequence-information and measurement-level comparisons.

## 9. Leakage-resistant experimental design

### 9.1 Distinguish two learning stages

Global method development: choose representation, hyperparameters, and scoring recipe using development blocks only.
Panel-local deployment: given a new panel, use its reference/support genomes to score assays or fit a lightweight local classifier. This use of support labels is permitted and should be explicit. Query genomes from that panel remain hidden until evaluation.

Thus an unseen-genus experiment can still have a labeled local support set: the genus was absent from global development, not necessarily absent from the deployment reference library. Do not call this zero-shot species identification.

### 9.2 Proposed split hierarchy

1. Group repeated assemblies/records of the same isolate; keep a group in one role.
2. Assign overlapping panels to connected taxonomic blocks.
3. Freeze global train/development/locked test blocks, e.g. approximately 60/20/20 by block, adjusting to available diversity.
4. Within each block, assign support versus hidden query genomes by isolate group, stratified by species if feasible.
5. For development panels, perform inner grouped resampling only within support data to select local thresholds and assay rankings.
6. Apply the frozen strategy to support data in locked test panels, then evaluate selected assays on hidden query genomes.

Report two regimes separately: (i) new isolates of familiar species and (ii) globally unseen taxonomic blocks with local support. Optional family-held-out analysis is a stronger transfer test.

### 9.3 Exact identity and biological collisions

Genome grouping and zero shared amplicon hashes cannot always both be satisfied: distinct genomes may genuinely share an identical 16S amplicon. Do not delete these collisions from the realistic primary evaluation. Report overlap counts, shared hashes, and their label memberships.

Add a strict sequence-novelty challenge set by holding out sequence-hash components or excluding overlap, but label the induced selection bias and loss of ambiguous taxa. For representation training in protocol 2, decontaminate fine-tuning/test overlap more strictly. Never claim zero overlap merely because accessions differ.

Exact/reverse-complement hashes, overlap audits, relatedness measures, and pretraining exposure are separate issues. If pretrained corpora are unknown, label pretraining contamination as uncertain; do not claim fully unseen genomes based on downstream splits alone.

## 10. Methods to benchmark

### 10.1 Mandatory classical baselines

B0: common fixed assay under a declared primer/read configuration; a pragmatic baseline, not a universal standard.
B1: uniform random choice among feasible assays; quantify repeated random choice uncertainty.
B2: exact sequence matching with explicit shared-label ambiguity and novel-sequence fallback/rejection.
B3: nearest-reference global alignment/edit-distance scoring; no forced resolution of tied species.
B4: canonical k-mer frequency features, k in a small prespecified set such as 4, 6, 8; cosine/Jensen–Shannon-type comparisons or a regularized classifier where mathematically appropriate.
B5: per-position diagnostic variation and entropy-based summaries using a valid alignment. Entropy alone is not species discrimination; compare target-versus-neighbor information and within-target variation.
B6: Qscore replication/adaptation where executable; report any difference in objectives, supported assays, cost model, and reference inputs. [web:300]

Do not compare a highly tuned embedding model only against GC content. Use equal split rules and explicit search budgets. Random forests/linear models on composition features are legitimate competitive baselines, not straw men.

### 10.2 Frozen learned representations

Minimum useful set: one general genomic encoder (e.g. DNABERT-2 checkpoint) and species-aware DNABERT-S; optionally a compact Nucleotide Transformer or Micro16S where code/checkpoints/license and input compatibility permit. Do not include a model solely to enlarge a model table. [web:353][web:331][web:318]

Record exact repository, checkpoint revision/hash, tokenizer revision, license, documented pretraining composition, maximum token length, strand policy, pooling, and special-token handling.

Primary pooling: one predeclared method, such as masked mean pooling over nucleotide tokens, excluding padding and special tokens. It is model-specific; verify how returned hidden states map to the actual sequence. Compare final-layer pooling with at most one development-selected alternative. Never silently truncate. Report token counts and coverage.

Primary strand policy: orient by extracted 16S strand. Reverse-complement consistency is a sensitivity test. If ensembling both orientations improves development performance, freeze that policy before testing.

Cache embeddings by sequence hash + model revision + pooling + orientation + tokenization. Identical observed inputs must produce identical deterministic representations under the same settings. Evaluate in inference mode, with no gradients for frozen models.

### 10.3 Common discrimination evaluator

Use a common nearest-neighbor or regularized linear-probe evaluation to compare representations, with all scaling/PCA fitted on support-training data only. Report classifier-specific sensitivity rather than pretending geometry and classifier choice are separable in all circumstances.

For target-versus-neighbor classification, weight the target as one half of the binary task and distribute the non-target half equally across neighbor species. This prevents the number of neighbors or a common species dominating the metric. Also report each neighbor-specific confusion rate.

### 10.4 Candidate ranking rule

Primary recipe: choose the feasible assay with the highest inner-support cross-validated species-balanced binary accuracy. Prespecify tie tolerance, e.g. 0.01 accuracy units, then prefer lower resource cost/shorter usable assay and greater matching coverage according to a locked order.

Alternative label-free geometric rankers may use target separation margins or distance distribution overlap, but must be validated against held-out performance. Their raw scores are not automatically calibrated probabilities.

Prefer a constraint-first/Pareto presentation over an arbitrary weighted sum. Show discrimination, target compatibility, neighbor coverage, ambiguity, and cost separately. If a weighted utility is used, declare its weights before test inspection and repeat sensitivities.

## 11. Endpoints and mathematical definitions

### 11.1 Primary utility

For panel p and assay a define U(p,a) as species-balanced binary accuracy on independent query genomes, with unsupported/failed products counted as incorrect in the end-to-end endpoint. For each genome average across its evaluated copies/read realizations, then weight genomes equally within species. Report TPR and each neighbor TNR separately:

U = 0.5 * mean_target_genomes(correct contribution) + 0.5 * mean_neighbor_species(mean_genomes(correct contribution)).

Conditional sequence classification among successfully observed inserts is a secondary endpoint; do not let it replace the end-to-end result. If coverage is unknown for a reference-only dataset, label U as conditional and do not imply amplification validity.

Each method selects a_hat(p) from support data only. Primary contrast is the block-level difference in U(p,a_hat) versus a prespecified strong classical selector. Multiple target panels per block are averaged before inference.

### 11.2 Recommendation regret

R(p) = max over feasible a of U(p,a) - U(p,a_hat).

The max is an evaluation-only finite-query oracle, not a trainable label derived from test data. It is optimistic/noisy, especially with small query sets. Report bootstrap stability and optionally oracle estimates based on a separate large simulation evaluation pool. Never use this oracle to choose hyperparameters on locked panels.

### 11.3 Ambiguity audit and reference-specific ceiling

For exact observed sequences s and class y, under explicitly chosen genome/species weights w(s,y), the best deterministic sequence-only classifier on that empirical population has accuracy ceiling:

A_ceiling = sum_s max_y w(s,y) / sum_s sum_y w(s,y).

This is an empirical, error-free, distribution-dependent ceiling, not a universal species-resolution bound. For a binary target task, combine all non-target labels only after defining species-balanced weights. Also report the target fraction sharing an observed sequence with any declared neighbor and cross-species collision tables.

No classifier can assign different deterministic labels to identical inputs without additional information. Priors, sample context, or species-specific dropout can change a different task; they do not create nucleotide evidence.

### 11.4 Secondary outcomes

Target sensitivity, species-specific specificity, macro-F1, multiclass accuracy, AUROC, AUPRC with stated prevalence, calibrated Brier score if probabilities exist, top-k assay hit rate with tie handling, rank correlation, coverage of recommendations, and recommendation stability. UMAP/silhouette scores are exploratory only.

Runtime: cold and cached latency, RAM, accelerator memory, disk, model size, dependency burden. Full-length information comparator and assay-specific matching rates must appear alongside classification metrics.

## 12. Simulation and independent validation

### 12.1 Sequence-first simulation

Use hidden query genomes as templates. Compare error-free inserts first, then substitutions/indels and empirical error profiles. ART and InSilicoSeq provide established simulation precedents, but default profiles may not match the selected current chemistry; validate or state the mismatch. [web:374][web:376]

Proposed stress grid: low/moderate/high quality profiles, single/paired read layouts, target abundance levels (including zero), balanced/staggered communities, primer compatibility scenarios, and multiple independent seeds. Select an achievable subset before running factorial expansion.

Simulation seeds and millions of reads do not increase the biological sample size. PCR mismatch matching is not a calibrated amplification-yield simulator. Model primer dropout, relative-yield sensitivity, chimera formation, and sequencing errors as separate mechanisms; document limitations.

### 12.2 Read-processing validation

Use a fixed trimming/merging/denoising workflow and reference mapping across methods. DADA2 resolves sample sequence variants but ASVs are not equivalent to species; preserve genome/copy mappings. Never train a denoiser on pooled development/test samples if it violates the intended deployment scenario; specify batch-level error estimation policy. [web:380][web:297]

### 12.3 Empirical tiers

Tier 1: independent isolate sequences or reads with genome-supported identity.
Tier 2: existing mock community with known membership and primer metadata; validate only assays actually measured.
Tier 3: prospective matched assays on the same DNA materials, including the recommended candidate, a common comparator, and full-length if feasible.

For new laboratory work: choose low-risk available strains/materials with appropriate institutional oversight; establish identity using independent genome or validated orthogonal evidence; randomize assay/sample processing order; blind analysts to recommendation identity; include extraction/PCR negatives and positive controls; distinguish independently grown/extracted biological materials from technical PCR replicates. Use a laboratory-approved PCR protocol and pilot primer conditions rather than inventing universal cycling parameters.

A small 3–5-panel prospective pilot can test feasibility, not broad generalization. Determine final replicate counts from pilot variance and the panel-level precision target. Do not validate a recommended V1–V3 assay using only V4 data or crop V4 reads to simulate an unmeasured region.

## 13. Statistical analysis plan

Primary estimand: mean paired difference in selected-assay U across independent taxonomic blocks on the locked evaluation set. Report absolute units, 95% cluster-bootstrap intervals, and the full panel distribution. A practical starting bootstrap count is 2,000, with more if interval stability is poor.

Resample taxonomic blocks first; optionally genomes within blocks. Never bootstrap sequence pairs as if independent. Average technical seeds within a block before the primary comparison. If few blocks exist, report wide intervals and descriptive results rather than manufacturing power with more reads.

Primary comparison must be named before testing, e.g. best development-selected frozen selector versus best development-selected classical selector. If multiple confirmatory contrasts are retained, use a prespecified multiplicity correction such as Holm; label all other analyses exploratory. Nested cross-validation folds are not independent replicates for a naive t-test.

Optional exploratory mixed model: per-panel or genome-level success as a function of method, assay, sequence similarity, reference density, and interactions, with appropriate block/genome grouping. Confirm distributional assumptions; a mixed model is not a cure for too few independent blocks. Do not regress on test-derived features to redefine the primary method.

Report both macro block means and performance by difficult subgroup. Predeclare treatment of unavailable assays: candidate failure in the intended universe, not automatic deletion. Provide a common-feasible-candidate sensitivity analysis without presenting it as the only result.

## 14. Precision and sample-size planning

No fixed sample size is defensible before surveying independent blocks. Conduct a development-only pilot, estimate block-level paired variability, and simulate the interval width/power for scientifically meaningful improvements (e.g. 0.02, 0.05, 0.10 in U as hypothetical effect sizes).

Prespecify an acceptable uncertainty target with the advisor, rather than taking all available records as sufficient. Examine realistic block counts such as 10, 20, 40, 60. Incorporate block dependence and heterogeneous effects. Keep the locked test untouched during planning. Report final eligible, excluded, and effective independent block counts.

## 15. Software implementation and release

Proposed commands, not existing installed functionality:

```text
regioneer references build --manifest accessions.tsv
regioneer panels validate --panels panels.tsv
regioneer assays extract --primers primers.tsv
regioneer benchmark run --config study.yaml
regioneer recommend --target TAXON_ID --neighbors neighbor_ids.txt --constraints constraints.yaml
```

Output: ranked table, provenance JSON, candidate failures, uncertainty, unresolved neighbors, primer/read configuration, reference counts, ambiguity audit, and HTML/Markdown report. Distinguish geometric score, measured validation utility, and predicted utility; do not call all of them confidence.

Fit deployment calibration using out-of-fold development predictions. Recommend only within supported organism/data regimes; indicate insufficient reference support. Use flags such as recommended, conditional, no candidate meets criterion, and insufficient evidence. A feasible default policy could require target compatibility >=0.95 and a development-calibrated lower confidence estimate of U >=0.90, but these are proposed research thresholds, not diagnostic standards. Exact acceptance rules must be frozen before release testing.

Reference bootstraps provide rank stability conditional on the database; they do not guarantee out-of-database correctness. Label the source of uncertainty explicitly.

Unit tests: strand orientation; IUPAC matching; primer clipping; coordinate mapping; no amplification; multiple copies; contradictory species labels; deterministic same-input embedding; long-token handling; no test fitting; unknown taxa; missing metadata; output provenance.

Release a CLI and small redistributable fixture first. Add an R Shiny interface only after reproducibility and scientific validity are established. Separate code license from reference sequence/model licenses.

## 16. Resources, timeline, and gates

Illustrative schedule, not a promise: weeks 1–3 curation/pilot; 4–6 expanded classical analysis; 7–9 frozen embeddings; 10–12 simulation/locked evaluation; 13–16 external validation, software hardening, writing. Laboratory sequencing procurement may extend this.

Start classical work on CPU; profile candidate models on accessible GPU hardware before committing. Do not budget a specific GPU-memory size without measuring the checkpoint, token length, precision, and batch size. Cache unique sequences to avoid paying repeatedly for duplicated references.

Gate 1: no release without accession/label audit and valid split manifest.
Gate 2: no embedding claim without competitive alignment/k-mer controls.
Gate 3: no assay recommendation without explicit primer/read constraints.
Gate 4: no universal claim from 1–2 panels.
Gate 5: no claim of experimental validation from reference-only simulation.

## 17. Planned figures, tables, and reporting

Figures: pipeline/split diagram; reference coverage; ambiguity by assay; paired selected-assay performance; regret distribution; class/neighbor confusion; reference-density learning curves; error/primer robustness; uncertainty calibration; accuracy-resource frontier.

Tables: panel taxonomy/provenance; exact primers; candidate/read configurations; model revisions; primary test contrasts; unavailable assays; computational costs; limitations and source licenses.

Required artifacts: locked protocol, accession snapshot, split files, exclusion log, dependency locks, model hashes, immutable predictions, analysis scripts, tests, and a runnable small example. Archive code with a persistent release identifier if licensing permits.

## 18. Interpretation rules

If classical methods win: publish the benchmark and deploy the simpler method.
If embeddings help only easy panels: do not advertise improved close-species resolution.
If full-length dominates: emphasize when short assays are inadequate and the evidence for escalation.
If all short assays fail a panel: report no suitable evaluated candidate.
If no independent experimental data exist: call outputs in-silico recommendations awaiting assay validation.

Source citations in this document refer to the shared annotated literature review. All formulas, threshold examples, workflow interfaces, and project-specific objectives are proposed protocol elements unless explicitly attributed.
