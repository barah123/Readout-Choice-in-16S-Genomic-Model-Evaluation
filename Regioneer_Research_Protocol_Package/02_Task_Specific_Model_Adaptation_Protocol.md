# Protocol 2: Task-specific adaptation of genomic representations for 16S amplicon selection

Version 1.0 | Prospective protocol | 7 November 2026

## 1. Title, scope, and relationship to protocol 1

Working title: Task-specific adaptation of genomic encoders for target-versus-neighbor discrimination and selection of informative 16S amplicons.

Design: controlled representation-learning and transfer study. It reuses a versioned reference/extraction/evaluation framework from protocol 1 but requires an independently protected evaluation set if protocol-1 results have already influenced development.

Primary intervention: adapt a pretrained encoder using a lightweight projection, parameter-efficient updates, or full fine-tuning. Training a foundation model from scratch is explicitly out of scope for the core study. A compact randomly initialized encoder can be a control, not a claim to have created a broadly pretrained genomic foundation model.

Terminology: use genomic language model or sequence encoder, not necessarily genomic LLM. An encoder that outputs embeddings is not a conversational model.

## 2. Literature rationale and prior-art limits

DNABERT-2 provides a multi-species pretrained sequence backbone. DNABERT-S adapts it to species-aware embeddings using MI-Mix and curriculum contrastive learning. Species-aware learning is therefore existing work, not the novelty claim. [web:353][web:331]

Micro16S uses triplet and pair objectives to encode phylogenetic relationships and promote region invariance. Its preprint reports limitations at lower taxonomic ranks and weaker downstream sample classification than classical baselines. This is relevant precedent but not evidence that task-specific amplicon ranking cannot work. [web:318]

Parameter-efficient approaches reduce the number of updated parameters, but benefits depend on architecture and task. A genomic study of LoRA-based adaptation and a separate efficiency investigation support evaluating resource/performance trade-offs rather than assuming full fine-tuning is essential. [web:313][web:319]

Scientific gap: whether a sequence representation adapted using short-16S target-versus-neighbor examples improves held-out assay recommendations beyond strong classical/frozen systems, especially for globally unseen organism groups with local support references.

## 3. Objectives and hypotheses

Primary objective: test whether task-specific adaptation reduces recommendation regret and/or improves the independent query performance of the selected feasible amplicon relative to the strongest prespecified non-adapted selector.

Primary hypothesis: a development-selected adaptation strategy improves selected-assay utility U on independent taxonomic blocks by a practically meaningful amount. Do not assume the result is positive.

Secondary hypotheses:
- Hard negatives from close species improve local discrimination more than taxonomically distant negatives.
- Same-region objectives better preserve useful assay differences than unconditional cross-region invariance.
- Improvement depends on reference density and label quality.
- Lightweight adaptation may recover much of the attainable improvement at lower operational cost.

A no-improvement result, especially under strict splits and strong baselines, is publishable evidence if the evaluation is sufficiently precise and the tested adaptations are credible.

## 4. Biological feasibility and unavoidable ambiguity

The model observes a sequence, not the true genome label. A target and a neighbor can share identical amplicons; no deterministic sequence-only embedding can make those identical inputs species-specific.

Compute the empirical collision/ceiling analysis from protocol 1 before training. Do not train a contrastive loss to push apart identical observed strings with different species labels: deterministic encoders cannot satisfy that demand. Taxonomy errors, shared 16S alleles, and intragenomic variants must remain visible.

Keep ambiguous examples in end-to-end evaluation. For training, use an explicit policy: mask impossible cross-label constraints; use set-valued species compatibility where warranted; and/or train at the coarsest defensible equivalence group. Report which policy was used. Do not convert unresolved examples into confidently separated training labels by deleting difficult neighbors.

Distinguish representation generalization from classification label support. A globally unseen species can still be classified using a labeled support reference embedded at deployment. Without such references, exact assignment to an arbitrary unseen species name is not the claimed task.

## 5. Study population and data construction

### 5.1 Shared data contract

Use genome-backed 16S sequences with stable accessions, genome/isolate grouping, copy mappings, reconciled species labels, assay boundaries, and model-ready observed inserts. Reference-only sequences without genome linkage form a separately labeled auxiliary corpus, not silently equivalent training rows.

Store taxonomy hierarchy, relatedness, sequence hashes, copy IDs, primer compatibility, and support/query roles. No manual label editing driven by model errors.

### 5.2 Data scale and diversity

Start with available audited development blocks, not an arbitrary million-sequence target. For feasibility, aim for enough species per block and independent genomes per species to form valid positive pairs and support/query evaluation. A preliminary target might be 30–60 independent development/evaluation blocks across multiple genera; final scope follows the eligible data and precision analysis.

Track counts of independent genomes and unique observed sequence types, not just region-expanded rows. Nine amplicons from one locus and many simulated reads are not nine biological replicates. Avoid training only on two favorite panels and calling the representation general-purpose.

### 5.3 Train/development/test hierarchy

Level A: global adaptation-training blocks.
Level B: global development blocks used for model selection and stopping.
Level C: locked independent test blocks, preferably with an additional family-held-out challenge.

Within each block assign genomes to local support and hidden query roles. The adapted encoder is trained only using global training blocks; model selection uses development blocks; test-panel local support may fit a fixed deployment classifier but may not update the global encoder in the primary transfer experiment.

Separate experiments:
1. Familiar-species/new-genome adaptation.
2. Unseen-genus/block transfer with local support.
3. Optional new-panel fine-tuning with a fixed support budget; this is few-shot adaptation and must not be conflated with frozen deployment.

If an independent test was used in protocol 1 to identify weak cases or choose architecture, it is now development data. Obtain new protected blocks or a prospective reference release and audit duplicate genomes across releases.

### 5.4 Training and pretraining contamination audit

Remove repeated isolate/genome records across global roles. Audit exact full-length, amplicon, and reverse-complement identity between fine-tuning and test data. For the primary clean-transfer analysis, exclude overlapping observed hashes from supervised adaptation while retaining their occurrence in the ecological test population; also publish realistic-overlap sensitivity results. Same-species novelty and unseen-genus generalization answer different questions.

Strict similarity-component splits can connect many species because 16S is conserved; do not apply an arbitrary low identity threshold until the entire corpus becomes one component. Report component sizes, relatedness distribution, and how the split changes the intended task.

Pretraining exposure may be unknown even after downstream decontamination. Record declared corpora; use release-date restrictions only as partial evidence, not proof. Distinguish pretraining exposure, supervised fine-tuning leakage, and deployment support overlap in reporting.

## 6. Experimental arms and escalation

A0: best classical selector from protocol 1, frozen recipe.
A1: frozen general genomic model + common downstream scorer.
A2: frozen species-aware DNABERT-S + identical scorer. [web:331]
A3: frozen backbone + trainable small projection head.
A4: parameter-efficient fine-tuning + the same projection/scorer.
A5: full fine-tuning only if A3/A4 and resources justify it.
A6: optional compact randomly initialized sequence encoder, to test whether external pretraining contributes.

Optional: reproducible Micro16S representation arm if available and compatible; label the supporting work as a preprint. [web:318]

Keep local classifier, calibration policy, candidate assays, and query observations fixed across arms where possible. If a method needs different scoring, include a matched-scoring ablation so changes are not misattributed to representation learning.

## 7. Proposed model architecture

Input x is a primer-trimmed observed insert from one candidate assay. The backbone f_theta yields hidden states, pooling P yields h, and an optional projection g_phi yields a normalized vector z:

h = P(f_theta(x)); z = g_phi(h) / ||g_phi(h)||.

A panel-conditioned score is computed later from target and neighbor support embeddings. Do not concatenate a species name or taxonomy ID into a sequence-only encoder and then claim it recovered species identity from DNA.

Starting projection: one or two linear layers with a modest hidden dimension and 64/128/256-dimensional output tested on development data only. These dimensions are prospective hyperparameters, not biological constants.

Record layer pooling, token limits, padding mask, orientation, numerical precision, and deterministic inference behavior. Use official checkpoint/tokenizer interfaces and pin revisions; avoid assuming every model exposes the same modules or supports the same LoRA targets.

## 8. Training examples and sampling

### 8.1 Primary same-assay pairs

Anchor: one observed insert from a training genome.
Positive: same species, same assay, preferably a different isolate/genome and, when available, a different copy.
Negative: a different species from the same declared near-neighbor block and same assay.

Same-assay matching prevents region identity from becoming the easy species shortcut. Restrict or weight positive constraints when a species contains strongly divergent copies; verify label validity before forcing compactness.

An example batch design is P species by K genomes/sequence instances per species, sampling within a taxonomic block and assay. P/K are tuned to memory and data availability; record the number of independent genomes actually represented. If a species lacks valid positives, do not fabricate them by counting identical technical read copies as independent examples.

### 8.2 Negative mining

Start with random within-block negatives, then development-selected semi-hard/hard negatives. Include distant negatives only as a smaller auxiliary component so the loss is not solved mainly by phylum-level separation.

Mine only from the adaptation-training partition using a frozen baseline or current training model. Never inspect test labels or query performance to select hard examples. Screen identical cross-label sequences and uncertain taxonomic cases before using them as negatives.

### 8.3 Permitted augmentation

Reverse complement if consistent with the representation/strand policy; realistic read truncation; empirical sequencing noise; limited ambiguity masking where supported. Do not assume arbitrary nucleotide substitutions preserve fine-grained identity: a single diagnostic base can determine discrimination. Keep clean and noise-augmented losses/results separate.

Do not augment with unsequenced template context. Paired reads require an explicit architecture/masking strategy, not invented filled-in sequence.

## 9. Training objectives

### 9.1 Primary simple objective: supervised contrastive learning

For anchor i, P(i) is the set of valid same-species positives from compatible observed-assay conditions; A(i) contains valid comparison examples excluding self and masked impossible pairs. With cosine similarity sim and temperature tau:

L_i = -(1/|P(i)|) * sum_{p in P(i)} log[ exp(sim(z_i,z_p)/tau) / sum_{a in A(i)} exp(sim(z_i,z_a)/tau) ].

Skip anchors with no valid positives; report their fraction. Use masks to exclude exact cross-label identity constraints and uncertain labels. This is a proposed baseline loss; DNABERT-S has its own published MI-Mix/C2LR procedure and should not be described as identical to this implementation. [web:331]

### 9.2 Alternative triplet objective

L_triplet = max(0, d(z_anchor,z_positive) - d(z_anchor,z_negative) + margin).

Use a limited development-selected margin grid, normalized embeddings, and an explicit valid-triplet policy. Micro16S provides relevant precedent for triplet and pair losses, but reproducing it requires its exact published formulation rather than substituting this simplified expression. [web:318]

### 9.3 Do not force unconditional cross-region invariance

A weakly informative V4 observation and a highly informative V1–V3 observation from the same genome need not support equally certain species predictions. Forcing identical embeddings can conceal that difference. Primary training is same-assay; cross-region alignment is a separately weighted ablation, not the default assumption.

If used, align only supported information or employ uncertainty-aware/masked compatibility objectives. Compare assay-ranking quality, not just an invariance metric. A model must not be rewarded for pretending that a shared short sequence contains species-specific information available only in another region.

### 9.4 Optional direct ranking objective

Only after the representation experiment is stable, learn an assay ranker from support-summary features: within-target variation, target-neighbor distances, shared sequence fraction, primer compatibility, length, and reference density.

Training ranking targets must come from held-out genomes internal to global training blocks, not from global development/test queries. For each training panel create an inner support/evaluation split, estimate assay utilities U_a and U_b, and form a preference only when the estimated difference exceeds a declared uncertainty/tie tolerance.

A pairwise logistic example is L_rank = log(1 + exp[-y_ab*(s_a-s_b)]), with y_ab in {-1,+1}. Ties can be omitted or modeled explicitly. Labels are noisy finite-sample performance estimates, not biological ground truth.

Use nested grouping to ensure neither the genome used to construct ranking labels nor a near-duplicate leaks into the features for that evaluation role. Treat this ranker as an additional arm, not an unreported change to the encoder-only hypothesis.

## 10. Optimization and compute plan

Begin with A3, then A4. Example development-only starting grids:
- backbone learning rate: 1e-5, 3e-5, 1e-4;
- projection learning rate: 1e-4, 3e-4, 1e-3;
- temperatures: 0.05, 0.1, 0.2;
- LoRA rank: 4, 8, 16;
- projection dimension: 64, 128, 256;
- effective batch size determined by valid positive/negative content and profiling, not GPU memory alone.

Do not run the full Cartesian product. Use a fixed staged search budget, e.g. a small pilot followed by at most 10–20 development runs per principal adaptation family if resources permit, documenting the actual budget. Give classical comparators a credible hyperparameter search too.

Use early stopping on development-panel recommendation utility or a predeclared proxy; never on locked test accuracy. Save last/best checkpoints, optimizer state, random seeds, and loss components. Inspect gradient norms, embedding variance, positive/negative distances, and valid constraint counts. Lower loss without improved development utility is not success.

Mixed precision, accumulation, and checkpointing are implementation options after numerical validation. Profile short/long amplicons and representative batches. Record GPU type, measured peak memory, elapsed time, tokens, optimizer steps, trainable/total parameters, and energy estimate only if actually measured.

LoRA module names depend on architecture. Inspect the model and verify that only intended weights require gradients; compare counts before training. A frozen backbone must remain unchanged under a checkpoint hash/tensor test. Avoid unsupported assertions that adaptation fits on a particular consumer GPU without a measured run.

Repeat the final selected setting across several seeds, e.g. 3–5 if feasible. Seeds describe optimization variability, not independent biology.

## 11. Evaluation and estimands

### 11.1 Primary

Use the protocol-1 support/query recommendation evaluation. The adapted strategy sees support references and feasible assays, selects an assay without query labels, and is scored on independent query genomes using the same observed-sequence/assay policy as its comparators.

Primary contrast: adapted selected-assay utility minus development-selected strongest non-adapted comparator, averaged over independent test blocks. Predeclare one adaptation family or a fixed development selection rule before opening the test.

Regret versus the evaluation-only best feasible assay is a complementary endpoint. Report whether the model changes recommendations, whether those changes help, and whether classification improves when the assay is fixed. These distinguish better representation from better selection.

### 11.2 Secondary

Species-balanced accuracy, neighbor-wise errors, macro-F1, AUROC/AUPRC with prevalence, Brier score/calibration, rank stability, ambiguity subgroup performance, resource cost, and abstention risk/coverage.

For transfer: report results by familiar species, unseen genera/blocks, and reference density. A classifier cannot predict an unseen species label absent from the allowed support/label mechanism; define this honestly.

UMAP, nearest-neighbor plots, phylogenetic association, and clustering metrics are diagnostic visualizations, not primary evidence of experimental utility.

### 11.3 Calibration and abstention

Fit temperature scaling or another prespecified calibration only on development predictions with proper support/query separation. Use an explicit unknown/out-of-panel test set for open-set rejection; do not equate a low maximum softmax probability with validated unknown-species detection.

Thresholds may target research-level sensitivity/specificity, but they must be selected before test evaluation. Plot accuracy among accepted queries against the accepted fraction. Report no-recommendation rates and confidence interval coverage. If a bootstrap interval is conditional on the reference library, label it as such.

### 11.4 Statistics

Paired block bootstrap, averaging nested panels and seeds appropriately; 95% intervals for the primary contrast and regret difference. Do not treat training seeds, sequence pairs, or overlapping panels as independent replicates. Use a prespecified multiplicity strategy for confirmatory ablations; otherwise label them exploratory.

Define a practical improvement threshold before test opening, e.g. a hypothetical 0.03 absolute U improvement with acceptable resource overhead. This is a decision criterion to discuss, not a universal scientific standard. Interpret interval width and cross-block consistency, not a p-value alone.

Precision/power planning uses development-only block-level variability and realistic effect scenarios. If independent test blocks are scarce, emphasize uncertainty and limit the claim rather than generating more synthetic reads to inflate n.

## 12. Required ablations and negative controls

Essential:
1. Classical versus frozen general versus frozen species-aware versus projection-only versus PEFT.
2. Same-assay loss versus cross-region alignment.
3. Random versus near-neighbor hard negatives.
4. Genome-grouped versus taxonomic-block transfer.
5. Clean versus empirical-error observations.
6. Reference density and species balancing.
7. Identical-input/shared-sequence ambiguity subset.
8. Equal downstream classifier/scoring where feasible.

Useful controls:
- Label permutation within the allowed training structure; no true-label test information used for selection.
- Diagnostic-base perturbation/counterfactuals chosen using training references, recognizing these are mechanistic probes rather than natural genomes.
- Length-only and GC-only classifiers to detect shortcuts.
- Same-input deterministic embedding check.
- Copy-aware versus one-copy-per-genome sensitivity.
- Model pretraining exposure uncertainty and temporal/source holdout.

Never select counterfactual bases using test labels, then present performance as independent validation.

## 13. Interpretation of outcomes

Adaptation helps only familiar species: limit to within-panel adaptation; do not claim universal transfer.
Adaptation improves fixed-assay classification but not recommendations: report a representation/classification result, not an assay-design breakthrough.
Adaptation improves geometric separation only: insufficient evidence for deployment.
Adaptation appears to resolve identical cross-label inputs: investigate leakage, nondeterminism, label-derived inputs, or context shortcuts immediately.
PEFT matches full fine-tuning: prefer the simpler operational option if other properties are similar.
Classical remains best: stop model escalation unless a new mechanistic hypothesis and independent evaluation justify it.

A universal model is not required for a useful paper. A clear negative result about close-species 16S and learned representations can prevent inappropriate adoption.

## 14. Ethics, reproducibility, and licensing

Use public or authorized data with accession provenance. Remove unintended host sequence where relevant. Avoid uploading controlled data to hosted embedding services. Do not distribute source sequences or checkpoints contrary to terms. Separate code, metadata, model, and database licenses.

No claims of diagnostic suitability. Select validation organisms/materials under institutional biosafety review; use purified material or appropriate surrogates where sufficient. Do not make unsupported pathogenicity predictions.

Release config, split hashes, model revisions, full test predictions, failed runs relevant to selection, compute logs, and evaluation scripts. Record human curation and model-selection decisions. Archive a model card documenting intended use, unsupported taxa, pretraining uncertainty, ambiguity behavior, and assay constraints.

## 15. Timeline and decision gates

Illustrative after the shared benchmark: weeks 1–2 split/data audit; 3–4 frozen/projection baselines; 5–7 PEFT experiments; 8–9 ablations and calibration; 10–12 locked test and reproducibility; additional time for laboratory validation.

Gate A: enough independent audited data and a protected test set.
Gate B: a reproducible failure of the frozen approach that adaptation could plausibly address.
Gate C: development gain versus strong classical/species-aware baseline, not just training loss reduction.
Gate D: improvement survives biological holdout, ambiguity checks, and cost assessment.
Gate E: release only with calibrated limits and an abstention policy.

## 16. Manuscript and software outputs

Main figures: data/split hierarchy; objective and valid-pair masks; training diagnostics; block-level selected-assay improvement; regret; taxonomy-transfer breakdown; ambiguity ceiling; calibration/risk-coverage; error and reference-density curves; compute frontier.

Main tables: all checkpoints and trainable parameter counts; datasets and provenance; ablation budget; primary independent comparisons; failure subgroups; model-card restrictions.

Provisional contribution statement: We evaluate whether task-specific adaptation of genomic representations improves target-conditioned 16S amplicon selection under realistic biological and assay constraints. Do not claim first species-aware DNA embeddings or first learned 16S representations. [web:331][web:318]
