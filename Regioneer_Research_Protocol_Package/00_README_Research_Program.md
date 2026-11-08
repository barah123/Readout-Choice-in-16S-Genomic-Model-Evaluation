# Regioneer research program: implementation guide

Version: 1.0 | Prepared: 7 November 2026 | Status: prospective protocol, not executed research

## Deliverables and reading order

1. `01_Regioneer_Benchmark_and_Tool_Protocol.md`: Research idea 1—benchmark classical and frozen learned representations for target-conditioned 16S amplicon selection; implement a validated recommender.
2. `02_Task_Specific_Model_Adaptation_Protocol.md`: Research idea 2—test whether task-specific representation learning or fine-tuning improves recommendations beyond the strongest frozen/classical system.
3. `03_Annotated_Literature_Review.md`: shared critical literature review, prior-art map, source identifiers, and reading priorities.
4. `04_Implementation_Schemas_and_Checklists.md`: proposed repository, data contracts, configuration templates, pseudocode, audit checklists, and reporting forms.

The documents are a guide for implementation, not a claim that regioneer or its commands currently exist. Commands prefixed `regioneer` and workflow filenames are proposed interfaces. Numeric design choices are prospective starting values, not established universal standards or a completed power calculation. Change them before locking the protocol; record subsequent deviations.

## Central scientific distinction

Research idea 1 asks: Can a target-specific scoring approach recommend experimentally feasible 16S amplicons that discriminate a target from specified near-neighbor species, and do frozen embeddings help?

Research idea 2 asks: Does adapting a pretrained sequence representation improve that recommendation task on independent organism groups, after strong classical and frozen-model baselines are established?

Neither project is primarily about making attractive UMAP clusters. The final endpoint is the performance of the recommended assay on independent biological sequences and, where available, independent experimental material.

## Evidence update that affects novelty

Qscore is closer prior art than a generic region-comparison paper: its publication describes an online system that takes taxonomy or OTU inputs, accepts expected microbial structures or taxa of interest, and recommends amplification configurations. Therefore, do not claim the first user-facing or taxa-informed 16S recommender. The proposed distinction is an explicitly target-versus-neighbor, genome-held-out, ambiguity-aware benchmark with realistic assay constraints and a direct test of learned representations. [web:300]

DNABERT-S already learns species-aware DNA embeddings. Micro16S is a 2026 preprint about phylogenetic 16S representations. Neither species-aware embeddings nor learned 16S representations are novel by themselves. Transfer to short amplicons and the quality of target-conditioned assay recommendations remain hypotheses to evaluate. [web:331][web:318]

## Claims you may and may not make

Appropriate before results:
- We will benchmark target-conditioned amplicon-selection strategies.
- We will evaluate whether learned representations add value to classical sequence methods.
- We will quantify conditions under which no evaluated amplicon meets the selected performance criteria.

Not appropriate before results:
- AI identifies the universally optimal 16S region.
- All target species can be distinguished by selecting the right region.
- An embedding can resolve species with identical observed sequences.
- This is the first software to recommend sequencing regions from taxa of interest.
- A handful of organism panels establishes general applicability.

## Recommended staged execution

Stage A: repository, accession provenance, taxonomy reconciliation, and a small pipeline-verification pilot.
Stage B: expanded classical benchmark and explicit sequence-ambiguity audit.
Stage C: frozen representations, common downstream evaluation, and assay simulation.
Stage D: locked independent evaluation; release a classical-only recommender if it is best.
Stage E: adaptation only if data, resources, and a testable failure mechanism justify it.
Stage F: independent laboratory validation where possible; publish the level of evidence actually obtained.

Research idea 1 can stand alone. Research idea 2 is contingent, and a no-improvement result is valid. Do not repeatedly inspect one test set to decide which model to release. If stage-1 test results guide stage-2 development, acquire a new locked stage-2 test set.

## Division of work

Python: accessions and sequence processing, in-silico PCR, feature extraction, embeddings, training, predictions, and software CLI.
R: panel-level summaries, cluster-aware uncertainty, mixed-effects exploratory analyses, publication figures, and reproducible reports.
Workflow engine: choose either Snakemake or Nextflow; avoid maintaining both for the first release.
Application layer: a CLI first; add R Shiny or another GUI only after the scientific core is stable.

## Scope and uncertainties

Default domain is bacteria. Archaea require separately validated extraction, taxonomy, and primers and should not be advertised as supported merely because a source database includes them. Default task is identification of a sequence from one organism in a prespecified panel. Environmental community detection and clinical diagnostics require additional endpoints and validation.

Public sequence data generally avoid new participant recruitment, but institutional determinations, source terms, human-read handling, controlled-access restrictions, and biosafety requirements remain applicable. This guide is not institutional ethics or laboratory approval.

The literature review is targeted, not a systematic or exhaustive review. Publication details supplied in the catalogue are those verified from the linked records; where only a title and source are established, do not fill in missing author or DOI fields by guessing. Record database/checkpoint versions at implementation rather than treating any release as permanently current.
