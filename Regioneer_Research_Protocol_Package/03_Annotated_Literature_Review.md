# Shared annotated literature review and prior-art map

Prepared 7 November 2026. This is a targeted critical review for implementation, not an exhaustive systematic review. Citations link to the source records used to prepare the protocols. Results reported by a publication are not automatically validated for the proposed regioneer population or assay task.

## 1. Main conclusions from the literature

1. The biological question is established: discriminatory information varies across 16S segments and taxa. Full-length measurements can improve resolution, but do not justify universal species/strain-identification claims. [web:297]
2. Assay selection is more than choosing a named variable region. Primer coverage, sequence observability, read configuration, reference taxonomy, and analysis can change the result. [web:360][web:300]
3. The closest tool prior art includes Qscore, not only region-extraction utilities. Its publication explicitly describes a taxa-/OTU-informed online recommendation system. A proposal that claims novelty simply because users submit taxa and receive a sequencing configuration would overlap this work. [web:300]
4. Species-aware DNA embeddings already exist. DNABERT-S must be considered when claiming representation-learning novelty. [web:331]
5. Learned 16S embeddings already exist as a research direction. Micro16S is especially pertinent but was a preprint in the version reviewed. It optimizes phylogenetic structure/region invariance rather than exactly the same assay-selection objective. [web:318]
6. Baseline strength and scope of validation matter. Large-sequence embedding benchmarks are not direct evidence of reliable short-16S discrimination. [web:347][web:348]

## 2. Biological resolution and the value of full-length controls

### L01. Johnson et al., 2019 — essential biological anchor

Title: Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome analysis.
Venue: Nature Communications 10, 5029.
DOI: 10.1038/s41467-019-13036-1.
Record: https://www.nature.com/articles/s41467-019-13036-1

The study combines reference-derived experiments, a bacterial mock community, and biological samples. Its central points are that subregions differ in resolution, full-length data retain more information, and intragenomic 16S copies can contain meaningful variants. It also shows why single-nucleotide resolution should not be equated automatically with one variant per species. [web:297]

Relevance: requires full-length controls and genome/copy-aware bookkeeping. It directly motivates showing which target-neighbor distinctions are lost when an amplicon is shortened.

Limitations for reuse: some evaluations classify reference-derived sequences against the originating database; historical taxonomy, clustering rules, and platform errors differ from the proposed independent-genome design. Do not import reported percentages or 99% clustering as universal expected performance.

Implementation action: reproduce one small qualitative subregion comparison with current curated references, then make held-out isolate performance—not same-database re-identification—the primary validation.

### L02. Comparison of the full-length sequence and sub-regions of 16S…, 2024

Record: https://pmc.ncbi.nlm.nih.gov/articles/PMC11264597/

The comparison reports superior taxonomic resolution from full-length sequences while noting remaining limitations even at full length. It supports treating short-region optimization as a constraint-dependent problem, not a substitute for considering full-length sequencing. [web:306]

Implementation action: include a full-length information comparator and report what is gained/lost by shorter assays. Retrieve the complete bibliographic export before citing in a manuscript; the title shown here is abbreviated from the source record.

### L03. Stoddard et al., 2015 — rrnDB

Title: rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development.
Venue: Nucleic Acids Research 43(D1): D593–D598.
DOI: 10.1093/nar/gku1201.
Record: https://pubmed.ncbi.nlm.nih.gov/25414355/

The database documents variation in ribosomal operon copy number and supplies curated resources for interpretation. Copy number affects the mapping of molecules to organisms, and divergent copies within a genome complicate apparent sequence diversity. [web:359]

Implementation action: retain genome-level membership and copy counts; weight genomes separately from copies. rrnDB can contextualize expected copy number but cannot repair a collapsed draft assembly or establish that all copies were recovered.

## 3. Region extraction, primers, and assay recommendation

### L04. V-Xtractor, 2010 — extraction rather than recommendation

Title: V-Xtractor: An open-source, high-throughput software tool to identify and extract hypervariable regions of small subunit (16S/18S) ribosomal RNA gene sequences.
Record: https://www.sciencedirect.com/science/article/abs/pii/S0167701210002745
Code record: https://github.com/carden24/V-Xtractor

The software identifies hypervariable regions using hidden Markov models. This provides clear precedent for automating extraction, but extraction does not by itself establish which region best separates a particular species panel. [web:282][web:293]

Implementation action: evaluate extraction boundaries and failure cases; preserve coordinate and strand provenance. If using a maintained fork, record that fork's commit and differences rather than citing only the original tool.

### L05. Klindworth et al., 2013 — essential primer baseline

Title: Evaluation of general 16S ribosomal RNA gene PCR primers for classical and next-generation sequencing-based diversity studies.
Venue: Nucleic Acids Research 41(1): e1.
DOI: 10.1093/nar/gks808.
Record: https://pubmed.ncbi.nlm.nih.gov/22933715/

The paper evaluates 175 primers and 512 primer pairs against a SILVA reference set and experimentally examines a selected bacterial pair. Primer choice changes coverage and taxonomic spectrum. [web:360]

Relevance: distinguishes sequence information from what is actually amplified. A discriminatory region that consistently misses target strains is a poor recommendation.

Implementation action: source exact primer sequences, including IUPAC variants, and report mismatch rules by position. Do not treat the historical best pair as universally optimal under modern reference diversity or different study objectives.

### L06. mopo16S, 2018 — direct optimization prior art

Title: Optimizing PCR primers targeting the bacterial 16S ribosomal RNA gene.
Venue: BMC Bioinformatics.
DOI: 10.1186/s12859-018-2360-6.
Record: https://link.springer.com/article/10.1186/s12859-018-2360-6
Code: https://github.com/sambofra/mopo16s

The method jointly optimizes primer efficiency/specificity, coverage, and coverage bias. It demonstrates that multi-objective primer selection is established work. [web:332][web:337]

Relevance: a baseline and design reference for constraint handling. It is not automatically equivalent to target-versus-neighbor assignment utility, so compare objectives explicitly.

Implementation action: inspect input assumptions, candidate primer generation, output objective scales, and software reproducibility. If a literal rerun is infeasible, report a faithful simplified objective benchmark as an adaptation, not an exact replication.

### L07. Zhang et al., 2023 — Qscore; closest recommendation-system precedent

Title: Comprehensive Assessment of 16S rRNA Gene Amplicon Sequencing for Microbiome Profiling across Multiple Habitats.
Venue: Microbiology Spectrum 11(3): e00563-23.
DOI: 10.1128/spectrum.00563-23.
Record: https://pmc.ncbi.nlm.nih.gov/articles/PMC10269731/
Code described by the paper: https://github.com/qdu-bioinfo/qscore

Qscore combines amplification sensitivity, taxonomic annotation, sequencing type/length, and a cost term. The study describes a global assessment of 35,889 species and habitat-conditioned analyses using 157,390 microbiomes. Its online service takes taxonomy or OTU inputs and returns recommended amplification configurations, with customizable objective weights. [web:300]

This materially narrows the novelty claim. User-facing, taxon-informed selection is already published. Regioneer needs to distinguish its target-negative panel definition, independent evaluation, ambiguity handling, assay-specific uncertainty, learned/classical comparison, or validated performance—not merely its interface.

Important limitations to consider: simulated amplification/error models simplify laboratory effects; database completeness and taxonomy affect assignment; cost proxies are platform/time specific. Reported community-level performance close to a shotgun comparator under particular simulations does not prove clinical species discrimination across arbitrary target complexes.

Implementation action: reproduce representative Qscore behavior where code permits; compare under matched references and assay sets; publish where objectives differ. A smaller targeted discrimination objective can be useful even if it does not dominate a community-profile objective.

### L08. HVRLocator, 2026 — recent extraction/annotation tool

Title: HVRLocator: a computationally efficient tool for identifying…
Record: https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giag040/8644286
DOI shown in record: 10.1093/gigascience/giag040.

The record describes locating amplicon endpoints and identifying corresponding hypervariable regions. This is relevant to boundary QA and assay annotation, not automatically evidence of a target-specific ranking system. [web:285]

Implementation action: inspect the complete paper, supported inputs, and code before adopting. Retain the abbreviated title until a full bibliographic export is verified; do not invent missing metadata.

## 4. Databases, taxonomy, and reproducibility

### L09. Quast et al., 2013 — SILVA foundation

Title: The SILVA ribosomal RNA gene database project: improved data processing and web-based tools.
Venue: Nucleic Acids Research 41: D590–D596.
DOI: 10.1093/nar/gks1219.
Bibliographic record verified through: https://pmc.ncbi.nlm.nih.gov/articles/PMC5210518/

SILVA is a key rRNA reference resource. The historical paper should be cited with the actual release used; citing a database article does not specify a reproducible snapshot. [web:367]

Implementation action: record release identifier, download date, filters, checksums, and taxonomic completeness. Do not treat all species labels or all records as equally verified.

### L10. SILVA in 2026 — resource update

Title: SILVA in 2026: a global core biodata resource for rRNA…
Record: https://pmc.ncbi.nlm.nih.gov/articles/PMC12807666/

Use this update to understand the resource state and choose a documented release, not to assume that the newest database automatically improves the target task. [web:368]

Implementation action: inspect release changes and licensing at execution. Preserve a primary frozen snapshot and analyze a newer snapshot as a temporal robustness test after duplicate auditing.

### L11. Robeson et al., 2021 — RESCRIPt

Title: RESCRIPt: Reproducible sequence taxonomy reference database management.
Venue: PLOS Computational Biology.
DOI: 10.1371/journal.pcbi.1009581.
Record: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009581

RESCRIPt supplies a Python/QIIME 2 framework for retrieving, filtering, comparing, and generating reference databases with provenance. Its focus on database composition and standardized processing is directly relevant to a reproducible benchmark. [web:373]

Implementation action: adopt the provenance principles even if the whole project is not implemented as a QIIME 2 plugin. Keep source sequences, taxonomic crosswalks, and filters inspectable.

### L12. GTDB release 10 — genome-derived taxonomy context

Title: GTDB release 10: a complete and systematic taxonomy for 715 230…
Record: https://espace.library.uq.edu.au/view/UQ:4cf482e

The record describes a large genome-based bacterial/archaeal taxonomy release. Such taxonomy can define organism blocks independently of the representation being tested. [web:356]

Implementation action: do not mix GTDB species assignments and NCBI names without a crosswalk. Genome-derived taxonomic boundaries need not be perfectly separable by one marker gene. Confirm full publication metadata and the release actually used before manuscript submission.

## 5. Classical methods and frozen genomic representations

### L13. 16S Classifier, 2015 — early ML baseline

Title: 16S Classifier: A Tool for Fast and Accurate Taxonomic…
Record: https://pmc.ncbi.nlm.nih.gov/articles/PMC4315456/

The study uses random forests and short-sequence features for 16S taxonomic classification and compares hypervariable-region behavior. It establishes that machine learning and region-dependent classification are not new contributions by themselves. [web:289]

Implementation action: ensure the benchmark includes strong composition/classification approaches, not only naive GC differences. Export the complete citation before manuscript use.

### L14. DNABERT-2 — general backbone

Title: DNABERT-2: Efficient Foundation Model and Benchmark For Multi-Species Genome.
Status: conference-paper version identified as ICLR 2024.
Record: https://arxiv.org/pdf/2306.15006v2.pdf

DNABERT-2 uses a multi-species sequence pretraining framework and an efficient tokenization/architecture design. It is an appropriate general baseline, not automatically specialized for short conserved bacterial marker discrimination. [web:353]

Implementation action: pin the checkpoint/tokenizer and inspect pooling/context requirements. Do not silently truncate sequence or compare differently pooled checkpoints without documenting that choice.

### L15. Zhou et al., 2025 — DNABERT-S; mandatory adaptation prior art

Title: DNABERT-S: pioneering species differentiation with species-aware DNA embeddings.
Venue: Bioinformatics 41(Supplement_1): i255–i264.
DOI: 10.1093/bioinformatics/btaf188.
Record: https://academic.oup.com/bioinformatics/article/41/Supplement_1/i255/8199363
Official code record: https://github.com/MAGICS-LAB/DNABERT_S

The published work adapts DNABERT-2 for species-aware representations, introducing Manifold Instance Mixup and curriculum contrastive learning. It evaluates species classification, clustering, and metagenomic tasks. Use the peer-reviewed paper rather than treating the older 2024 preprint as the final evidence. [web:331][web:338]

Relevance: your proposed objective overlaps established species-aware learning. The new question is transfer to short 16S inserts and reliable assay recommendations, not the general possibility of pulling same-species DNA together.

Implementation action: benchmark pretrained DNABERT-S before attempting custom adaptation. Inspect training sequence lengths, organisms, corpus overlap, and inference assumptions. Do not assume performance on genomic fragments transfers to nearly identical amplicons.

### L16. Benchmarking DNA Foundation Models for Genomic Sequence Classification, 2024

Records: https://pubmed.ncbi.nlm.nih.gov/39185205/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC11343214/

The source describes comparison of DNABERT-2, Nucleotide Transformer v2, and HyenaDNA embeddings across 57 datasets. It supplies a general evaluation framework and model-selection context. The version linked in PMC is a preprint; verify final publication status before manuscript citation. [web:347][web:351]

Implementation action: borrow reproducible extraction and comparison ideas, not cross-task superiority claims. Regulatory-genomics accuracy is not evidence of target-neighbor amplicon discrimination.

### L17. Unsupervised evaluation of pre-trained DNA language model…, 2025

Record: https://link.springer.com/article/10.1186/s12864-025-11913-2
DOI: 10.1186/s12864-025-11913-2.

This work provides model/representation comparison context, including architectures and embedding dimensions. Unsupervised separation or clustering should remain secondary to held-out recommendation utility. [web:348]

Implementation action: evaluate pooling, normalization, sequence-length effects, and dependence on the downstream classifier. Verify the full title/author export before citing formally.

### L18. Genome language modeling (GLM): a beginner’s cheat sheet, 2025

Record: https://academic.oup.com/biomethods/advance-article/doi/10.1093/biomethods/bpaf022/8093260
DOI: 10.1093/biomethods/bpaf022.

This is orientation material explaining tokenization, representations, and genomic language modeling. It is useful for methodological terminology, not primary evidence of performance on the proposed task. [web:267]

Implementation action: use for team onboarding; use the original model/task papers to justify experiments.

## 6. Learned 16S representations and adaptation

### L19. Bishop, Herbold, Dobson, and Ogilvie, 2026 — Micro16S

Title: Micro16S: Universal Phylogenetic 16S rRNA Gene Representations for Deep Learning of the Microbiome.
Status: bioRxiv preprint; version reviewed was the March 2026 v1 record.
DOI: 10.64898/2026.03.21.713432.
Record: https://www.biorxiv.org/content/10.64898/2026.03.21.713432v1.full

Micro16S combines triplet and pair losses to map 16S sequences into a phylogenetically structured space, promoting region invariance. It reports better region invariance than k-mer representations but limitations at lower taxonomic ranks. Its downstream microbiome sample classifier was generally outperformed by classical ML across the six tasks described. [web:318]

Important distinction: downstream disease/demographic sample classification is not the same endpoint as taxonomic classification, and neither is identical to amplicon selection. Do not collapse these result categories into 'embeddings do not work.'

Implementation action: examine whether region invariance is useful or counterproductive for preserving assay-specific resolution; include a same-assay versus cross-region-loss ablation. Treat published performance as preliminary until peer review and independent replication.

### L20. seqLens: Optimizing Language Models for Genomic Predictions, 2026

Record: https://academic.oup.com/mbe/article/43/7/msag139/8715817

The record reports LoRA-based adaptation of genomic models with substantially fewer trainable parameters than complete updates. This supports a staged adaptation strategy, not a guarantee that PEFT will outperform frozen features or fit an arbitrary GPU. [web:313]

Implementation action: profile actual checkpoint, token length, precision, and batch composition. Record trainable and total parameters separately; measure downstream assay utility and resource costs.

### L21. Optimizing Genomic Language Models for Efficient Training…, 2026

Record: https://openreview.net/forum?id=80aL5fuIzJ

The submission discusses efficient adaptation/evaluation of genomic encoders. Its record is useful for adaptation trade-offs; venue/acceptance status should be checked before describing it as peer reviewed. [web:319]

Implementation action: compare the efficiency/performance frontier. Do not use reduced trainable-parameter count as a substitute for measured memory, time, and independent accuracy.

## 7. Sequencing errors and empirical validation

### L22. Callahan et al., 2016 — DADA2

Title: DADA2: High-resolution sample inference from Illumina amplicon data.
Venue: Nature Methods 13: 581–583.
DOI: 10.1038/nmeth.3869.
Record: https://pmc.ncbi.nlm.nih.gov/articles/PMC4927377/pdf/nihms782534.pdf

DADA2 models/corrects amplicon errors and infers sequence variants. Its resolution makes it suitable for read-level validation, but an ASV is not automatically a species and does not remove marker ambiguity. [web:380]

Implementation action: freeze primer trimming, quality processing, error estimation, merging, and chimera policy. Retain expected template/copy identity for simulated reads and report sequence inference separately from taxonomic assignment.

### L23. Huang et al., 2012 — ART

Title: ART: a next-generation sequencing read simulator.
Venue: Bioinformatics 28(4): 593–594.
DOI: 10.1093/bioinformatics/btr708.
Record: https://academic.oup.com/bioinformatics/article/28/4/593/213322

ART simulates sequencing reads using empirical platform profiles. It is appropriate for controlled error stress tests, but defaults reflect particular technologies and do not automatically reproduce a current run or PCR bias. [web:374]

Implementation action: document profile origin and fit to the intended chemistry; keep clean-sequence, PCR, and sequencing-error experiments distinct.

### L24. Gourlé et al., 2019 — InSilicoSeq

Title: Simulating Illumina metagenomic data with InSilicoSeq.
Venue: Bioinformatics 35(3): 521–522.
DOI: 10.1093/bioinformatics/bty630.
Record: https://academic.oup.com/bioinformatics/article/35/3/521/5055123

InSilicoSeq supplies parallelized Illumina read simulation. It is a candidate error-generation component, not an end-to-end validated amplicon-PCR simulator. [web:376]

Implementation action: test whether the simulation tool handles short amplicon templates as intended and whether paired-end geometry is appropriate. Validate synthetic read length/coverage and the expected species/allele counts before benchmarking models.

### L25. Performance evaluation of six popular short-read simulators, 2022

Record: https://www.nature.com/articles/s41437-022-00577-3
DOI: 10.1038/s41437-022-00577-3.

The comparison discusses ART, DWGSIM, InSilicoSeq, Mason, NEAT, and wgsim and highlights simulator-specific assumptions and profiles. This supports assessing simulator limitations rather than treating simulated data as perfect experimental ground truth. [web:375]

Implementation action: compare at least one empirical validation subset and one alternative error scenario. Robustness across two simulators is useful but does not replace laboratory evidence.

## 8. Synthesis: what the proposed studies add

Protocol 1 could contribute a transparent target-versus-neighbor evaluation, shared-sequence/copy-aware ambiguity reporting, assay-specific constraints, independent recommendation utility, and a usable implementation. The contribution is strongest when it demonstrates better decisions or reliable abstention relative to Qscore-inspired and classical baselines.

Protocol 2 could contribute evidence on whether learned species-aware representations can be adapted to short conserved-marker discrimination and whether that improves actual region choice. New losses are not required for a useful study; a well-controlled transfer result or negative result is scientifically informative. A genuinely new objective must be motivated by a measured failure, not merely by the desire to claim algorithmic novelty.

## 9. Literature review maintenance plan

Before preregistration, update searches for: 16S region/amplicon selection; target-specific primer optimization; Qscore alternatives; 16S species-aware representation learning; short-amplicon genomic model adaptation; benchmark leakage and taxonomic transfer. Search PubMed, scholarly indexes, software repositories, and preprints separately.

For each candidate record capture title, authors, DOI/accession, year, publication status, code/checkpoints, input/target type, reference source, split policy, whether independent experimental validation exists, and direct overlap with the proposed claim. Include negative results and unavailable/failed software in the audit.

A bibliographic export is required before submission. Do not cite a secondary discussion as though it were the original model paper. Recheck abbreviated entries and preprint status. The source citations in this package are traceability aids, not a substitute for a final reference-manager library.
