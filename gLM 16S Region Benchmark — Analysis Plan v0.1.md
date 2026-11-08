# Analysis Plan: A Genomic Language Model Benchmark for Selecting 16S rRNA Hypervariable Regions for Target Organisms

**Version:** 0.1 (draft protocol)
**Positioning:** Extends the systematic review's method-dependence thesis (breastfeeding vs formula feeding; composition robust, function method-dependent) with a primary-data study: the profiling region itself is a method choice, and we aim to predict its fitness for a target organism before sequencing.

---

## 1. Study Overview

### 1.1 Primary hypothesis (H2 from prior discussion)

For a given target organism (or set of closely related organisms), the discriminative quality of each 16S rRNA hypervariable region (V1–V9) can be predicted from pretrained genomic language model (gLM) embeddings of that region, and these predictions generalize to empirical taxonomic resolution observed in published region-comparison studies and mock communities.

### 1.2 Secondary hypothesis (H1, embedded as the benchmark framework)

gLM-based region scoring outperforms (or matches, with better transfer across taxa) simple sequence-based discriminability metrics — pairwise identity, k-mer divergence, and Shannon entropy — which serve as the classical baselines. The benchmark itself is the evaluation framework; a null result (gLM does not beat baselines) is a reportable finding consistent with the Micro16S precedent (classical ML outperforming a genomic LM on 16S benchmarks; https://www.biorxiv.org/content/10.64898/2026.03.21.713432v1.full.pdf).

### 1.3 Outcome

An R package/workflow (`regioneer` working title): input = target organism(s) + optional nearest-neighbor set; output = ranked hypervariable regions with (a) embedding-separability score, (b) baseline identity/entropy scores, (c) in silico primer coverage, and (d) recommended region + primer pair.

---

## 2. Reference Data

### 2.1 Full-length 16S reference sequences

| Source | Use | Notes |
|---|---|---|
| SILVA SSU 138.3 (https://www.arb-silva.de/) | Primary. Aligned reference sequences; slicing aligned columns gives consistent region boundaries | Non-redundant 99% and 100% OTU files |
| GTDB (https://gtdb.ecogenomic.org/) | Taxonomy ground truth; representative genomes | Extract rRNA from representative genomes with barrnap |
| Living Tree Project LTP (https://imedea.uib-csic.es/microbiology/ltp/) | Curated type-strain sequences for species-level panels | High-quality, type-strain anchored |
| EZBioCloud 16S (https://www.ezbiocloud.net/) | Optional cross-check for species-level identity thresholds | |

**Key technical advantage of SILVA:** sequences are aligned to a fixed reference alignment. Region boundaries (V1–V9, in *E. coli* coordinates) correspond to alignment column ranges, so regions can be extracted consistently across all sequences — no per-sequence primer matching needed for the discriminability analysis (primer matching is a separate step, §7).

### 2.2 Target organism panels

**Panel A — Classical discriminability benchmark (from Chakravorty et al. 2007, https://pmc.ncbi.nlm.nih.gov/articles/PMC2562909/):**
- *Bacillus anthracis* vs *B. cereus* group (V6 known to separate by a single polymorphism)
- *Yersinia pestis* vs *Y. enterocolitica* (V1V2 and V6V7V8 separate; V4/V5/V7/V8 weak)
- *Staphylococcus aureus* vs coagulase-negative staphylococci (V1 best)
- *Mycobacterium* spp. (V2 best per Chakravorty)
- *Haemophilus* spp. (V3 best per Chakravorty)
- Enterobacteriaceae close pairs (known hard case for all regions)

**Panel B — Microbiome-relevant (ties to your review):**
- *Bifidobacterium longum* subsp. *infantis* vs *longum* vs *suis* (infant gut; the organism-level story of breastfeeding vs formula feeding)
- *Lactobacillus*/*Limosilactobacillus* vaginal species group
- *Escherichia*/*Shigella* (gut)
- *Streptococcus* oral group (e.g., *S. mitis*/*pneumoniae* complex — a canonical unresolvable-by-16S case)
- *Cutibacterium acnes* phylotypes (skin; V1–V3 shown comparable to full-length, https://journals.asm.org/doi/10.1128/msystems.00399-24)

Panel A provides the validation anchor (published expectations), Panel B provides novelty and direct relevance to your review's population (healthy term infants).

### 2.3 Empirical validation datasets (mock communities and published region comparisons)

- ZymoBIOMICS microbial community standard (even and log-distributed; known composition)
- QIIME2 mockrobiota (https://docs.qiime2.org/ — mock communities with full-length expected taxonomy)
- Published full-length 16S datasets (PacBio/ONT) from which sub-regions can be sliced in silico: skin microbiome full-length vs sub-region study (https://journals.asm.org/doi/10.1128/msystems.00399-24); long vs short amplicon comparison (https://pmc.ncbi.nlm.nih.gov/articles/PMC13712307/)
- V1–V2 sputum region study (https://www.nature.com/articles/s41598-023-30764-z) and respiratory taxa (V1–V2 highest resolving power)
- Plant-associated genera region comparison (V1–V3 best for half of analyzed genera; https://academic.oup.com/ismecommun/article/4/1/ycae034/7624640)

---

## 3. Region Extraction (in silico)

1. Download SILVA aligned SSU (e.g., `SILVA_138.3_SSURef_Nr99_tax_silva.fasta` + alignment).
2. Define V-region column ranges using *E. coli* J01695 coordinates as in Chakravorty 2007 (V1: 69–99, V2: 137–242, V3: 433–497, V4: 576–682, V5: 822–879, V6: 986–1043, V7: 1117–1173, V8: 1243–1294, V9: 1435–1465; verify against current alignment).
3. Slice each region (plus each common two-region amplicon: V1–V2, V2–V3, V3–V4, V4–V5, V1–V3, V6–V9) for every reference sequence.
4. QA: distribution of region lengths, gap fraction, N per region per taxon.

Deliverable: a tidy data object (one row per sequence × region) consumed by both the baseline and gLM arms.

---

## 4. Models and Baselines

### 4.1 gLMs to test

| Model | Params | Rationale | Practical notes |
|---|---|---|---|
| DNABERT-2 (https://huggingface.co/zhihan1996/DNABERT-2-117M) | 117M | Multi-species pretraining; BPE tokens; modest GPU needs; strong on GUE benchmark tasks | Best first candidate; easy `transformers` inference |
| HyenaDNA (medium-160k, https://huggingface.co/crabfishable...) | ~7M–60M | Very long context; single-nucleotide resolution | Trained on human reference genome — test transfer to bacterial 16S explicitly |
| Nucleotide Transformer v2 (https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-ref) | 500M | Pretrained on 3,202 human + 850 diverse species genomes | Context ~6k nt, fine for single regions |
| DNABERT-S | — | Contrastive species-aware embeddings; designed for metagenomic binning | https://arxiv.org/html/2407.11435v2 |

**Excluded with justification (write this in the paper):** gLM (Hwang et al., https://pmc.ncbi.nlm.nih.gov/articles/PMC10991518/) and gLM2 operate at contig scale with protein-level tokens (operon structure, co-regulation), not short amplicon slices — not directly applicable. Evo/Evo 2 are generative whole-genome models; computationally disproportionate for a ~500 nt scoring task.

### 4.2 Baselines (must-have; reviewers will demand them)

- B1: Pairwise identity — for each region, within-target mean identity vs min target-to-nearest-neighbor identity; discriminability margin = (min between-group identity − max within-group identity). This is the classical in silico approach.
- B2: k-mer frequency (4–6-mer) features + random forest / logistic classification per region.
- B3: Per-site Shannon entropy per region (region-wide and target-specific).
- B4: Simple heuristic — region length and hypervariable-site density.

If the gLM does not beat B1, report it. That result is consistent with Micro16S (classical ML beating a genomic LM) and strengthens the benchmark's credibility.

---

## 5. Embedding-Based Region Scoring

### 5.1 Extraction

- For each sequence × region, compute embedding via mean pooling over token embeddings (last layer; also record layer-wise scores — gLM papers show lower layers sometimes carry more functional signal, https://pmc.ncbi.nlm.nih.gov/articles/PMC10991518/).
- Python side (`transformers`, GPU batch inference), outputs parquet/feather.

### 5.2 Scoring metrics (computed in R)

For each (target group, region, model):
1. **Separation margin:** silhouette width using target-vs-neighbor labels.
2. **Nearest-neighbor classification accuracy** (leave-one-out) — does the nearest embedding neighbor of each target sequence belong to the target group?
3. **Margin ratio:** distance to nearest non-target / distance to nearest target.
4. **Embedding rank:** rank regions by metrics 1–3; compare with baseline ranks (Spearman).
5. **Phylogenetic congruence:** correlate embedding distance matrices with GTDB phylogenetic distance per region (Mantel test).

### 5.3 Head-to-head validation logic

- **Gold standard A (published expectations):** For Panel A organisms, our region ranking should recover published results (V6 separates *B. anthracis*/*B. cereus*; V1 separates *S. aureus*; V4/V5/V7/V8 weak; enterobacteriaceae unresolvable). Metric: top-k agreement, weighted Kendall tau.
- **Gold standard B (mock communities):** Slice full-length mock reads into regions, run taxonomy assignment (DADA2 + SILVA GTDB classifier), measure per-region assignment accuracy vs known composition; test whether our gLM score predicts empirical accuracy across regions and organisms (regression, Spearman).
- **Gold standard C (habitat studies):** Predict best region for plant-associated and sputum taxa from embeddings; compare with published conclusions (V1–V3 plant; V1–V2 sputum).

---

## 6. Statistical Analysis Plan (R)

- Primary analysis: Spearman correlation between gLM region scores and empirical taxonomic-resolution ranks, per panel and pooled; bootstrap 95% CIs over sequences.
- Model comparison: gLM vs B1–B4 via paired comparisons over (organism × region) cells; report effect sizes, not just p-values.
- Multiple testing: Benjamini-Hochberg across region × organism cells for "region X resolves target Y" claims.
- Pre-specify a null-result interpretation: gLM matching but not exceeding B1 is acceptable if it generalizes across taxa where B1 requires per-taxon threshold tuning.
- Power: reference panels are deterministic (no sampling), but mock-community validation uses bootstrap over reads; report uncertainty.

---

## 7. Primer Coverage Layer (the "can you actually sequence it" step)

Embedding discriminability ignores amplification bias. Integrate:

1. TestPrime (SILVA; https://www.arb-silva.de/search/testprime/) or PrimerProspector (https://github.com/capsidsgit/primerprospector) coverage of standard primer pairs (27F/338F/515F/806R/1492R etc.) for each target group.
2. Report a combined recommendation matrix: `discriminability score × primer coverage = effective score`.
3. Known failure cases to include as tests: *Bifidobacterium* with 27F (limited amplification); EMP 515F/806R covering ~83% of bacterial/archaeal sequences (https://www.technologynetworks.com/genomics/articles/16s-rrna-sequencing-principles-design-and-limitations-416582).
4. Compare recommendations against KuafuPrimer (few-shot ML primer design; https://www.biorxiv.org/content/10.64898/2026.03.29.714677v1) where their habitat-level conclusions overlap our organism-level ones.

---

## 8. Workflow Architecture

```
├── 00_setup/            # environments: python (transformers) + R (tidyverse)
├── 01_data/             # SILVA, GTDB, LTP downloads; panel definitions (YAML)
├── 02_regions/          # slicing script (Python or R + Biostrings)
├── 03_baselines/        # B1–B4 (R)
├── 04_embeddings/       # DNABERT-2 / HyenaDNA / NTv2 inference (Python)
├── 05_scoring/          # separability metrics (R)
├── 06_validation/       # mock communities, published-study recovery (R)
├── 07_primers/          # coverage layer (external tools, parsed in R)
└── 08_package/          # R package wrapper: target → ranked regions report
```

R does: region slicing QA, all statistics, visualization, package. Python does: model inference only. Glue via feather/parquet.

---

## 9. Phased Timeline (suggested)

| Phase | Work | Milestone |
|---|---|---|
| 0 | Environments, SILVA download, region extraction, QA | Reproducible region table |
| 1 | Panel A baselines (B1–B3) | Recover Chakravorty expectations with B1 (sanity check) |
| 2 | DNABERT-2 embeddings + scoring on Panel A | First gLM vs baseline comparison |
| 3 | Add HyenaDNA + NTv2; Panel B (incl. *B. longum* subspp.) | Full model × panel matrix |
| 4 | Mock-community empirical validation | Score-predicts-accuracy regression |
| 5 | Primer coverage layer; recommendation matrix | Tool output v1 |
| 6 | R package + manuscript | Submission draft |

**Kill-switch criteria (pre-specify):** if Phase 1's B1 baseline cannot recover published expectations, fix region extraction before any gLM work; if Phase 2 shows gLM embeddings fail to separate Panel A groups that B1 separates, stop model expansion and pivot the paper to the benchmark + negative finding.

---

## 10. Known Pitfalls and Mitigations

| Pitfall | Mitigation |
|---|---|
| gLM embeddings reflect pretraining taxonomy bias (e.g., model memorizes *E. coli*) | Hold out one Panel A group entirely; test transfer on Panel B |
| Mean pooling washes out informative positions | Add max pooling and attention-weighted pooling variants |
| SILVA alignment gaps make region slices uneven | Trim columns with >50% gaps; QA length distributions |
| Embedding reconstruction privacy issues (Evo 2/NTv2 invertible) | Not a data risk here (public references), but cite as consideration (https://juleskreuer.eu/research/dna-embedding-inversion) |
| Reviewers: "why not just use alignment?" | B1 is that argument — we test it head-to-head; frame gLM value as transferability across taxa without threshold tuning |
| Single-region studies can't resolve species regardless of region (identical V3–V4/V4 sequences between species) | Treat as expected nulls; the tool should output "no region resolves this target — use shotgun/whole-genome" |

---

## 11. Manuscript Framing

- Title candidate: "Predicting the informative 16S rRNA hypervariable region for target organisms: a genomic language model benchmark"
- Tie-in paragraph: profiling region choice is a method variable, exactly the class of method-dependence your systematic review identifies (16S-inferred function vs metagenomics). This study makes that dependence predictable a priori.
- Target journals: *Microbiome*, *mSystems*, *Bioinformatics*, or *ISME Communications*.

---

## 12. Key References

- Chakravorty et al. 2007, region discriminatory power: https://pmc.ncbi.nlm.nih.gov/articles/PMC2562909/
- Hwang et al. 2024, gLM (contig-scale): https://pmc.ncbi.nlm.nih.gov/articles/PMC10991518/
- DNABERT-S / gLM2 survey context: https://arxiv.org/html/2407.11435v2
- Micro16S (classical ML beats genomic LM on 16S): https://www.biorxiv.org/content/10.64898/2026.03.21.713432v1.full.pdf
- KuafuPrimer (ML primer design): https://www.biorxiv.org/content/10.64898/2026.03.29.714677v1
- Skin full-length vs sub-regions: https://journals.asm.org/doi/10.1128/msystems.00399-24
- Plant-associated region comparison: https://academic.oup.com/ismecommun/article/4/1/ycae034/7624640
- 16S sequencing principles/limitations (primer coverage): https://www.technologynetworks.com/genomics/articles/16s-rrna-sequencing-principles-design-and-limitations-416582
- Foundation models for microbiome research (evaluation framework): https://pmc.ncbi.nlm.nih.gov/articles/PMC13573924/
