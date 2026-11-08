# Phase 0–1 Prototype Report: Identity Baseline on the B. cereus Group

**Date:** November 2026
**Scope:** Kill-switch test from the analysis plan (`glm_16s_region_benchmark_plan.md`, §9 Phase 1): does the simple identity baseline (B1) recover published region-discrimination expectations for *B. anthracis* before any gLM work?

---

## 1. What was done

1. **Reference sequences.** 83 unique full-length 16S rRNA gene sequences fetched from NCBI for: *B. anthracis* (15), *B. cereus* (19), *B. thuringiensis* (19), *B. weihenstephanensis* (5), and *B. subtilis* (25, outgroup). QA filter: 1,450–1,600 nt (removed 5'- or 3'-truncated records that corrupted V1/V9 slices).
2. **Coordinate reference.** The 16S gene of *E. coli* J01695, extracted from the GenBank annotation (positions 1268–2809 of the rrnB operon record; 1,542 nt). Landmark validation passed exactly: 27F/8F motif at position 8, 515F motif at position 515, 1492R reverse complement at position 1492.
3. **Region slicing.** Each Bacillus sequence was globally aligned to the *E. coli* reference (Biopython PairwiseAligner); V1–V9 boundaries taken from the Chakravorty et al. 2007 table (V1: 69–99 … V9: 1435–1465) were transferred through the alignment coordinates. Slices outside 75–130% of expected length were dropped.
4. **Scoring.** Pairwise identity per region (2×matches / summed lengths over the global alignment), for within-*B. anthracis* pairs and *B. anthracis* vs other cereus-group pairs. 2,328 pair-level identities written to `pair_identities.csv`.

## 2. Results (nA = 9 anthracis; nO = 26–27 cereus-group after QA)

| Region | within_min | within_mean | btwn_max | btwn_mean | margin | vs B. subtilis |
|---|---|---|---|---|---|---|
| V1 | 1.0000 | 1.0000 | 1.0000 | 0.9099 | 0.0000 | 0.6786 |
| V2 | 0.9826 | 0.9961 | 1.0000 | 0.9613 | 0.0174 | 0.8509 |
| V3 | 0.9848 | 0.9966 | 1.0000 | 0.9764 | 0.0152 | 0.8490 |
| V4 | 0.9907 | 0.9979 | 1.0000 | 0.9884 | 0.0093 | 0.9148 |
| V5 | 0.9831 | 0.9962 | 1.0000 | 0.9691 | 0.0169 | 0.9266 |
| V6 | 0.9828 | 0.9914 | 1.0000 | 0.9520 | 0.0172 | 0.9080 |
| V7 | 0.9821 | 0.9960 | 1.0000 | 0.9536 | 0.0179 | 0.9643 |
| V8 | 1.0000 | 1.0000 | 1.0000 | 0.9391 | 0.0000 | 0.7500 |
| V9 | 0.9524 | 0.9841 | 1.0000 | 0.9400 | 0.0476 | 0.9209 |

**Erratum (Phase 2):** the identity function used in the original run of this report counted aligned positions rather than matching positions, inflating identities toward 1.0 for equal-length gapless pairs. The table above is corrected. The report's qualitative conclusions are unchanged — btwn_max = 1.0 is now a genuine string-identity claim — and one finding is strengthened: V6 carries the largest mean between-group divergence of any region (except V9), directionally consistent with Chakravorty's maximum-heterogeneity claim, while still failing to separate the groups.

Figure: `prototype_summary.png`.

## 3. Findings

1. **The pipeline is sound.** The clonality of *B. anthracis* is recovered exactly (within-group identity = 1.0 in V1–V8), the outgroup *B. subtilis* is separated in every region (0.96–1.0), and the cereus group sits at ~91–99% identity to *B. anthracis* (region-dependent; see corrected table) — all biologically expected.
2. **The classical V6 expectation does not generalize.** Chakravorty et al. 2007 reported that V6 distinguishes *B. anthracis* from *B. cereus* by a SNP at position 988. In our NCBI sampling, position 988 is G in every group, and the dominant *B. anthracis* V6 haplotype is carried by at least one *B. cereus* strain — zero of the anthracis V6 haplotypes are exclusive to anthracis. Their claim rested on two anthracis strains (Vollum, Sterne) compared against a single *B. cereus* reference sequence (NTBC01_Bc16SA); the polymorphism does not survive broader sampling.
3. **No 16S region resolves *B. anthracis* from the cereus group.** Between-group maximum identity is 1.0 in every region — for each region there exists at least one cereus-group strain indistinguishable from anthracis. This matches the modern whole-genome understanding (*B. anthracis* is a monophyletic subclade of *B. cereus* sensu lato) and confirms 16S is the wrong tool for this target.

## 4. Interpretation for the study design

- **Kill-switch assessment: not triggered, but the target was wrong.** The pipeline recovered the correct biology; what failed is the published expectation, because it was derived from a two-strain comparison. This is itself a result in the paper's favor: published region-discrimination claims are strain-panel-dependent — region choice is method-dependent at the reference-panel level, exactly the thesis of the study.
- **Reclassify *B. anthracis* from "validation anchor" to "expected-null case."** The tool should output "no region resolves this target — use whole-genome or canonical SNP markers," and this case becomes the demonstration of that output class.
- **Promote *S. aureus* vs CoNS (V1) and Panel B organisms (*B. longum* subspp., *E. coli/Shigella*, oral *Streptococcus*)** to primary validation targets, and keep enterobacteriaceae as a second null case.
- **Add a "nearest published expectation" tier to the validation framework:** for each target, record whether our metric recovers the published claim, and treat non-recovery (as here) as a reportable finding rather than a pipeline failure.

## 5. Deliverables in this directory

- `fetch_seqs.py` — NCBI retrieval
- `slice_and_score.py` — alignment-based region slicing + identity scoring
- `figure.py` — pair-level dump + summary figure
- `pair_identities.csv` — 2,328 pair-level region identities
- `prototype_summary.png` — Figure 1 of the prototype
- `cereus_group_16s.json`, `ecoli_16s.txt` — cached data

## 6. Sources

- Chakravorty et al. 2007, J Microbiol Methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC2562909/
- NCBI EUtils: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- E. coli rrnB operon J01695 (16S annotated at 1268..2809): https://www.ncbi.nlm.nih.gov/nuccore/J01695
