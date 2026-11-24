# Paper plan: why 16S region recommendations do not transfer between settings

**Status:** plan, not a draft. Every claim below is mapped to evidence already in hand.
**Date:** 2026-11-23.

---

## The claim

Published recommendations about which 16S region to sequence disagree with each other, and
the disagreement is not a failure of any individual study. It follows from how the task
works: whether two organisms can be told apart depends on whether they happen to share an
identical sequence in the region being used, and which organisms are present differs between
settings. A region recommendation is therefore a statement about a particular set of
organisms, not a property of the region, and it cannot be transferred to a different set.

Two further results bound how confidently anyone can make such a recommendation at all:
standard evaluations overstate how well methods work, and reference databases are often too
thin to support the calculation.

---

## What this paper is and is not

**It is** an explanation of a disagreement that already exists in the literature, with
measurements that supply the mechanism, plus a tool that computes the answer for a given set
of organisms.

**It is not** a demonstration that region recommendations disagree. The literature already
shows that, and three papers are cited for it. Our data covers 20 bacterial groups of
clinical and environmental interest, not a cross-habitat design, so we explain the
disagreement rather than reproduce it. **This distinction must be explicit in the abstract.**
Claiming to have shown non-transport ourselves would be a claim our design does not support.

---

## Structure, with evidence per section

### 1. The disagreement exists

Three studies, three different recommendations:

| Setting | Recommendation | Source |
|---|---|---|
| Plant-associated genera, 16 genera, 6,821 genomes | V1–V3 over V3–V4 | Costa et al. 2024, ISME Commun, `10.1093/ismeco/ycae034` |
| Respiratory, 33 sputum samples + mock | V1–V2 highest resolving power | Sci Rep 2023, `10.1038/s41598-023-30764-z` |
| Vaginal, in silico + experimental | V4 optimal at species level | J Infect Dis 2018, `10.1093/infdis/jiy508` |

Costa et al. also report that the best region differs *between genera within their own study*,
and describe their method as a blueprint for choosing a region per genus. Cite this as
agreeing with us, not competing — they observe the pattern; we give the mechanism.

*Caveat to carry:* the respiratory result is weak (AUC 0.736, interval 0.566–0.906, 10-species
mock). Say so rather than stacking it as equal evidence.

### 2. The mechanism: shared sequences, not region quality

| Claim | Evidence | Source file |
|---|---|---|
| Performance is a property of the organism pair, not the method | block identity explains 98.41% of variance, method identity 0.30% | `step1c`, variance decomposition |
| What governs it is the ambiguity ceiling | ceiling predicts block accuracy at R² = 0.96; a raw collision indicator gives 0.23 | §2A analysis |
| Every method sits just below that ceiling | accuracy = ceiling − 2.11 pp (SD 2.48), method-independent | `fold_matched_ceiling_audit.csv` |
| Difficulty is near-binary, not graded | 63 target–neighbour pairs: those sharing a haplotype score 85.9%, those differing at all score 97–99%; the method gap stays ≤0.25 pp in every stratum | `divergence_stratified_results.json` |

This is the core of the paper. The last row is the figure.

### 3. Why the region matters, and why "longer" is the wrong rule

| Claim | Evidence |
|---|---|
| Region content decides collisions, not length | V1–V2 at 294 nt resolves better than V3–V4 at 429 nt (collision 0.067 vs 0.218; ceiling 98.55 vs 96.56) |
| The field's two default amplicons are the weaker choices | V4 and V3–V4 are the worst of five tested |
| Longer spans do split some shared sequences, but not the hard ones | of 8 collision-heavy blocks, 2 resolved, 1 improved, 5 unchanged |

Present the V1 result **as a replication** of Costa et al. in a different organism set, not as
a new finding. Agreement across two very different designs is the point.

### 4. Standard evaluations overstate performance

| Claim | Evidence |
|---|---|
| Grouping by isolate does not make folds sequence-disjoint | 92.0% of held-out isolates share their exact sequence with training |
| Removing that overlap costs about 11 points | k-mer 96.36 → 85.42; DNABERT-2 96.16 → 84.92 |
| And it removes the ability to distinguish methods | contrast −0.50 pp, interval [−3.20, +2.07]; block SD rises 0.37 → 5.16 |
| A literal lookup table scores 90.3% under the standard design and 50.0% without reuse | `C_exact` |
| Verified, not assumed | zero train/test overlap and exact matching scoring exactly 0.500 in every fold |

### 5. Reference data is often too thin to answer

| Claim | Evidence |
|---|---|
| A naive calculation over-promises resolvability | wrong 20.5% of the time at 2 organisms per species, 7.9% at 20 |
| The error is predictable from depth but not correctable | log-log fit R² = 0.98; coverage does not predict individual error (R² = 0.006); a depth-keyed band does not hold out of sample (0.70–0.83 against nominal 0.95) |
| A 10% tolerance needs about 12 organisms per species | from the five-point curve |
| Complete genomes rarely reach that | median 9 per species; 39.4% at or above 12 |
| 16S references do | median 64; 85.1% at or above 12 |
| Draft genomes do not help | median 59, but drafts lose 16S copies to repeat collapse |

*Caveat:* these species passed an eligibility screen requiring ≥10 isolates, so all counts are
upper bounds.

### 6. What follows

- A region recommendation is only meaningful relative to a declared set of organisms.
- It should come with a refusal when reference depth is inadequate, at a stated tolerance.
- The copy-aware version of the calculation is not supported by available data for most
  species, and should be restricted to the ~40% that have enough complete genomes.

---

## Figures

1. The three published recommendations side by side, with the organism sets that produced them.
2. Accuracy against ambiguity ceiling, both arms, fold-matched. *(exists: Fig 7)*
3. Difficulty by target–neighbour divergence, showing the near-binary structure and the flat
   method gap. *(exists: Table 10 data)*
4. Collision fraction and ceiling by span, showing V1-containing spans against V4/V3–V4.
   *(exists: `step1_amplicon_spans.json`)*
5. Error rate of a naive resolvability call against reference depth, with the tolerance
   thresholds marked. *(exists: `step2_reference_incompleteness.json`)*
6. Available reference depth per species by source, against the threshold from figure 5.
   *(regenerating: `reference_depth_by_species.json`)*

Figures 2–6 are already computed. Figure 1 is a literature summary.

---

## What is deliberately left out

- The language-model comparison. It is a negative result at seven spans and belongs in the
  gLM manuscript, which already carries it. Mentioning it here invites the paper to be read
  as a model-benchmarking paper.
- The ranking endpoint. 37–44% hit rate with 0.57–0.65 pp regret says the ranking is both
  unreliable and unimportant. It is a sentence in the discussion, not a result.
- Any calibrated prediction of resolvability. It failed its test; §5 reports the failure as a
  limit, not a method.

---

## Honest weaknesses to state in the paper

1. Twenty bacterial groups, one organism set. We explain non-transport; we do not demonstrate
   it across habitats.
2. The subsampling treats the full reference as truth, but the full reference is itself
   incomplete. The error rates are floors.
3. Species labels come from submitter metadata, not verified provenance.
4. Conditional endpoint throughout: no primers, no amplification, no read layout. Everything
   assumes the sequence was observed.

---

## Next step

Decide the target journal before drafting, because it sets the length and whether the tool is
a supplementary artifact or a co-equal contribution. The content supports either a methods
paper or a resource paper, and those are different documents.
