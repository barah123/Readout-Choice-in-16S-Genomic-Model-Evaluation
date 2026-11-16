# Protocol amendment 1: regioneer after the gLM benchmark

**Status:** proposed amendment to `01_Regioneer_Benchmark_and_Tool_Protocol.md` v1.0.
**Date:** 2026-11-16.
**Trigger:** completion of the gLM 16S region benchmark, including its post hoc
sequence-disjoint, similarity-disjoint, tuned-readout and divergence-stratified analyses.

This document does not rewrite the protocol. It records which protocol questions the gLM
study has already answered, which it has made harder, and which it leaves untouched, and
proposes the specific section edits that follow. Every number cited is from a deposited
result file in the gLM repository, named at the point of use.

---

## 1. The headline

**The protocol anticipated this outcome and said what to do about it.** H3 already declined
to assume that frozen representations beat classical methods; H1 and H2 were written to be
recoverable "even if H3 fails"; Gate 2 required competitive k-mer controls; §9.3 anticipated
exact-sequence collisions; and §18 already commits: *"If classical methods win: publish the
benchmark and deploy the simpler method."* gLM has run H3 ahead of schedule and returned a
clear answer.

What changes is not whether regioneer is viable but **what it is for**. The ranking function
is weaker than assumed and the resolvability function is stronger. The evidence supports
moving the second to the centre of the tool and demoting the first.

---

## 2. Evidence base

| Finding | Value | Source |
|---|---|---|
| Frozen DNABERT-2 vs k-mer, locked test | −0.21 pp, 95% CI [−0.59, +0.16], TOST p = 4.1 × 10⁻⁴ | `locked_test_results.json` |
| Random 256-dim projection of k-mer counts | 96.2%, above all four models | `locked_test_scores.json` |
| Exact-sequence overlap, held-out vs training | **92.0%** of evaluated isolates | `sequence_disjoint_results.json` |
| Accuracy cost of removing that overlap | **≈ 11 pp**, both arms | same |
| Contrast under sequence-disjoint folds | −0.50 pp [−3.20, +2.07], **inconclusive** | same |
| Do embeddings generalise better to unseen sequence? | −0.30 pp [−2.92, +2.22], p = 0.837; better in 3 of 14 blocks | same |
| Block-level SD, isolate-grouped → sequence-disjoint | **0.37 pp → 5.16 pp (14×)** | same |
| 99%-identity similarity-disjoint | −1.47 pp, same verdict, audit max cross-fold identity 0.9894 | `similarity_disjoint_results.json` |
| Region-ranking regret, locked test | k-mer 0.65 pp, DNABERT-2 0.57 pp; paired −0.09 pp [−0.46, +0.26] | `region_ranking_summary.json` |
| Best-region hit rate | **37% (k-mer), 44% (DNABERT-2)** | same |
| Difficulty vs divergence | near-binary on shared haplotypes: 13-pt swing between strata, method gap ≤ 0.25 pp in every stratum | `divergence_stratified_results.json` |
| Evaluated span length | **31–107 nt**, median 58 | frozen region coordinates |

---

## 3. Disposition of the hypotheses

**H3 — frozen embeddings may improve some panels.** *Answered: no, for frozen
representations on short regions.* Equivalence on the locked test; no advantage on
difficulty-selected panels; and, decisively, no slower degradation when exact matches are
removed, which was the one mechanism that would have justified them. A random projection of
k-mer counts matched every learned representation. **Fine-tuning (protocol 2) is now the only
live route for learned representations**; frozen embeddings should not be on regioneer's
critical path.

**H2 — estimates deteriorate under held-out evaluation.** *Confirmed, and it is the dominant
effect in the whole study.* Not a caveat but the main result: ~11 points, with variance
rising 14-fold.

**H1 — the preferred assay varies across panels.** *Untouched and still live.* gLM fixed one
region per block for most analyses and did not test assay-level variation across panels.

**H4 — amplification constraints change the preferred assay.** *Untouched.* gLM modelled no
primers, no amplification, no read layout.

---

## 4. Proposed section amendments

### §3 Hypotheses — restate H3 as closed
Record that H3 has been tested externally for frozen representations on single short regions
and was not supported, cite the gLM benchmark, and narrow the remaining H3 to adapted or
fine-tuned representations under protocol 2. Do not re-run the frozen-embedding comparison as
a primary arm.

### §9.2 Split hierarchy — promote sequence-disjointness from add-on to primary
§9.3 currently treats a sequence-novelty challenge set as an *addition* to a realistic
primary evaluation. gLM shows the two regimes give different verdicts: equivalence holds under
isolate grouping and becomes inconclusive once exact sequences are separated. A design in
which the primary analysis retains 92% sequence overlap cannot distinguish a method that
generalises from one that looks up.

Proposed: report **both regimes as co-primary**, pre-declared, with the sequence-disjoint
regime carrying the generalisation claim and the isolate-grouped regime carrying the
operational claim. Single-linkage clustering on the identity graph is the construction that
guarantees every cross-fold pair falls below threshold; greedy or centroid clustering does
not. Audit the realised split directly — gLM's audit reports maximum cross-fold identity and
verifies that an exact-match classifier scores exactly 0.500 under a correct disjoint split.

### §10.2 Frozen learned representations — demote to a reference arm
Keep one frozen encoder as a reference, drop the rest. The full checkpoint/pooling/strand
apparatus in §10.2 is disproportionate to an arm now known to be equivalent to k-mer counting
and to a random projection of it. **Retain the random-projection control**: it is the cheapest
and most informative single comparator gLM identified, and it bounds any claim that a learned
space is doing representational work.

### §10.4 Candidate ranking rule — the tie tolerance was right, the ranking is weak
§10.4's prespecified tie tolerance of 0.01 accuracy units is vindicated: regions tie often
enough that a method can record near-zero regret while matching the single best region well
under half the time. But that same fact undercuts the ranking as a product. Report **regret
and hit rate together, always**, and state plainly that a low hit rate with low regret means
the choice rarely matters rather than that the recommender is accurate.

### §11 Endpoints — make resolvability primary, ranking secondary
§11.3's ambiguity ceiling and collision tables are currently supporting analysis. The
divergence stratification says this is where the information is: performance is governed by
whether target and neighbour share sequence at the evaluated assay, not by graded divergence.

Proposed reordering:
- **Primary:** per-assay resolvability — empirical ceiling, target-neighbour collision table,
  and an explicit abstention when no candidate clears a prespecified criterion.
- **Secondary:** U(p, â) and recommendation regret, retained as now.

Note also that the ceiling must be aggregated the same way as the accuracy it bounds. A
whole-population ceiling does not bound a mean of per-fold accuracies; in gLM 3 of 66
block-by-method cells exceeded one, while none exceeded the fold-matched version.

### §14 Precision planning — the block counts are low by two orders of magnitude
§14 proposes examining "realistic block counts such as 10, 20, 40, 60". Those are adequate
under isolate grouping and not remotely adequate under sequence-disjoint evaluation. Using
gLM's measured block-level SDs, n = (1.96·SD/h)²:

| 95% CI half-width | Isolate-grouped (SD 0.37 pp) | Sequence-disjoint (SD 5.16 pp) |
|---|---|---|
| ±1.00 pp | 1 | **102** |
| ±0.75 pp | 1 | **182** |
| ±0.50 pp | 2 | **409** |
| ±0.25 pp | 9 | **1,637** |

**If the sequence-disjoint regime is co-primary, regioneer needs on the order of 100 blocks
for a ±1 pp interval, not 40–60.** This is the single most consequential planning consequence
and it should be settled before curation scope is fixed, not after. Note the SD is itself
estimated from 14 blocks and is therefore uncertain; treat 100 as a floor.

A second feasibility constraint: in gLM, only 14 of 20 blocks could support exact-disjoint
folds at all, and 9 of 20 at 99% identity, because the excluded blocks had ≤4 distinct target
haplotypes. Panel curation must therefore select for **haplotype diversity**, not only isolate
count, or blocks will be unusable in the regime that matters.

### §16 Gates — add one
Proposed Gate 6: *no generalisation claim from an evaluation whose train/test exact-sequence
overlap has not been measured and reported.* gLM's experience is that this quantity is large,
invisible by default, and decisive.

### §18 Interpretation rules — the first rule now applies
"If classical methods win: publish the benchmark and deploy the simpler method." That branch
is now taken for frozen representations. The practical consequence is a **pure-classical R
package**: no torch, no GPU, no model checkpoints, no `trust_remote_code`, far lighter
dependencies and a far easier release. This is a gain for a translational deliverable, not a
loss.

---

## 5. What gLM cannot settle

**Region length is the big one.** gLM evaluated strict hypervariable spans of 31–107 nt
(median 58). regioneer's unit is the assay — primer pair, trimming, read layout — i.e. real
amplicons of roughly 253 nt (515F-806R V4) to 460 nt (V3–V4). At 31–107 nt a canonical
4/5/6-mer profile approaches a lossless encoding of the sequence, so the classical baseline
sits near the information ceiling **by construction**. That argument weakens as amplicons
lengthen, and the gLM equivalence result should not be assumed to transfer. This is a genuine
open question and regioneer's design is the right one to answer it.

Also untouched: amplification and primer-coverage constraints (H4); full-length versus short
assays (secondary Q3); intragenomic copy variation with genome-linked copies, which gLM could
not address because standalone 16S deposits carry no parent-genome link and Tier A complete
genomes fix properly; and open-panel rejection.

---

## 6. Net effect on the deliverable

| | v1.0 | After this amendment |
|---|---|---|
| Headline function | rank candidate assays | report resolvability and abstain when unresolvable; rank secondarily |
| Learned representations | main arm | one frozen reference arm; adaptation deferred to protocol 2 |
| Evaluation regime | isolate-grouped primary, novelty as challenge set | both co-primary, sequence-disjoint carries generalisation |
| Block target | 10–60 | ~100+ for ±1 pp, curated for haplotype diversity |
| Runtime dependencies | torch, checkpoints, GPU profiling | pure classical |
| Scientific contribution | a better recommender | a leakage-resistant benchmark that establishes when 16S *cannot* resolve, plus the evidence for not using a language model |

---

## 7. Decisions required

1. Accept or reject sequence-disjoint evaluation as co-primary. **Everything in §14 depends on
   this**, and the honest reading of gLM is that a 92%-overlap primary cannot support a
   generalisation claim.
2. Confirm the block-count target against curation capacity. If ~100 blocks is not reachable,
   the sequence-disjoint arm should be declared descriptive from the outset rather than
   reported as underpowered after the fact.
3. Confirm the shift of resolvability to primary and ranking to secondary.
4. Decide whether a frozen encoder is retained at all, or whether the gLM benchmark is simply
   cited and the arm dropped.
5. Decide whether to run a short amplicon-length pilot before committing, since that is the
   one assumption gLM leaves genuinely open and it conditions the whole rationale.
