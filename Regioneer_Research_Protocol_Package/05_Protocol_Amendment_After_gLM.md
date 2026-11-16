# Protocol amendment 1 (rev 5): re-ranking Protocol 1 after the gLM pilot

**Status:** proposed amendment to `01_Regioneer_Benchmark_and_Tool_Protocol.md` v1.0.
**Date:** 2026-11-16. **Supersedes revs 1–4** of this file.

**New in rev 5:** Step 1c has been run — the k-mer/encoder comparison AT amplicon length.
Equivalence transfers; §2C records it and §9's demotion is now supported at regioneer's own
unit rather than by extrapolation. Rev 4's Finding 3 is correspondingly narrowed.

**Rev 4 added:** the amplicon-span check has been RUN (§2B), and it changes three things —
the best candidate assay is not the one the field uses, the ICC at assay spacing is worse
than at region spacing, and headroom grows with length so the gLM equivalence cannot be
assumed to transfer. §7 and the execution plan (§12) are updated accordingly.

**What changed in rev 3.** A viability appraisal proposed moving the study's centre of
gravity from the method contrast to resolvability, and supplied a decisive check. That check
has now been run on deposited gLM results (§2A). It confirms the reframing by a wider margin
than the appraisal claimed, corrects which quantity carries it, and — against the appraisal's
expectation — shows the reframing does **not** dissolve the sample-size problem. §6 and §7 are
rewritten accordingly.

**What changed in rev 2, and why.** Rev 1 read as a corrective to v1.0. That framing was
wrong on three counts, all identified in external appraisal and verified here:

1. **It used a secondary-endpoint result to restructure the primary endpoint.** v1.0 §11.1
   defines U as *end-to-end*, with failed products counted incorrect, and states in the same
   paragraph that conditional sequence classification "is a secondary endpoint; do not let it
   replace the end-to-end result." The gLM benchmark modelled no primers, no amplification and
   no read layout. It measured that secondary endpoint, on 31–107 nt spans that are not assays
   under §5.
2. **Most of its proposals already exist in v1.0**, as secondary or supporting clauses. The
   honest operation is promotion, not correction.
3. **It misstated v1.0's block target** as "10–60". §14's 10/20/40/60 are counts to *examine
   in a precision simulation*; the curation target is §7.1's "approximately 30–60
   non-overlapping taxonomic blocks."

Rev 2 is therefore a **re-ranking memo**. For each item: what v1.0 already says, what the gLM
pilot shows, and what should move. Two items require real debate — the block-count estimate
and the learned arm — and they are marked as such.

---

## 1. What the gLM pilot does and does not license

**Licensed, endpoint-independently.** A 256-dimension Gaussian random projection of k-mer
counts reached 96.2% on the locked test, above all four learned models, and frozen embeddings
degraded no more slowly than k-mer counting when exact sequence matches were removed
(−0.30 pp, 95% CI [−2.92, +2.22], p = 0.837; better in 3 of 14 blocks). That argument is about
whether the learned space does representational work. It does not depend on which endpoint is
scored, and it supports a **research-investment** decision.

**Not licensed.** Any claim about U. U contains amplification and read-feasibility terms that
gLM set to zero by construction.

**H4 is the reason, not a side note.** Rev 1 listed H4 — amplification constraints can change
the preferred assay — as "untouched". H4 is precisely the mechanism by which end-to-end U
diverges from conditional accuracy. If H4 holds, a method's conditional discrimination need
not determine its U-ranking. H4 being open is *why* the gLM result cannot be transported to
the primary endpoint.

**One asymmetry worth preserving.** The ambiguity ceiling is a property of the reference
library, not of a classifier or an endpoint. If a target and a neighbour share an identical
observed amplicon, no sequence-only method separates them, whatever the amplification does.
So **collision ⟹ unresolvable end-to-end** holds one-way and endpoint-independently, and an
abstention rule built on collisions is valid now. The converse — no collision ⟹ resolvable —
requires the amplification term and must wait for H4.

**Throughout this memo**, any statement sourced to gLM is qualified as *conditional-endpoint,
sequence-only, sub-assay-length*.

---

## 2. What v1.0 already contains

| Proposal | Already in v1.0 as | Operation |
|---|---|---|
| Exact-match leakage diagnostic | **B2**, a mandatory baseline (§10.1) | print adjacent to every headline accuracy, not buried in a baselines table |
| Report regret with hit rate | §11.4 top-k hit rate with tie handling | pair them in the same display as §11.2 regret |
| Overlap reporting | §9.3 already *requires* overlap counts and shared hashes | promote to blocking (Gate 6) |
| Resolvability / abstention | §3 Q1 and Q5; §5 *Unresolvable*; §11.3; §15 flag set | promote from secondary to primary |
| Cheap non-learned control | §10.1's anti-straw-man rule and B4 | add the random projection explicitly |
| Haplotype-aware curation | §7.2 eligibility, stated in isolate groups | change the unit (see §4) |
| Pilot to estimate block-level variability | §14, already prescribed | execute it, at assay length |
| Tie tolerance in ranking | §10.4, 0.01 accuracy units | vindicated; keep |

v1.0 anticipated almost all of this. The genuinely new content is two items: an external
estimate of block-level SD, and the proposal to change the learned arm's status.

---

## 2A. The empirical case for the reframing, and one correction to it

A variance decomposition over the 20 locked-test blocks × 6 methods, on deposited results:

| Source | Share of total variance in balanced accuracy |
|---|---|
| Between **block** | **98.41%** |
| Between **method** | **0.30%** |
| Residual | 1.28% |

**Block identity explains 324 times what method identity explains.** Whatever determines the
answer, it is a property of the panel, not of the comparison procedure. This is the strongest
single argument for demoting the method contrast, and it is stronger than the ≈50× the
appraisal estimated.

**But the quantity is the ceiling, not the collision indicator.** The appraisal proposed
regressing accuracy on the fraction of target genomes sharing an observed sequence with a
neighbour. Both predictors were tested against block-mean accuracy:

| Predictor | R² | p |
|---|---|---|
| Collision fraction (share a sequence with any neighbour) | **0.234** | 0.030 |
| Empirical ambiguity ceiling at the evaluated region | **0.964** | 2 × 10⁻¹⁴ |

The collision indicator is weak because a target sequence colliding with *some* neighbour
sequence does not make the *class* unresolvable. `proteus_morganellaceae` has a collision
fraction of 0.967 and 98.1% accuracy; `rhizobium_agrobacterium` 0.914 and 98.0%; `rothia`
0.875 and 99.3%. The ceiling weights how many isolates sit on each side of each shared
haplotype; the indicator does not. **§11.3 already defines the right quantity.**

**The mechanism.** Across all 33 blocks, accuracy sits a near-constant distance below the
fold-matched ceiling — **2.11 pp (SD 2.48)** — and that distance is method-independent:
0.92 pp for k-mer counting against 1.13 pp for DNABERT-2 on the locked test. The shortfall is
the same order as the entire between-method spread (0.78 pp). Accuracy is approximately
*ceiling minus a small constant*, which is why the ceiling predicts it and the method does
not.

One caveat to state whenever this is reported: the ceiling bounds accuracy by construction, so
a high R² is partly mechanical. The non-trivial content is the **uniformity** of the
shortfall. That is what converts a bound into a usable predictor.

**This clears the appraisal's own kill criterion in advance.** Its criterion 1 was that the
project fails if a calibrated resolution probability is no better than the raw exact-collision
indicator. On deposited data the raw indicator explains 23% of between-block variance and the
properly weighted quantity explains 96%. The gap between naive and calibrated is already
demonstrated, before any block is curated.

---

## 2B. The amplicon-span check, run

The cached full-length locked-test records (6,761 across 20 blocks) were re-sliced at
primer-trimmed amplicon inserts alongside the strict regions, and haplotypes, collisions,
ceilings and k-mer accuracy recomputed per span
(`scripts/step1_amplicon_span_check.py` → `results/step1_amplicon_spans.json`).

| Span | nt | Haplotypes | Collision fraction | Ceiling | k-mer | Headroom |
|---|---|---|---|---|---|---|
| V4 strict | 107 | 32 | 0.284 | 95.53 | 94.22 | 1.31 |
| **A_V4** 515F–806R | 253 | 51 | 0.246 | 96.04 | 94.31 | 1.73 |
| **A_V1V2** 27F–338R | 294 | 171 | **0.067** | **98.55** | 95.76 | 2.78 |
| **A_V6V8** 939F–1378R | 401 | 69 | 0.207 | 97.60 | 96.06 | 1.55 |
| **A_V3V4** 341F–806R | 429 | 74 | 0.218 | 96.56 | 94.53 | 2.03 |
| **A_V1V3** 27F–534R | 488 | 189 | **0.058** | **98.67** | 95.81 | 2.86 |

Aggregate: strict regions 27 haplotypes, collision 0.449, ceiling 90.73; amplicons 111
haplotypes, collision 0.159, ceiling 97.49.

**Finding 1 — region content beats length, and the field's default assays are the weaker
choice.** V1–V2 at 294 nt resolves better than V3–V4 at 429 nt (collision 0.067 against
0.218; ceiling 98.55 against 96.56). The two most widely used amplicons, V4 and V3–V4, are the
*worst* of the five tested for target-versus-neighbour discrimination. This is a direct,
practically useful result of exactly the kind regioneer exists to produce, and it argues
against any framing in which "longer amplicon" is the recommendation.

**Finding 2 — length partly splits shared haplotypes, but a hard core persists.** Of the 8
blocks with >20% collisions at V4 strict, moving to V3–V4 resolved 2 (`rothia` 0.25→0.00,
`stenotrophomonas` 0.31→0.01), improved 1 (`veillonella` 0.63→0.10) and left 5 unchanged
(`brucella` 0.77→0.75, `francisella` 0.96→0.96, `methylobacterium` 0.85→0.81, `proteus`
0.93→0.89, `rhizobium` 0.74→0.70). For protocol 2 this means the addressable set grows but
its hard core does not.

**Finding 3 — headroom GROWS with length.** The gap between the ceiling and k-mer accuracy is
**1.58 pp across strict regions and 2.19 pp across amplicons**. Rev 1's near-lossless argument
was asserted by analogy to 31–107 nt spans; at 253–488 nt it weakens, and there is measurably
more room for a better method. Whether anything *claims* that room is a separate question,
answered in §2C.

---

## 2C. Step 1c — the equivalence transfers to amplicon length

§2B left one question open: more headroom exists at amplicon length, but ceilings bound what
*any* method can do and say nothing about whether methods differ. Step 1c measures it
directly, like for like with gLM — same probe, same isolate-grouped folds, same pooling,
conditional endpoint, 12,566 unique amplicon slices embedded with DNABERT-2
(`scripts/step1c_amplicon_methods.py` → `results/step1c_amplicon_methods.json`).

| Span | nt | k-mer | DNABERT-2 | Δ (pp) | 95% CI | Verdict at ±1 pp |
|---|---|---|---|---|---|---|
| V2 strict | 106 | 95.17 | 94.84 | −0.33 | [−0.93, +0.15] | equivalent |
| V4 strict | 107 | 94.22 | 94.21 | −0.00 | [−0.43, +0.43] | equivalent |
| A_V4 | 253 | 94.31 | 94.25 | −0.06 | [−0.38, +0.23] | equivalent |
| A_V1V2 | 294 | 95.76 | 95.57 | −0.19 | [−0.48, +0.08] | equivalent |
| A_V6V8 | 401 | 96.06 | 95.75 | −0.31 | [−0.62, −0.02] | equivalent |
| A_V3V4 | 429 | 94.53 | 94.18 | −0.35 | [−1.03, +0.14] | inconclusive |
| A_V1V3 | 488 | 95.81 | 96.04 | **+0.23** | [−0.17, +0.69] | equivalent |

**The equivalence transfers.** Six of seven spans are equivalent at the ±1 pp margin and the
seventh is inconclusive, not favourable. No span shows a DNABERT-2 advantage whose interval
excludes zero. The largest positive difference anywhere is +0.23 pp at the longest amplicon,
with an interval spanning zero.

**And the extra headroom goes unclaimed by either arm.** At A_V1V3 the ceiling is 98.67 and
the arms reach 95.81 and 96.04 — DNABERT-2 closes 0.23 pp of a 2.86 pp gap. At A_V1V2 it is
*worse* than k-mer. So §2B's larger headroom at amplicon length is not room a frozen encoder
exploits; it is room neither composition counting nor a frozen representation reaches. What
occupies it is a real open question — plausibly within-species and intragenomic variation
rather than anything a better sequence representation would recover — and it is not an
argument for the learned arm.

**Consequence.** §9's demotion of the learned arm is now supported **at regioneer's own unit**
rather than extrapolated from 31–107 nt spans. The pure-classical package in §18 follows
without the length caveat that rev 4 had to carry.

**Scope, unchanged.** This remains the conditional endpoint — no primers, no amplification, no
read layout. It speaks to §11.1's *secondary* endpoint. H4 is still why nothing here
transports to U, and the assay-length pilot with in-silico PCR is still the step that changes
that.

---

## 3. §11.3 — the ceiling must be fold-matched, with a demonstration

v1.0's A_ceiling is a whole-population quantity over weights w(s,y). U in §11.1 is a nested
mean evaluated per fold. These are not comparable, and the pooled quantity can sit **below** a
legitimately attainable accuracy. Minimal counterexample, verified:

- fold 1: target `a`, neighbour `b` → fold ceiling **1.000**
- fold 2: target `c`, neighbour `a` → fold ceiling **1.000**
- pooled over both folds → A_ceiling **0.750**

Pooling manufactures a collision between `a`-as-target and `a`-as-neighbour that exists in no
evaluation fold. This is the mechanism behind the 3 of 66 block-by-method cells that exceeded
a whole-block ceiling in gLM, while none exceeded the fold-matched version.

**Edit:** compute A_ceiling per evaluation fold using that fold's weights, then aggregate by
the identical nested-mean scheme used for U. State that the whole-population ceiling is a
separate descriptive quantity about the reference library and bounds nothing. Two lines, no
cost, and it protects the endpoint this memo proposes to make primary.

---

## 4. §7.2 — eligibility is specified in the wrong unit, and conflicts with §7.1

§7.2 asks for "at least 5 distinct genome/isolate groups per species, preferably 10 or more."
Under heavy amplicon sharing, isolate-group count and distinct-haplotype count diverge, and it
is the second that determines whether a block can support a sequence-disjoint split. In gLM,
only **14 of 20** blocks supported exact-disjoint folds and **9 of 20** at 99% identity; the
rest had ≤4 distinct target haplotypes.

**Edit:** add a second eligibility criterion in §7.2 — *distinct observed-amplicon haplotypes
per species, at the evaluated assay*. This is **assay-dependent**: a block may be eligible at
V1–V3 and ineligible at V4, which §7.2's current single-number form cannot express. Rev 1 put
this in §14; it belongs in eligibility.

**The conflict with §7.1, and its resolution.** §7.1 requires "easy, intermediate, and
intrinsically ambiguous panels" and forbids excluding hard panels because methods perform
poorly. A haplotype-diversity rule excludes exactly the intrinsically ambiguous blocks — which
are the evidence the resolvability endpoint needs. Declare **two frames**:

- **Diversity-eligible frame** — blocks meeting the haplotype criterion; carries the
  method-contrast endpoint.
- **Full frame** — every curated block including the ineligible ones; carries resolvability
  and abstention.

Low-haplotype blocks are not curation failures. They are results, and they are the ones the
tool exists to flag.

---

## 5. §9 — promote the sequence-disjoint regime to co-primary

§9.3 currently frames sequence novelty as an *addition* to a realistic primary evaluation and
warns about its induced selection bias. That caution is well placed, which is why the
proposal is co-primary status rather than replacement.

gLM's grounds: equivalence held under isolate grouping (−0.21 pp) and became **inconclusive**
once exact sequences were separated (−0.50 pp, 95% CI [−3.20, +2.07]), with both arms losing
≈11 points and the block-level SD rising from 0.37 to 5.16 pp. A primary analysis in which
92% of evaluated units share a sequence with training cannot distinguish generalisation from
lookup. *(conditional-endpoint, sequence-only, sub-assay-length.)*

**Edit:** report **both regimes, pre-declared and co-primary** — isolate-grouped carrying the
operational claim, sequence-disjoint carrying the generalisation claim. Use single-linkage
clustering on the identity graph, which is the construction that guarantees every cross-fold
pair falls below threshold; greedy and centroid methods do not. **Audit the realised split**:
report maximum cross-fold identity, and verify that B2 scores exactly 0.500 under a correct
disjoint split, which is a necessary condition.

**Gate 6 (new):** no generalisation claim from an evaluation whose train/test exact-sequence
overlap has not been measured and reported. This promotes an existing §9.3 reporting
requirement to blocking.

---

## 6. §11 — resolvability primary, ranking secondary

gLM's region-ranking results: regret 0.65 pp (k-mer) against 0.57 pp (DNABERT-2), paired
difference −0.09 pp [−0.46, +0.26], with best-region hit rates of **37% and 44%**. Regions tie
often enough that a selector misses the argmax more than half the time at almost no cost.
§10.4's prespecified 0.01 tie tolerance is vindicated; the ranking as a *product* is not.

Separately, difficulty is near-binary on shared sequence rather than graded with divergence:
a 13-point swing between strata, with the method gap ≤0.25 pp in every stratum.

**Edit:**
- **Primary:** per-assay resolvability, carried by the **fold-matched ambiguity ceiling**
  (§2A: R² = 0.96 against block accuracy, where a raw collision indicator gives 0.23), with
  the target–neighbour collision table as a supporting display and explicit abstention when
  no candidate clears a prespecified criterion. The open research problem is estimating that
  ceiling under **reference incompleteness** — a target with 3 sequenced genomes and a
  neighbour with 500 gives a badly biased sample of true haplotype diversity — which is what
  converts a database lookup into a calibrated statement.
- **Secondary:** U(p, â) and regret, retained as now, with regret and hit rate always
  displayed together and a note that low regret with a low hit rate means the choice rarely
  matters, not that the selector is accurate.

The abstention half of this is licensed now (§1, one-way asymmetry). The positive
resolvability claim awaits H4.

---

## 7. §14 — recomputed against both candidate estimands

v1.0 §14 already says no fixed sample size is defensible before surveying blocks and
prescribes a development-only pilot to estimate block-level paired variability. Rev 1
performed that procedure with external data; it was **executing §14, not correcting it**.

What needs changing is the illustrative grid, and the answer now depends on which estimand
§11 makes primary.

### 7.1 If the method contrast stays primary

Using gLM's block-level SD of 5.16 pp under sequence-disjoint folds, t-based:

| 95% CI half-width | Isolate-grouped (SD 0.37) | Sequence-disjoint (SD 5.16) |
|---|---|---|
| ±1.00 pp | 1 | **105** |
| ±0.50 pp | 2 | **415** |

The SD is estimated from 14 blocks, so its own 95% interval is [3.74, 8.31] pp; since n scales
as SD², the requirement spans **≈55 to ≈265 blocks**. Those 14 blocks are also a *selected*
subsample — the more haplotype-diverse ones, i.e. the frame §4 proposes to curate for — so the
direction of bias is not obvious. Plan against the upper end, or declare the arm powered for
±2 pp.

### 7.2 If resolvability becomes primary — and it does not dissolve the problem

The appraisal argued that moving to a calibration deliverable makes blocks test cases rather
than replicates, and that 20–40 blocks would suffice. **That is too optimistic, by a factor of
roughly 1.5 to 3**, because predictions within a block are not independent.

Each block contributes one prediction per candidate assay, so cells accumulate faster than
blocks — but only to the extent that those predictions are independent, and they are not.
**ICC has now been measured at both spacings** (§2B data):

| Spacing | m | ρ | Calibration slope | Three-bin curve |
|---|---|---|---|---|
| gLM's 9 strict regions | 8.2 | 0.496 | 56 | 112 |
| **regioneer's 5 candidate amplicons** | 5.0 | **0.649** | **72** | **144** |

**ρ is worse at assay spacing, not better**, and the reason is structural: candidate amplicons
**nest** — V1–V2 ⊂ V1–V3, V4 ⊂ V3–V4 — so they are more correlated than adjacent distinct
regions. An appraisal suggested ρ ≈ 0.2 might bring this inside §7.1's existing 30–60 target.
Measured, it does not: ρ = 0.649.

Nor does adding candidate assays rescue it. At the measured ρ, going from 4 to 12 assays per
panel moves the slope requirement only from 63 to 54 blocks, because each extra nested assay
contributes little independent information.

So: the reframing still buys a cheaper **minimum** — **72 blocks against the method contrast's
105** — but a proper calibration curve costs **more** (144). §7.1's 30–60 target is adequate
for neither.

Three caveats, all of which argue for the conservative end: the ICC is itself estimated from
20 blocks; "resolvable" here is a ceiling-threshold proxy rather than an assay-level outcome;
and the effective-cell targets are conventional rules of thumb from prediction-model
validation, not derived for this design.

### 7.3 What follows

Whichever estimand is chosen, **every number here exceeds §7.1's 30–60 planning target**, and
the honest options are: raise the curation target toward ~60 and accept the minimal
calibration deliverable; or declare the study a calibrated case series under §7.1's own
instruction and drop the population claim. Decide before curation scope is fixed, not after.

---

## 8. §16 — a de-scoping ladder, because no gate catches scope

v1.0 commits to a four-way taxonomy crosswalk, dual extraction with boundary auditing, 30–60
curated blocks, validated in-silico PCR over a degenerate primer table, seven mandatory
baselines including a Qscore replication, two to four learned encoders with full provenance, a
multi-factor simulation grid, DADA2 read-processing validation, three empirical tiers ending
in wet-lab work, a CLI with twelve unit-test classes, and ten figure classes — on a 13–16 week
schedule. Gates 1–5 check **validity, not scope**. Nothing forces "this is too big, cut X."

This memo adds pressure: a co-primary regime, a larger block target, and — if both regimes and
both endpoints go co-primary — a compounding design and multiplicity burden.

**Edit — declare the ladder now, with date triggers, not judgement calls:**

- **Tier 1, defensible core.** Blocks under corrected §7.2; Track R regions *and* Track A
  assays for a small primer set (V4, V3–V4, V1–V3, full-length); B0–B5; both evaluation
  regimes; U, resolvability, regret and hit rate; error-free plus one error profile. This is a
  complete paper on its own.
- **Tier 2, if time allows.** B6 Qscore; one frozen encoder as reference arm; wider simulation
  grid; empirical Tiers 1–2.
- **Tier 3, defer explicitly.** Prospective wet-lab validation; Shiny interface; additional
  encoders; full factorial stress grid.

---

## 9. §10.2 — demote the learned arm, do not delete it

**Rev 1 proposed dropping it. That was wrong and is withdrawn.** The contribution this memo
claims — a benchmark establishing when 16S cannot resolve, plus evidence against reaching for
a language model — requires evidence in *this study's* endpoint. Citing gLM supplies evidence
at a different endpoint and sub-assay length.

**Edit:** retain **one** frozen encoder as a reference arm under the same splits. It converts
a citation into a result for roughly one line in a table. Drop the rest, and compress §10.2's
checkpoint/pooling/strand apparatus accordingly. Add the **random-projection control**, which
is the cheapest and most informative single comparator gLM identified and bounds any claim
that a learned space is doing representational work.

§18's dependency argument is unaffected: ship the package pure-classical and keep the encoder
arm in the analysis repository. §15 already separates the CLI from optional components.

**H3 disposition:** restate as *deprioritised on cost-effectiveness grounds for frozen
representations at sub-assay length*, resting on the random-projection control, **not** as
answered. Adaptation remains live under protocol 2.

---

## 10. Execution plan

Revised after external review of an earlier ordering, which correctly identified that the
estimand decision had been made to gate on a pilot that was itself gated on it.

**0a — choose the primary estimand. Today, no new data needed.** A strategic choice, not a
statistical one. Separating it from 0b removes the circularity.

**Step 1 — DONE (§2B).** Amplicon re-slice, haplotypes, collisions, ceilings, k-mer accuracy,
and ρ at assay spacing. Outputs deposited.

**0b — set the block target. Settled by Step 1, not by the pilot.** At the measured ρ = 0.649:
72 blocks for a calibration slope, 144 for a curve, against 105 for the method contrast. All
exceed §7.1's 30–60.

**Step 1c — DONE (§2C).** Equivalence transfers: 6 of 7 spans equivalent, none favouring the
encoder, and the extra headroom claimed by neither arm.

**Step 1b — in parallel, off the critical path.** Secure Tier 2 public mock-community data
with usable primer metadata. This has been carried as an aside through three revisions, which
is how such items die. It is the difference between a tool and a proposal.

**Step 2 — protocol edits needing no new data.** Decisions 1, 2, 3, 6 below, plus three items
raised in review and not previously captured:

- **Freeze §15's abstention thresholds, with a stated basis, before the locked set is
  touched.** Once resolvability is primary these *are* the primary decision rule. This is the
  item most likely to be fudged under pressure, because a small adjustment visibly improves
  the headline.
- **State the 2 × 2 explicitly.** Evaluation *regime* (isolate-grouped, sequence-disjoint) and
  *endpoint* (resolvability, method contrast) are orthogonal. Name which cells are
  confirmatory and which descriptive; "both regimes co-primary, resolvability primary" is not
  well formed.
- **Write the de-scoping ladder into §16** before curation starts.

**Step 3 — assay-length pilot (§7.1's 2–4 groups).** Its job is now to exercise the in-silico
PCR path that H4 depends on and to confirm the block target, not to establish it.

**Step 4 — curation**, two frames, with an **adequacy gate at 4 → 5**: fewer than ~20 blocks
having both adequate haplotype diversity and adequate reference support → report case studies
under §7.1's own instruction and drop the population claim.

**Step 5 — freeze and locked evaluation**, single use, against the 2 × 2 declared in Step 2.

**Step 6 — release.** Documentation, fixtures and packaging only: **build Steps 3–5 inside the
package from the start**, or the analysis gets written as scripts and rewritten later.

---

## 11. Decisions required

| # | Decision | Can be taken |
|---|---|---|
| 1 | Sequence-disjoint regime co-primary (§5) | **now** |
| 2 | §11.3 fold-matching and §7.2 unit/two-frame edits (§3, §4) | **now** |
| 3 | De-scoping ladder into §16 (§8) | **now**, before curation |
| 4 | Resolvability primary, ranking secondary (§6) | abstention half now; positive half after step 2 |
| 5 | Block-count target and which endpoint it is powered for (§7) | after step 2 supplies an assay-length SD; note every option exceeds §7.1's 30–60 |
| 6 | Learned arm reduced to one reference encoder (§9) | **now** |

Decisions 1, 2, 3 and 6 are independent of the open scientific questions. Decisions 4 and 5
should wait for the assay-length pilot, which is the step that converts gLM's
conditional-endpoint finding into something that speaks to U.
