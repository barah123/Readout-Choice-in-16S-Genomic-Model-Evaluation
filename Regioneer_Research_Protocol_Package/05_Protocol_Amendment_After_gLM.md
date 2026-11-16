# Protocol amendment 1 (rev 2): re-ranking Protocol 1 after the gLM pilot

**Status:** proposed amendment to `01_Regioneer_Benchmark_and_Tool_Protocol.md` v1.0.
**Date:** 2026-11-16. **Supersedes rev 1** of this file.

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
- **Primary:** per-assay resolvability — fold-matched ceiling, target–neighbour collision
  table, and explicit abstention when no candidate clears a prespecified criterion.
- **Secondary:** U(p, â) and regret, retained as now, with regret and hit rate always
  displayed together and a note that low regret with a low hit rate means the choice rarely
  matters, not that the selector is accurate.

The abstention half of this is licensed now (§1, one-way asymmetry). The positive
resolvability claim awaits H4.

---

## 7. §14 — execute the pilot v1.0 already prescribes, at assay length

v1.0 §14 already says no fixed sample size is defensible before surveying blocks and
prescribes a development-only pilot to estimate block-level paired variability. Rev 1
performed that procedure with external data; it was **executing §14, not correcting it**.

What needs changing is the illustrative grid. Using gLM's SD of 5.16 pp, t-based:

| 95% CI half-width | Isolate-grouped (SD 0.37) | Sequence-disjoint (SD 5.16) |
|---|---|---|
| ±1.00 pp | 1 | **105** |
| ±0.50 pp | 2 | **415** |

**This estimate is soft and must be presented as such.** SD is estimated from 14 blocks, so
its own 95% interval is [3.74, 8.31] pp, and since n scales as SD² the requirement spans
**≈55 to ≈265 blocks**. Those 14 blocks are also a *selected* subsample — the more
haplotype-diverse ones, i.e. the frame §4 proposes to curate for — so the direction of bias is
not obvious. Plan against the upper end, or declare the arm powered for ±2 pp.

**And the grid may answer the wrong estimand.** §14's n is for a between-method contrast. If
§6 is adopted and resolvability becomes primary, the planning quantity is a block-level
proportion: **≈96 blocks for ±10 pp, ≈384 for ±5 pp** at worst case. Recompute for whichever
endpoint ends up primary.

Against §7.1's target of 30–60 blocks, every one of these numbers is larger. That is the
substantive planning consequence and it should be settled before curation scope is fixed.

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

## 10. Next steps, in order

1. **Two-hour check, before anything else.** Count distinct observed-amplicon haplotypes per
   species, and exact k-mer-profile collisions at k = 4, 6, 8, as a function of span length,
   across blocks already in hand. This tests directly whether the classical baseline is near
   the information ceiling at amplicon length — which rev 1 asserted by analogy to gLM's
   31–107 nt spans rather than measuring. It also calibrates the §7.2 threshold empirically.
2. **Run the §7.1 pilot (2–4 groups) at assay length**, with Track A amplicons through §8.2's
   in-silico PCR. One pilot tests whether the gLM equivalence transports to 253–460 nt,
   supplies the §14 SD in the right units on the right endpoint, and exercises the
   feasibility code H4 depends on.
3. **Make the §11.3 fold-matching edit and the §7.2 unit change.** Both small, both
   unambiguous, neither depends on the H3 debate.
4. **Write the §8 ladder into §16** before curation starts.

---

## 11. Decisions required

| # | Decision | Can be taken |
|---|---|---|
| 1 | Sequence-disjoint regime co-primary (§5) | **now** |
| 2 | §11.3 fold-matching and §7.2 unit/two-frame edits (§3, §4) | **now** |
| 3 | De-scoping ladder into §16 (§8) | **now**, before curation |
| 4 | Resolvability primary, ranking secondary (§6) | abstention half now; positive half after step 2 |
| 5 | Block-count target and which endpoint it is powered for (§7) | after step 2 supplies an assay-length SD |
| 6 | Learned arm reduced to one reference encoder (§9) | **now** |

Decisions 1, 2, 3 and 6 are independent of the open scientific questions. Decisions 4 and 5
should wait for the assay-length pilot, which is the step that converts gLM's
conditional-endpoint finding into something that speaks to U.
