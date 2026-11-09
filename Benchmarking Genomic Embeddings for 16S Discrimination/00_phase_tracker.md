# Protocol 3 — phase tracker

**Protocol:** `05_Third_Idea_Critical_Review_and_Benchmark_Protocol.md` (this directory)
**Scope:** idea 3 only (frozen-representation benchmark). Ideas 1 (`regioneer`) and 2
(adaptation) are downstream and tracked in `../Regioneer_Research_Protocol_Package/` —
not started.
**Git baseline:** tag `pre-audit-baseline` (commit `57f427f`) — everything as it stood
before this protocol's Phase 1 audit. Per protocol §5.1, pre-fix results from
`glm_prototype_bundle (1)/` are preserved there as engineering diagnostics, not relied
on as empirical findings going forward.

This file tracks the six phases in protocol §18. Update at the end of each phase.

| Phase | Work | Status |
|---|---|---|
| 1 | Audit: archive code, fix ceiling/tie/weight issues, synthetic fixtures, current-panel summaries | **Done** — see `phase1_audit/` and summary below |
| 2 | Pilot: DNABERT-2, HyenaDNA, NTv2, DNABERT-S if runnable; classical controls; 3–5 dev blocks; primary + linear-probe readouts | **Done** — `phase2_pilot/phase2_pilot_report.md` |
| 2-full | Scaled to the §6.2 block target: 31 declared blocks (26 eligible), 4 models, full readout hierarchy | **Done** — `phase2_pilot/phase2_full_report.md` |
| 3 | Freeze: final model/checkpoint list, panel eligibility, primary metric, readout rule, multiplicity plan, margins, compute budget | **Done** — `phase3_freeze/FREEZE.md` (binding) |
| 4 | Scale: locked model × panel × region matrix on audited data | **Done** — `phase4_locked_test/`, 20/22 blocks eligible and scored |
| 5 | Evaluate: freeze predictions, unseal labels, prespecified summaries | **Done** — `phase5_evaluate/phase5_results.md` |
| 6 | Diagnose/report: exploratory analyses, reproducible release, conditional conclusions | **Done** — `phase6_report/FINAL_REPORT.md`, `README.md` |

## Phase 1 deliverables (done)

- [x] `phase1_audit/audit_report.md` — systematic walkthrough of protocol §5.2's checklist
  against both existing panels. **Two previously-undocumented gaps found:** (1) genome/
  isolate/copy identity is not tracked anywhere in the existing fetch pipeline — confirmed
  concretely by parsing raw NCBI headers: at least 1 (anthracis/cereus) and 7 (staph/CoNS)
  strain identifiers appear on multiple separate records, invisible to the current
  exact-string dedup; this is not retroactively fixable and is now a Phase 2/4 data-collection
  requirement, not a Phase 1 patch. (2) `kmer_matrix()` fits its vocabulary on pooled
  train+test sequence content rather than train-only — a transductive (not label) leak, to be
  rebuilt fit-per-fold when Phase 2's classical control is implemented.
- [x] `phase1_audit/synthetic_fixture.py` — 6-part known-ground-truth suite (alignment/
  extraction under a known indel, oracle-ceiling formula on 3 hand-computed cases, LOHO-CV
  recovery, permutation-null collapse, tie-aware NN, identity self-check). All pass. One
  instructive failure along the way, kept as a documented finding rather than silently
  patched: an initial haplotype design shared k-mer-level support across the class boundary
  despite looking zero-overlap at the string level, and with only 1 haplotype/class/fold
  remaining, k-mer+RF inverted every prediction (0.0 balanced accuracy, reproducible across
  5 seeds) — independent mechanistic confirmation of why ≥5 haplotypes/class is the
  preferred minimum, not just ≥2.
- [x] `phase1_audit/env_smoke_tests.md` — re-ran real (not import-only) inference smoke
  tests for HyenaDNA and NTv2 rather than inheriting the prior session's claim. Both confirmed:
  deterministic, finite, correct output shapes across 29/107/1500nt. **New:** got DNABERT-S
  running too (flagged as essential prior art by the critical review, not previously in either
  env) — same MosaicBERT codebase as DNABERT-2, the same 4 documented patches transferred
  directly, loaded and ran cleanly on the first attempt. All four Phase 2 core models
  (DNABERT-2, DNABERT-S, HyenaDNA, NTv2) are now confirmed runnable.

## Phase 2 headline (full detail in `phase2_pilot/phase2_pilot_report.md`)

**The Phase 3a result this project treated as its cleanest finding does not survive.** Phase 3a
reported DNABERT-2 at staph V1 as "exactly chance, 0.500" against k-mer's 0.997. Re-run with a
linear probe (R2) instead of RF-on-raw-embeddings, on 57 target isolates instead of 19:
DNABERT-2 = **0.964**, HyenaDNA = **0.979**, against a ceiling of 0.982. The representation did
contain the information; the prior readout could not reach it. This is exactly the error the
protocol's claim ladder (§4) and decision table (§22) exist to prevent, and it is a correction
to our own prior claim rather than a new result about any model.

Three further pilot findings:
- **R0 (unsupervised geometry) misranks representations** and must be diagnostic-only: DNABERT-S
  scores 0.542 at staph V1 by geometry but 0.945 by probe; several cereus regions score R0 ≈ 0.8
  where the ceiling is only ≈ 0.55, i.e. the geometry is separating something that is not the
  class label.
- **Frozen models vs. classical is close to a tie** under a matched probe — mean Δ of +0.007
  (cereus) and +0.002 (staph), no region exceeding +0.025. And `C_random` (random projection of
  the k-mer space) nearly matches `C_kmer`, so small deltas between high-dimensional feature
  sets should not be over-read.
- **The cereus panel is not uniformly null.** At 92–114 target isolates rather than 9, V6 (0.80)
  and V7 (0.70) support real discrimination, while V3/V4/V5/V8 remain near their low ceilings.
  V6 is where Chakravorty et al.'s original anthracis/cereus claim sat — at this sample size it
  partially replicates, where `prototype_report.md` had recorded it as failing to.

## Phase 2-full headline (detail in `phase2_pilot/phase2_full_report.md`)

Scaled from 2 panels to **31 pre-declared taxonomic blocks** (committed before fetching,
commit `1918e6f`), 154 species, 12,356 records, 277 scored cells. **26 blocks pass the
pre-declared eligibility gate**; the 5 failures are thin-reference organisms (*Y. pestis*,
*B. pertussis*, Borrelia, Leptospira, Legionella neighbours) retained as descriptive-only.

**Primary contrast, 25 independent species-rank blocks:** best linear probe on frozen
embeddings minus best classical control = **+0.0014, 95% bootstrap CI [−0.0032, +0.0057]**;
family-grouped stringent version −0.0001. 18 of 25 blocks are ties within ±0.01. Under §15's
vocabulary this is *evidence of practically negligible improvement* — a narrow interval
excluding gains above ~0.006 — **not** "inconclusive" and **not** evidence that the
representations lack information. The flexibility reference matters: the models beat a random
projection of the k-mer space by only +0.0036.

**The dominant effect is the readout, not the model.** Spread across R0→R1→R2 within a model
is ~0.21–0.27; spread across the four models within R2 is 0.007 — about thirty-fold. Any
benchmark reporting a single pooled-embedding-distance number measures its readout far more
than its representation, which is exactly how this project's own Phase 3a went wrong.
DNABERT-S, the one model contrastively trained for species-aware embeddings, does **not**
beat DNABERT-2 here (0.942 vs 0.946) — its 10 kbp training does not automatically transfer to
short conserved 16S windows, answering an open question from the protocol's §2.2.

**Our a-priori ambiguity predictions were systematically pessimistic.** Expected-class
ordering held (A 0.994 > B 0.975 > C 0.912), but textbook "16S-unresolvable" groups —
*E. coli*/*Shigella* 0.931, *Enterobacter* 0.934, *Acinetobacter* 0.990 — are substantially
separable given ~100 isolates per class and the right region. Only campylobacter (0.782)
behaves as the literature predicts.

## Phase 3 freeze (binding) — `phase3_freeze/FREEZE.md`

**The decision that mattered most: all 31 Phase 2 blocks are now development data and are
burned for confirmatory use.** Readout choice, comparator set, pre-named arms, margin and
precision plan were all derived by looking at them, so re-using them as a test set would be
the "repeated test-guided model selection" failure protocol §18 names. Phase 2's results keep
their full value as *development* findings — the frozen≈classical headline is now a
hypothesis to be tested, not a confirmed result.

Frozen: four checkpoints pinned by sha256 / hub revision; **R2 linear probe** as primary
readout (development showed readout dominates model ~30-fold); primary estimand
**DNABERT-2:R2 − C_kmer**, both arms pre-named from development (20/25 and 17/25) so no
selection sits inside the statistic; isolate-grouped 5-fold CV; region chosen support-only
and applied identically to both arms; **equivalence margin ±0.010** via bootstrap TOST;
confirmatory family of exactly 5 tests (1 primary + 4 Holm-corrected), everything else
explicitly exploratory; ≤3 h CPU budget; single-use unsealing checklist.

Precision plan from development SD (0.01138): ≥14 eligible blocks gives an expected CI
half-width ≈0.006, inside the ±0.010 margin. **22 locked-test blocks declared** across 22
families with **zero family overlap** with development (94 species queries), assuming ~35%
eligibility loss. If fewer than 10 prove eligible the primary analysis will be reported as
underpowered, not spun.

Three deviations from the protocol document are logged with justification in
`phase3_freeze/DEVIATIONS.md` (margin tightened from the illustrative ±0.03; single pre-named
arms; the development/test re-split).

## Phase 5 — CONFIRMATORY RESULT (locked test, run once)

Unsealing checklist verified by git, not asserted: `FREEZE.md` committed 07:20:42 and
unmodified since; first locked-test record fetched 08:54:19, **94 minutes later**.

**Primary endpoint — DNABERT-2:R2 − C_kmer over 20 independent blocks:
mean −0.00211, 95% bootstrap CI [−0.00595, +0.00163].
Verdict at the frozen ±0.010 margin: EQUIVALENT (practically negligible difference).**
17 of 20 blocks tie within ±0.01. Achieved CI half-width 0.0038 vs 0.005 planned.

This **replicates the development finding on entirely independent taxa** — 20 families
disjoint from all 31 development blocks — with the opposite sign and a comparably narrow
interval. Two block sets, 45 blocks total, same conclusion.

Secondary (Holm-corrected): DNABERT-S and NTv2 also equivalent; **HyenaDNA is reliably
slightly worse than k-mer** (−0.0088, Holm p=0.004) though below the practical margin;
best-of-4 does not survive correction (p=0.075).

**The most uncomfortable exploratory number:** `C_random`, a random projection of the k-mer
space, tops the table at 0.9616 — nominally above every learned representation. Best-of-4
minus C_random is +0.0021 [−0.0005, +0.0052]. At this difficulty most measured performance
reflects classifier flexibility over a high-dimensional space, not representation content.

Explicitly **not** supported: that gLMs contain no 16S information (all reach 0.95–0.96
against a ~0.98 ceiling), that they are useless for 16S generally (one task, one readout
family, clean sequences), or that classical methods are superior (the finding is
equivalence, not superiority).

## Phase 6 — diagnostics and final report

**Diagnostic A (difficulty stratification) is the most important follow-up.** Half the locked
test blocks (10/20) were saturated — classical error under 1%, so no method could differ.
Restricted to the 5 genuinely hard blocks the difference is +0.0026 [−0.0046, +0.0120]:
centred near zero but wide. And the gLM advantage **increases with task difficulty**
(Spearman rho = +0.595, p = 0.006 for best-of-4). Exploratory, found after the fact, n=5 —
but it is the strongest lead for future work and the reason the report does not conclude
that these models are useless here. It also means the equivalence verdict is tight on easy
problems and genuinely uncertain on hard ones, which is now stated explicitly.

**Diagnostic B (pooling sensitivity)** on the same hard blocks: CLS token, max pooling and
mean+max concatenation all land within ±0.003 of k-mer. No pooling choice rescues the gLM arm.
Intermediate *layers* were deliberately not tested — MosaicBERT stores them unpadded
internally, and extracting them correctly would mean reimplementing part of the model, so the
check was abandoned rather than done wrongly. Recorded as an explicit open limitation.

**Release:** `README.md` (how to reproduce, ~3 h on a laptop) and
`phase6_report/FINAL_REPORT.md` (findings in plain language, no jargon).

## Phase 7 — two pre-registered extensions, both negative

**Extension A (hard panels): hypothesis NOT supported.** Phase 6 suggested the models might
gain ground as problems get harder (rho=+0.595, p=0.006, n=5). Tested properly on 13 panels
screened for difficulty using classical methods only (mean k-mer accuracy 0.818 vs 0.961 in
Phase 5; hardest 0.552): mean DNABERT-2 − k-mer = **−0.0019 [−0.0128, +0.0066]**. 4 panels
favour the model, 3 favour k-mer, 6 tie. Leave-one-out shows no hidden single-panel effect.
**The Phase 6 correlation is reported as not replicated**, with no regrouping to rescue it,
as pre-committed. Equivalence now holds across **58 independent groups** spanning easy to
genuinely hard.

One loose end: best-of-4 minus k-mer = +0.0070 [−0.0010, +0.0135]. Most likely upward bias
from taking the maximum of four noisy arms; distinguishing that from real model
complementarity needs a pre-registered test with the selection rule fixed in advance.

**Extension B (imperfect sequences): equivalence robust.** Six degradation conditions, both
arms given identical damaged sequences. No condition produced a difference both meaningful
and supported. The practical finding is unrelated to the models: **read truncation costs
~4 points, more than a tenfold increase in sequencing error (~2–3 points)** — primer
placement matters more than classifier choice.

Extension B also required a self-correction. Three conditions appeared to show the model
tolerating noise better, with a clean mechanism (one base change alters up to 15 k-mer
counts). It was driven **entirely by one small panel** (*francisella*, 63 records, +0.21
shift); excluding it moved the estimate from +0.0099 to −0.0006, and leave-one-out flipped
the sign. In the worst-case condition the mean pointed one way while 10 of 16 panels pointed
the other.

**Incidental finding, model-independent and arguably the most useful for practitioners:**
groups picked specifically for reputedly conserved 16S (*Rickettsia*, *Chlamydia*, *Xylella*,
*Synechococcus*/*Prochlorococcus*, *Desulfovibrio*, *Thermus*) all screened at ceiling 1.000
and k-mer 1.000. Genuine difficulty sits almost entirely in **subspecies comparisons and
tight species complexes** (*M. abscessus* subsp. 0.577, *X. campestris* group 0.647). Third
independent time this project has found the "16S cannot separate these" claim to fail given
adequate isolates and the right region.

## What Phase 1 did *not* and could not fix

The genome/isolate-identity gap (finding 1 above) is the one real blocker carried into Phase 2:
any new panel fetched for the pilot needs isolate/BioSample metadata attached at fetch time
(`efetch -format docsum` or GenBank flat-file, not bare FASTA), and splits must hold out
genome/isolate groups, not just haplotype strings, per protocol §12.1.

## Framing correction carried in from the review (applies to everything from here on)

Prior sessions' language ("gLM does not capture the signal") is retired. The defensible
question, per protocol §1 and the external review: *do frozen genomic representations
preserve information useful for fine-scale 16S discrimination, and how does accessibility
depend on model, region, panel, and readout* — not a yes/no verdict on "AI vs. classical."
Use the claim ladder (protocol §4) and decision table (§22) before writing any conclusion
sentence. The *S. aureus*/CoNS V1 result (k-mer+RF = 0.997, DNABERT-2 geometry+RF = 0.500)
stands as a real, audited data point — but it is one model, one readout (R3-ish bounded
probe on raw embeddings), one panel. It supports claim-ladder rung 2–3 ("no practically
useful discrimination under this tested frozen readout on this panel"), not rung 5
("absence of signal anywhere in the model").
