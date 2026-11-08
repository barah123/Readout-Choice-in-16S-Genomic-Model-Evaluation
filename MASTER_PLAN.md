# Master Plan — gLM 16S Region Benchmark

**Last updated:** 2026-11-06
**Protocol:** `gLM 16S Region Benchmark — Analysis Plan v0.1.md` (project root)
**Working title (paper):** "Predicting the informative 16S rRNA hypervariable region for
target organisms: a genomic language model benchmark"
**Working title (tool):** `regioneer` (R package, not yet started)

This file is the living status record. Update it at the end of every phase — don't let
`prototype_report.md` / `phase2_report.md` / etc. be the only place a decision is recorded.

---

## 1. Where this came from

The project started as a Perplexity planning conversation (`literature.pdf`) that produced
the v0.1 analysis plan, extending the user's systematic review on breastfeeding vs. formula
feeding (infant gut microbiome; compositional differences robust, functional differences
method-dependent). The idea: profiling *region* choice is itself a method variable, and its
effect on taxonomic resolution should be predictable a priori from gLM embeddings — a
related but independent instance of the SR's "method choice shapes conclusions" thesis (not
a test of the same mechanism — see analysis plan §11 framing note, already corrected).

A prototype (Phases 0, 1, 2) was then built and run in a separate session/sandbox and handed
over as `glm_prototype_bundle (1)/`, along with an external review of the analysis plan. This
Claude Code session picked up from there, ran the review's requested diagnostic (Phase 2b, §4
below), and is now the place future phases should continue from.

---

## 2. Phase status

| Phase | Work | Status | Headline result |
|---|---|---|---|
| 0 | Environments, SILVA/NCBI download, region extraction, QA | **Done** (prototype-scale) | 83 sequences fetched, QA'd to 43; V1–V9 sliced via alignment to *E. coli* J01695 |
| 1 | Panel A baselines (B1 identity) — kill-switch check | **Done** | Pipeline correct; *B. anthracis* reclassified from validation anchor to **expected-null case** (see §3) |
| 2 | DNABERT-2 embeddings + scoring on Panel A | **Done** | Embeddings track true identity (Spearman ≈ −0.78 to −0.80); agree with B1's region ranking (0.73–0.76); caught and fixed a B1 identity-function bug |
| 2b | Reviewer-requested k-mer/entropy composition diagnostic | **Done** | Mixed: pair-level evidence says embedding ≈ composition; region-level LOO-NN accuracy says they diverge — **superseded by Phase 2c, see below** |
| — | **Environment setup**: `glm-benchmark` (torch/transformers, DNABERT-2 + HyenaDNA verified) and `glm-ntv2` (pinned transformers 4.34.1, NTv2 verified) | **Done** | See §3.5 below and `CLAUDE.md` Environment section |
| 2c | Proper B2 baseline (k-mer + **RF classifier**, LOHO-CV, permutation null, BH correction) + B3/B4 as formal baselines | **Done** (this session) | Resolved the Phase 2b divergence: it was two artifacts, not two feature spaces (see §3.6). **Zero (region × arm) cells survive BH correction** |
| 2d | Fixed a batch-composition label-leakage bug in the DNABERT-2 embedding cache discovered while building 2c | **Done** (this session) | New cache `dnabert2_embeddings_fixed.npz`; old cache must not be used for scoring (see §3.6) |
| 3a | Full B1-B4 + gLM baseline on *S. aureus* vs. CoNS (Panel A partial-separation target) | **Done** (this session) | **First real (non-null) head-to-head in the project.** At V1 (published-best region, oracle ceiling = 1.000): k-mer+RF = 0.997 (near-perfect), DNABERT-2 embedding = 0.500 (exact chance). See §3.7 below |
| 3 | HyenaDNA + NTv2 on the same panel; Panel B (incl. *B. longum* subspp.); *Mycobacterium*, *Haemophilus* | **Environments ready; recommended to run HyenaDNA on the Phase 3a panel next — see §4 item 1 and §5** | — |
| 4 | Mock-community empirical validation | Not started | — |
| 5 | Primer coverage layer | Not started (plan says move earlier — see §6) | — |
| 6 | R package + manuscript | Not started | — |

**Kill-switch status (plan §9):** **partially triggered, pending a second model — see §5.**
Phase 1's B1 baseline on the anthracis panel recovered the underlying biology correctly
(clonal within-anthracis identity, outgroup separation) — it just revealed the *published* V6
expectation doesn't generalize past a 2-strain comparison, a reportable finding not a pipeline
failure (`prototype_report.md` §4). Phase 3a is the first panel where the kill-switch's actual
condition ("gLM fails to separate groups that B1 separates") is testable with real signal
present, and at V1 specifically **B1 itself fails** (0.523, does not beat permutation null) so
the literal trigger isn't met there — but k-mer (the plan's other classical baseline) solves it
almost perfectly while gLM does not, which is the substance the kill-switch cares about even
though worded around B1 specifically. See `phase3a_report.md` §5 and §5 below before deciding
whether to invoke the "stop model expansion" half of the kill-switch.

---

## 3. What's been done — narrative

### Phase 0–1: Identity baseline on the *B. cereus* group (`prototype_report.md`)
Fetched 83 full-length 16S sequences (NCBI) for *B. anthracis*, *B. cereus*, *B.
thuringiensis*, *B. weihenstephanensis* (target + near-neighbors), and *B. subtilis*
(outgroup). QA'd to 43 unique, near-full-length sequences. Sliced V1–V9 by aligning each
sequence to *E. coli* J01695 and transferring Chakravorty et al. 2007 coordinates. Computed
pairwise identity per region for within-anthracis, anthracis-vs-cereus-group, and
anthracis-vs-outgroup comparisons (2,328 pairs).

**Result:** the pipeline is correct (within-anthracis identity = 1.0 in most regions,
outgroup cleanly separated everywhere), but Chakravorty's claim that V6 separates *B.
anthracis* from *B. cereus* by a single SNP **does not replicate** on broader NCBI sampling —
their claim rested on two anthracis strains vs. one cereus reference. Between-group maximum
identity is 1.0 in *every* region: no 16S region can give a perfect anthracis/cereus assay,
consistent with modern whole-genome phylogenetics (anthracis is a monophyletic subclade of
cereus sensu lato). **Action taken:** reclassified *B. anthracis* from "validation anchor" to
"expected-null case" for the rest of the benchmark.

### Phase 2: DNABERT-2 embeddings on the same slices (`phase2_report.md`)
Embedded the same 387 QA'd slices with DNABERT-2-117M (CPU, mean-pooled last-layer, three
compatibility patches documented in `dnabert2_patched/README_PATCHES.md`). Scored LOO
nearest-neighbor accuracy, silhouette, and pairwise embedding-vs-identity correlation per
region.

**Results:** embedding cosine distance tracks true sequence identity strongly (pooled
Spearman ≈ −0.80; the two anthracis V6 haplotypes, one SNP apart, separate cleanly in
embedding space — the needed validity check). Region rankings agree between the identity arm
and the embedding arm (Spearman 0.73). V6 separability is at chance (LOO 0.528), confirming
the expected-null call from Phase 1. V1/V4 carry the most partial signal (LOO ≈ 0.80) but no
region supports a perfect assay. **Bonus finding:** the embedding arm's identical-sequence
sanity check exposed a real bug in the Phase 0–1 identity function (it counted aligned, not
matching, positions — inflated same-length gapless pairs to 1.0). Fixed; Phase 0–1 numbers
recomputed; conclusions survived and one finding strengthened (V6 has the largest mean
between-group divergence of any region, consistent with Chakravorty's heterogeneity claim
even though it still doesn't separate the groups). **Protocol change adopted:** the
identical-sequence check is now a standing QA step for any future baseline/embedding
function, not a one-off.

### Phase 2b: Is the embedding signal just composition in disguise? (`phase2b_kmer_diagnostic_report.md`)
An external review of the analysis plan (full text in `initial_idea`) flagged the central risk
of the whole study: 16S is short, highly conserved, and non-coding, so it's not obvious a
priori that a gLM's embeddings carry separability signal beyond simple k-mer composition —
exactly the pattern the plan's own cited precedent (Micro16S) describes. The review asked for
one explicit diagnostic before trusting any future "gLM beats baseline" claim: regress the
gLM separability score on composition-based baseline scores and report the residual. This
session ran it on the existing Phase 2 data (no new embedding needed).

**Results (full numbers in the report):**
- **Pair-level (N=2,328, robust):** composition (k-mer+GC+length) explains 86% of the
  variance in embedding distance. Predicting true identity, composition alone (R²=0.711) and
  embedding alone (R²=0.714) are **statistically tied**; combining both adds only +0.034 R².
  The part of the embedding distance *not* explained by composition barely tracks true
  identity at all (Spearman drops from −0.78 to +0.25 and flips sign). **This arm of the
  evidence supports the reviewer's concern.**
- **Region-level (N=9, exploratory):** a continuous separability measure (silhouette) agrees
  almost perfectly between the k-mer arm and the embedding arm (Spearman 0.95) — consistent
  with the pair-level story. But a decision-relevant measure (LOO nearest-neighbor
  classification accuracy) does **not** agree (Spearman −0.07), and the embedding beats raw
  k-mer nearest-neighbor classification in 7 of 9 regions, sometimes by a wide margin (V1:
  0.49 vs 0.80; V4: 0.33 vs 0.81). **This arm of the evidence complicates the reviewer's
  concern** — but the k-mer "baseline" used here is a cheap raw-distance nearest-neighbor
  proxy, not the random-forest classifier the plan's actual B2 specifies, and N=9 is thin.

**Verdict carried forward:** the reviewer's concern is *largely* supported but not fully
settled. The unresolved part (region-level classification divergence) is exactly the thing a
proper B2 (k-mer + RF, as specified) would settle, and that hasn't been built yet. This
doesn't change the plan's existing priority of moving to partial-separation targets next —
if anything it sharpens the case for doing so (see §4).

### Phase 2c + 2d: the real B2 baseline, and two errata (`phase2cd_report.md`)

Built the baseline Phase 2b deferred — k-mer + random forest under leave-one-haplotype-out
(LOHO) CV with a 200-draw label-permutation null and Benjamini-Hochberg correction across all
(region × arm) cells — to settle the region-level LOO-NN divergence §5 flagged as the most
load-bearing open question. Doing so required looking at per-region class structure, which
surfaced two defects that invalidate every region-level number in Phases 2 and 2b:

- **Erratum A:** the target class isn't 9 independent sequences per region, it's **1–3
  unique slice haplotypes** (V1 and V8: exactly 1, shared with the outgroup class). 72–90% of
  samples have tied nearest neighbors at distance zero, so Phase 2/2b's LOO-NN accuracy was
  measuring `argmin`'s tie-break order (which correlates with the species-grouped input order),
  not separability. Proof: three regions (V1, V4, V7) reported accuracy *above* the
  Bayes-optimal ceiling for that region — a literal impossibility — and switching only the
  tie-break rule moves individual regions by up to 0.39 accuracy, reordering the region
  ranking that is the study's entire output.
- **Erratum B:** the Phase 2 DNABERT-2 embedding cache (`dnabert2_embeddings.npz`) leaks the
  class label through batch composition. Its mean-pooling mask dropped the batch's last
  *column* as `[SEP]`, correct only for the batch's longest sequence; shorter, padded rows
  kept their real `[SEP]` in the pooled mean instead. Sequences were batched in species-grouped
  order, so the resulting distortion correlated with class label — directly producing the two
  impossible-looking "significant" cells in the corrected analysis below (V2, V4 on the old
  cache). Fixed in `phase2d_reembed_fixed.py` (mask built from the tokenizer's own
  `attention_mask`; each unique slice embedded once and broadcast, making the embedding a pure
  function of the sequence by construction) → `dnabert2_embeddings_fixed.npz`.

**Corrected results, full numbers in `phase2cd_report.md`:**
- **Pair-level (N=2,328) Phase 2b conclusions are unaffected** — composition still explains
  86% of embedding distance; embedding's incremental R² over composition is +0.031 (was
  +0.034). Not a materially different finding.
- **Region-level, leakage-controlled: nothing survives.** Under LOHO-CV + permutation null +
  BH correction across all 14 estimable (region × arm) cells, **zero cells beat chance**
  (smallest q = 0.14). The embedding arm is flat (0.48–0.50 balanced accuracy, every p > 0.29)
  in every region; on the pre-fix cache the same pipeline had shown V2 and V4 "beating" the
  null (p = 0.010, 0.005) — both vanish once the batch leak is closed. Where the k-mer RF
  shows any nominal signal at all (V2, V5) it is nominally *better* than the gLM arm, the
  opposite of what Phase 2b's raw-distance proxy had suggested.
- V1 and V8 are **not estimable** (1 anthracis haplotype — LOHO-CV removes the class from
  training entirely); Phase 2/2b had reported numbers there anyway (0.800, 0.528).

**Verdict:** the Phase 2b region-level divergence is resolved, but not by either baseline
winning — both arms of that comparison were artifacts. On this (documented expected-null)
target, no method of the three (B1 identity, B2 k-mer+RF, gLM embedding) demonstrates
separability beyond the Bayes-optimal ceiling imposed by haplotype collapse (0.63–0.82
balanced accuracy, not 1.0). This reinforces rather than changes the existing priority of
moving to a partial-separation target (§4), and adds a new hard gate: **check unique slice
haplotypes per class per region before committing compute to any candidate target** (§5).

### Phase 3a: *S. aureus* vs. CoNS — the first real (non-null) head-to-head (`phase3a_report.md`)

Per §4 item 2's new hard gate, ran the haplotype pre-check (`haplotype_gate_staph.py`) on
*S. aureus* vs. pooled coagulase-negative staphylococci (CoNS) *before* committing fetch/embed
compute. Unlike anthracis/cereus, it passed with real structure: **V1 (Chakravorty's predicted-
best region) has zero cross-class haplotype ambiguity and an oracle ceiling of 1.000** — perfect
separation is structurally achievable here, which the anthracis panel could never show by
construction. Fetched 787 raw NCBI records (first attempt at `retmax=40` yielded only 5 unique
target sequences — most "16S ribosomal RNA"-titled NCBI records for these species are short
clinical-typing amplicons, not full-length genes; `retmax=300` fixed this, yielding 195 unique
near-full-length sequences, 19 *S. aureus* / 176 CoNS) and ran the full B1-B4 + gLM pipeline
(`phase3a_staph_baseline.py`) with every Phase 2c/2d lesson applied from the start (LOHO-CV,
per-arm permutation null, BH correction, embedding purity check) rather than retrofitted.

**Headline result:** at V1, **B2 (k-mer + RF) scores 0.997 balanced accuracy — 0.003 off the
oracle ceiling — while the gLM arm (DNABERT-2 embedding + RF) scores exactly 0.500 (chance,
p=0.662).** Verified this is not an artifact by checking the raw embedding geometry directly: the
two unique *S. aureus* V1 haplotypes sit 0.155 cosine distance apart, while one of them sits only
0.078 from its nearest CoNS neighbor — closer to an off-target sequence than to its own species'
other haplotype. The embedding space itself isn't degenerate (CoNS-to-CoNS distances span
0.036–0.825), it just doesn't encode "same species" at this specific 29-nt region. V9 shows the
same pattern more mildly (k-mer hits the ceiling exactly at 1.000, gLM recovers partial signal at
0.912). V2 is the one region where gLM does not lose (all three arms succeed, gLM nominally
best at 0.971 vs. k-mer's 0.962, though both are near the 0.997 ceiling). Full per-region, per-arm
table with BH q-values in `phase3a_report.md` §3.

**Verdict:** this is the Micro16S precedent (plan §1.2) reproduced in its cleanest form yet — not
a null result from an uninformative panel, but a genuine loss on a target where the classical
baseline essentially solves the problem outright. It also partially engages the plan's §9
kill-switch language (see the Phase status table above and §5 below) — but it's one gLM model
(DNABERT-2) on one panel, and V1 is the shortest region in the set (29nt), so whether this
reflects a general gLM limitation or something specific to DNABERT-2's tokenization of very short
sequences is still open. **This is a decision point for the user, not something resolved
unilaterally this session** — see §5.

### Environment setup

Before starting Phase 3 work, set up and fully verified (actual forward passes + determinism
checks, not just package installs) the two conda environments needed for every model in the
plan's §4.1 table:

- **`glm-benchmark`** — Biopython/numpy/pandas/scipy/matplotlib/scikit-learn + torch 2.15.0 +
  transformers 5.22.1. **DNABERT-2** required a 4th CPU-compatibility patch beyond the three
  already documented (`triton` has no macOS/Apple Silicon wheel at all — a platform gap the
  original three patches, made on a different sandbox, never hit); fixed by replacing the
  triton kernel file with a stub `ImportError`, loaded from a local patched copy
  (`dnabert2_local/`) rather than the bare HF repo id. **HyenaDNA** needed zero patches.
- **`glm-ntv2`** — a separate environment, transformers pinned to **4.34.1**. Nucleotide
  Transformer v2's custom code is old enough (saved under transformers `4.32.0.dev0`) that
  patching it forward to transformers 5.x hit a cascading series of config/API breaks, the
  last of which originated *inside* `transformers`' own internal weight-tying finalization
  code rather than the model repo's files — a signal to stop patching (risk of silently wrong
  tied weights) and pin an older transformers instead. Zero patches needed once pinned;
  determinism check confirms exact reproducibility (cosine similarity = 1.0).

Full technical detail, including the "why two environments" reasoning to apply to any future
model, is in `CLAUDE.md`'s Environment section — read that before running any embedding step.
This also answers the "local vs. HPC" question from this session: all three planned models run
fine CPU-only on this machine (M4, 16GB RAM) at the panel sizes used so far; HPC would only
become relevant if the reference panels were scaled up to the full SILVA/GTDB databases rather
than the curated target + nearest-neighbor panels the plan specifies.

---

## 4. What's next — recommended order

These are concrete, in priority order, synthesizing the plan's own kill-switch logic, the
external review's three "fixes before Phase 0," and what Phase 2b surfaced.

1. ~~Build the real B2 baseline (k-mer + random forest, plan §4.2) before anything else.~~
   **Done this session (Phase 2c/2d)** — see narrative above and `phase2cd_report.md`. Answer:
   zero (region × arm) cells beat a permutation null after BH correction; the Phase 2b
   region-level divergence this was meant to resolve turned out to be two artifacts (LOO-NN
   tie-breaking + an embedding batch-leakage bug), not a real gLM-vs-k-mer difference.

2. ~~Move off the expected-null target — now with a mandatory pre-check.~~ **Done this session
   (Phase 3a)** — *S. aureus* vs. CoNS passed the haplotype gate and produced the project's
   first real head-to-head result: k-mer+RF ≈ oracle ceiling at V1, gLM at exact chance. See
   narrative above and `phase3a_report.md`. The haplotype gate itself (count unique slice
   haplotypes per class per region before committing compute; ≥2 minimum to be estimable, ≥5
   suggested for a well-powered cell) is now a standing step for any further target —
   *B. longum* subspp. (Panel B) still needs it applied before any fetch/embed compute there.

3. **Decide how to respond to Phase 3a's kill-switch-adjacent result before spending more
   compute on new panels.** Two options, not mutually exclusive:
   - **(a) Run HyenaDNA (zero patches needed, already verified in `glm-benchmark`) on the exact
     same *S. aureus*/CoNS panel** (`staph_group_16s.json`, same QA'd sequence list) to test
     whether the V1 chance-level failure is general to gLMs or specific to DNABERT-2's
     tokenization of very short (29nt) sequences. Cheap relative to compute already spent —
     the panel, k-mer/B1 results, and oracle ceilings are already built and don't need rerunning.
   - **(b) Treat Phase 3a as sufficient and invoke §9's kill-switch now** — pivot the paper
     toward the benchmark-plus-negative-finding framing the plan already pre-authorizes,
     skipping further model expansion.
   - Recommendation: (a) first, since it's cheap and changes the conclusion's generality either
     way — a DNABERT-2-specific failure and a general-gLM failure are different papers. This
     is flagged as a user decision point in `phase3a_report.md` §5.2, not resolved here.

4. **Elevate the held-out-taxon transfer test to a headline, pre-registered hypothesis**, not
   a line item buried in the plan's §10 Known Pitfalls. The actual question anyone adopting
   `regioneer` cares about is whether embeddings generalize to a taxon absent from a gLM's
   pretraining — everything run so far is a sanity/comparison check, not that test. Design
   this explicitly: hold out one Panel A group entirely, test transfer on Panel B.

5. **Compute primer coverage (plan §7) in parallel starting now, not at Phase 5 — and recommend
   actual primer-bounded amplicons, not abstract V-region labels.** A region that wins on
   discriminability but can't be amplified for the target clade isn't a usable recommendation
   (known case already flagged in the plan: *Bifidobacterium* + 27F). This also reframes what
   the tool's final output unit should be: real assays are built from primer pairs, which often
   span multiple V-regions or only part of one (e.g. mopo16S's validated primers ran V2→V6/V7,
   not a clean single-region window) — "sequence V3" is a useful intermediate research result,
   but "this primer pair produces a resolving amplicon" is the actionable one. Doesn't block
   other items — can run concurrently since it only needs the target panels' taxonomy, not new
   sequence data.

6. ~~Set up a dedicated torch+transformers environment~~ **Done** — `glm-benchmark`
   (DNABERT-2 + HyenaDNA) and `glm-ntv2` (Nucleotide Transformer v2, pinned transformers
   4.34.1) are both built and verified. See §3.5 above and `CLAUDE.md` §Environment.

---

## 5. Open risks / things to watch

- **R vs. Python workflow deviation.** The plan specifies R for all stats/viz/packaging;
  everything built so far is Python. Decide deliberately whether to port before the package
  stage (Phase 6) or formally revise the plan — don't let this drift unexamined.
- **Reviewer's framing correction for §11 (manuscript tie-in)** has already been applied in
  the plan read so far — confirm the actual manuscript draft (when started) says "a further,
  independent instance of" the SR's method-dependence thesis, not "extends" or "resolves" it.
- **Identical-sequence QA check** must be run on any new baseline/scoring function before its
  output is trusted — this is how the Phase 0–1 bug was caught, and it won't catch itself a
  second time if skipped.
- ~~The B2 proxy-vs-real-classifier gap~~ **Resolved (Phase 2c, this session)** — see
  narrative above and `phase2cd_report.md`. Do not cite Phase 2/2b's region-level LOO-NN
  numbers in any manuscript text; they are superseded (impossible on their own terms — three
  regions exceeded the Bayes-optimal ceiling) and the corrected analysis found no surviving
  signal in either arm.
- **New: embedding-cache purity check, standing requirement for every future embedding run**
  (HyenaDNA, NTv2, DNABERT-S if ever added). Phase 2's DNABERT-2 cache leaked the class label
  through batch composition (a mean-pooling mask bug + species-grouped batch order) and
  silently inflated/distorted every region-level embedding score in Phases 2 and 2b — caught
  only because Phase 2c happened to check string-identical slices for identical embeddings.
  **Before scoring any embedding cache, assert that string-identical input sequences produce
  bit-identical (or near-bit-identical) embedding vectors.** The dedupe-then-broadcast pattern
  in `phase2d_reembed_fixed.py` (embed each unique string once, look up by string) satisfies
  this by construction and is cheaper than per-row embedding — reuse it directly rather than
  re-deriving the pooling fix.
- ~~Haplotype collapse is panel-dependent — re-check it for every new target~~ **Validated
  this session**: applied to *S. aureus*/CoNS before any fetch/embed compute (Phase 3a) and it
  correctly predicted which regions could show real separation (V1, V9: zero ambiguity) before
  any model was run. Still needs applying to *B. longum* subspp. (Panel B) before that target's
  compute is committed.
- **New (Phase 3a): decide the kill-switch response before running more panels.** DNABERT-2
  embeddings scored exact chance (0.500) at the one region (V1) with a real, non-trivially-
  achievable, oracle-consistent signal that k-mer+RF nearly solved (0.997). This engages plan
  §9's kill-switch language but on one model/one panel — see MASTER_PLAN §4 item 3 and
  `phase3a_report.md` §5.2 for the recommended next check (HyenaDNA on the same panel) before
  deciding whether to pivot the paper or continue model expansion. Don't spend further panel
  compute (Panel B, *Mycobacterium*, *Haemophilus*) until this is resolved, since the answer
  changes whether that compute is well spent.
- **New (Phase 3a): NCBI `esearch` retmax matters a lot for full-length yield.** A first fetch
  attempt at `retmax=40` per species yielded only 5 unique *S. aureus* sequences after QA (most
  "16S ribosomal RNA"-titled NCBI records for clinical species are short typing amplicons, not
  full-length genes) — looked like a genuine scarcity until checked against NCBI's true total
  hit counts (5,534 for *S. aureus*). Fixed with `retmax=300` → 195 unique sequences. Check true
  NCBI counts (`esearch` with `retmax=0`, read `.esearchresult.count`) before concluding a
  candidate panel is too thin to pursue, and don't reuse `retmax=40` from `fetch_seqs.py` as a
  template without reconsidering it for a new species.
- **New: no result so far has been tested for robustness to realistic assay perturbation.**
  Every score to date (B1/B2/B3/B4/gLM, all phases) is computed on clean reference sequences.
  None of it has been checked against primer mismatches, sequencing error, intragenomic 16S
  copy-number variation within a single genome, or uncertain/incomplete reference coverage —
  all factors known in the literature to affect real-world region/primer performance (mopo16S's
  entire premise; Johnson et al. 2019 on intragenomic copy variation). A region or amplicon that
  looks like a clean winner on curated NCBI reference sequences could still fail in practice.
  This should be an explicit, named gap in any manuscript claim of a "validated recommendation,"
  not just an assumed limitation — and is a candidate Phase 4/5 workstream alongside
  mock-community validation, not something the current B1–B4/gLM scripts attempt.

---

## 6. File index (bundle contents, in phase order)

See `CLAUDE.md` for the full annotated directory tree. Reports, in the order to read them for
full context: `prototype_report.md` → `phase2_report.md` → `phase2b_kmer_diagnostic_report.md`
→ `phase2cd_report.md` (B2 RF baseline + both errata — supersedes the region-level numbers in
the two reports before it) → `phase3a_report.md` (first real non-null head-to-head, *S. aureus*
vs. CoNS) → this file.
