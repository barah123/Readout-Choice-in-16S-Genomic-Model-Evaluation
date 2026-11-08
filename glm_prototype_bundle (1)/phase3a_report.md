# Phase 3a Report — S. aureus vs. coagulase-negative staphylococci: the first real (non-null) head-to-head

**Date:** 2026-11-06
**Scripts:** `fetch_staph.py`, `haplotype_gate_staph.py`, `phase3a_staph_baseline.py`
**Outputs:** `staph_group_16s.json`, `staph_dnabert2_embeddings.npz`, `phase3a_staph_summary.json`
**Panel:** target = *Staphylococcus aureus* (n=19 unique near-full-length sequences); near-neighbor
class = pooled coagulase-negative staphylococci (CoNS: *S. epidermidis*, *S. haemolyticus*,
*S. saprophyticus*, *S. hominis*, *S. lugdunensis*, *S. capitis*; n=176), from 787 raw NCBI
records fetched this session.

---

## 1. Why this panel, and what makes it different from the anthracis panel

MASTER_PLAN §4 item 2 called for moving off *B. anthracis* vs. *B. cereus* group — Panel A's
documented expected-null case — to a target with real published partial separation. Chakravorty
et al. 2007 name V1 as the best-separating region for *S. aureus* vs. CoNS. Per the new gate
Phase 2c established (count unique slice haplotypes per class per region *before* committing
compute), this panel was checked first (`haplotype_gate_staph.py`): V1 has **zero cross-class
haplotype ambiguity** and an oracle ceiling of **1.000** — unlike every region of the anthracis
panel, perfect separation is structurally achievable here. This is the first panel in the project
where a method comparison can produce a real winner rather than two methods both capped at chance.

## 2. Method

Four arms, all scored under **identical leave-one-haplotype-out (LOHO) CV** (so a slice
haplotype present in both classes never straddles a train/test split), each with its own
200-draw sequence-level label-permutation null and a shared Benjamini-Hochberg correction
across all 18 estimable (region × arm) cells:

- **B1** — identity: tie-aware 1-NN on pairwise sequence identity (the plan's simplest baseline).
- **B2** — k-mer (4/5/6-mer) + random forest.
- **gLM** — DNABERT-2 embedding + random forest. The embedding is built the Phase 2d way from
  the start (each unique slice embedded once, broadcast to every sequence sharing it) — verified
  pure function of sequence (purity check: max within-haplotype deviation = exactly `0.000e+00`
  in all 9 regions) rather than discovered broken after the fact.
- **B3/B4** (Shannon entropy, length/variable-site density) reported as diagnostic rankers, not
  scored classifiers, matching their role in Phase 2c.

Slicing, identity, k-mer, entropy, and RF-CV functions are copied verbatim from
`phase2c_b2_rf_baseline.py` per the CLAUDE.md convention.

## 3. Results

### 3.1 Oracle ceiling

| region | target haplotypes | CoNS haplotypes | ambiguous | oracle acc | oracle balanced acc |
|---|---|---|---|---|---|
| **V1** | 2 | 27 | **0** | **1.000** | **1.000** |
| V2 | 4 | 28 | 1 | 0.995 | 0.997 |
| V3 | 2 | 21 | 1 | 0.985 | 0.991 |
| V4 | 2 | 23 | 1 | 0.908 | 0.710 |
| V5 | 1 | 16 | 1 | 0.903 | 0.543 |
| V6 | 3 | 26 | 2 | 0.908 | 0.940 |
| V7 | 1 | 13 | 1 | 0.903 | 0.733 |
| V8 | 1 | 11 | 1 | 0.903 | 0.756 |
| **V9** | 3 | 14 | **0** | **1.000** | **1.000** |

V5, V7, V8 have only 1 target haplotype and are not estimable for any scored arm.

### 3.2 Full result, all arms, BH-corrected

| region | arm | balanced acc | p | q (BH) | survives | vs. oracle |
|---|---|---|---|---|---|---|
| **V1** | B1 identity | 0.523 | 0.065 | 0.116 | No | −0.477 |
| **V1** | **B2 k-mer** | **0.997** | 0.005 | 0.015 | **Yes** | **−0.003** |
| **V1** | **gLM embed** | **0.500** | 0.662 | 0.662 | **No** | **−0.500** |
| V2 | B1 identity | 0.902 | 0.005 | 0.015 | Yes | −0.095 |
| V2 | B2 k-mer | 0.962 | 0.005 | 0.015 | Yes | −0.035 |
| V2 | gLM embed | 0.971 | 0.005 | 0.015 | Yes | −0.026 |
| V3 | B1 identity | 0.546 | 0.020 | 0.045 | Yes | −0.445 |
| V3 | B2 k-mer | 0.526 | 0.174 | 0.261 | No | −0.465 |
| V3 | gLM embed | 0.500 | 0.488 | 0.576 | No | −0.491 |
| V4 | B1 identity | 0.500 | 0.348 | 0.448 | No | −0.211 |
| V4 | B2 k-mer | 0.494 | 0.512 | 0.576 | No | −0.216 |
| V4 | gLM embed | 0.500 | 0.572 | 0.606 | No | −0.210 |
| V6 | B1 identity | 0.531 | 0.055 | 0.109 | No | −0.410 |
| V6 | **B2 k-mer** | 0.638 | 0.020 | 0.045 | Yes | −0.302 |
| V6 | gLM embed | 0.526 | 0.080 | 0.130 | No | −0.414 |
| V9 | B1 identity | 0.507 | 0.229 | 0.317 | No | −0.493 |
| **V9** | **B2 k-mer** | **1.000** | 0.005 | 0.015 | **Yes** | **−0.000** |
| **V9** | **gLM embed** | **0.912** | 0.005 | 0.015 | **Yes** | **−0.088** |

(permutation p floor = 1/(200+1) = 0.005; `n/e` rows for V5/V7/V8 omitted.)

### 3.3 B3/B4 (diagnostic, unscored)

| region | B3 entropy (all) | B3 entropy (target only) | B4 mean length | B4 variable-site density |
|---|---|---|---|---|
| V1 | 0.484 | 0.031 | 29.0 | 1.000 |
| V2 | 0.265 | 0.031 | 113.9 | 0.958 |
| V3 | 0.185 | 0.005 | 65.0 | 0.898 |
| V4 | 0.176 | 0.075 | 107.0 | 0.914 |
| V5 | 0.167 | 0.000 | 59.0 | 0.946 |
| V6 | 0.242 | 0.020 | 60.0 | 0.982 |
| V7 | 0.107 | 0.000 | 56.0 | 0.929 |
| V8 | 0.120 | 0.000 | 52.0 | 0.745 |
| V9 | 0.334 | 0.012 | 30.4 | 0.800 |

## 4. The headline: at the one region literature predicts, k-mer nails it and gLM is at exact chance

**V1 is the cleanest test case the project has produced.** It is structurally perfect-separable
(oracle = 1.000, zero ambiguous haplotypes) and not trivially solved by raw identity (B1 fails,
0.523, p=0.065) — so it is a fair, non-trivial test of whether a trained classifier can recover
real but non-obvious structure. **B2 (k-mer + RF) comes within 0.003 of the oracle ceiling
(0.997, p=0.005). gLM (DNABERT-2 embedding + RF) is exactly at chance (0.500, p=0.662).**

This is not an artifact. Checked directly on the raw embedding vectors
(`staph_dnabert2_embeddings.npz`, region V1):

- The two unique *S. aureus* V1 haplotypes sit **0.155 cosine distance apart**.
- One of those two haplotypes sits only **0.078 cosine distance from its nearest CoNS
  neighbor** — closer to an off-target sequence than to the other *S. aureus* haplotype.
- This isn't a degenerate, information-free embedding space: CoNS-to-CoNS pairwise distances
  span 0.036–0.825, real compositional structure is present and resolved.

DNABERT-2's pretrained representation of this 29-nt region genuinely does not place the two
*S. aureus* haplotypes closer to each other than to CoNS sequences, while 4-6-mer composition
and a random forest separate them almost perfectly. This is the Micro16S precedent the analysis
plan pre-registers (§1.2: "a null result \[for gLM\] is a reportable finding") reproduced in its
cleanest form yet — not a null result from an uninformative panel (as anthracis/cereus was), but
a genuine loss on a target where the classical baseline essentially solves the problem.

**V9 shows a softer version of the same pattern:** k-mer hits the oracle ceiling exactly (1.000),
gLM recovers real but weaker signal (0.912, still BH-significant). **V2 is the one region where
gLM does not lose** — all three arms succeed, and gLM is nominally best (0.971 vs. k-mer's 0.962),
though both are close enough to the 0.997 ceiling that this is not a strong win.

**V3 is a thin, borderline case worth flagging honestly rather than folding into the headline.**
B1 technically survives BH correction (0.546, q=0.045) while k-mer (0.526) and gLM (0.500) do
not — by the letter of this test, B1 "separates" a group gLM fails to. But 0.546 is barely above
the 0.495 permutation-null mean; with 195 LOHO folds the test has power to detect a small, real
effect that is nonetheless practically weak (compare to V1 or V9, where the surviving arm is
within 0.09 of a ceiling of 1.0). This is a different kind of "B1 separates it" than Chakravorty's
original single-SNP claims, and shouldn't be read as equivalent evidence to V1.

## 5. Interpretation for the study design

**5.1 This triggers the plan's own pre-registered kill-switch language, in its narrow form, at
V1.** Plan §9: *"if Phase 2 shows gLM embeddings fail to separate Panel A groups that B1
separates, stop model expansion and pivot the paper to the benchmark + negative finding."* At
V1 specifically, B1 does *not* separate the groups (0.523, fails) — so the kill-switch's literal
trigger condition (gLM fails where B1 succeeds) is not met there. It comes closer to being met at
V3 (B1 survives at 0.546, gLM does not), but as §4 above argues, that's a thin statistical
margin, not a Chakravorty-style clean separation. **The more important fact is independent of
the kill-switch's exact wording: on the one region with real achievable structure and a
non-trivial solution, the classical baseline the plan already asks you to compare against
essentially solves it and the gLM arm does not solve it at all.** That is the substantive
negative result the benchmark is designed to be able to report, and this panel is the first time
the project has produced it with real statistical standing (oracle ceiling, LOHO-CV, permutation
null, BH correction, and a directly-verified embedding geometry) rather than the earlier panel's
confound of methodological artifacts.

**5.2 This is one gLM model on one panel — not yet grounds to invoke the "stop model expansion"
half of the kill-switch.** DNABERT-2's failure at V1 could be model-specific (e.g. its BPE
tokenizer may under-resolve a 29-nt sequence into very few tokens — V1 is the shortest region in
the panel) rather than a general property of genomic language models. HyenaDNA is already
verified working in `glm-benchmark` with zero patches; running it on this exact panel
(`staph_group_16s.json` / the same QA'd `seqs` list) is cheap relative to the compute already
spent here, and would directly test whether this is a DNABERT-2 artifact or reproduces. **This is
a genuine decision point for the user, not something to resolve unilaterally** — the plan's own
language says "stop and pivot," but a one-model, one-panel result is thin evidence for a full
paper-direction change.

**5.3 Construct validity is good.** V1 being the predicted-best region and showing real,
non-trivial, oracle-consistent structure — while B1 alone can't solve it — is itself a
confirmation that this panel and pipeline behave as expected (the "Gold standard A" check the
plan specifies in §... region ranking recovering published expectations). That's worth noting
independent of the gLM-vs-k-mer result: the pipeline finds real biology when real biology is
there, unlike the anthracis panel where it couldn't by construction.

**5.4 Protocol notes carried forward correctly this time.** Every lesson from Phase 2c/2d was
applied from the start here: the haplotype gate before any compute commitment (§5.3 of
`phase2cd_report.md`), LOHO-CV instead of plain LOO-NN, a permutation null per arm, BH correction,
oracle-ceiling reporting alongside every score, and an embedding purity check built into the
embedding step rather than retrofitted. No new errata this phase.

## 6. Recommended next step

1. **User decision needed:** given §5.2, decide whether to (a) treat V1 as sufficient evidence
   to pivot per §9's literal kill-switch language, or (b) run HyenaDNA (and/or NTv2) on this same
   panel first, to see whether the gLM failure at V1 is DNABERT-2-specific before deciding on the
   paper's framing. Recommendation: (b) — it's cheap (environment already verified, panel and
   embeddings cache already built) and would change the conclusion's generality either way.
2. If proceeding with (b): re-run the embedding + RF + permutation-null steps of
   `phase3a_staph_baseline.py` with HyenaDNA in place of DNABERT-2 on the same QA'd sequence list,
   reusing the k-mer/B1/oracle results (they don't depend on the embedding model).
3. Independent of (1)/(2): this panel is now a validated, real-signal test case. Future model
   comparisons (NTv2, HyenaDNA) and future panels (Panel B, *B. longum* subspp.) should be
   compared against this result, not the anthracis null panel, when assessing whether a gLM arm
   adds value to `regioneer`.
