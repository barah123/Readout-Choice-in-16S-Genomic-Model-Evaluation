# Phase 2 (pilot) report — four frozen models × readout hierarchy × classical controls

**Protocol:** `../05_Third_Idea_Critical_Review_and_Benchmark_Protocol.md` §18 Phase 2
**Status:** development-phase pilot. **No confirmatory inference here by design** — protocol
§18 reserves p-values, multiplicity correction and equivalence testing for Phase 5 on locked
data. Everything below is descriptive and may inform Phase 3's freeze decisions.

## What was run

| component | choice |
|---|---|
| panels | *S. aureus* vs. 6 CoNS species; *B. anthracis* vs. *B. cereus* group |
| models | DNABERT-2, DNABERT-S, HyenaDNA, NTv2 (all four of the protocol §9 core) |
| readouts | R0 geometry (unsupervised pair AUC), R1 nearest-reference, R2 linear probe |
| classical | C_exact (exact lookup), C_align (identity 1-NN), C_kmer (4/5/6-mer + same probe as R2), C_random (matched-dim random projection) |
| CV | 5-fold `StratifiedGroupKFold` grouped on **isolate**, not haplotype string |
| metric | species-balanced accuracy, isolate-weighted, against matched isolate-weighted ceiling |

## Data: rebuilt with isolate metadata (Phase 1 audit findings 1 & 4)

`fetch_with_metadata.py` replaces the old bare-FASTA fetch, pairing `efetch` with `esummary`
so each record carries a structured `strain` field. Coverage is essentially complete: **794/796
staph records and 469/469 cereus records have a usable strain**, so isolate grouping is real,
not a regex guess.

Honest sizing of the leak Phase 1 flagged: **3 isolates / 6 records** in staph and **0** in
cereus carry more than one distinct sequence string. So the specific leak — same isolate,
different string, straddling a split — was real but small, and would not by itself have
changed Phase 3a's numbers. It is now closed structurally regardless.

The far bigger change is **panel size**, because the old fetch used `retmax=40`:

| panel | old (Phase 3a) | now |
|---|---|---|
| *S. aureus* target isolates | 19 | **57–61** |
| *B. anthracis* target isolates | 9 | **92–114** |
| cereus panel total records | 85 | **469** |

Note 49 sequence strings in staph (19 in cereus) are shared across >1 isolate. These are
biologically real and are kept, with the dependence handled by grouping rather than by
deletion (protocol §8.3).

## Result 1 — the Phase 3a headline does not survive, and the cause is the readout

Phase 3a reported staph V1 as the project's cleanest negative: k-mer+RF = 0.997, DNABERT-2 =
**0.500, exactly chance**, and a direct inspection of embedding geometry appeared to confirm
the representation simply did not encode species there. Re-run under this protocol:

| staph V1 | score | ceiling 0.982 |
|---|---|---|
| C_kmer (classical) | 0.981 | at ceiling |
| C_align | 0.966 | |
| DNABERT-2 **R0 geometry** | 0.923 | |
| DNABERT-2 **R1 nearest-reference** | 0.950 | |
| DNABERT-2 **R2 linear probe** | **0.964** | near ceiling |
| HyenaDNA R2 linear probe | **0.979** | at ceiling |

**DNABERT-2 at staph V1 is not at chance. It is at 0.964 against a 0.982 ceiling.** This is
precisely the failure mode the protocol's §4 claim ladder and §22 decision table were written
to prevent: Phase 3a measured one readout (RF on raw embeddings, on a 19-isolate panel) and
reported it in language that implied the representation lacked the information. Both the
readout and the sample size were doing the work.

This is the single most important finding of the pilot, and it is a correction to our own
prior claim, not a new negative result about any model.

## Result 2 — R0 geometry badly misranks representations

R0 (unsupervised cosine separation) disagrees sharply with what a probe can extract:

- staph V1: DNABERT-S R0 = 0.542 (near chance) → but R2 = 0.945.
- cereus V3/V5/V8: R0 = 0.79–0.83 for every model → but R2 ≈ 0.50 and the **ceiling is only
  0.549–0.576**. R0 here is "separating" something that cannot be the class label, since no
  classifier can exceed the ceiling. It is reading genuine structure (likely species-group
  substructure within the pooled neighbour class) that is not the target contrast.

Practical consequence for Phase 3: **R0 must not be used to rank models or regions.** It is a
diagnostic only. This directly supports the protocol's insistence on a readout hierarchy, and
it is a concrete reason DART-Eval's geometry-vs-probing distinction matters here.

## Result 3 — the frozen-vs-classical comparison is close to a tie

Best model readout minus best classical control, per region:

| panel | model-best wins (>+0.01) | loses (<−0.01) | within ±0.01 | mean Δ |
|---|---|---|---|---|
| cereus_group | 3/9 | 1/9 | 5/9 | **+0.007** |
| staph_CoNS | 1/9 | 0/9 | 8/9 | **+0.002** |

No region in either panel shows a frozen model beating the best classical control by more than
**+0.025**. The honest descriptive summary is that under a matched probe, frozen genomic
representations and 4/5/6-mer composition extract **nearly the same information** from these
16S windows — neither a clean embedding win nor the clean embedding failure Phase 3a implied.

Two caveats that cut against over-reading this, both pointing the same direction:

1. **C_random is nearly as good as C_kmer** (staph V1: 0.980 vs 0.981; cereus V2: 0.664 vs
   0.662). A random projection of the k-mer space retaining almost all performance means much
   of what is being measured reflects classifier flexibility over a high-dimensional
   composition space, not a specifically meaningful feature set. Protocol §11 anticipates
   exactly this; it argues for treating small Δ between any two high-dimensional feature sets
   as uninformative.
2. Several "wins" sit where the ceiling is low (cereus V3/V5/V8, ceiling 0.55–0.58). A +0.02
   delta against a 0.55 ceiling is not evidence of useful discrimination.

## Result 4 — the cereus panel is not uniformly null

Phase 0–2d treated *B. anthracis*/cereus as the expected-null panel. With 92–114 target
isolates instead of 9, that is too coarse:

| region | ceiling | best classical | best model |
|---|---|---|---|
| V6 | 0.840 | 0.804 | 0.800 |
| V7 | 0.757 | 0.702 | 0.709 |
| V2 | 0.737 | 0.662 | 0.668 |
| V3/V4/V5/V8 | 0.549–0.576 | ~0.50–0.51 | ~0.50–0.53 |

V6 and V7 support real, substantial discrimination (0.70–0.80) that the 9-isolate version of
this panel could not reveal. The panel is better described as **region-dependent**: largely
uninformative at V3/V4/V5/V8, genuinely informative at V6/V7. Notably V6 is where Chakravorty
et al. originally claimed anthracis/cereus separation — a claim `prototype_report.md` recorded
as failing to replicate. At this sample size it partially does replicate.

## Implications for Phase 3 (freeze)

1. **Primary readout should be R2 (linear probe).** R0 is diagnostic-only (Result 2); R1 sits
   between. Fixing this before locked testing is exactly what §18's freeze step is for.
2. **C_random must stay in the confirmatory comparison set**, not just as a footnote — Result 3
   shows it changes how a small embedding-vs-k-mer delta should be read.
3. **Report the ceiling beside every cell, and suppress cells with low ceilings** from headline
   claims. A margin against a 0.55 ceiling is not comparable to one against 0.98.
4. **Both current panels are now better-powered but are still only 2 blocks.** Protocol §6.2
   asks for ~15–30 independent taxonomic blocks. Phase 3 should decide the real block list;
   nothing here supports a general claim about genomic language models.
5. **Prior claims to retract or soften in any manuscript text:** the Phase 3a "gLM at exact
   chance on V1" result and the framing built on it. The underlying numbers were correctly
   computed for what they measured; the error was treating one readout at n=19 isolates as a
   statement about the representation.

## Files

`fetch_with_metadata.py` → `staph_panel_meta.json`, `cereus_panel_meta.json` ·
`isolate_audit.py` · `build_blocks.py` → `blocks.json`, `block_summary.json` ·
`embed_all_models.py` → `emb_{dnabert2,dnabert_s,hyenadna,ntv2}.npz` ·
`run_readouts.py` → `readout_results.json` · `summarize.py`
