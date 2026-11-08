# Phase 2 (full) — 31 taxonomic blocks, 4 frozen models, readout hierarchy vs. classical controls

**Protocol:** `../05_Third_Idea_Critical_Review_and_Benchmark_Protocol.md` §18 Phase 2, scaled
to the §6.2 block target.
**Status:** development-phase. **No confirmatory inference** — §18 reserves p-values,
multiplicity correction and equivalence testing for Phase 5 on locked data. Bootstrap CIs
below describe sampling variability across blocks; they are not confirmatory tests.
**Panel set was committed before any fetch** (`blocks_definition.py`, commit `1918e6f`).

## 1. Scale

| | pilot (2 panels) | full |
|---|---|---|
| blocks declared | 2 | **31** (30 species-rank + 1 subspecies stratum) |
| families | 2 | 27 |
| species queried | 11 | 154 |
| records | 1,265 | **12,356** |
| isolate × region rows | 10,142 | **97,621** |
| unique slices embedded / model | ~900 | **7,878** |
| block × region cells scored | 18 | **277** |

NCBI fetch: 154/154 species, **0 failures**. Strain metadata coverage ~100%, so isolate
grouping is a real key rather than a regex guess.

## 2. Eligibility — 26 of 31 blocks pass the pre-declared gate

Gate (declared before fetching): target ≥10 independent isolates **and** ≥2 neighbour species
with ≥5 isolates, in at least one region.

Region-cell accounting across all 277 cells, each failure mode coded distinctly per §19 —
never collapsed into a single "no signal":

| status | cells |
|---|---|
| eligible | **218** |
| insufficient_target_isolates | 44 |
| insufficient_neighbour_species | 9 |
| no_coverage_one_class (slicing QA removed a class) | 6 |

The 5 ineligible blocks are retained as **descriptive-only**, not deleted:

| block | reason |
|---|---|
| legionella | 55 target isolates but 0/4 neighbour species reach 5 |
| yersinia | *Y. pestis*: 5 isolates |
| bordetella | *B. pertussis*: 6 |
| borrelia | 6 |
| leptospira | 4 |

These are organisms now characterised almost entirely by whole-genome sequencing rather than
16S deposits (*Y. pestis* is additionally a select agent). That is a **reference-availability
artifact, not biology** — and it is the empirical realisation of the protocol's "class D
low-reference" stratum, which could not be declared a priori.

## 3. Primary contrast — frozen representations vs. classical, across independent blocks

Prespecified estimand (§13): mean paired difference, best R2 linear probe minus best classical
control, one value per independent block. Region chosen by the **classical** arm's strongest
cell and applied identically to both arms, so region choice cannot favour the models.

**n = 25 species-rank eligible blocks**

| quantity | value |
|---|---|
| mean Δ (best-R2 − best-classical) | **+0.0014** |
| 95% bootstrap CI over blocks | **[−0.0032, +0.0057]** |
| blocks favouring model / classical / tie (±0.01) | 5 / 2 / **18** |
| stringent, family-grouped (n = 22) | **−0.0001** [−0.0048, +0.0042] |

**The interval excludes any improvement larger than ~0.006 balanced-accuracy units.** Under
the protocol's §15 vocabulary this is *evidence of practically negligible improvement* rather
than "inconclusive" — the CI is narrow and straddles zero, on 25 independent blocks spanning
6 phyla. It is **not** evidence that the representations lack information (§22: that extension
is unsupported); R2 probes reach 0.94 mean balanced accuracy, close to the classical 0.947.

Largest single-block gains are small and sit where the ceiling leaves room: aeromonas V3
+0.022, klebsiella V3 +0.020, escherichia/shigella V1 +0.018. Largest losses: corynebacterium
V2 −0.033, staphylococcus V1 −0.017.

**The flexibility reference matters here.** Best-R2 minus `C_random` (a random projection of
the k-mer space) is **+0.0036 [+0.0000, +0.0071]** — the frozen models beat a random
projection by about as little as they differ from a tuned k-mer baseline. Per §11 this argues
against reading any of these small deltas as evidence about representation quality.

## 4. The readout hierarchy is the dominant effect — far larger than model choice

Averaged over the 25 eligible species-rank blocks:

| model | R0 geometry | R1 nearest-reference | R2 linear probe |
|---|---|---|---|
| DNABERT-2 | 0.741 | 0.918 | **0.946** |
| DNABERT-S | 0.722 | 0.913 | 0.942 |
| HyenaDNA | 0.670 | 0.909 | 0.939 |
| NTv2 | 0.724 | 0.912 | 0.940 |
| *(classical)* | – | – | *0.947* |

**Spread across readouts within a model is ~0.21–0.27. Spread across models within R2 is
0.007.** Which readout you use dominates which model you use by roughly thirty-fold. Any
benchmark reporting a single pooled-embedding-distance number — as our own Phase 3a did — is
measuring its readout far more than its representation.

Notably **DNABERT-S does not beat DNABERT-2** at R2 (0.942 vs 0.946) despite being the one
model of the four explicitly contrastively trained for species-aware embeddings. Its training
used 10 kbp genomic fragments; transfer to short conserved 16S windows is evidently not
automatic. This was flagged as an open question in the protocol's §2.2 and is now answered
for this task.

## 5. Calibration of our own a-priori predictions

`blocks_definition.py` recorded an expected separability class per block before any data was
seen. Measured ceilings:

| expected class | n | mean measured ceiling | range |
|---|---|---|---|
| A (informative) | 7 | 0.994 | 0.982–1.000 |
| B (partial) | 10 | 0.975 | 0.782–1.000 |
| C (highly ambiguous) | 9 | 0.912 | 0.848–0.990 |

Ordering is correct (A > B > C) but **we systematically over-predicted ambiguity**: blocks
expected to be 16S-unresolvable still reach a mean ceiling of 0.912. Textbook
"16S-indistinguishable" groups — *E. coli*/*Shigella* (0.931), *Klebsiella* (0.854),
*Enterobacter* (0.934), *Acinetobacter* (0.990), *Listeria* (0.960) — are substantially
separable once ~100 isolates per class are available and the best region is used. The
received wisdom appears to be partly a sample-size and region-choice artifact.

The one block that behaves as the literature predicts is **campylobacter (ceiling 0.782)**,
the lowest in the eligible set.

## 6. Subspecies stratum (reported separately, never pooled — §3)

`bifidobacterium_longum_subsp`: ceiling 0.930, classical 0.910, best R2 0.910. Even at
subspecies level the ceiling is not near 0.5 — but n is small and this is one block, so it
supports no general subspecies claim.

## 7. What this changes

1. **Phase 3a's headline is now doubly retracted.** The pilot already showed DNABERT-2 at
   staph V1 is 0.964, not 0.500. At full scale staphylococcus V1 is one of only two blocks
   where classical beats the models at all, by 0.017 — an ordinary small margin, not a
   qualitative failure.
2. **The honest overall finding is a near-exact tie with a narrow interval**, not a win for
   either side. On this evidence a `regioneer`-style tool has no accuracy reason to carry a
   genomic language model: classical k-mer composition matches it at ~1/1000th the compute.
   That is a defensible, publishable result *because* the interval is narrow and the block
   count is adequate.
3. **The strongest methodological contribution is the readout-vs-model decomposition** (§4),
   which quantifies a confound that has produced at least one wrong conclusion in this very
   project.
4. **Phase 3 freeze inputs:** primary readout R2; `C_random` stays in the comparator set;
   report ceilings beside every cell and suppress low-ceiling cells from headline claims;
   25 species-rank blocks + 1 subspecies stratum is the locked block list; 5 descriptive-only
   blocks stay in the accounting.

## 8. Limits

Development-phase, no confirmatory inference. Single pooling rule (masked mean) and single
layer (last) per model — §17.2 diagnostics not yet run, so "R2 is the best accessible readout"
is not established, only "R2 ≫ R0/R1". Blocks are genus-level; the family-grouped sensitivity
(n=22) is reported but the Enterobacteriaceae collapse is the only material one. Clean
reference sequences only — no primer, sequencing-error, or intragenomic-copy perturbation.

## 9. Files

`blocks_definition.py` (pre-committed) · `fetch_blocks.py` → `raw_species/` ·
`build_blocks_full.py` → `blocks_full.json`, `block_summary_full.json` ·
`embed_full.py` → `emb_full_*.npz` · `run_readouts_full.py` → `readout_results_full.json` ·
`summarize_full.py` → `block_level_results.json`
