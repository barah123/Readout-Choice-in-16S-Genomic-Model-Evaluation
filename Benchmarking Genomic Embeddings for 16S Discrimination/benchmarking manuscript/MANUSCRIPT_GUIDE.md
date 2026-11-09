# Manuscript writing guide

Everything needed to write the paper. Every number below names the file it comes from, so
nothing has to be retyped from memory or recomputed.

---

## Suggested title

> **Readout choice, not model choice, determines whether genomic language models help:
> a pre-registered benchmark on fine-scale 16S discrimination**

Alternative, if the venue prefers the result foregrounded:

> Pretrained genomic language models match k-mer counting for 16S species discrimination:
> a pre-registered benchmark across 58 bacterial groups

**Do not** title it "genomic language models fail at 16S." That claim is not supported — the
models reach 95–96% against a 97% ceiling. The finding is equivalence, and the contribution
is the methodology that makes that equivalence credible.

## One-paragraph summary

Pretrained genomic language models are increasingly proposed for sequence-identification
tasks, but it is unclear whether their representations carry information beyond simple
composition statistics. We benchmarked four published models (DNABERT-2, DNABERT-S,
HyenaDNA, Nucleotide Transformer v2) against 4/5/6-mer counting on target-versus-neighbour
16S discrimination across 58 bacterial groups, using a pre-registered analysis plan committed
before the confirmatory data were collected. Under a matched linear probe the two approaches
were statistically equivalent (−0.2 pts, 95% CI −0.6 to +0.2), and this held on an
independent locked test set, on groups selected specifically for difficulty, and under
simulated sequencing degradation. Two methodological findings explain why such comparisons
often mislead: on a like-for-like endpoint the choice of how a representation is read from a model
affected accuracy about four times more than the choice of model, and a random projection of
k-mer counts matched every learned representation while permuted labels collapsed every
method to chance. Within these panels the lowest accuracies occurred in subspecies and
species-complex comparisons rather than in genera reputed to have conserved 16S.

---

## Figures

| file | shows | key numbers |
|---|---|---|
| `fig1_main_result` | **Main result.** Per-group differences on the locked test + the three replications | −0.2 pts, CI −0.6 to +0.2; 17/20 within ±1 pt |
| `fig2_readout_dominates` | Readout vs model vs layer, **metrics kept separate** (R0 is pair AUC; R1/R2 are balanced accuracy) | R1→R2 +2.9 vs model 0.73 pts (~4x) |
| `fig3_methods_and_random` | All methods, bounded by the permuted-label floor and the empirical ceiling | chance 50.1%, random proj. 96.2%, k-mer 96.1%, ceiling 97.0% |
| `fig4_difficulty` | Saturation on the locked test; the hard-panel replication | 10/20 saturated; hard-panel mean −0.2, CI −1.3 to +0.7 |
| `fig5_robustness` | Degradation: both arms fall together | truncation −4 pts vs 1% error −2.8 pts |
| `fig6_region_ranking` | Task T2: regret from choosing the wrong region | locked test ~0.6 pts; hard panels ~2 pts |
| `fig7_ceiling` | Achieved accuracy against the information ceiling | methods sit close to the ceiling throughout |

All figures are 300 dpi PNG plus vector PDF, Okabe-Ito colourblind-safe palette.

**Figure 2 keeps metrics strictly separate.** R0 is a pair-level AUC and sits on its own
axis, explicitly labelled as not comparable; every quantitative statement uses R1 vs R2,
which share an endpoint. Panel d puts only balanced-accuracy quantities on one axis, with
the layer bar shaded differently because it was measured on the hard panels.

## Tables

| file | content |
|---|---|
| `table1_study_design` | three experiments, groups declared/eligible/analysed, record counts |
| `table2_endpoints` | primary + four Holm-corrected secondary contrasts with CIs and verdicts |
| `table3_method_performance` | all eight methods, mean accuracy and SD |
| `table4_readout_hierarchy` | R0/R1/R2 per model, metrics labelled; R1→R2 paired effect with CI |
| `table5_robustness` | six degradation conditions, both arms |
| `table6_layer_sweep` | HyenaDNA and NTv2 across five depths |
| `table7_region_ranking` | task T2: chosen-region accuracy, regret, hit rate |

Each exists as `.csv` (for typesetting) and `.md` (for reading).

---

## The numbers, with sources

### Scale
- 3 experiments, **58 groups analysed** (25 development, 20 locked test, 13 hard).
- 68 groups declared; 30,797 records fetched, **24,009 analysed**.
- Source: `results/tables/table1_study_design.csv`.

### Primary result (locked test, n = 20)
- DNABERT-2 − k-mer = **−0.21 pts**, 95% CI **[−0.60, +0.16]**, p = 0.253.
- Verdict: **equivalent** at the ±1 pt margin fixed in advance.
- Source: `results/locked_test_results.json` → `primary`.

### Replication
- Development (25 groups): **+0.14 pts**
- Locked test (20 groups): **−0.21 pts**
- Hard panels (13 groups): **−0.19 pts**, CI [−1.28, +0.66]
- Sources: `results/block_level_results.json`, `results/locked_test_results.json`,
  `results/extension_A_results.json`.

### Secondary contrasts (Holm-corrected)
- DNABERT-S −0.11 pts (equivalent); NTv2 −0.28 (equivalent);
  **HyenaDNA −0.88, Holm p = 0.004** — reliably worse but below the practical margin;
  best-of-4 +0.29, does not survive correction (p = 0.075).
- Source: `results/tables/table2_endpoints.csv`.

### Methodological finding 1 — the two controls
- Random projection of k-mer counts: **96.2%**, nominally above every learned representation.
  It is a distance-preserving compression (Johnson–Lindenstrauss), so it shows **redundancy
  in the k-mer feature space**, NOT that accuracy comes from classifier flexibility.
- **Label permutation** through the full pipeline: every arm falls to **50.1%** (chance),
  so the measured accuracy is genuine sequence-to-label signal.
- Sources: `results/tables/table3_method_performance.csv`,
  `results/label_permutation_control.json`.

### Methodological finding 2 — readout matters more than model
- **R0 is a pair-level AUC; R1 and R2 are isolate-level balanced accuracy. Never subtract
  across them.** An earlier draft did, producing an invalid "30x" claim.
- Same-endpoint effect, R1 → R2: **+2.8 to +3.0 pts** (DNABERT-2 +2.83 [+1.37, +4.57]).
- Spread across models at R2: **0.73 pts**. Ratio **~4x**.
- Layer depth: 0.53 pts, but measured on the hard panels, not the development set.
- Source: `results/tables/table4_readout_hierarchy.csv`, `results/extension_C_*.json`.

### The self-correction (worth a paragraph in Discussion)
- An earlier analysis in this project reported DNABERT-2 at **0.500 (chance)** on
  *S. aureus* V1. With a linear probe and 57 isolates instead of 19 the same model on the
  same bacteria scores **0.964**. Documented in `reports/06_pilot_and_self_correction.md`.

### Robustness
- Truncating 10 bases each end: **−4.2 pts** (k-mer), **−4.4** (DNABERT-2).
- 1% sequencing error: −2.8 / −1.8. Combined: −5.5 / −5.0.
- The apparent noise-robustness advantage was **one panel** (*francisella*); excluding it
  moves +0.0099 → −0.0006. Source: `reports/03_extension_B_degraded_sequences.md`.

### Region ranking (T2)
- Locked test: k-mer regret 0.65 pts, DNABERT-2 0.57; difference −0.09 [−0.46, +0.26].
- Hard panels: k-mer 1.95, DNABERT-2 2.02; difference +0.07 [−0.61, +0.65].
- **No method reliably identifies the single best region** (37–44% on the locked test), but
  the cost of a wrong choice is small there (~0.6 pts) and three times larger on hard groups
  (~2 pts).
- Source: `results/region_ranking_summary.json`, `results/tables/table7_region_ranking.csv`.

### Where difficulty lives
- Genera reputed to have conserved 16S (*Rickettsia*, *Chlamydia*, *Xylella*,
  *Synechococcus*/*Prochlorococcus*, *Desulfovibrio*, *Thermus*) screened at **100%**.
- Hardest groups are subspecies/complexes: *M. abscessus* subsp. **55%**,
  *X. campestris* group **63%**.
- Source: `reports/02_extension_A_hard_panels.md` §6.

---

## Suggested structure

1. **Introduction** — region choice matters; genomic LMs proposed; prior art (DNABERT-S,
   Micro16S, Qscore, Chakravorty 2007, Johnson 2019); the gap is a controlled comparison with
   adequate baselines.
2. **Methods** — panels and eligibility; isolate-level grouping; the readout hierarchy;
   classical controls including the random projection; information ceilings; pre-registration
   and the locked test.
3. **Results** — Fig 1 (main), Fig 3 (methods + random), Fig 2 (readout dominates),
   Fig 4 (difficulty), Fig 5 (robustness), Fig 6 (region ranking), Fig 7 (ceilings).
4. **Discussion** — why equivalence rather than failure; the two methodological traps and the
   documented self-correction; practical guidance; limitations.
5. **Data and code availability** — manifest with SHA-256 in `supplementary/S7`, `S8`.

## Claims to make, and claims to avoid

**Make:**
- Equivalence within ±1 pt, replicated three times, robust to degradation.
- Readout choice matters ~4x more than model choice, on a like-for-like comparison.
- A random projection of k-mer counts matches learned representations (redundancy), while permuted labels collapse to chance (the signal is real).
- Difficulty tracks taxonomic rank, not genus reputation.
- No accuracy case for a language model in a 16S region-recommendation tool.

**Avoid:**
- "Genomic language models fail / contain no information" — contradicted by 95–96% accuracy.
- Any claim about tasks, models or data types not tested.
- "Classical methods are superior" — the finding is equivalence, not superiority.
- Treating the best-of-4 signal (+0.29 pts) as real; it is most likely selection bias.

## Known reviewer objections, and where each is answered

| objection | answer |
|---|---|
| "Your panels were too easy" | Experiment 3, 13 difficulty-screened groups — `reports/02` |
| "Only clean sequences" | Experiment 4, six degradation conditions — `reports/03` |
| "You only tested the final layer" | Extension C for 2 of 4 models — `reports/04`; remains a stated limitation for the MosaicBERT pair |
| "Only one task" | Region ranking (T2) — `results/tables/table7` |
| "You tuned until it worked" | Plan committed 94 min before the first confirmatory sequence; `supplementary/S2`, `S3` |
| "Negative results are unreliable" | Equivalence testing with a pre-set margin, not a non-significant p-value |
| "Newer models might differ" | Open limitation; mitigated by model choice mattering ~4x less than readout |

## Still to do before submission

1. Write the manuscript text.
2. Decide the venue (benchmark/methods journals fit best).
3. Optionally extend the layer sweep to DNABERT-2/DNABERT-S — requires implementing and
   validating their unpad/repad logic.
4. Deposit code and the accession manifest in a public repository.
