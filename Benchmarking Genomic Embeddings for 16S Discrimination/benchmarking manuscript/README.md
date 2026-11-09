# Benchmarking manuscript — complete materials

Everything needed to write the paper. Start with **`MANUSCRIPT_GUIDE.md`**: it gives the
suggested title and structure, every number with the file it came from, the claims to make
and avoid, and where each predictable reviewer objection is already answered.

## Contents

```
MANUSCRIPT_GUIDE.md     write from this
figures/                7 figures, 300 dpi PNG + vector PDF
results/tables/         7 tables, CSV (typesetting) + Markdown (reading)
results/                source JSON for every number in the tables and figures
reports/                the full narrative reports, in reading order
supplementary/          protocol, frozen plan, deviations, panel definitions, data manifest
scripts/                every analysis script, runnable end to end
```

## Reports, in reading order

| file | content |
|---|---|
| `00_FINAL_REPORT.md` | the findings in plain language — the best single overview |
| `01_confirmatory_result.md` | the pre-registered locked-test result |
| `02_extension_A_hard_panels.md` | difficulty-selected replication; where 16S difficulty lives |
| `03_extension_B_degraded_sequences.md` | robustness to imperfect sequences |
| `04_extension_C_layers.md` | intermediate-layer sweep |
| `05_development_benchmark.md` | the 31-group development benchmark |
| `06_pilot_and_self_correction.md` | how an earlier wrong conclusion was found and fixed |
| `07_audit_of_prior_work.md` | audit of the project's own earlier pipeline |

## Supplementary

`S1` protocol · `S2` frozen analysis plan · `S3` protocol deviations · `S4` extension
pre-registrations · `S5`/`S6` panel definitions (committed before data collection) ·
`S7`/`S8` data manifest with SHA-256 checksums · `S9` synthetic validation tests

## Reproducing

Scripts in `scripts/` expect the phase directories one level up and the two conda
environments described in `../../CLAUDE.md`. Figures and tables regenerate with:

```bash
cd scripts
python make_tables.py
python make_figures.py
```

Raw downloads and model weights are excluded from version control. The exact sequences used
are listed in `supplementary/S8_accession_manifest.tsv`, with checksums in `S7`.

## The headline

Pretrained genomic language models and k-mer counting perform **equivalently** on
target-versus-neighbour 16S discrimination: −0.2 percentage points, 95% CI −0.6 to +0.2,
replicated across 58 bacterial groups in three experiments and robust to sequence
degradation.

The more reusable findings are methodological: **how a representation is read from a model
matters ~30× more than which model is used**, and **a random projection of k-mer counts
matches every learned representation** at this difficulty. Both are easy to overlook, and
overlooking either can produce a confident but wrong conclusion — as happened three times
during this work, each documented rather than removed.
