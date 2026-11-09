# Can pretrained DNA language models pick the best 16S region? A benchmark

This folder contains a complete, pre-registered benchmark testing whether pretrained genomic
language models carry more useful information than simple k-mer counting, for the task of
telling a target bacterial species apart from its closest relatives using one 16S region.

**Short answer: no. They perform the same.** Details in
`phase6_report/FINAL_REPORT.md`.

## How to read this folder

Read in this order:

| file | what it is |
|---|---|
| `phase6_report/FINAL_REPORT.md` | **start here** — the findings, in plain language |
| `05_Third_Idea_Critical_Review_and_Benchmark_Protocol.md` | the protocol this followed |
| `phase3_freeze/FREEZE.md` | the analysis plan, fixed before the final data was collected |
| `phase5_evaluate/phase5_results.md` | the confirmatory result |
| `phase7_extensions/extension_A_results.md` | follow-up: hard groups |
| `phase7_extensions/extension_B_results.md` | follow-up: imperfect sequences |
| `phase7_extensions/extension_C_results.md` | follow-up: intermediate model layers |
| `00_phase_tracker.md` | running log of every phase |

## Why the results can be trusted

The analysis plan was written and committed to version control **before** the final dataset
existed. Git timestamps show the plan was fixed at 07:20:42 and the first final-dataset
record was downloaded at 08:54:19 the same day — 94 minutes later. The plan was not edited
afterwards. This matters because it rules out the most common way benchmarks mislead:
adjusting the method until the numbers look good.

The study used three separate sets of bacteria, 58 groups in total. The first (25 groups)
made every design decision. The second (20 groups, from entirely different bacterial
families) tested those decisions once. The third (13 groups) was chosen specifically because
the task is hard there, to check that the answer was not an artifact of easy cases. All three
gave the same answer, as did a fourth run on deliberately degraded sequences.

## Reproducing the work

Two conda environments are required (see `../CLAUDE.md` for how they were built):
`glm-benchmark` for DNABERT-2, DNABERT-S and HyenaDNA, and `glm-ntv2` for Nucleotide
Transformer v2, which needs an older `transformers`.

```bash
# Development set (31 groups) — all design decisions came from here
cd phase2_pilot
python fetch_blocks.py                    # downloads from NCBI; resumable
python build_blocks_full.py               # aligns and slices V1–V9
python embed_full.py dnabert2 dnabert_s hyenadna
conda run -n glm-ntv2 python embed_full.py ntv2
python run_readouts_full.py
python summarize_full.py

# Locked test set (22 groups) — used once
cd ../phase4_locked_test
python fetch_test.py
python build_test.py
python embed_test.py dnabert2 dnabert_s hyenadna
conda run -n glm-ntv2 python embed_test.py ntv2
python score_test.py
python make_manifest.py

# The pre-registered analysis
cd ../phase5_evaluate
python analyze_locked_test.py

# Exploratory follow-ups (clearly labelled as such)
cd ../phase6_report
python diag_difficulty.py
python diag_layers.py

# Phase 7 extensions, each pre-registered before its data was collected
cd ../phase7_extensions
python screen_hard_panels.py          # difficulty screen, classical methods only
python screen_round2.py               # deeper screen (Amendment A-1)
python run_extension_A.py all         # fetch + build the selected hard panels
python score_extension_A.py embed dnabert2 dnabert_s hyenadna
conda run -n glm-ntv2 python score_extension_A.py embed ntv2
python score_extension_A.py score     # Extension A: hard panels
python robustness_perturbation.py     # Extension B: degraded sequences
python extension_C_layers.py hyenadna # Extension C: layer sweep
conda run -n glm-ntv2 python extension_C_layers.py ntv2 embed
python extension_C_layers.py ntv2 score
```

Expected total runtime: roughly 3 hours on a laptop (Apple M4, 16 GB, no GPU). Downloading
from NCBI is the slowest step. Model weight files and raw downloads are excluded from version
control; the exact sequences used are recorded in
`phase4_locked_test/accession_manifest.tsv`, with SHA-256 checksums in `manifest_hashes.json`.

## What is deliberately not included

- No fine-tuning or model training. Every model is used exactly as published.
- No claim about tasks other than the one tested.
- No wet-lab validation.

## Honest record of mistakes

`phase1_audit/audit_report.md` and `phase2_pilot/phase2_pilot_report.md` document two errors
found in this project's own earlier work, both of which had produced a published-looking but
wrong conclusion. They are kept rather than quietly removed, because the way they were caught
is part of what the benchmark is for. `phase3_freeze/DEVIATIONS.md` records every departure
from the original protocol, with dates and reasons.
