# Readout choice in 16S genomic model evaluation

Analysis code, pre-registered plans, derived data and manuscript materials for a benchmark
comparing four pretrained genomic language models against k-mer counting on fine-scale 16S
rRNA discrimination.

**Headline result: the two approaches perform the same.** Mean difference −0.21 percentage
points of balanced accuracy, 95% confidence interval −0.60 to +0.16, against a 1-point
equivalence margin fixed before the confirmatory data were collected (paired-*t* TOST,
p = 4.1 × 10⁻⁴). Replicated across 58 bacterial groups in three experiments and robust to
simulated sequence degradation.

Two findings are more reusable than the headline:

- **How a representation is read from a model matters about four times more than which model
  is used** (+2.9 against 0.73 accuracy points, on a like-for-like endpoint).
- **A random projection of k-mer counts matches every learned representation**, while
  permuting labels reduces every method to chance. Either control alone leaves the
  interpretation underdetermined.

---

## Start here

| file | what it is |
|---|---|
| [`Benchmarking.../benchmarking manuscript/reports/00_FINAL_REPORT.md`](Benchmarking%20Genomic%20Embeddings%20for%2016S%20Discrimination/benchmarking%20manuscript/reports/00_FINAL_REPORT.md) | the findings in plain language |
| `.../benchmarking manuscript/MANUSCRIPT_GUIDE.md` | every number paired with the file it came from |
| `.../benchmarking manuscript/Manuscript_and_Author_Notes_Package/` | the PLOS-formatted manuscript and author notes |
| `.../phase3_freeze/FREEZE.md` | the analysis plan, committed before the confirmatory data existed |
| `.../00_phase_tracker.md` | running log of every phase |

## Repository layout

```
Benchmarking Genomic Embeddings for 16S Discrimination/
  phase1_audit/           audit of prior work, synthetic validation fixtures
  phase2_pilot/           development benchmark, 31 declared blocks
  phase3_freeze/          frozen analysis plan, locked-test panel definitions
  phase4_locked_test/     confirmatory data pipeline and manifests
  phase5_evaluate/        the pre-registered confirmatory analysis
  phase6_report/          diagnostics and the final report
  phase7_extensions/      hard panels, degradation, intermediate layers
  benchmarking manuscript/
    figures/              working PNG and vector PDF
    figures/plos_tiff/    300 dpi submission TIFFs
    results/tables/       manuscript tables as CSV and Markdown
    reports/              narrative reports in reading order
    supplementary/        protocol, frozen plan, deviations, manifests
    scripts/              every analysis script
glm_prototype_bundle (1)/ earlier prototype work the audit reports refer to
```

## Reproducing

Two conda environments are needed; Nucleotide Transformer v2 requires an older
`transformers` and cannot share the main one. See `requirements.txt` and the Environment
section of `CLAUDE.md`.

```bash
cd "Benchmarking Genomic Embeddings for 16S Discrimination/phase2_pilot"
python fetch_blocks.py && python build_blocks_full.py
python embed_full.py dnabert2 dnabert_s hyenadna
conda run -n glm-ntv2 python embed_full.py ntv2
python run_readouts_full.py && python summarize_full.py
```

The locked test, extensions, figures and tables follow the same pattern; full commands are
in the manuscript folder's README. Expect roughly three hours on a laptop with no GPU, most
of it downloading from NCBI.

**Not included:** model weight files (identified by SHA-256 or hub revision instead) and raw
NCBI downloads. The exact sequences used are listed in
`supplementary/S8_accession_manifest.tsv` with checksums in `S7`.

## On the pre-registration claim

Analysis plans and panel definitions were committed before the corresponding data were
retrieved. The freeze commit is dated 94 minutes before the first locked-test record.

**These are author-recorded git timestamps in a repository that had no remote until this
deposit, so they document the intended sequence rather than provide independent
certification.** Git commit dates are set by the committer and are not externally
verifiable. `S10_preregistration_evidence.md` gives the commit identifiers, the commands to
check them, and an explicit statement of the limits. Circumstantial support: 2 of 22
pre-declared locked-test blocks and 5 of 31 development blocks later failed the authors' own
eligibility criteria, and a pre-registered directional hypothesis did not replicate and is
reported as such.

## Honest record of corrections

Several confident-looking conclusions in this project turned out to be artifacts, and each
is documented rather than removed:

- an early analysis put DNABERT-2 at chance (50%) on a staphylococcal panel; with a linear
  probe and more isolates the same model reached 96.4%
- a claimed "thirty-fold" readout effect compared a pair-level AUC with balanced accuracy;
  the like-for-like figure is about fourfold
- an apparent noise-robustness advantage was driven entirely by one small panel
- a reported TOST p-value came from an under-resolved bootstrap tail count rather than a
  proper test

`phase3_freeze/DEVIATIONS.md` logs every protocol deviation with its date and justification.

## Licence

Code is MIT; data, text and figures are CC BY 4.0. Third-party material, including patched
model files, retains its original licence. See `LICENSE`.

## Citation

See `CITATION.cff`. Author metadata must be completed before the Zenodo release is created.
