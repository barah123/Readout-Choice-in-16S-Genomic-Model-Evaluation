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
| 2 | Pilot: DNABERT-2, HyenaDNA, NTv2, DNABERT-S if runnable; classical controls; 3–5 dev blocks; primary + linear-probe readouts | Not started |
| 3 | Freeze: final model/checkpoint list, panel eligibility, primary metric, readout rule, multiplicity plan, margins, compute budget | Not started |
| 4 | Scale: locked model × panel × region matrix on audited data | Not started |
| 5 | Evaluate: freeze predictions, unseal labels, prespecified summaries | Not started |
| 6 | Diagnose/report: exploratory analyses, reproducible release, conditional conclusions | Not started |

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
