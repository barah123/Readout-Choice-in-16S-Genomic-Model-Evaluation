# Protocol deviations log

Every change to `FREEZE.md` or `test_blocks_definition.py` after their commit must be
recorded here using the form from `04_Implementation_Schemas_and_Checklists.md` §12, and
must state whether it was made **before or after** locked-test inspection. A deviation made
after inspection invalidates the confirmatory status of the affected endpoint.

Template:

```
Date:
Protocol/version:
Change:
Reason:
Before or after protected test inspection:
Who approved:
Affected hypotheses/endpoints:
Whether a new independent test is required:
Sensitivity analysis retaining original plan:
```

---

## Deviations recorded against the original Protocol 3 document

### D-001 — equivalence margin tightened from the protocol's illustrative ±0.03 to ±0.010

- **Date:** 2026-11-08
- **Protocol/version:** Protocol 3 v1.0 §15
- **Change:** confirmatory equivalence margin set to ±0.010 absolute balanced-accuracy units,
  not the ±0.03 given as a planning example in §15.
- **Reason:** §15 explicitly calls ±0.03 "not a universal accepted threshold" and requires a
  biological/application justification. Development showed a block-level SD of 0.01138 and an
  observed CI half-width of 0.0045 at n=25, so ±0.03 would be so wide that equivalence would
  be declared almost automatically and the test would carry no information. ±0.010 is the
  smallest difference that could plausibly change a region recommendation, and is attainable
  at the planned n.
- **Before or after protected test inspection:** **before** — locked test not yet fetched.
- **Affected endpoints:** primary equivalence test only.
- **New independent test required:** no.
- **Sensitivity analysis:** the CI will also be reported against ±0.03 so the result can be
  read under the protocol's original illustrative margin.

### D-002 — primary contrast uses single pre-named arms rather than a best-of selection

- **Date:** 2026-11-08
- **Protocol/version:** Protocol 3 v1.0 §13
- **Change:** primary estimand fixed as **DNABERT-2:R2_linear − C_kmer**, instead of
  "best frozen-model strategy vs best classical strategy".
- **Reason:** §13 permits either a named model or a development-selected best-model rule, but
  §14.2 requires any permutation/bootstrap reference to repeat every selection step contained
  in the statistic. Naming both arms from development (DNABERT-2 best at R2 in 20/25 blocks;
  C_kmer best classical in 17/25) removes two selection steps and makes the primary estimand
  unambiguous. The best-of-4 variant is retained as a Holm-corrected secondary.
- **Before or after protected test inspection:** **before**.
- **Affected endpoints:** primary; adds one secondary.
- **New independent test required:** no.
- **Sensitivity analysis:** best-of-4-R2 vs C_kmer reported as secondary test #5.

### D-003 — all 31 Phase 2 blocks reclassified as development; new locked test set created

- **Date:** 2026-11-08
- **Protocol/version:** Protocol 3 v1.0 §12.1, §18
- **Change:** the 31 blocks in `blocks_definition.py` are declared development-only and
  burned for confirmatory use. 22 new blocks, from 22 families disjoint from development,
  are declared as the single-use locked test set.
- **Reason:** readout choice, comparator set, pre-named arms, margin and precision plan were
  all derived from results on those 31 blocks. Re-using them as a test set is the
  "repeated test-guided model selection" failure named in §18's stopping rules.
- **Before or after protected test inspection:** **before** — the test taxa have never been
  fetched, sliced, embedded or scored.
- **Affected endpoints:** all confirmatory endpoints.
- **New independent test required:** this *is* that new test.
- **Sensitivity analysis:** Phase 2 development results remain reported in full, labelled
  development/exploratory throughout.
