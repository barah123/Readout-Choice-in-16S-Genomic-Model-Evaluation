# Phase 1 audit report

**Scope:** every claim made in `glm_prototype_bundle (1)` (Phases 0–2d, 3a) against protocol
§5.2's checklist, with concrete evidence re-derived this session — not the prior reports'
own self-description taken at face value, per the external review's instruction ("The quoted
claims... are user-provided descriptions, not independently audited findings").
**Baseline:** git tag `pre-audit-baseline` (commit `57f427f`).
**Panels audited:** *B. anthracis* vs. *B. cereus* group (43 sequences); *S. aureus* vs.
pooled CoNS (195 sequences).

---

## 1. "How many independent genomes, isolates, full-length genes, copies, unique region
sequences, and labels exist?"

**Finding: only the last two (unique region sequences, labels) are actually known. Genome,
isolate, and copy counts are not tracked anywhere in the existing pipeline, and cannot be
retrofit onto the already-fetched data without new NCBI metadata lookups.**

Evidence: both fetch scripts (`fetch_seqs.py`, `fetch_staph.py`) pull NCBI `nuccore` records
via `esearch`/`efetch` and keep only `(species_group, accession, sequence)` — the rest of
each FASTA header is discarded except the first whitespace-delimited token (the accession).
Checked the raw headers directly this session:

```
B_anthracis | QB222994.1 Bacillus anthracis strain RTB13 16S ribosomal RNA gene, partial sequence
S_aureus    | QB243438.1 Staphylococcus aureus strain NISCIUTQ30 16S ribosomal RNA gene, partial sequence
```

These are standalone Sanger 16S submissions (typical of isolate-ID studies), not
genome-assembly-derived excerpts — there is no assembly accession or BioSample linkage to
recover. A `strain` identifier *is* present in the free-text header and is currently thrown
away. A crude regex extraction of `strain (\S+)` found:
- cereus/anthracis panel: 84/85 records had a parsable strain token; **1 strain id
  (`B_cereus`/`3`) appears on 2 separate records** — i.e. at least one isolate was deposited
  more than once and is currently counted as 2 independent sequences.
- staph/CoNS panel: 778/787 records had a parsable strain token; **7 strain ids appear on 2–4
  records each** (e.g. `S_saprophyticus`/`ATCC` ×4, `S_hominis`/`AGATA` ×4) — at least 7
  isolates are represented by multiple records each, invisible to the current
  exact-sequence-string dedup.

This is a lower bound, not a precise count — the regex is crude (e.g. "ATCC" is a truncated
match of a longer catalog string, not necessarily one real isolate) and NCBI's own isolate
metadata (`BioSample`, `strain` qualifier) was never queried via `esummary`/`efetch` XML, only
the FASTA defline. **Verdict: gap, not fixable retroactively.** Closing it requires re-fetching
with isolate/BioSample metadata attached (`efetch -format docsum` or GenBank-flat-file mode,
not FASTA), which is Phase 2/4 scope, not something to patch onto the archived data.

## 2. "Are raw sequences or haplotypes collapsed? Which weights survive collapse?"

**Compliant, with one caveat.** Every scored pipeline since Phase 2c collapses to unique
per-region slice strings ("haplotypes") and explicitly carries `(n_target_seqs,
n_neighbor_seqs)` per haplotype through to the oracle-ceiling and balanced-accuracy
calculations (`defaultdict(lambda: [0, 0])` in both `phase2c_b2_rf_baseline.py` and
`phase3a_staph_baseline.py`). Raw sequence counts are never silently discarded.

**Caveat:** the collapse unit is the *haplotype string*, not the genome/isolate. Per finding
1, a handful of same-isolate pairs exist that are not string-identical (if they differ by
even one sequencing-error base) and would NOT be caught by this collapse — they would be
treated as two independent haplotypes, understating the true ambiguity in a way the current
reports do not surface. This is the same gap as finding 1, restated at the modeling level.

## 3. "Can the same observed sequence map to multiple species/subspecies?"

**Compliant, directly measured.** `hap_ambig` (count of haplotypes with members in both
classes) is computed and reported for every region in both panels — e.g. anthracis/cereus V1
and V8 each have exactly 1 such haplotype (all 9 anthracis sequences identical to a cereus-
group sequence); staph/CoNS V1 and V9 have 0. This is the single best-audited item in the
existing pipeline — it directly drives the oracle-ceiling computation, not a side note.

## 4. "Does nearest-neighbor evaluation include the query itself, another copy from its
parent genome, or an identical technical duplicate?"

**Partially compliant.** `loho_cv` and `loho_nn_cv` both hold out *all* rows sharing the
query's haplotype string together (`train = ~test` / `train = np.where(hap != h)[0]`), so a
query never sees itself or a string-identical duplicate during training. This correctly
closes the "identical technical duplicate" sub-case.

**It does not close the "other copy from the same parent genome" sub-case**, because genome
identity isn't tracked (finding 1). If the same isolate contributed two records differing by
one base (plausible — sequencing/submission noise), they'd land in two different haplotype
groups and could cross the train/test split as if independent. Given finding 1 found 1–7
same-strain record clusters per panel, this is a real, nonzero-probability leak, not a
theoretical one. **Verdict: gap, carries forward to Phase 2's split design (protocol §12.1:
"copies/regions from the same genome cannot cross roles").**

## 5. "What is the tie policy when multiple species are equally near?"

**Compliant.** Fixed in Phase 2c specifically *because* the original Phase 2/2b LOO-NN
accuracy was measuring `argmin` tie-break order, not separability (three regions reported
accuracy above their own Bayes-optimal ceiling). The fix — tie-aware fractional credit
averaged over all tied neighbors (`nn_expected_acc`, `loho_nn_cv`) — is documented in
`phase2cd_report.md` §2 and is the one issue protocol §18's "fix... tie... issues" instruction
most directly maps to. Confirmed already fixed; no further action this phase.

## 6. "Do batches, normalization, pooling, file order, or model mode change embeddings?"

**Compliant, independently re-verified this session, not just cited from the prior report.**
Re-ran the direct check rather than trusting `phase2cd_report.md`'s numbers:

```
reg   OLD cache maxdev   NEW cache maxdev   (max L2 spread among string-identical slices)
V1         7.126e-01          0.000e+00
V2         2.291e-01          0.000e+00
V4         2.136e-01          0.000e+00
V5         3.594e-01          0.000e+00
V6         3.841e-01          0.000e+00
V7         4.090e-01          0.000e+00
regions with batch-dependent OLD embeddings: 6/9
```

The original `dnabert2_embeddings.npz` cache is batch-composition-dependent in 6 of 9
regions (mechanism: `mask[:, -1] = False` drops the batch's last *column*, correct only for
the batch's longest row — see `phase2cd_report.md` §3 for the full mechanism). The corrected
cache (`dnabert2_embeddings_fixed.npz`, and the staph panel's `staph_dnabert2_embeddings.npz`,
built pure-function-of-sequence from the start) show exactly `0.000e+00` deviation in all 9
regions. **Verdict: fixed and confirmed, for DNABERT-2.** This check has not yet been run for
HyenaDNA or NTv2, since neither has been used for a real embedding pass yet — it must be
repeated for every new model in Phase 2, not assumed to carry over.

## 7. "Does any feature derive from labels, target-specific batching, filenames, or
test-set summary statistics?"

**One real gap found, previously undocumented.** `kmer_matrix()` builds its k-mer vocabulary
from *all* slices passed to it in one call:

```python
vocab = sorted({(k, s[i:i + k]) for k in ks for s in slices for i in range(len(s) - k + 1) ...})
```

Both `phase2c_b2_rf_baseline.py` and `phase3a_staph_baseline.py` call this once per region on
the *full* target+neighbor population (`slices = [seqs[i][2][name] for i in idx]`, `idx = tgt
+ cns`), **before** the LOHO-CV fold split happens. This means the feature space (which k-mers
exist as columns at all) is defined using the unlabeled sequence content of what will later be
held-out "test" rows in every fold — a transductive feature-construction leak. It is *not* a
label leak (held-out labels are never used, only held-out sequence content, and only to decide
which k-mers are representable as columns, not their values), and it's a common, usually
low-impact shortcut in k-mer pipelines generally — but it does not meet the stricter standard
protocol §10/04-schema asks for ("fit any scaler/PCA/normalization on support-training only").

**Verdict: real finding, not yet fixed.** Low expected impact (it affects feature-space
definition, not feature values or labels) but should be rebuilt as fit-per-fold when Phase 2's
classical control is implemented under the new protocol, not patched into the archived script.

No other label-, batch-, or filename-derived feature was found. `per_site_entropy` and
`variable_site_density` are order-invariant (column-wise aggregation over an unordered
multiset of slices, confirmed by reading — no indexing by position in the input list that
could correlate with species-grouped fetch order).

## 8. "Do all models observe precisely the same nucleotide bases?"

**Trivially compliant so far, not yet a real test.** Only one model (DNABERT-2) has been used
for actual scoring to date. Confirmed both arms draw from the same `seqs[i][2][name]` slice
per sequence/region (`Xk = kmer_matrix(slices)`, `Xe = embeddings[name][idx]`, same `slices`/
`idx`) — so k-mer and DNABERT-2 never saw different inputs for the same row. This check
becomes substantive only once a second model (HyenaDNA/NTv2) is added in Phase 2 — must be
re-verified then, not assumed from this single-model result.

## 9. "Are the purported environment checks only imports, or successful tokenization and
inference on realistic inputs?"

**Not re-verified this session — scheduled separately (tracked as a distinct Phase 1 task).**
`CLAUDE.md` and `MASTER_PLAN.md` both describe HyenaDNA as "confirmed working" and NTv2 as
passing "a determinism check (cosine similarity = 1.0)" from an earlier session. Per the audit
principle, these are claims to re-verify, not facts to inherit. Results will be written to
`env_smoke_tests.md` in this same directory once run.

---

## Summary verdict

| checklist item | verdict |
|---|---|
| 1. genome/isolate/copy counts | **gap** — not tracked; not retroactively fixable; defer to Phase 2/4 data collection |
| 2. haplotype collapse + weights | compliant (caveat: collapse unit is string, not genome — same root cause as #1) |
| 3. cross-label identical sequences | compliant, directly measured |
| 4. NN self/parent-genome/duplicate exclusion | partial — string-identical case closed; same-genome-different-string case open (root cause #1) |
| 5. tie policy | compliant, fixed in Phase 2c |
| 6. batch/pooling invariance | compliant for DNABERT-2, re-verified this session; must repeat per new model |
| 7. label/batch/filename-derived features | **one gap found** — k-mer vocab fit on pooled train+test content, not train-only |
| 8. identical inputs across models | compliant so far (only one model used); real test deferred to Phase 2 |
| 9. real inference smoke tests vs. import-only | pending — separate re-verification this phase |

Two concrete, previously-undocumented findings came out of this audit (#1/#4's genome-identity
gap, #7's k-mer vocab leak) that neither `phase2cd_report.md` nor `phase3a_report.md` flagged.
Both are addressed as Phase 2 design requirements (genome-aware fetch + splits; per-fold k-mer
vocab fitting), not patches to the archived Phase 2c/3a scripts.

## Synthetic fixture (protocol §8.4)

`synthetic_fixture.py` validates implementation mechanics against hand-computed ground truth,
independent of biological ambiguity — alignment-based region extraction under a known indel,
the oracle-ceiling formula (zero-overlap, fully-shared, and hand-constructed partial-ambiguity
cases), LOHO-CV recovery on a verified-zero-feature-overlap region, permutation-null collapse,
and tie-aware nearest-neighbour credit on a hand-constructed tie. All six parts pass.

One part failed on the first attempt and is kept as a documented finding rather than quietly
fixed: an initial haplotype design that looked zero-overlap at the *string* level actually
shared k-mer-level support across the class boundary (both classes included a haplotype ending
in "CCCCC" and one ending in "GGGGG"). With only 1 haplotype per class surviving in a given
LOHO fold, k-mer+RF locked onto the crossed suffix instead of the fully-diagnostic prefix and
inverted every held-out prediction — balanced accuracy = 0.0, reproducible across 5 seeds and
both with/without `class_weight='balanced'` (i.e. a feature-design confound, not stochastic
noise). This is independent, mechanistic confirmation of why `phase2cd_report.md` §5.3's
"≥5 haplotypes/class preferred, ≥2 is only a bare minimum" guidance exists — it is not merely
a sample-size heuristic, it is the point at which a single held-out point per class becomes
vulnerable to any partial feature confound, demonstrated here on a case with zero biological
ambiguity at all.
