# Final review — response to the twelve required pre-submission checks

**Date:** 2026-11-13
**Manuscript:** *Pretrained DNA models and k-mer counting on 16S hypervariable-region
discrimination* (`Manuscript_Draft_Pretrained_DNA_Models_16S.md`,
`Manuscript_PLOS_CompBiol.docx`)
**Repository state:** `main` at commit `46d33a9`, 37 commits, 10 phase tags

Every number below was recomputed from the frozen result files in this repository during this
pass, or read directly off the environment and checkpoint files. Nothing was taken from the
earlier notes (`Query_Resolution.md`, `Scientific_Confirmations.md`) without re-deriving it,
and nothing is marked resolved unless the supporting file exists and is named here. Where a
check could not be satisfied, it is marked **NOT RESOLVED** or **AUTHOR ACTION** rather than
narrated as closed.

Units throughout: **pp** = percentage points of isolate-weighted, species-balanced accuracy.
The unit of analysis is the block (bacterial group); n = 20 eligible locked-test blocks.

---

## Summary of the twelve checks

| # | Check | Status |
|---|---|---|
| 1 | Multiplicity plan | **Verified against the frozen plan.** Table 2 is correct as structured. Methods wording must change; two families now specified separately, adjusted results supplied. |
| 2 | Secondary equivalence | **Resolved.** Component one-sided outputs now computed and deposited for all four secondary contrasts. HyenaDNA does **not** meet formal equivalence. |
| 3 | Best-of-four rule | **Resolved by documentation.** Selection uses held-out outcomes. It is an oracle comparison; selection gap quantified at +0.66 pp. |
| 4 | Bootstrap difference test | **Resolved.** Fully specified, including that it is a percentile-interval inversion, not a shifted-null test, and that HyenaDNA's p sits at the estimator's resolution floor. |
| 5 | Panel accounting | **Verified in full**, including the 58 ↔ 59 relation. One framing inconsistency found and reported. |
| 6 | Difficult-panel interval | **Resolved.** Exact frozen output is −1.27 pp. The −1.28 is a rounding error with a traceable origin. Placeholder "≈ −1.3 to +0.7" must be removed. |
| 7 | Noise sensitivity estimate | **Resolved.** +0.99 and +0.85 are two different estimators of the same condition on the same 20 panels. Both collapse without one panel. |
| 8 | Ambiguity ceiling | **Resolved; Figure 7 rebuilt.** A fold-matched ceiling is a strict bound on all 66 cells and changes the mean by 0.001 pp. Methods edit still outstanding. |
| 9 | Orientation and grouping | **NOT RESOLVED.** Orientation defect disclosed but unrepaired. Sequence-disjoint analyses not performed. Grouping exposure newly measured: 92.0% of held-out isolates share their exact slice with the training folds. |
| 10 | Software and patched models | **Resolved.** Environment files did not previously exist; exported and deposited this pass. All four checkpoint identities verify exactly. Patches remain numerically unvalidated — must be stated. |
| 11 | Reproducibility and figures | **Partially resolved.** URL, licence, commit and tags verified. **All seven figures audited, rebuilt where needed, and regenerated to PLOS spec; both sets of legends rewritten (§11a–11b).** Release DOI, two of three manifests, unit-level predictions and figure insertion remain open. |

Three cross-cutting corrections affecting more than one check are consolidated in
[§13](#13-consolidated-corrections-required-before-submission).

---

## 1. Multiplicity plan

**Concern.** The source Methods place five comparisons in one confirmatory family, but the
table adjusts four secondary difference tests and leaves the primary unadjusted. Verify the
frozen plan, specify separate difference and equivalence families, and supply the correct
adjusted results.

**What was checked.** `phase3_freeze/FREEZE.md` §7 (= `S2 Text` §7), committed 2026-11-08
before the locked test was fetched, read verbatim against `Table 2` and against
`scripts/analyze_locked_test.py`.

**Finding — the table is right and the prose is wrong.** The frozen plan declares a
confirmatory family of exactly five tests with an explicitly *hierarchical* correction
structure:

> | # | test | correction |
> |---|---|---|
> | 1 | PRIMARY: DNABERT-2:R2 − C_kmer, equivalence at ±0.010 | none (single primary) |
> | 2–5 | SECONDARY: DNABERT-S, HyenaDNA, NTv2, best-of-4-R2, each vs C_kmer | Holm across these 4 |

Table 2 implements this exactly, and `analyze_locked_test.holm()` is applied to the four
secondary p-values only. **No silent change was made to the multiplicity plan.** What is
wrong is the Methods sentence, which reads as though one correction were applied across all
five. It must say: *a hierarchical confirmatory family of five — one prespecified primary
tested unadjusted, four secondary contrasts Holm-adjusted among themselves.*

**Two families, specified separately.** The plan tests two distinct null hypotheses and the
manuscript conflates them in one table:

- **Difference family** — H₀: δ = 0. Inference by the two-sided percentile bootstrap
  (§4 below). Primary unadjusted; Holm across the four secondary.
- **Equivalence family** — H₀: |δ| ≥ 1.00 pp. The frozen decision rule (FREEZE.md §4) is
  *interval containment*: equivalence is declared only if the entire 95% bootstrap CI lies
  inside ±1.00 pp. **No multiplicity adjustment was prespecified for equivalence.** Component
  one-sided tests are supplied in §2 below, and Holm across the four secondary equivalence
  tests is reported there as a post hoc sensitivity, labelled as such.

**Adjusted results — difference family** (recomputed, not retyped):

| Contrast | Mean (pp) | 95% CI (pp) | p (bootstrap) | p (Holm) | p (paired t) | p (Holm, t) |
|---|---|---|---|---|---|---|
| DNABERT-2 − k-mer (**primary**) | −0.211 | [−0.595, +0.163] | 0.253 | *not adjusted* | 0.302 | *not adjusted* |
| DNABERT-S − k-mer | −0.106 | [−0.468, +0.207] | 0.580 | 0.580 | 0.566 | 0.566 |
| HyenaDNA − k-mer | −0.885 | [−1.483, −0.377] | ≤0.001 | ≤0.004 | 0.0077 | 0.031 |
| NTv2 − k-mer | −0.281 | [−0.683, +0.113] | 0.168 | 0.336 | 0.197 | 0.394 |
| best-of-4 − k-mer | +0.286 | [+0.026, +0.571] | 0.025 | 0.075 | 0.061 | 0.182 |

Every value reproduces the frozen `phase5_evaluate/locked_test_results.json`. The paired-t
columns are added only because the bootstrap p for HyenaDNA is resolution-limited (§4); they
are a sensitivity, not a replacement. **No adjusted p-value was invented** — all are Holm
output from the raw values above, produced by the deposited script
`scripts/audit_multiplicity_and_equivalence.py`.

**Status: verified; Methods wording correction required.**

---

## 2. Secondary equivalence

**Concern.** Provide component equivalence-test outputs for secondary models if formal
equivalence is being claimed. Table 2 retains the difference-test numbers but should no
longer present every small interval as a reported formal equivalence test.

**Finding.** Component outputs did not exist for the secondary arms; only the primary had a
TOST. They are now computed on the same paired block-level estimand, margin ±1.00 pp,
n = 20 blocks, 19 df:

| Contrast | Mean (pp) | SD (pp) | t (H₀: δ ≤ −1) | p_lower | t (H₀: δ ≥ +1) | p_upper | **p_TOST** | 90% CI (pp) |
|---|---|---|---|---|---|---|---|---|
| DNABERT-2 − k-mer (**primary**) | −0.211 | 0.889 | 3.969 | 4.11 × 10⁻⁴ | −6.089 | 3.72 × 10⁻⁶ | **4.11 × 10⁻⁴** | [−0.555, +0.133] |
| DNABERT-S − k-mer | −0.106 | 0.814 | 4.908 | 4.89 × 10⁻⁵ | −6.077 | 3.82 × 10⁻⁶ | **4.89 × 10⁻⁵** | [−0.421, +0.208] |
| HyenaDNA − k-mer | −0.885 | 1.328 | 0.389 | 0.351 | −6.347 | 2.17 × 10⁻⁶ | **0.351** | [−1.398, −0.371] |
| NTv2 − k-mer | −0.281 | 0.940 | 3.423 | 1.43 × 10⁻³ | −6.096 | 3.66 × 10⁻⁶ | **1.43 × 10⁻³** | [−0.644, +0.082] |
| best-of-4 − k-mer | +0.286 | 0.641 | 8.967 | 1.48 × 10⁻⁸ | −4.977 | 4.19 × 10⁻⁵ | **4.19 × 10⁻⁵** | [+0.038, +0.534] |

Holm across the four **secondary** equivalence tests (post hoc; no adjustment was
prespecified for this family): DNABERT-S 1.68 × 10⁻⁴; NTv2 2.85 × 10⁻³; best-of-4
1.68 × 10⁻⁴; HyenaDNA 0.351. The three equivalence declarations survive adjustment.

Distribution-free corroboration (bootstrap TOST, B = 20,000, seed 1729; tail counts have a
floor of 1/B = 5 × 10⁻⁵, so a zero count must be reported as *< 5 × 10⁻⁵*, never as zero):
primary p = 1.0 × 10⁻⁴, 90% CI [−0.529, +0.105] pp; DNABERT-S < 5 × 10⁻⁵, [−0.407, +0.171];
NTv2 3.0 × 10⁻⁴, [−0.627, +0.052]; best-of-4 < 5 × 10⁻⁵, [+0.066, +0.524]; HyenaDNA 0.332,
[−1.379, −0.433]. Both routes agree on every verdict.

**Three consequences for the manuscript.**

1. **HyenaDNA must not appear with an equivalence verdict of any kind.** p_TOST = 0.351 and
   both its 90% and 95% intervals fall outside the margin. The current Table 2 entry
   ("INCONCLUSIVE at this margin") is correct and should be retained, but the Results text
   must say that formal equivalence was tested and *not* established, not merely that the
   interval was wide.
2. **DNABERT-S and NTv2 can now be stated as formal equivalence results**, which the abstract
   currently asserts on interval containment alone ("also met the reported equivalence
   criterion"). That claim is now supported: p_TOST = 4.89 × 10⁻⁵ and 1.43 × 10⁻³, Holm-adjusted
   1.68 × 10⁻⁴ and 2.85 × 10⁻³.
3. **The frozen containment rule is conservative, and that should be said.** A TOST at
   α = 0.05 corresponds to the **90%** interval; the frozen plan requires the **95%** interval
   to lie inside the margin. The plan therefore makes equivalence harder to declare than the
   nominal α implies. Both intervals are tabulated above so a reader can check either.
4. **best-of-4 should carry no equivalence verdict**, for the separate reason given in §3: it
   is an oracle contrast, so applying the containment rule to it mechanically produces a
   verdict about a quantity nobody can obtain in practice.

**Status: resolved.** Outputs deposited in `results/multiplicity_equivalence_audit.json`.

---

## 3. Best-of-four rule

**Concern.** State exactly how the best model was chosen, and whether held-out labels or
scores influenced selection. If selection used test outcomes, identify it as an oracle
comparison rather than an operational rule.

**Finding — selection used held-out outcomes. It is an oracle.** In
`scripts/analyze_locked_test.py` the contrast is

```python
best = np.array([max(blocks[b][f"{m}:R2"] for m in MODELS) - blocks[b]["C_kmer"]
                 for b in names])
```

where `blocks[b][f"{m}:R2"]` is model *m*'s per-block mean held-out balanced accuracy across
the five outer folds — **the same quantity the contrast itself is computed from**. The maximum
is taken per block, after the held-out scores are known. Held-out labels enter through those
scores. There is no training-set or inner-fold criterion, and no prespecified rule for
choosing among the four.

**It is also not a stable choice.** The per-block argmax over the 20 locked-test blocks is
DNABERT-2 in 12 blocks, DNABERT-S in 5, NTv2 in 2, HyenaDNA in 1 — so no single model could
have been named in advance to reproduce it.

**The selection gap, quantified.** Taking a per-block maximum over four noisy measurements is
upward-biased by construction, and the size of that bias is directly readable:

| Quantity | Value (pp) |
|---|---|
| Mean of the four individual contrasts | **−0.371** |
| Per-block oracle maximum (best-of-4) | **+0.286** |
| **Selection gap** (attributable to per-block maximisation) | **+0.657** |
| Operational pre-named rule (always DNABERT-2, named from development) | **−0.211** |

The entire sign reversal — from −0.37 pp to +0.29 pp — is produced by the selection step. The
same pattern appears on the hard panels, where best-of-4 is +0.70 pp
[−0.10, +1.35] against a primary of −0.19 pp, and `extension_A_results.md` §4 already names
selection bias as the default explanation.

**Required wording.** Report best-of-4 as an **oracle (hindsight) comparison** conditioned on
held-out outcomes; report the +0.66 pp selection gap beside it; state that no operational rule
achieves it and that the deployable pre-named rule gives −0.21 pp; remove the equivalence
verdict from its Table 2 row (§2). The existing Table 2 note — "Best-of-four is a hindsight
comparison and is not a deployable model-choice rule" — is the right content but is doing this
work in a footnote while the row itself reads like a result.

**Status: resolved by documentation.**

---

## 4. Bootstrap difference test

**Concern.** Document how the two-sided p = 0.253 and the secondary p-values were calculated.
The earlier attachment explains the paired-t TOST equivalence procedure but not this separate
test.

**Full specification** (`analyze_locked_test.boot_ci`, reproduced verbatim in the deposited
audit script):

- **Estimand:** mean paired per-block difference in balanced accuracy.
- **Resampling unit:** the block. Folds, pairwise distances and seeds never count as n.
- **Draws:** B = 2,000 (frozen). **Seeds:** 1729 for the primary; **1736** (`SEED + 7`) for
  all four secondary contrasts — a different seed from the primary, which should be stated.
- **Interval:** percentile, 2.5th and 97.5th of the bootstrap means.
- **p-value:** `p = 2 × min(Pr[mean* ≤ 0], Pr[mean* ≥ 0])`, computed on the **unshifted**
  resampling distribution.

**The critical disclosure.** Because the resampling distribution is not shifted to the null,
this p is an **inversion of the percentile interval**, not a shifted-null bootstrap
hypothesis test. It must be described that way. The same construction was already criticised
and corrected once in this project for the equivalence test
(`Scientific_Confirmations.md` §6, where a percentile tail count was wrongly reported as a
TOST p); the difference test uses the identical construction and needs the identical caveat.

**Resolution limit, with the tail counts.** The smallest non-zero two-sided p attainable at
B = 2,000 is 2/2000 = **0.001**, and a zero tail count returns exactly 0.0:

| Contrast | Lower tail | Upper tail | Reported p |
|---|---|---|---|
| DNABERT-2 − k-mer (primary) | 1747 | 253 | 0.253 |
| DNABERT-S − k-mer | 1420 | 580 | 0.580 |
| **HyenaDNA − k-mer** | 1999 | **1** | **0.001 — at the floor** |
| NTv2 − k-mer | 1832 | 168 | 0.168 |
| best-of-4 − k-mer | 25 | 1975 | 0.025 |

HyenaDNA's p = 0.001 rests on a **single** bootstrap draw out of 2,000. It is the smallest
non-zero value the estimator can return and must be reported as **≤ 0.001**, with its
Holm-adjusted 0.004 likewise as **≤ 0.004**. The paired-t equivalent is better resolved
(p = 0.0077, Holm 0.031) and supports the same qualitative conclusion, so the finding is
unaffected — only its stated precision is.

The primary p = 0.253 is well away from the floor (253 of 2,000 draws) and needs no caveat
beyond the construction note.

**Status: resolved.** Specification and all tail counts deposited in
`results/multiplicity_equivalence_audit.json`.

---

## 5. Panel accounting

**Concern.** Confirm 31 declared development blocks, 25 analysed species-level blocks, the
separately reported eligible subspecies stratum, and the five eligibility failures. Explain
how the total of 58 relates to that separate stratum.

**All four confirmed**, from `phase2_pilot/block_summary_full.json`,
`phase4_locked_test/test_block_summary.json` and
`phase7_extensions/hard_block_summary.json`:

| Set | Declared | Eligible | Analysed | Eligibility failures |
|---|---|---|---|---|
| Development | 31 | 26 | **25** | `legionella`, `yersinia`, `bordetella`, `borrelia`, `leptospira` |
| Locked test | 22 | 20 | 20 | `actinomyces`, `bartonella` |
| Hard panels | 15 | 13 | 13 | `salmonella_subspecies`, `campylobacter_jejuni_subsp` |
| **Total** | **68** | **59** | **58** | 9 |

- **31 declared development blocks:** confirmed.
- **Five eligibility failures:** confirmed and named above. All five failed the frozen
  eligibility gate (FREEZE.md §6: ≥10 independent target isolates and ≥2 neighbour species
  with ≥5 each) and are retained as descriptive-only, not deleted.
- **The separately reported eligible subspecies stratum** is `bifidobacterium_longum_subsp`.
  It is the only non-species-rank block in the development set, is counted among the 31
  declared and the 26 eligible and in the 12,356 records, and is excluded from the 25 analysed
  and from every pooled estimate.
- **How 58 relates to it:** **58 = 59 eligible − 1.** The single eligible unit not in the 58
  is the subspecies stratum, reported on its own in the development report. Equivalently
  58 = 25 + 20 + 13.

**One inconsistency found, and it is in the framing rather than in any statistic.** The stated
reason for holding `bifidobacterium_longum_subsp` out of pooled estimates is that species-rank
and subspecies-rank panels are separate strata. That rule is not applied to the hard panels:
`hard_block_summary.json` records `rank: "mixed"` for all 15 — the repository carries **no
per-panel rank for this set** — and 7 of the 13 eligible panels are infraspecific by their own
declared names: six at subspecies or pathovar rank (`bacillus_subtilis_subsp`,
`lactobacillus_delbrueckii_subsp`, `lactococcus_lactis_subsp`,
`mycobacterium_abscessus_subsp`, `streptococcus_dysgalactiae_subsp`,
`xanthomonas_campestris_pv`) plus `clostridium_botulinum_groups`, an intraspecific
physiological-group comparison. "58 analysed groups" is therefore not a rank-homogeneous
population.

No reported statistic is affected: the hard panels never enter the confirmatory estimate and
carry their own interval. What must change is the abstract and Methods framing of 58 as a
single population — and `extension_A_results.md` §5's "the equivalence now holds across 58
independent bacterial groups", which pools three heterogeneous strata into one claim.
`extension_A_results.md` §6 makes the opposite point itself: difficulty concentrates in
subspecies comparisons, so rank is exactly the variable that should not be pooled over
silently.

**Status: verified in full; one framing correction required.**

---

## 6. Difficult-panel interval

**Concern.** The source reports the lower limit as −1.28 in one place and −1.27 in another.
The revision carries "approximately −1.3 to +0.7" pending confirmation of exact output.

**Exact frozen output**, `phase7_extensions/extension_A_results.json`:

```
mean : -0.0019091379940167842
ci   : [-0.012748114717514806, +0.00660273317171494]
n    : 13
```

In percentage points to two decimals: **−0.19 pp, 95% CI [−1.27, +0.66]**.

**−1.27 is correct. −1.28 is a rounding error with a traceable origin.** −0.0127481… rounds
to −0.0127 at four decimals, but `extension_A_results.md` §3 records it as **−0.0128**, and
that value was carried into the manuscript Results (line 119). The abstract (line 21) has
−1.27 and is correct. The two manuscript figures therefore disagree because of one
transcription step, not because of two analyses.

**Required.** Use **−1.27** in both places; correct `extension_A_results.md` §3 to
[−0.0127, +0.0066]; delete the "approximately −1.3 to +0.7" placeholder — the exact output is
confirmed and no longer needs hedging.

**Status: resolved.**

---

## 7. Noise sensitivity estimate

**Concern.** The narrative sensitivity estimate starts at +0.99 pp, whereas the 1%
substitution row of Table 5 gives +0.85 pp. Identify the condition, subset and estimator
behind +0.99 before retaining that comparison.

**Finding — same condition, same subset, two different estimators.** Both come from the
single `robustness_perturbation.py` run (`phase7_extensions/robustness_results.json`), both
on **all 20 locked-test panels**, both under the **1% substitution** condition. They differ
only in estimand:

| Estimator | Definition | Value |
|---|---|---|
| **Level** (Table 5) | mean over panels of (DNABERT-2 − k-mer) *under 1% substitution* | **+0.853 pp** |
| **Shift** (narrative) | mean over panels of the *per-panel change* in (DNABERT-2 − k-mer) from clean to 1% substitution | **+0.990 pp** |

The two are related exactly by the clean-condition baseline:
+0.990 = +0.853 − (−0.137), where −0.137 pp is the clean-condition level difference. The shift
estimator is the differences-in-differences quantity printed by the script's own "delta shift"
block (`robustness_perturbation.py`, final summary table); the level estimator is its
`delta_mean`. There is no inconsistency once each is labelled.

**Both collapse when one panel is removed**, which is the point the narrative is making:

| Condition | Level, all 20 | Level, excl. *francisella* | Shift, all 20 | Shift, excl. *francisella* | Panels +/−/tie |
|---|---|---|---|---|---|
| 1% substitution | +0.853 | **−0.204** | +0.990 | **−0.060** | 9 / 5 / 6 |
| copy variation | +0.369 | −0.095 | +0.506 | +0.049 | 9 / 5 / 6 |
| 1% error + trimming | +0.343 | −0.398 | +0.480 | −0.254 | 6 / 10 / 4 |

*francisella* is the smallest panel in the set (63 records, 55 isolates) and shifts by +0.21
under 1% error. Under **either** estimator the apparent noise-robustness advantage is
attributable to it, and under both the sign reverses on its removal. The conclusion in
`extension_B_results.md` Finding 3 stands; only the estimator labelling needs to change.

**Required.** Keep the comparison, and label it: state that Table 5 reports the level
difference at each condition and that the leave-one-panel-out narrative uses the paired shift
from clean, give both leave-out values (−0.20 pp level, −0.06 pp shift), and note that only
six of 20 panels moved at all under the 1% condition. Either estimator may lead, but the
manuscript must not place +0.99 and +0.85 in adjacent sentences without saying they are
different quantities.

**Status: resolved.**

---

## 8. Ambiguity ceiling

**Concern.** Prefer compatible aggregation for accuracy and ceiling. The existing full-block
ceiling is not a strict upper bound on mean fold-level accuracy. Confirm Figure 7 reflects
this limitation.

**The incompatibility, stated precisely.** The reported accuracy is the mean over five outer
folds of per-fold balanced accuracy, each computed on a held-out fold with its own class
composition. The plotted ceiling is `build_blocks_full.ceilings` evaluated once on the **whole
block** at the modal chosen region. A single fold can exceed a whole-population bound, so
their mean can too. The two quantities are not aggregated over the same thing, and no
relabelling fixes that.

**Measured exceedance.** Across both sets plotted in Figure 7 — 20 locked-test and 13 hard
panels, two arms each, 66 cells — **3 cells exceed the full-block ceiling**:

| Set | Block | Arm | Accuracy | Full-block ceiling | Excess |
|---|---|---|---|---|---|
| locked test | `veillonella` | k-mer | 0.8524 | 0.8494 | +0.302 pp |
| locked test | `veillonella` | DNABERT-2 | 0.8524 | 0.8494 | +0.302 pp |
| hard panels | `bradyrhizobium` | k-mer | 0.9176 | 0.9157 | +0.187 pp |

**Figure 7 does not currently reflect this.** `scripts/make_figures.py:413–451` plots the
full-block ceiling on x against achieved accuracy on y, draws the y = x line, and annotates
it *"ceiling (no method can exceed this)"*. The figure plots the k-mer series only (33
points), **2 of which lie above that line**. The annotation is contradicted by the data it
sits on.

**The compatible aggregation, computed.** Applying the identical ceiling formula to each
held-out fold at that fold's own chosen region, then averaging over the five folds — matching
the accuracy's aggregation exactly:

| Quantity | Value |
|---|---|
| Cells above the **full-block** ceiling (66 cells) | **3** |
| Cells above the **fold-matched** ceiling (66 cells) | **0** |
| Mean ceiling, full-block / fold-matched | 0.9251 / **0.9252** |
| Per-block difference (fold-matched − full-block) | mean **+0.001 pp**, range −1.84 to +1.78 pp |

The fold-matched ceiling is a strict upper bound on every reported cell, and it changes the
mean by one thousandth of a percentage point. For `veillonella` it is 0.8524 — exactly the
reported accuracy, so that cell sits on the bound rather than above it. The substantive
message of Figure 7 is unchanged.

One honest caveat to carry: like the full-block ceiling, the fold-matched ceiling is computed
on held-out composition and so bounds what any sequence-based rule *could* achieve on that
fold. It is not a bound achievable by a rule fitted on training data alone. Both versions
share this character; the fold-matched one merely bounds the right quantity.

**Applied.** Figure 7 has been rebuilt on the fold-matched ceiling, the annotation now reads
"fold-matched ambiguity ceiling (attainable maximum)" and is positioned above the diagonal
where no data can fall, and both arms are plotted. The aggregation is stated in the legend.
Per-block values for both versions are deposited in `results/fold_matched_ceiling_audit.csv`
(33 blocks × both ceilings × both arms), generated by
`scripts/audit_fold_matched_ceiling.py`. The **Methods** edit stating the aggregation is still
required — see §13 item 10.

**Status: resolved; figure rebuilt, Methods edit outstanding.**

---

## 9. Orientation and grouping

**Concern.** The reported exclusion sensitivity is retained, but the pipeline still needs an
orientation fix. Additional exact-sequence-disjoint or sequence-cluster-disjoint analyses
remain recommended. They were not performed for this revision.

**Both halves are accepted. This item is NOT resolved.**

**Orientation — defect real, disclosed, unrepaired.** `audit_strand_sensitivity.py` found
22 of 24,009 records (0.09%) deposited reverse-complemented. A reverse-complemented 16S still
aligns globally to the forward *E. coli* J01695 reference with contiguous query spans that
fall inside the 0.75–1.30× retention window, so these records are **retained and silently
mis-sliced** rather than rejected — the opposite of what the earlier note asserted. All 22
were retained for 8 or 9 of the 9 regions. The leave-out sensitivity (re-executing
`score_test.score_block` on affected locked-test blocks with those isolates removed, reusing
the existing embedding cache, after confirming the re-execution reproduced every frozen
per-block value exactly) moves the primary estimate from −0.2108 to −0.2074 pp, the 95% CI
from [−0.5949, +0.1625] to [−0.5916, +0.1697], p from 0.253 to 0.257, with the verdict
unchanged and a largest single-block shift of 0.067 pp. The defect is immaterial to every
reported result. **The fix has not been applied**, deliberately: the locked test is scored
once under a frozen plan, so repairing the extractor now would replace a frozen confirmatory
result with an unfrozen one. It must remain disclosed in Methods and Limitations as an
unrepaired pipeline defect, with the orientation check named as required work for any
follow-on analysis.

**Grouping — exposure now measured, and it is large.** Cross-validation is grouped on isolate
(`StratifiedKFold` over an isolate-level table, which is grouped-by-isolate by construction
since each isolate contributes one modal slice per region; this is the implementation of
FREEZE.md §3's `StratifiedGroupKFold` and should be noted as an equivalence, not a
deviation). Isolate grouping does **not** make the folds sequence-disjoint. Measured at each
locked-test block's modal chosen region:

| Quantity | Value |
|---|---|
| Isolate groups sharing their exact slice with ≥1 other isolate group | **5,552 / 5,977 (92.9%)** |
| Isolate groups in a haplotype carrying **both** classes | **1,348 / 5,977 (22.6%)** |
| Held-out isolates whose exact slice is present in the training folds | **5,431 / 5,902 (92.0%)** |

Per-block the last figure ranges from 71.2% (`paenibacillus`) to 99.0% (`veillonella`).

This is substantially biology rather than a data-handling error — identical 16S sequences
across distinct strains are expected, and the cross-class component is already charged as
irreducible error by the ambiguity ceiling (`err = Σ min(a, b)`), which is why the ceiling is
reported beside every cell. But the consequence for interpretation is sharp and currently
unstated: **the reported accuracies measure performance within a known haplotype population,
not generalisation to unseen haplotypes.** With 92% of held-out isolates sharing a string with
training data, an exact-sequence-disjoint evaluation would be scored on roughly 8% of the
current units and would answer a different question.

**Not performed.** Neither exact-sequence-disjoint nor sequence-cluster-disjoint
cross-validation was run, here or previously. They should not be claimed, approximated, or
described as supported by the ceiling. The honest position is: the exposure is now quantified
(table above, reproducible from `phase4_locked_test/test_blocks.json` plus the frozen fold
definition), the limitation is stated with its direction, and the disjoint analyses are named
as required future work.

**Status: NOT RESOLVED — two items remain open by design, both must be stated as limitations.**

---

## 10. Software and patched models

**Concern.** Verify all software versions and model revisions against saved environment files.
Document the compatibility patches and any validation. No external version or
numerical-equivalence verification was performed.

**First finding: there were no saved environment files.** The repository contained only a
hand-maintained `requirements.txt`, and the versions in `Scientific_Confirmations.md` §9 were
queried from *live* conda environments — which is not a saved artefact and cannot be checked
by a reader. Four environment files have been exported during this pass and deposited:

```
supplementary/S11_environment_glm-benchmark.yml   (79 lines, conda env export --no-builds)
supplementary/S11_environment_glm-ntv2.yml        (55 lines)
supplementary/S11_piplock_glm-benchmark.txt       (56 packages, pip freeze)
supplementary/S11_piplock_glm-ntv2.txt            (32 packages)
```

**Versions verified against those files and against `requirements.txt` — all match exactly:**

| Environment | Verified |
|---|---|
| `glm-benchmark` | Python 3.11.18, NumPy 2.4.6, SciPy 1.17.1, scikit-learn 1.9.1, PyTorch 2.15.0, Transformers 5.22.1, Biopython 1.88, joblib 1.6.0 |
| `glm-ntv2` | Python 3.11.18, PyTorch 2.15.0, Transformers 4.34.1 |

**Model identity verified against FREEZE.md §1 — all four match exactly:**

| Model | Declared | Verified this pass |
|---|---|---|
| DNABERT-2 | sha256 `7ff39ec77a48…`, 468,354,983 B | `7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba`, 468,354,983 B ✓ |
| DNABERT-S | sha256 `f3cfc3d05418…`, 468,320,053 B | `f3cfc3d0541859df64759e758cc4bc40fe1d5760799a7d098a9f9f09501c1e5f`, 468,320,053 B ✓ |
| HyenaDNA | hub rev `7ebf71773d22c0ede2cc55cb2be15ee8c289e1ce` | cached snapshot `7ebf71773d22c0ede2cc55cb2be15ee8c289e1ce` ✓ |
| NTv2 | hub rev `06615c1660c892fc199840c18123f8385b3542a8` | cached snapshot `06615c1660c892fc199840c18123f8385b3542a8` ✓ |

**Patches — documented, and still not numerically validated.** Four patches were applied to
the two MosaicBERT repositories (DNABERT-2, DNABERT-S) for CPU and `transformers`-5.x
compatibility: a config fallback, a meta-device initialisation ordering fix, a CUDA-only
kernel branch, and removal of a `triton` dependency that has no macOS/arm64 wheel. All four
are listed in `dnabert2_patched/README_PATCHES.md`. They are *asserted* to make no numerical
change, but **no numerical-equivalence test against an unpatched reference was run, because
the unpatched code does not execute in this environment at all.** This must be stated
explicitly in Methods rather than implied by the phrase "compatibility-only". Nucleotide
Transformer v2 required **zero** patches — it runs unmodified under the pinned
`transformers==4.34.1` in `glm-ntv2`, with a determinism-check cosine similarity of exactly
1.0 — which is the right contrast to draw: where patching would have meant reimplementing
framework internals, the dependency was pinned instead.

**Not performed:** no verification against external sources (PyPI release records, Hugging
Face model cards, upstream commit history). All checks above are internal consistency checks
against this machine's environments and caches, and should be described as such.

**Status: resolved, with the unvalidated-patch statement required in Methods.**

---

## 11. Reproducibility and figures

**Concern.** Add the repository URL, release DOI, commit identifier, licences, manifests,
predictions, test fixtures and supporting files. Insert the actual figures and verify every
legend, axis and numerical annotation against final outputs. Computational cost remains a
useful missing addition.

**Verified present:**

| Item | Status |
|---|---|
| Repository URL | `https://github.com/barah123/Readout-Choice-in-16S-Genomic-Model-Evaluation` — in `CITATION.cff` and `.zenodo.json`; `origin` configured, `main` tracks `origin/main` |
| Commit identifier | `46d33a9` (HEAD); 37 commits; 10 phase tags (`phase1-audit-complete` … `phase7-extensions`, `manuscript-package`) |
| Licence | MIT — `LICENSE`, `CITATION.cff:37`, `.zenodo.json:6` |
| Test fixtures | `supplementary/S9_synthetic_validation_tests.py` — known-ground-truth tests of alignment/slicing, the oracle-ceiling formula, grouped CV and the permutation null |
| Environment files | `supplementary/S11_*` — deposited this pass (§10) |
| Frozen plans and panel definitions | `S1`–`S6`, `S10` |
| Figures as files | 7 figures × PDF + PNG, plus `figures/plos_tiff/Fig1–7.tiff` |

**Gaps that remain open:**

1. **Release DOI — AUTHOR ACTION.** `.zenodo.json` is prepared but neither it nor
   `CITATION.cff` carries a `doi` field; Zenodo mints one on release. The manuscript cannot be submitted with AQ02/AQ20 open.
2. **Manifests cover one set of three.** `S7_data_manifest_hashes.json` and
   `S8_accession_manifest.tsv` (6,933 rows) cover the **locked test only**. There is no
   accession manifest or hash record for the development set (12,356 records) or the hard
   panels (4,720 records). The Methods and the Supporting Information list imply manifests for
   all analysed data.
3. **Unit-level predictions do not exist.** `Query_Resolution.md` AQ16 states "unit-level
   predictions in `S3 Data`". **No such file exists anywhere in the repository**, and the
   Supporting Information list promises "immutable predictions". Either generate the file from
   the frozen scores or withdraw the claim; it is currently unsupported, and
   `Query_Resolution.md` AQ16 must be corrected either way.
4. **Figures are not inserted.** Neither `Manuscript_PLOS_CompBiol.docx` nor
   `Manuscript_Draft_Pretrained_DNA_Models_16S.docx` contains a single embedded image
   (`word/media` is empty in both). The figures exist as files but the manuscript does not
   display them. Each was therefore audited standalone against its source data and its
   rendered PNG — see §11a — but none has been checked in page context, at final column
   width, beside its caption.
5. **Two table cells render as a literal `;`.** Table 2's `p (Holm)` cell for the primary row
   and Table 4's `R0 pair AUC` cell for k-mer both print `;` where a not-applicable marker
   belongs — a `build_plos_manuscript.py` artifact. Both are cells where the value is
   legitimately absent (the primary is not Holm-adjusted; k-mer has no R0), so the fix is a
   proper em-dash with a note, not a number.
6. **Figure 7's annotation is contradicted by its own data** (§8) — the one legend error found
   by checking annotations against frozen outputs.

**Computational cost — the data already exist and should simply be added.** `FREEZE.md` §10
records measured development cost: fetch ≈ 65 min; slicing 12,356 records ≈ 25 min on 8 cores;
embedding 7,878 unique slices in 45 s (DNABERT-2), 62 s (DNABERT-S), 474 s (HyenaDNA), 187 s
(NTv2); readouts, 277 cells ≈ 1 min on 6 workers — all CPU, MacBook Air M4, 16 GB, against a
pre-committed locked-test budget of ≤ 3 h wall-clock. This is directly relevant to the
manuscript's own argument that a method requiring ~1000× the compute of k-mer counting needs
more than 1 pp to justify itself, and it belongs in the Results rather than only in the frozen
plan.

**Status: partially resolved.** Items 1–4 above are blockers.

### 11a. Figure-by-figure audit

Each figure was checked by re-deriving its plotted values from the frozen files
(`scripts/make_figures.py` read line by line) and by inspecting the rendered PNG.

| Fig | Numbers correct? | Defect found | Severity | Now |
|---|---|---|---|---|
| 1 | panel a yes; **panel b plotted the wrong estimand** | double-oracle development row, no interval, invisible zeros | **high** | **rebuilt** |
| 2 | yes | panel d axis named the aggregation, not the quantity | low | **relabelled** |
| 3 | yes | palette collision, unlabelled control line, whole-block ceiling | **medium** | **rebuilt** |
| 4 | yes | no pooled estimate or interval | low | **rebuilt** |
| 5 | yes (matches Table 5 exactly) | no intervals anywhere | low | **rebuilt** |
| 6 | yes | regret only; truncated-baseline bars | low | **rebuilt** |
| 7 | **annotation contradicted by its own data** | wrong ceiling aggregation; model arm computed then discarded | **high** | **rebuilt** |

All seven were changed and regenerated as 300 dpi RGB LZW TIFF.

**An eighth defect surfaced while verifying them: the script's own PLOS geometry check was
measuring the wrong thing.** `save()` validated `fig.get_size_inches()`, the *declared* size,
but `savefig` runs with `bbox_inches="tight"`, so the saved raster is larger than the
declaration. Fig 1 was declared 7.20 in and reported `width_ok=True`, while the TIFF actually
on disk was **7.61 in — over the PLOS 7.5 in maximum**, and would have been caught at
production, not here. The check now measures the saved raster, and Fig 1 was narrowed to fit.
Measured on the files themselves:

| | Fig1 | Fig2 | Fig3 | Fig4 | Fig5 | Fig6 | Fig7 |
|---|---|---|---|---|---|---|---|
| width (in) | 7.35 | 7.13 | 4.93 | 7.12 | 5.53 | 7.13 | 5.13 |
| height (in) | 4.53 | 6.13 | 3.53 | 3.13 | 5.13 | 5.32 | 4.33 |

All within 2.63–7.5 in wide and ≤ 8.75 in tall, all 300 dpi, all LZW, all under 0.4 MB.
What follows is what was found and what was done.

**Figure 1 — panel b plots a different contrast from the one it is compared against.**
The locked-test and hard-panel rows are the pre-named DNABERT-2 − C_kmer contrast
(−0.21 and −0.19 pp). The development row is hardcoded as `0.14` and is the stored
`model − classical` quantity from `block_level_results.json`: a **double oracle**, taking the
best of four models per block (DNABERT-2 in 20/25, DNABERT-S 3, NTv2 1, HyenaDNA 1) *and* the
best of two classical comparators per block (`C_kmer` 17/25, `C_align` 8/25). Three rows on a
shared axis under one heading, one of which is selected twice over.

The matched contrast, recomputed for the same 25 development blocks at each block's selected
region, is **+0.12 pp, 95% CI [−0.23, +0.48]** (10 blocks favour the model, 7 the k-mer arm,
8 tie). The plotted point is therefore nearly right by coincidence — 0.02 pp from the matched
value — but it is the wrong quantity, and it is the one positive row in the panel.

Also in Figure 1:

- **The development row has no interval** while the other two do, so the least reliable
  estimate reads as the most precise. The CI above should be drawn.
- **The ±1 pt shaded band is unexplained.** Its on-plot definition was removed in commit
  `4cc712a` on the stated basis that the captions would carry it; they do not (see below).
- **Seven of 20 blocks have a difference of exactly zero and render as no bar at all**
  (`veillonella`, `rothia_micrococcaceae`, `nocardia`, `leuconostoc`, `fusobacterium`,
  `flavobacteriaceae`, `cutibacterium`). A reader cannot distinguish an exact zero from a
  missing block. (`comamonadaceae` at −0.002 and `proteus_morganellaceae` at +0.008 pp are
  near-zero but not exact, and keep hairline bars.)
- Panel a says "percentage points", panel b says "pts", for the same unit.

A related pre-registration note, not an error: `FREEZE.md` §5 sized the locked test using
SD = 0.01138, which is the SD of this same double-oracle development delta, not of the primary
contrast. The direction is favourable — the achieved locked-test SD is 0.00889, so the study
was planned more conservatively than it needed to be — but the basis should be stated.

**Figure 2 — numbers verified; panel d relabelled on review.** Panel d recomputes to readout
(R1→R2) 2.89 pts, model choice 0.73 pts, layer depth 0.53 pts, ratio 3.97×, matching the
manuscript's "+2.8 to +3.0 points against a 0.73-point model spread, roughly fourfold". R0 is
on its own axis and labelled as a pair AUC; the earlier invalid AUC-to-accuracy subtraction is
gone. `paired_ci()` is computed in panel b and not drawn, but it is **not** dead code: it is
the source of the PLOS Fig 2 legend's quoted paired effect, which recomputes to +2.83 pts,
95% CI [+1.37, +4.57] — exactly as the legend states. It is kept, with a comment saying so.

The panel d axis read "spread (balanced-accuracy pts)", which named the aggregation rather
than the quantity and invited the reader to treat three differently-constructed numbers as one
kind of thing. It now reads **"difference in balanced accuracy (percentage points)"**, which is
literally true of all three bars — R2 minus R1, best minus worst model, best minus worst depth
— and the legend names which difference each bar is.

**The layer bar is now defined precisely in Methods**, because its construction is the least
obvious of the three. It aggregates in this order: for each model and depth, balanced accuracy
is averaged across the 13 hard panels, giving five depth means; the range of those five is that
model's depth effect; the two model values are averaged, giving 0.53 pp. The Methods text also
now records the five depths (0.25/0.50/0.75/0.90/1.00 of encoder depth — HyenaDNA layers
2, 4, 7, 8, 9 of 10; NTv2 layers 7, 14, 22, 26, 29 of 30), the correctness gate that ran
before any result was computed (extracted final layer must be `allclose` to the standard
forward output at 1e-6; both passed, the run aborts otherwise), why DNABERT-2 and DNABERT-S
were not swept (internal unpadding returns `(total_tokens, hidden)`), and two directional
caveats: averaging panels before taking the range understates per-panel depth effects, and
depth was chosen on the same panels used for the comparison, so any best-depth figure is
optimistic rather than a bound.

**Figure 3 — two defects, both cosmetic in origin and substantive in effect.**

1. **Palette collision.** `C["ntv2"]` and `C["exact"]` are both `#CC79A7`. NTv2 and
   Exact match are drawn in the identical pink. In a figure whose docstring claims a
   colourblind-safe palette, two series are indistinguishable by colour.
2. **The dark-red vertical line at 50.1% carries no label and the figure has no legend.**
   It is the label-permutation floor (mean 50.12% over 20 blocks × 3 arms;
   `label_permutation_control.json`, gap 46.0 pts — the number is sound). Its on-plot text was
   removed in `4cc712a`. As rendered, a reader sees an unexplained red line.

The dashed ceiling line has now been **recalculated on the fold-matched aggregation**, for
consistency with Figure 7 — two figures drawing the same quantity under two different
definitions would be worse than either. The verified value is unchanged at the displayed
precision: the locked-test mean is **97.0071% fold-matched against 97.0380% whole-block**, both
97.0% to one decimal. The line is now labelled "fold-matched ambiguity ceiling, 97.0%" in the
legend and the caption says so.

The figure's old caption asked the author to "specify whether error bars describe between-block
standard deviation or uncertainty in the mean" — **there are no error bars**, though the
between-block SDs are in Table 3. That remains the case; the two reference lines carry the
context instead.

**Figure 4 — numbers verified.** The "10 of 20 groups saturated" annotation is exactly right
(10 blocks have k-mer error below 1 pt). Panel b shows the 13 per-panel differences but no
pooled estimate or interval, which the caption asks for; add −0.19 pp [−1.27, +0.66] and mark
that it crosses the margin.

**Figure 5 — numbers verified against Table 5 to the decimal**, all six conditions, both arms.
The only gap is that the CIs exist in `robustness_results.json` and are not drawn; the caption
asks for them.

**Figure 6 — numbers verified.** Locked-test regret: k-mer 0.654, DNABERT-2 0.566,
DNABERT-S 0.673, NTv2 1.196, HyenaDNA 1.612 pts; hard panels 1.946 / 2.017 / 1.887 / 1.609 /
1.372 pts. The figure plots regret only, not selected-region accuracy, and does not show the
hit rate — defensible, since §2 of `Scientific_Confirmations.md` establishes that regret is
the better-behaved metric under ties, but the caption should say that is why.

**Figure 7 — the §8 problem, plus a production defect.** The annotation "ceiling (no method
can exceed this)" is **overprinted on the data points** in the top-right cluster and is barely
legible as rendered, and it is contradicted by two of the 33 plotted points. Additionally,
`fig7_ceiling()` computes `v["dnabert2:R2"]` into `pts` and **never plots it** — the figure
contains no model data at all, only the k-mer series. In a paper comparing models to k-mer
counting, Figure 7 shows only the baseline. Either plot both arms or retitle the figure to
what it actually shows. The axis label "information ceiling" also differs from the
manuscript's term "empirical ambiguity ceiling"; pick one.

### 11b. Figure legends — correction to an earlier finding, and what was wrong

An earlier pass of this document stated that "the figure captions are not captions". **That
was half right and is corrected here.**

- **`build_plos_manuscript.py` — the submission artifact — carries real legends.** They are
  prose, they define the ±1 pt shaded band in Fig 1, they identify the red permutation line in
  Fig 3 at 50.1%, and their quoted statistics check out (Fig 2's "+2.83 points, 95% CI +1.37
  to +4.57" recomputes exactly). Commit `4cc712a`'s claim that the captions carry the
  information removed from the plots **was** honoured here.
- **`Manuscript_Draft_Pretrained_DNA_Models_16S.md` — the author-notes draft — did not.** Its
  seven caption blocks were author-facing specifications written in the imperative
  ("Display…", "Remove ratio labels…", "Resolve the Veillonella discrepancy before describing
  every marker as a strict test-set bound"). A reader of that file had no legends at all, and
  the two artifacts disagreed about what the figures showed.

Both are now rewritten against the rebuilt figures, and they agree. The draft's placeholder
blocks are replaced with the submission legends, and the PLOS legends are updated to describe
what each figure now shows: the exact-zero markers and the matched development contrast
(Fig 1), the truncated-baseline caveat (Fig 3), the pooled hard-panel interval and its
relation to the margin (Fig 4), the paired-difference panel (Fig 5), why regret rather than
hit rate is plotted (Fig 6), and the fold-matched ceiling with the exceedance count it
resolves (Fig 7).

Six of the seven changes the draft's old specifications requested are now implemented:
label development as exploratory (Fig 1), add the hard-panel interval (Fig 4), add paired
error bars (Fig 5), show selected-region accuracy alongside regret (Fig 6), state where
ceiling and accuracy use different distributions (Fig 7), and resolve the *Veillonella*
discrepancy (Fig 7). The seventh — splitting the permutation null into its own panel
(Fig 3) — was **not** done: the floor is a single scalar bounding all eight methods, and a
labelled reference line carries it without spending a panel. That is a judgement call, not an
oversight, and it is recorded here as one.

---

## 12. What this pass changed in the repository

| Path | Change |
|---|---|
| `scripts/audit_multiplicity_and_equivalence.py` | **new** — reproduces the bootstrap difference test verbatim, adds component TOST outputs and both Holm families, and records the best-of-four selection accounting |
| `results/multiplicity_equivalence_audit.json` | **new** — output of the above (§1–§4) |
| `results/fold_matched_ceiling_audit.csv` | **new** — per-block full-block and fold-matched ceilings for all 33 plotted blocks (§8) |
| `supplementary/S11_environment_glm-benchmark.yml` | **new** — conda export (§10) |
| `supplementary/S11_environment_glm-ntv2.yml` | **new** — conda export (§10) |
| `supplementary/S11_piplock_glm-benchmark.txt` | **new** — pip freeze (§10) |
| `supplementary/S11_piplock_glm-ntv2.txt` | **new** — pip freeze (§10) |
| `scripts/audit_isolate_denominators.py` | **new** — block-by-block, fold-by-fold denominator reconciliation (§14.3) |
| `results/isolate_denominator_audit.csv` | **new** — output of the above |
| `scripts/analysis_sequence_disjoint.py` | **new** — exact-sequence-disjoint evaluation (§15) |
| `results/sequence_disjoint_results.json` | **new** — per-block values, split audit and summary statistics |
| `scripts/analysis_similarity_disjoint.py` | **new** — 99%-identity single-linkage disjoint folds + cross-fold audit (§16.1) |
| `results/similarity_disjoint_results.json` | **new** — output of the above |
| `scripts/analysis_tuned_readout.py` | **new** — per-representation `C` search, both splitters (§16.2) |
| `results/tuned_readout_results.json` | **new** — output of the above |
| `scripts/analysis_divergence_stratified.py` | **new** — per-neighbour decomposition + divergence bins (§17.1) |
| `results/divergence_stratified_results.json` | **new** — 63 target–neighbour pairs |
| `scripts/analysis_r1_locked_test.py` | **new** — R1 → R2 on the locked test, with fidelity check (§17.2) |
| `results/r1_locked_test_results.json` | **new** — output of the above |
| `scripts/audit_fold_matched_ceiling.py` | **new** — generates the fold-matched ceilings Figure 7 now plots (§8) |
| `scripts/make_figures.py` | **modified** — Figures 1, 3, 4, 5, 6, 7 rebuilt (§11a); palette collision fixed |
| `figures/fig{1..7}.png`, `.pdf`, `plos_tiff/Fig{1..7}.tiff` | **regenerated** — all seven, PLOS spec verified |
| `Manuscript_and_Author_Notes_Package/build_plos_manuscript.py` | **modified** — seven figure legends rewritten (§11b) |
| `Manuscript_and_Author_Notes_Package/Manuscript_PLOS_CompBiol.docx` | **rebuilt** from the above |
| `Manuscript_and_Author_Notes_Package/Manuscript_Draft_Pretrained_DNA_Models_16S.md` | **modified** — submission legends; ceiling-aggregation and layer-summary Methods |
| `build_plos_manuscript.py` (Methods) | **modified** — ceiling aggregation, intermediate-layer sweep and Fig 2d layer summary now defined precisely |

No frozen artefact was modified, no model was re-run and no analysis was re-scored. The
figures are redrawn from the same frozen values they always plotted, with two exceptions that
are themselves corrections: Figure 1b now plots the matched contrast instead of the
double-oracle one, and Figure 7 plots the fold-matched ceiling instead of the whole-block one.
No Results or Methods prose has been edited — those corrections are still listed in §13.

---

## 13. Consolidated corrections required before submission

**Status 2026-11-14: items 1–28 are applied. Items 29–34 are author-only and remain open.**
All manuscript edits were made in `build_plos_manuscript.py` (the source; the `.docx` is a
build artifact) and mirrored in `Manuscript_Draft_Pretrained_DNA_Models_16S.md`. The docx was
rebuilt and verified; a stale-value sweep returns zero hits for every corrected number.

**Numerical corrections — applied**

1. **TOST p-value.** `p = 0.0001` replaced by **4.1 × 10⁻⁴** in all three places, described as
   a paired-t TOST with the bootstrap interval as *corroboration, not independent
   replication*, since both analyse the same block-level differences. Component one-sided
   values given.
2. **Hard-panel lower limit** −1.28 → **−1.27** in the manuscript and in
   `extension_A_results.md`; the "approximately −1.3 to +0.7" placeholder is gone.
3. **HyenaDNA p-values** reported as bounds, **≤ 0.001** raw and **≤ 0.004** Holm, with the
   better-resolved paired-t (0.0077 / Holm 0.031) given alongside.
4. **CI rounding** harmonised: `.zenodo.json` −0.60 → **−0.59**, matching the frozen
   [−0.5949, +0.1625].

**Methods — applied**

5. Multiplicity described as **hierarchical**: primary unadjusted, four secondary Holm-adjusted
   among themselves, with an explicit statement that one correction was *not* applied across
   all five.
6. **Difference and equivalence families specified separately**, with their distinct nulls and
   decision rules, and the Holm-adjusted equivalence p-values labelled post hoc.
7. **Bootstrap difference test fully specified** — estimator, block as unit, B = 2,000, seeds
   1729 / 1736, percentile interval — and named as an *inversion of the percentile interval on
   an unshifted resampling distribution*, with its 0.001 resolution floor.
8. The **95% containment rule is stated as conservative** relative to the 90% interval that
   TOST at α = 0.05 implies.
9. **Best-of-four defined as an oracle**, with the per-block winner counts, the **+0.657 pp**
   selection gap and the −0.211 pp operational alternative.
10. **Ceiling aggregation** stated, figures on the fold-matched version.
11. **Splitter equivalence** noted (`StratifiedKFold` on an isolate-level table *is* the
    isolate-grouped split) — an equivalence, not a deviation.
12. **Patches not numerically validated** stated plainly; version checks described as internal,
    not against external release records.
13. **Computational cost added** from the frozen plan, tied to the margin rationale.

**Results and Abstract — applied**

14. DNABERT-S and NTv2 equivalence reported as **formal TOST results**.
15. HyenaDNA equivalence **tested and not established** (p_TOST = 0.351) — distinguished from
    an interval that merely happens to be wide.
16. **Equivalence verdict removed** from the best-of-four row; it reads "oracle contrast; no
    equivalence verdict".
17. **Both noise estimators labelled** — Table 5's level difference (+0.85) and the narrative's
    paired shift from clean (+0.99), related by the −0.14 clean baseline, with both leave-out
    values (−0.20 / −0.06) and the 9/5/6 block split.
18. **"58" no longer framed as one population**: three strata, never pooled, 7 of 13 hard
    panels infraspecific; `extension_A_results.md` corrected too. The **6,599 split
    population** is now stated.
19. **Sequence reuse written into Results and Limitations**, not left under Limitations alone —
    92.0% exposure, the disjoint verdicts, the negative degradation test and the divergence
    stratification (§15–§17).

**Figures and tables — applied**

20–24. Figures 1b, 3, 4b, 5, 6, 7 rebuilt; all seven captions rewritten. *(§11a–11b)*
25. **Precision-plan SD provenance stated**: 0.01138 belongs to the double-oracle development
    contrast, the achieved primary SD was 0.00889, so the plan was conservative.
26. **Figures embedded in the draft** for review. The PLOS build deliberately carries legends
    only — PLOS Computational Biology requires figures as separate TIFF/EPS files, which are in
    `figures/plos_tiff/`. Annotations were verified against frozen outputs; verification *in
    page context at final column width* still needs doing once the journal template is applied.
27. **The two `;` cells fixed.** They were an em-dash corruption in the draft only; the source
    CSVs were always correct. Replaced with "— (unadjusted)" and "— (no embedding)".
28. **The development sign flip is now stated**: the matched contrast is +0.12 pp against the
    locked test's −0.21 pp, with a note that the earlier +0.14 was selected twice over.

**Author-only, and blocking — still open**

29. Release DOI (Zenodo). *(§11)*
30. Accession manifests and hash records for the development and hard-panel sets. *(§11)*
31. Unit-level predictions file, or withdraw the `S3 Data` claim and correct
    `Query_Resolution.md` AQ16. *(§11)*
32. Author, affiliation, funding, competing-interest and ethics statements
    (AQ02, AQ20–AQ23); verification of all 13 references against publisher records.
33. Decide whether to rename the `copy_var` condition to *within-isolate sequence variation*.
34. Fix `PYTHONHASHSEED`, or replace the `hash((bid, cond))`-based perturbation seed, for
    exact reproducibility of the robustness run.

**Two caveats raised by this pass's own new analyses**

35. The tuned-readout grid selects its **lower boundary** (`C = 0.01`) in most fold-arm
    decisions. Extend it to 1 × 10⁻⁴ before the tuned comparison is presented as final. *(§16.2)*
36. The similarity-disjoint threshold is **not a fixed mismatch budget** across 31–107 nt
    regions; 6 of 9 evaluable blocks tolerate zero mismatches at 99%. Reported per block, and
    must stay reported. *(§16.1)*

---

## 14. Adjudication of the external review (2026-11-14)

An external review was received arguing that the manuscript measures "inside a lookup table".
The author's assessment accepts the sequence-overlap criticism as the paper's most important
weakness while rejecting several of its statistical arguments. This section records what was
**verified against the repository** on 2026-11-14, separating claims that hold from claims
that do not. It is an adjudication, not an endorsement of either document.

### 14.1 The central criticism is correct and is already this document's §9

The external review and §9 of this document independently reach the same place: 92.0% of
held-out isolates share their exact regional slice with the training folds, so the benchmark
establishes performance **within a known haplotype population**, not generalisation to unseen
variants. That is the binding limitation on interpretation, and the author's position — that
it constrains rather than invalidates a narrowly stated comparison — is the defensible one.
The frozen primary estimand is a statement about this population, and it remains valid as such.

What the external review adds, and what should be adopted: exact-sequence-disjoint evaluation,
a stricter similarity-disjoint evaluation where class support permits, results split by whether
a test sequence appears in training, and per-target-neighbour rather than pooled-neighbour
reporting. The divergence-stratified analysis is the cheapest of these and is genuinely
missing — percent identity is already computed for `C_align`.

Three corrections to how the external review frames it, all of which the author's assessment
gets right:

- **"You are measuring inside a lookup table" overstates what the 90.3 → 96.2 gap shows.**
  That gap does not identify which observations produce the improvement. `C_exact` assigns the
  majority label among identical training sequences and **abstains at 0.5 when the test
  sequence is absent**, with abstentions scored rather than dropped. Exact-match coverage,
  accuracy conditional on a match, and the abstention contribution have to be separated before
  the 6 points can be attributed.
- **Cluster-disjoint folds are not "the ~8% unseen sequences".** Reassigning whole clusters to
  folds creates a new evaluation in which previously repeated sequences become unseen at
  training time. The evaluable population is not the current 8%.
- **"Under a day on the M4" should not be promised.** Feasibility depends on cluster sizes,
  class support after reassignment, nested selection, and how many blocks remain evaluable.

A cluster-disjoint split also needs an operational definition, because clustering at a
threshold does not guarantee every cross-fold pair falls below it. The split must be **audited
directly for cross-fold similarity**, with identity definition, alignment coverage, strand
handling and ambiguous-base treatment all specified — the last two being exactly the pipeline
areas §9 already flags as unrepaired.

### 14.2 Claims checked against the repository

| External review claim | Verdict |
|---|---|
| Table 5 contradicts Methods; standalone 1% row missing | **Not reproducible.** Table 5 currently carries all six conditions including "1% substitution error" at +0.85. Already corrected in an earlier revision. |
| Fig 3 draws a 92.52% ceiling across a locked-test-only chart | **Wrong, and traceable.** "92.52" appears nowhere in any manuscript source. 92.5% is the **pooled 33-block** mean from this document's own §8 audit; the reviewer appears to have read it here and mistaken it for a manuscript claim. Fig 3's line has always been locked-test-only, and is now the fold-matched 97.0%. |
| Counts 6,760 → 5,977 → 5,902 don't reconcile; 75 isolates lost | **The invariant holds; the inference from it does not.** See §14.3. Nothing is lost. |
| Table 7 regret: 0.57 − 0.65 = −0.08, text says −0.09 | **Text is right.** Unrounded means give −0.0879 pp → −0.09. The stored paired contrast is −0.0879 pp, 95% CI [−0.458, +0.259]. −0.08 is the difference of already-rounded values. Add a note rather than changing the number. |
| HyenaDNA gains +0.60 at layer 8, closing most of its deficit | **Misreads Table 6.** +0.60 is layer 8 **versus k-mer on the hard panels**; the final layer is already +0.44 there. The layer-selection gain is **+0.16**, which is what Results states. It also cannot be subtracted from the −0.885 pp locked-test contrast: different datasets, and depth was chosen on the same panels it was evaluated on. |
| Cite Rahimi & Recht instead of Johnson–Lindenstrauss | **Decline.** The arm is a linear Gaussian projection of existing k-mer features. Random Fourier features add a nonlinear map to approximate a kernel — a different construction. JL stays, provided it is not presented as proving sufficiency for classification; report the empirical result plainly and keep the theory modest. |
| Hash-seeded perturbations are irreproducible | **Correct and unrepaired.** Already §13 item 34. |
| Manifests, predictions, DOI, declarations missing | **Correct.** Already §13 items 29–32; the `S3 Data` predictions file still does not exist. |

### 14.3 The denominator chain, reconciled

The author identified the right remedy — a block-by-block, fold-by-fold accounting — and it
now exists (`scripts/audit_isolate_denominators.py` →
`results/isolate_denominator_audit.csv`). The chain closes completely:

| Step | Count |
|---|---|
| Isolate groups, all 22 **declared** locked-test blocks | **6,760** |
| less the 2 blocks failing the eligibility gate (`actinomyces`, `bartonella`) | −161 |
| Isolate groups, 20 **eligible** blocks — **the split population** | **6,599** |
| Held out across the five folds | **6,599** ✓ |
| Groups present at each block's **modal** region | **5,977** (−622 absent at that region) |
| Held out **and** present at that fold's **own** chosen region — the evaluated population | **5,902** (−697) |

The reviewer's invariant is correct and it **holds**: every group is held out exactly once,
6,599 = 6,599. The error is in assuming 5,977 and 5,902 must be equal. They are two different
projections of the same 6,599, differing by which region's availability is applied — the modal
region for one, each fold's own selected region for the other. The "75 missing isolates" is
exactly 697 − 622.

The mechanism is that **12 of 20 blocks had folds select more than one region**
(`cutibacterium` chose V9, V2, V4, V4, V2; `comamonadaceae` chose five distinct regions), and
per-region availability varies across the corpus from 100% at V4 to 48.7% at V9. Decisively,
the direction is not uniform: `fusobacterium` evaluates **132** against a modal-region count of
129, and `streptomyces` **366** against 364. A non-modal region can carry *more* isolates. Loss
cannot produce that.

**What was genuinely missing is the 6,599.** The split population was never stated in the
manuscript or in this document, which is what made the chain look broken. It should be added
to Methods, and the overlap figures in §9 should be labelled with which denominator they use —
92.9% is of 5,977 at the modal region, 92.0% is of 5,902 evaluated across folds.

### 14.4 Where the external review's statistics need correction

- **Saturation does not make equivalence automatic.** The author's correction is right: if
  k-mer sits 0.5 pp below the ceiling, that bounds how much *better* another method can be, not
  how much *worse*. Both interval bounds still carry information, and all four model contrasts
  are in fact negative. Saturation makes the *superiority* question less informative; it does
  not make the two-sided test vacuous.
- **Do not adopt `accuracy / ceiling` as a co-primary outcome.** It does not measure remaining
  headroom and ignores the 0.5 chance floor of balanced accuracy. The defensible descriptives
  are the **ceiling gap** `A_ceiling − A` and, where the denominator is safely above zero,
  **chance-corrected attainment** `(A − 0.5) / (A_ceiling − 0.5)`. Neither supplies an
  equivalence margin on its own, and a relative improvement over `A_ceiling − A_kmer` is
  unstable as k-mer approaches the ceiling. Keep the frozen absolute-scale primary; add ceiling
  gaps and a non-saturated-block restriction as **labelled post hoc** analyses. A newly chosen
  endpoint cannot become retrospectively co-primary.
- **The TOST criticism is half right.** The paired-t and the bootstrap interval are not
  independent corroboration — same block differences, different inferential procedure. "Mutual
  sensitivity analysis" is the accurate phrase and §2 should say so. But the frozen decision
  rule was already **bootstrap CI containment**, so "foreground the bootstrap" largely restates
  the existing hierarchy; and seven exact zeros justify examining robustness, not discarding
  the paired-t. **Do not swap in Wilcoxon as a drop-in**: it targets a different location
  parameter, so it changes the scientific question rather than testing the same one more
  robustly.
- **"Matched classifier" needs qualifying, not retracting.** A fixed `C = 1.0` across ~2,700
  sparse k-mer dimensions, 768, 1,024 and 256 dense ones is a legitimate *fixed-readout*
  comparison and answers a real question. It does not answer how each representation performs
  after comparable tuning. Add a nested, per-representation `C` search on a common grid,
  criterion and budget — and keep the fixed-readout result visible rather than recanting it.
  The "injectivity threshold" argument for sweeping `k` is speculative and should not be used
  as the stated rationale; sweep `k` as a sensitivity analysis on its own terms.
- **Four negative contrasts are not four replications.** Same blocks, same comparator, related
  procedures, and DNABERT-S is built on DNABERT-2. Report the uniform sign descriptively; do
  not convert it into an informal independent-evidence argument.
- **R0 versus R2 is a dissociation, not a readout gain.** The comparable readout effect stays
  R1 → R2. The manuscript already draws this line and should keep it.
- **The label-noise floor is asserted, not established.** A 1–2% mislabelling rate has not been
  measured here, and an assumed label-error rate does not translate directly into a floor on
  the *paired difference between two methods* — both arms see the same labels. The defensible
  programme is: check the conflicting accessions, screen for chimeras by a documented
  procedure, assess agreement against a curated reference treated as supporting evidence rather
  than ground truth, and report how excluding questionable records moves the contrast.

### 14.5 Two wording changes this adjudication forces

1. **"Prespecified" needs moderation**, in the title and abstract as well as Methods. What
   exists is an analysis plan recorded in local commits before locked-test retrieval, with no
   remote and no independent certification — already disclosed in `S10`. The defensible phrase
   is **"an analysis plan recorded before locked-test retrieval"**. This is a narrower claim
   than asserting what readers or PLOS must infer from the word.
2. **The headline claim should not be replaced wholesale.** The external review's proposed
   headline asserts that layer choice moves performance several times more than model identity.
   The repository does not support it: the layer effect is **0.53 pp** against a model spread of
   **0.73 pp** (Fig 2d), so layer moves performance *less*, not several times more. The readout
   effect (2.89 pp) is the one that dominates, and it is development-only. Adopt the external
   review's framing about saturation and sequence reuse; do not adopt its arithmetic.

### 14.6 Revision order

Following the author's ordering, with this document's items mapped onto it:

1. Reconcile denominators (**done**, §14.3) and repair the reproducibility defects
   (§13 items 34, 29–31).
2. ~~Exact-sequence-disjoint evaluation~~ **DONE — see §15.** Equivalence does not survive:
   −0.50 pp [−3.20, +2.07], inconclusive at ±1 pp, SD up 14-fold. Both arms lose ~11 pp.
   Similarity-disjoint (99%) also **DONE — see §16.1**: −1.47 pp, inconclusive, audit passed
   at max cross-fold identity 0.9894.
3. ~~Tuned-readout and per-method region-selection sensitivity analyses.~~ **DONE — §16.2,
   §16.3.** Tuning does not favour the models under disjoint folds and moves the
   isolate-disjoint verdict to inconclusive; per-method region selection leaves the verdict
   unchanged at −0.025 pp. Grid needs extending below 0.01 before the tuned result is final.
4. ~~Target–neighbour and divergence-stratified reporting.~~ **DONE — §17.1.** Difficulty is
   near-binary on shared sequence (13-point swing); the method gap is flat at ≤0.25 pp in
   every divergence bin.
5. Preserve the frozen primary analysis and separate every new analysis from it explicitly.
6. ~~Run R1 on the locked test, or retitle.~~ **DONE — §17.2.** The readout claim is
   confirmed: +1.98 pp [+0.69, +3.61], 2.5× the model spread. The title stands; the quoted
   "roughly fourfold" is development-only and must be labelled or replaced with 2.5×.

---

## 15. Exact-sequence-disjoint evaluation — result (2026-11-14)

**Post hoc. The frozen confirmatory analysis is unchanged and was not recomputed.**
Script: `scripts/analysis_sequence_disjoint.py`. Output:
`results/sequence_disjoint_results.json`. No re-embedding was needed — the locked-test cache
is keyed by slice string, so every vector required already existed.

### 15.1 Design

One variable changes and everything else is held fixed. The region is each block's modal
chosen region from the frozen run, held fixed across folds; the frozen pipeline re-selects a
region per fold, and fixing it removes selection variance so that the **splitter is the only
difference** between the two arms. Same probe, same k-mer featuriser, same embedding cache,
same seed.

- **Control — isolate-disjoint:** `StratifiedKFold(5)` over isolates, the frozen design.
- **Treatment — sequence-disjoint:** `StratifiedGroupKFold(5)` grouped on the exact slice
  string, so no sequence appears in both training and test.

Cross-class haplotypes are kept and assigned whole to one fold rather than dropped, which
preserves the population; the ambiguity ceiling already charges them as irreducible error.

### 15.2 Feasibility — 14 of 20 blocks, and the 6 losses are informative

Haplotype grouping needs enough distinct haplotypes per class to form folds. **14 of 20**
blocks have ≥5 in each class. The six that do not are excluded as not evaluable, with reasons
recorded rather than dropped silently:

| Block | Region | Isolates | Haplotypes (target / neighbour) |
|---|---|---|---|
| `francisella` | V4 | 55 | 3 (2 / 2) |
| `veillonella` | V8 | 98 | 4 (2 / 4) |
| `granulicatella` | V3 | 171 | 9 (2 / 7) |
| `nocardia` | V5 | 405 | 13 (2 / 11) |
| `leuconostoc` | V1 | 409 | 24 (3 / 21) |
| `rothia_micrococcaceae` | V7 | 402 | 37 (4 / 34) |

Every one has **≤4 target haplotypes**. `nocardia` has 405 isolates and 2 target haplotypes;
`leuconostoc` has 409 isolates and 3. These blocks cannot support a disjoint evaluation at any
useful fold count — which is itself a direct measurement of how little independent sequence
variation the panel design actually contains. All 14 evaluable blocks produced 5 clean folds
with no degenerate folds and full test coverage.

### 15.3 The split is verified disjoint, not assumed

Under a correct exact-sequence-disjoint split, `C_exact` must abstain on every test sequence
and score exactly 0.500. Asserted per fold:

| Check | Isolate-disjoint | Sequence-disjoint |
|---|---|---|
| max train/test exact overlap | — | **0.000000** |
| `C_exact` | 90.62% | **exactly 0.500 on all 14 blocks** |

### 15.4 Result

Mean balanced accuracy over the 14 evaluable blocks:

| Arm | Isolate-disjoint | Sequence-disjoint | Change |
|---|---|---|---|
| k-mer counting | 96.36 | 85.42 | **−10.94** |
| Random projection of k-mers | 96.43 | 85.76 | −10.67 |
| DNABERT-2 | 96.16 | 84.92 | **−11.24** |
| DNABERT-S | 96.14 | 86.46 | −9.68 |
| HyenaDNA | 95.77 | 86.13 | −9.64 |
| NTv2 | 95.98 | 82.03 | −13.96 |
| Exact match | 90.62 | 50.00 | −40.62 |

The primary contrast, DNABERT-2:R2 − C_kmer, n = 14:

| Split | Mean (pp) | 95% CI (pp) | SD (pp) | p_TOST | Verdict at ±1 pp |
|---|---|---|---|---|---|
| isolate-disjoint (control) | −0.20 | [−0.39, −0.02] | 0.37 | 1.04 × 10⁻⁶ | **EQUIVALENT** |
| **sequence-disjoint** | **−0.50** | **[−3.20, +2.07]** | **5.16** | **0.361** | **INCONCLUSIVE at this margin** |

**The control validates the design**: −0.20 pp against the frozen −0.21 pp. Per-block values
differ where fixing the region changes the answer, but the mean reproduces.

### 15.5 What this establishes, stated precisely

1. **The equivalence verdict does not survive sequence-disjoint evaluation.** It does not
   reverse — the point estimate stays near zero and the direction is unchanged — but the
   block-level SD rises **14-fold**, from 0.37 pp to 5.16 pp, and the interval becomes roughly
   ±2.6 pp against a ±1 pp margin. At n = 14 the data cannot support an equivalence claim on
   unseen haplotypes. Neither can they support a difference claim. **Inconclusive is the
   honest verdict, and it should be reported as such.**
2. **The scale of the sequence-reuse effect is now measured: ~11 points.** Both primary arms
   lose about 11 pp when identical sequences are removed from training. That is an order of
   magnitude larger than anything separating the methods, and it quantifies the §9 limitation
   that was previously only described.
3. **The upside hypothesis is not supported.** Frozen embeddings do **not** degrade more slowly
   than k-mer counts on unseen variants: DNABERT-2 loses 11.24 pp against k-mer's 10.94 pp,
   difference **−0.30 pp [−2.92, +2.22], p = 0.837**, with the model degrading less in only
   **3 of 14** blocks. The external review hoped this would be a positive novel result. It is
   not. It should be reported as a negative result, because it is the direct test of the
   mechanism the models are supposed to provide.
4. **`C_exact` collapsing 90.62 → 50.00 is the cleanest statement of what the benchmark was
   measuring.** A literal lookup table scored 90.6% under the published design and scores
   chance once it cannot look anything up.

One secondary observation, not a claim: DNABERT-S (−9.68) and HyenaDNA (−9.64) degrade less
than k-mer while NTv2 (−13.96) degrades more. These are not the pre-named arm, the differences
are small against a 5 pp block-level SD, and nothing should be built on them without a
pre-registered test.

### 15.6 Consequences for the manuscript

- The frozen primary analysis stands as a statement about its population and must be labelled
  that way throughout — title, abstract, Results, Discussion — not only under Limitations.
- This analysis belongs in **Results**, labelled post hoc, with the feasibility table. It is
  currently the paper's most generalisable finding.
- The honest framing of the whole paper is now supported by measurement rather than argument:
  *fine-scale 16S discrimination as commonly constructed is dominated by exact-sequence reuse;
  removing it costs ~11 points and leaves no method distinguishable at the resolution this
  design affords.*
- The similarity-disjoint (99%) analysis is now also done — see §16.1. **One prediction made
  here was wrong**: clustering at 99% does not restore any block to evaluability. It only
  merges haplotypes, so the group count falls and the evaluable set is a strict *subset* of
  the exact-disjoint 14.

---

## 16. Similarity-disjoint folds, tuned readout, region selection (2026-11-14)

**All post hoc. The frozen confirmatory analysis is unchanged and was not recomputed.**
Scripts: `analysis_similarity_disjoint.py`, `analysis_tuned_readout.py`. Outputs:
`results/similarity_disjoint_results.json`, `results/tuned_readout_results.json`.

### 16.1 Similarity-disjoint folds at 99% identity

**Operational definition, because the guarantee depends on it.** Identity is the project's
canonical measure — `2 × matches / (len_a + len_b)`, matches over aligned blocks only, from
Biopython global alignment at open −5 / extend −2 / match +2 / mismatch −1 — reused verbatim
from `run_readouts_full.identity_dist`, the same function `C_align` uses. All compared
sequences are slices of the same region, so alignment spans the full length of both. Strand:
slices are forward-oriented as extracted, and the 22 disclosed minus-strand records are **not**
reoriented, so this analysis inherits that defect. Ambiguity: N scores as a mismatch, making
an N-bearing slice *less* similar to its neighbours, so the direction is conservative for
disjointness.

**Clustering is single linkage** — connected components of the graph whose edges are pairs at
identity ≥ 0.99. This is the choice that delivers the guarantee: two sequences in different
components are by construction joined by no edge. Greedy and centroid methods (CD-HIT, UCLUST)
do not ensure every cross-cluster pair falls below the threshold.

**The audit passes.** Maximum identity between any training and any test sequence, computed
per fold across all evaluable blocks: **0.9894 < 0.99**, below threshold on every block.

**A threshold caveat that matters and is reported per block.** A fixed identity threshold is
not a fixed mismatch budget across regions of 31–107 nt. At a 31 nt V1 slice one mismatch
gives 96.8% identity — below 99% — so 99%-clustering collapses to *exact* clustering; at a
107 nt V2 slice one mismatch gives 99.1% and such pairs merge. **6 of the 9 evaluable blocks
tolerate zero mismatches at this threshold**, and `paenibacillus` (V1, 31 nt) has
99%-clustering identical to exact clustering. The effective tolerance is recorded for every
block.

**Feasibility falls, it does not rise.** Merging haplotypes reduces the group count, so only
**9 of 20** blocks remain evaluable, a strict subset of the exact-disjoint 14.

Matched comparison on the 9 blocks evaluable under both designs:

| Split | k-mer | DNABERT-2 | Contrast (pp) | 95% CI (pp) | SD | Verdict at ±1 pp |
|---|---|---|---|---|---|---|
| isolate-disjoint | 97.04 | 96.68 | −0.35 | [−0.58, −0.15] | 0.36 | **EQUIVALENT** |
| exact-sequence-disjoint | 87.21 | 85.63 | −1.59 | [−5.33, +1.89] | 5.85 | INCONCLUSIVE |
| 99%-similarity-disjoint | 87.01 | 85.54 | −1.47 | [−4.41, +0.92] | 4.35 | INCONCLUSIVE |

**The stricter threshold does not change the conclusion.** 99%-disjoint and exact-disjoint
agree closely (−1.47 against −1.59 pp), both inconclusive, both with block-level SD an order
of magnitude above the isolate-disjoint control. §15's finding is robust to tightening the
criterion.

### 16.2 Tuned readout — the fixed comparison was not biased against the models

Same grid, criterion and budget for every arm: `C ∈ {0.01, 0.1, 1, 10, 100}`, inner 3-fold on
training folds only, criterion balanced accuracy, inner folds grouped whenever the outer split
is grouped. Run under both splitters on the 14 blocks that support haplotype-disjoint folds.
The fixed-`C` result is reported alongside, not replaced.

| Split | Readout | Contrast (pp) | 95% CI (pp) | SD | Verdict at ±1 pp |
|---|---|---|---|---|---|
| isolate-disjoint | fixed C = 1.0 | −0.20 | [−0.39, −0.02] | 0.37 | **EQUIVALENT** |
| isolate-disjoint | **tuned** | **+0.12** | **[−0.49, +1.00]** | **1.52** | **INCONCLUSIVE** |
| sequence-disjoint | fixed C = 1.0 | −0.50 | [−3.20, +2.07] | 5.16 | INCONCLUSIVE |
| sequence-disjoint | **tuned** | **−0.60** | **[−4.79, +3.07]** | **7.76** | INCONCLUSIVE |

Gain from tuning, per arm (tuned minus fixed, pp):

| Split | k-mer | DNABERT-2 | DNABERT-S | HyenaDNA | NTv2 |
|---|---|---|---|---|---|
| isolate-disjoint | +0.03 | **+0.35** | +0.02 | +0.35 | +0.12 |
| sequence-disjoint | **+2.86** | +2.76 | −1.48 | −1.49 | −1.13 |

Three readings, in order of how much they matter:

1. **The reviewer's direction is right but the consequence is not what was claimed.** Under
   isolate-disjoint folds DNABERT-2 does gain more from tuning than k-mer (+0.35 against
   +0.03 pp), so the fixed readout was mildly conservative toward the model. But the effect of
   tuning is to **weaken the equivalence claim, not to make a model win**: the point estimate
   moves to +0.12 pp while the block-level SD quadruples, and the verdict goes
   EQUIVALENT → INCONCLUSIVE. Tuning adds selection variance faster than it adds signal.
2. **Under sequence-disjoint folds, tuning helps k-mer most** (+2.86 pp) and actively *hurts*
   the three non-primary models (−1.1 to −1.5 pp). The §15 conclusion is unchanged and is not
   an artefact of under-tuning the baseline.
3. **The verdict never becomes equivalence in any disjoint condition.** Four cells above, one
   equivalence verdict, and it is the fixed-readout isolate-disjoint control.

**A limitation that must be stated with this result.** `C = 0.01`, the *lower boundary* of the
grid, is selected in the large majority of fold-arm decisions — 47 to 50 of 70 for most arms
under isolate-disjoint. A boundary selection means the optimum may lie outside the grid, so
these tuned numbers are a lower bound on what tuning could achieve, not a converged optimum.
The grid should be extended downward (to 1 × 10⁻⁴) before the tuned comparison is presented as
final. That stronger regularisation is preferred throughout is itself consistent with small
per-block n against high feature dimension.

### 16.3 Per-method region selection — block-level, as requested

The frozen design selects one region by inner CV on the **classical arm only** and applies it
to every arm. Table 7 reports both arms at 96.1% under per-method selection, and the concern
was correctly raised that two rounded values do not establish a zero paired difference. The
block-level paired comparison already exists in `results/region_ranking_locked_test.json` and
needed no re-run:

| Selection rule | k-mer | DNABERT-2 | Paired contrast (pp) | 95% CI (pp) | p_TOST | Verdict |
|---|---|---|---|---|---|---|
| common region, chosen by k-mer (frozen) | — | — | −0.211 | [−0.601, +0.164] | 4.11 × 10⁻⁴ | EQUIVALENT |
| **each arm picks its own region** | 96.091 | 96.066 | **−0.025** | **[−0.339, +0.300]** | **6.36 × 10⁻⁶** | **EQUIVALENT** |

The paired difference is **−0.025 pp, not zero** — 7 of 20 blocks differ by exactly zero and
the rest do not. The headline verdict is **unchanged** under per-method selection, and the
shift from the frozen rule is +0.186 pp. This should be promoted from exploratory to a named
sensitivity analysis on the primary contrast, reported with the paired figure rather than the
two rounded means.

### 16.4 Where items 2 and 3 leave the paper

| Condition | Contrast (pp) | Verdict |
|---|---|---|
| Frozen primary (per-fold region, isolate-grouped) | −0.21 | EQUIVALENT |
| Per-method region selection | −0.03 | EQUIVALENT |
| Fixed region, isolate-disjoint | −0.20 | EQUIVALENT |
| Fixed region, isolate-disjoint, **tuned** | +0.12 | INCONCLUSIVE |
| **Exact-sequence-disjoint** | **−0.50** | **INCONCLUSIVE** |
| **99%-similarity-disjoint** | **−1.47** | **INCONCLUSIVE** |
| Sequence-disjoint, tuned | −0.60 | INCONCLUSIVE |

The equivalence verdict is robust to *region selection* and survives every variation that
keeps sequence reuse intact. It does not survive removing sequence reuse, under either
threshold, tuned or untuned. No condition reverses the direction, and no condition shows a
model advantage. That is a coherent and reportable result, and it is narrower than the claim
the manuscript currently makes.

---

## 17. Divergence stratification and the readout effect on the locked test (2026-11-14)

**Both post hoc. The frozen confirmatory analysis is unchanged.** Scripts:
`analysis_divergence_stratified.py`, `analysis_r1_locked_test.py`. Outputs:
`results/divergence_stratified_results.json`, `results/r1_locked_test_results.json`.

### 17.1 Per-neighbour decomposition and divergence stratification (item 4)

The pooled metric was decomposed, not replaced. The classifier is exactly the one the paper
reports — pooled target-vs-rest, frozen isolate-grouped folds, modal chosen region, nothing
retrained per neighbour. At test time the predictions are partitioned by neighbour species and
balanced accuracy is recomputed over {test targets} ∪ {test isolates of that species}, so this
decomposes the reported metric's own errors. Divergence is the project's canonical identity,
the same measure `C_align` and §16.1 use; the **maximum** cross-pair identity is used as the
stratifier because the hardest pair is what should drive difficulty.

63 target–neighbour pairs across all 20 blocks. Pooled means: k-mer 94.28%,
DNABERT-2 94.23%.

| Max target–neighbour identity | Pairs | Blocks | k-mer | DNABERT-2 | Difference (pp) |
|---|---|---|---|---|---|
| < 95% | 28 | 12 | 99.21 | 99.23 | **+0.03** |
| 95–98% | 7 | 5 | 97.37 | 97.30 | **−0.07** |
| 98–99% | 4 | 3 | 97.48 | 97.25 | **−0.23** |
| 99–99.9% | 2 | 1 | 100.00 | 100.00 | **+0.00** |
| ≥ 99.9% | 22 | 12 | 85.93 | 85.81 | **−0.12** |

Spearman accuracy against max identity: k-mer ρ = **−0.621** (p = 5.6 × 10⁻⁸), DNABERT-2
ρ = **−0.575** (p = 8.2 × 10⁻⁷).

**Three findings, and the third is the one that matters.**

1. **Difficulty is close to binary, not graded.** Pairs whose closest sequences differ at all
   (< 99% identity) are solved at 97–99%. Pairs at ≥ 99.9% — **22 of 63, and in every case the
   neighbour shares a haplotype with the target outright** — collapse to 85.9%. The
   intermediate bins hold 13 of 63 pairs. The determinant of performance is whether the two
   taxa share sequence at that region, not how far apart they are once they differ.
2. **The pooled metric does hide this.** A single 94.3% conceals a 13-point spread between
   easy and hard neighbours within the same blocks. The criticism that pooling discards the
   structure the premise is about is correct.
3. **The method gap is flat across every bin: +0.03, −0.07, −0.23, +0.00, −0.12 pp.** The
   hoped-for result — that one method falls off more slowly as divergence shrinks — does not
   occur. Difficulty swings **13 points** across bins while the method difference stays within
   **0.25 points** in all of them. Target–neighbour sequence identity determines performance
   by roughly **two orders of magnitude more** than method choice.

This strengthens the paper's central claim on a new axis and costs nothing: it is the cheapest
high-value analysis available, as the external review suggested, and it should be a figure.

### 17.2 The readout effect, now measured on the locked test (item 6)

The title foregrounds readout choice, but the R1 → R2 effect was development-only and
unconfirmed. Embeddings were cached, so the confirmatory measurement cost nothing.

**Fidelity first.** The outer folds were rebuilt with the frozen splitter and each fold reuses
the region the frozen run actually selected, read from `_chosen_regions` rather than
re-derived. The rebuilt R2 reproduces the frozen per-block values **exactly**, maximum absolute
deviation **0.00 × 10⁰**. R1 is then interpretable on the same footing.

| Model | R1 | R2 | R1 → R2 (pp) | 95% CI (pp) | Blocks improved |
|---|---|---|---|---|---|
| DNABERT-2 | 94.07 | 95.88 | +1.81 | [+0.47, +3.49] | 11 / 20 |
| DNABERT-S | 93.98 | 95.98 | +2.01 | [+0.75, +3.59] | 14 / 20 |
| HyenaDNA | 93.30 | 95.21 | +1.90 | [+0.46, +3.61] | 13 / 20 |
| NTv2 | 93.60 | 95.81 | +2.21 | [+0.87, +3.85] | 15 / 20 |
| **Mean over models** | | | **+1.98** | **[+0.69, +3.61]** | |

**The readout claim is confirmed, at a smaller magnitude than development suggested.**

| | Readout (R1 → R2) | Model spread at R2 | Ratio |
|---|---|---|---|
| Development (Fig 2d, manuscript) | +2.89 pp | 0.73 pp | 4.0× |
| **Locked test (this analysis)** | **+1.98 pp** | **0.78 pp** | **2.5×** |

Every per-model interval excludes zero, so the effect replicates on confirmatory data. It
shrinks by about a third, from +2.89 to +1.98 pp — the ordinary direction for a
development-to-confirmation move — and the ratio to model identity falls from 4.0× to 2.5×.

**Consequence for the title: it does not need changing, it needs a confirmatory number.** The
claim that readout choice moves performance more than model identity is now supported on the
locked test. What must change is the magnitude quoted alongside it: the manuscript's "roughly
fourfold" is the development figure, and the confirmatory figure is **2.5×**. Report both,
labelled by tier.

This also settles the external review's proposed headline, which asserted that *layer* choice
moves performance several times more than model identity. It does not: the layer effect is
0.53 pp against a model spread of 0.73–0.78 pp, so layer moves performance **less**. Readout is
the factor that dominates, and it is now confirmed.

---

## 18. Standing assessment

Of the twelve checks, **eight are now closed** (1, 2, 3, 4, 5, 6, 7, 10), **one is closed
analytically with an unapplied figure edit** (8), **one is partially closed with four
blockers plus a figure programme** (11), and **one is not resolved and should not be presented
as resolved** (9).

On the figures: six of seven were rebuilt and all seven regenerated to PLOS specification.
Two of the defects were substantive rather than cosmetic. Figure 1b plotted a double-oracle
development contrast beside two pre-named ones — the only positive row in the panel, produced
by the same selection mechanism check 3 flags in the best-of-four rule; it now plots the
matched contrast, +0.12 pp [−0.23, +0.48]. Figure 7 annotated a bound its own data crossed,
on an aggregation that could not bound it, while discarding the model series it computed; it
now plots the fold-matched ceiling and both arms, with no point above the diagonal. Both sets
of legends were rewritten against the rebuilt figures and now agree with each other.

Two findings change what the manuscript can claim. First, the secondary equivalence tests now
exist as formal tests, which **strengthens** the DNABERT-S and NTv2 statements and
**constrains** the HyenaDNA one — equivalence there was tested and failed. Second, the
grouping exposure (92% exact-slice overlap between held-out and training folds) means the
central result should be read as equivalence *within a known haplotype population*, which is
narrower than the current framing and is the most consequential limitation in the paper.

Neither finding touches the primary estimate. It remains −0.21 pp, 95% CI [−0.59, +0.16],
inside the prespecified ±1-point margin, with TOST p = 4.1 × 10⁻⁴, scored once under a plan
frozen before the test data were fetched.
