# Extension A — do the models help on genuinely hard problems?

**Pre-registered in `PREREGISTRATION.md` before any panel was screened or fetched.**
Hypothesis tested: on panels where k-mer counting is genuinely challenged, frozen genomic
language model representations achieve **higher** balanced accuracy.

**Answer: no. The hypothesis is not supported.**

---

## 1. Panel selection

30 candidates screened in round 1, 23 more in round 2 (Amendment A-1, after 11 round-1
candidates failed on insufficient data rather than on difficulty). Selection used the
ambiguity ceiling and k-mer accuracy only — **no genomic language model was involved**, so
selection could not favour the models.

15 panels passed the difficulty screen; 13 passed the eligibility gate and entered the
analysis, meeting the pre-set target of ≥12.

Two selected panels were lost to the gate, both for the same pre-registered reason — fewer
than two neighbour species with ≥5 isolates:

| panel | why |
|---|---|
| `salmonella_subspecies` | only 1 qualifying neighbour; non-*enterica* subspecies rarely deposited |
| `campylobacter_jejuni_subsp` | only 1 qualifying neighbour; *C. jejuni* subsp. *doylei* rarely deposited |

The gate was not relaxed after seeing the data.

## 2. These panels are genuinely hard

Unlike Phase 5, where half the panels were solved at over 99% by k-mer counting alone, every
panel here leaves real room:

| panel | ceiling | k-mer |
|---|---|---|
| mycobacterium_abscessus_subsp | 0.577 | 0.552 |
| xanthomonas_campestris_pv | 0.647 | 0.634 |
| pseudomonas_fluorescens_group | 0.841 | 0.766 |
| burkholderia_cepacia_complex | 0.840 | 0.806 |
| clostridium_botulinum_groups | 0.881 | 0.801 |
| lactobacillus_casei_group | 0.846 | 0.817 |
| bacillus_subtilis_subsp | 0.862 | 0.838 |
| nostoc_anabaena | 0.959 | 0.851 |
| lactococcus_lactis_subsp | 0.898 | 0.883 |
| pectobacterium_dickeya | 0.941 | 0.879 |
| bradyrhizobium | 0.916 | 0.918 |
| streptococcus_dysgalactiae_subsp | 0.955 | 0.944 |
| lactobacillus_delbrueckii_subsp | 0.960 | 0.929 |

Mean k-mer accuracy is 0.818, versus 0.961 on the Phase 5 panels. This is the difficulty
regime Phase 6 said was missing.

## 3. Primary result

| | |
|---|---|
| mean DNABERT-2 − k-mer | **−0.0019** |
| 95% bootstrap CI | **[−0.0128, +0.0066]** |
| panels favouring model / k-mer / tie (±0.01) | 4 / 3 / 6 |
| **verdict** | **H1 not supported — interval contains zero** |

**The Phase 6 correlation is reported as not replicated.** As pre-committed, no post-hoc
regrouping has been used to rescue it. On panels selected specifically for difficulty, the
two approaches remain indistinguishable.

The leave-one-out check confirms there is no hidden single-panel effect in either direction:
dropping any one panel moves the mean only between −0.0031 and +0.0021, and the sign is not
stable — exactly what a true null looks like.

## 4. One result that is not null, and should be followed up carefully

The secondary contrast — **best-of-4 models** minus k-mer — is **+0.0070
[−0.0010, +0.0135]**. The interval only just includes zero, and the point estimate is three
times larger than anything seen on easy panels.

Two readings, and the second is more likely:

1. Different models genuinely suit different hard panels, so picking the best per panel
   captures real complementary information.
2. **Selection bias.** Taking the maximum of four noisy measurements is upward-biased by
   construction. With a true difference of zero and four arms, "best-of-4" is expected to be
   positive even when no model is better.

Reading 2 is the default explanation and the one to assume. Distinguishing them requires a
rule for choosing the model *in advance* of seeing panel results — which is exactly why the
frozen plan named DNABERT-2 as the single primary arm. **This is a hypothesis for a future
pre-registered test, not a finding.**

Also notable: the single largest per-panel difference is `bradyrhizobium` at **−0.050**, k-mer
ahead. The largest effect on hard panels favours the classical method.

## 5. What this adds

Phase 5 established equivalence but could be criticised on the grounds that half its panels
were too easy to distinguish any method. **That criticism is now answered.** On 13 panels
chosen specifically because k-mer counting struggles with them — mean accuracy 0.818, one as
low as 0.552 — the two approaches are still indistinguishable.

Combined with Phase 5 (20 panels) and Phase 2-full (25 panels), the equivalence now holds
across **58 independent bacterial groups** spanning easy, intermediate, and genuinely hard
problems.

## 6. Where difficulty actually lives — an incidental finding

The screening produced a result independent of any model, and arguably more useful to
practitioners than the model comparison.

Groups chosen *specifically because* their 16S is reputed to be highly conserved —
*Rickettsia*, *Chlamydia*, *Xylella*, *Synechococcus*/*Prochlorococcus*, *Desulfovibrio*,
*Thermus*, *Arcobacter*, *Tannerella* — screened at a ceiling of **1.000** with k-mer accuracy
**1.000**. Perfectly separable.

Genuine 16S difficulty was concentrated almost entirely in **subspecies comparisons and tight
species complexes**: *M. abscessus* subspecies (0.577), the *Xanthomonas campestris* group
(0.647), the *Pseudomonas fluorescens* group (0.766), the *Burkholderia cepacia* complex
(0.806).

This is the third independent time in this project that the textbook claim "16S cannot
separate these organisms" has failed to hold up given adequate isolate numbers and the right
region. The earlier two were Phase 2-full (*E. coli*/*Shigella*, *Enterobacter*,
*Acinetobacter* all largely separable) and Phase 5 (the a-priori easy/partial distinction
collapsed entirely).

**Practical implication:** when assessing whether 16S can resolve a group, the rank of the
comparison matters more than the genus's reputation. Species-level comparisons are usually
tractable; subspecies-level ones usually are not.

## 7. Limitations

13 panels, so differences smaller than roughly 1 percentage point remain undetectable.
Selection required ≥10 target isolates and ≥2 neighbour species with ≥5 — the very hardest
comparisons, between organisms with few public deposits, are systematically excluded by that
requirement. Clean reference sequences only. The best-of-4 signal is unresolved and would
need its own pre-registered test with a model-selection rule fixed in advance.
