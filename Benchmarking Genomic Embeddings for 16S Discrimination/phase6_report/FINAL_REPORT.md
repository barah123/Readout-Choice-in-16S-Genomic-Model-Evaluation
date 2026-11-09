# Do pretrained DNA language models help identify bacteria from 16S sequences?

**Final report — 8 November 2026** (revised after the Phase 7 extension studies)

---

## 1. The question

To identify a bacterium from its 16S rRNA gene, laboratories usually sequence one short
section of that gene rather than all of it. The gene has nine variable sections, labelled V1
to V9. Which section you choose matters: a section that looks identical in your target
species and its close relatives cannot tell them apart.

Recently, "DNA language models" have been proposed for problems like this. These are large
neural networks trained on genome sequences, in the same way text language models are trained
on text. The hope is that they capture something about DNA that simpler methods miss.

**This study asked a single question: for telling one bacterial species apart from its
closest relatives using a 16S section, do these models work better than simply counting
short DNA patterns?**

Counting short patterns — called k-mer counting — is a technique from the 1990s. It involves
no neural network and runs in a fraction of a second.

## 2. What was done

Four published DNA language models were tested: DNABERT-2, DNABERT-S, HyenaDNA, and
Nucleotide Transformer v2. None was retrained or modified; each was used exactly as released.

For each model, the 16S sections were converted into numerical representations, and a simple
classifier was trained on top to separate the target species from its relatives. The same
classifier was trained on k-mer counts for comparison.

**Three separate experiments**, each with its plan written down before the data was collected:

| experiment | groups analysed | purpose |
|---|---|---|
| 1. Development | 25 | make every analysis decision |
| 2. Locked test | 20 | test those decisions once, on unrelated bacteria |
| 3. Hard panels | 13 | test specifically where the task is difficult |
| **total** | **58** | |

Each group is one target species against three to six of its closest relatives — for example
*Staphylococcus aureus* against six related staphylococci. The groups span six major
bacterial phyla. 30,797 sequence records were downloaded from NCBI, of which 24,009 entered
the analyses.

The second experiment is the one that carries the main claim. Its analysis plan was committed
to version control **before any of its data existed** — the version history shows the plan
fixed 94 minutes before the first sequence was downloaded, and unmodified afterwards. The
bacterial families in it do not overlap with the development set at all.

## 3. The main result

On the locked test set of 20 groups:

| | |
|---|---|
| Average difference (DNA language model minus k-mer counting) | **−0.2 percentage points** |
| 95% confidence interval | **−0.6 to +0.2 percentage points** |
| Groups where the two methods differed by more than 1 point | **3 of 20** |

**The two approaches perform the same.** The interval rules out any advantage for the
language models larger than about 0.2 percentage points, and any disadvantage larger than
about 0.6. The threshold agreed in advance, below which a difference would not change any
practical decision, was 1 percentage point — so the whole interval sits comfortably inside it.

A formal two one-sided tests (TOST) procedure against the pre-set ±1-point margin gives
**p = 0.0001**, so equivalence is supported by a test, not only by inspecting an interval.
The 90% interval conventionally used for TOST at α = 0.05 is [−0.53, +0.11], narrower still.

This is not "we couldn't find a difference because the study was too small." The study was
designed in advance to detect a 1-point difference and achieved better precision than
planned. There is nothing meaningful there to find.

**It replicated three times.** Development (25 groups): +0.1 points. Locked test (20
unrelated groups): −0.2 points. Hard panels (13 difficult groups): −0.2 points. Different
bacteria, different difficulty levels, same conclusion.

### Accuracy of each method (locked test)

| method | average accuracy |
|---|---|
| Random projection of k-mer counts | 96.2% |
| k-mer counting | 96.1% |
| DNABERT-S | 96.0% |
| DNABERT-2 | 95.9% |
| Nucleotide Transformer v2 | 95.8% |
| HyenaDNA | 95.2% |
| Direct sequence comparison (alignment) | 94.3% |
| Exact sequence matching | 90.3% |

## 4. Findings that matter more than the headline

### 4.1 A random projection of k-mer counts matches every learned representation

The top row is a control: k-mer counts passed through a random Gaussian projection to 256
dimensions. It scored 96.2% — nominally higher than every language model.

**A correction to an earlier version of this report.** That version described the projection
as "destroying biological meaning" and concluded that most measured accuracy therefore came
from classifier flexibility rather than from biological content. Both statements were wrong.
A random projection approximately preserves pairwise distances (the Johnson–Lindenstrauss
property), so it is a distance-preserving compression of the k-mer features, not a
destruction of them. It retains k-mer information, and nothing about its score speaks to
flexibility versus biology.

The control that does test that question is a **label permutation** through the complete
supervised pipeline — same features, same probe, same cross-validation, labels shuffled
within each group. Run on the locked test (20 groups × 20 permutations):

| arm | observed | permuted labels |
|---|---|---|
| k-mer counting | 96.3% | **50.1%** |
| Random projection of k-mer counts | 96.2% | **50.1%** |
| DNABERT-2 | 96.2% | **50.2%** |

Permuted accuracy is exactly chance for every arm. **The measured accuracy reflects genuine
sequence-to-label signal, not a classifier fitting noise.** The earlier flexibility claim is
withdrawn.

What the random-projection control does show is narrower but still useful: **the specific
identity and weighting of individual k-mer features is not what carries the signal.** A
random linear mixture of those counts performs identically, which indicates substantial
redundancy in the k-mer feature space. Any claim that one representation beats another by a
small margin should be checked against both controls — a random projection for redundancy,
and a label permutation for whether the task is being learned at all.

### 4.2 How you read the model matters more than which model you use

A neural network's output can be turned into a usable representation in several ways. Three
were tested: R0 (unsupervised geometry), R1 (nearest-reference lookup), and R2 (a trained
linear probe).

**A correction to an earlier version of this report.** That version reported a spread of
"21 to 27 percentage points" across readouts and a thirty-fold ratio over model choice. That
comparison was invalid. **R0 is a pair-level AUC, while R1 and R2 are isolate-level balanced
accuracy.** These are different quantities on different units and cannot be subtracted. The
large apparent gap was substantially an artifact of mixing metrics.

The defensible comparison uses R1 and R2, which share an endpoint:

| | effect |
|---|---|
| R1 → R2 (same metric, paired, 25 groups) | **+2.8 to +3.0 pts**; DNABERT-2 +2.83 [+1.37, +4.57] |
| Spread across the four models at R2 | **0.73 pts** |
| Ratio | **~4×**, not 30× |

**Readout choice still matters about four times more than model choice, and the readout
effect is statistically supported while the model effect is not.** That is a weaker claim
than the original but it rests on a like-for-like comparison.

R0 remains worth reporting for a different reason: geometric separation and supervised
performance disagree. DNABERT-S scores 0.542 by pair AUC at staphylococcal V1 yet 0.945 by
probe, and several *Bacillus* regions show high R0 where the ceiling and supervised accuracy
are both low. **A model that looks poor under raw embedding geometry can be fully competitive
under a trained probe** — which is the substantive point, and it does not require any
cross-metric subtraction.

The practical implication is unchanged: a study reporting that a DNA language model "failed"
a task, based on one simple reading method, may be measuring its own method rather than the
model.

This project made exactly that error. An earlier analysis concluded DNABERT-2 performed at
chance (50%) on *Staphylococcus aureus*. Re-examined with a better reading method and more
data, the same model on the same bacteria scored 96.4%. The information was there all along;
the earlier method could not reach it. The record is kept in this repository.

### 4.3 Difficulty does not rescue the models — tested and not supported

*This section previously reported a promising lead. It was tested directly and did not hold.*

Half the locked-test groups turned out easy: k-mer counting was already above 99%, leaving no
room for any method to show an advantage. Within the five genuinely hard groups there was a
suggestive pattern — the harder the problem, the better the models did relative to k-mer
counting (correlation 0.60, p = 0.006).

That was a hypothesis generated after seeing the data, on five groups. A third experiment
tested it properly, on 13 groups selected specifically for difficulty. Selection used only
sequence properties and the classical method, so it could not favour the models. These groups
are genuinely hard: average k-mer accuracy 0.82 versus 0.96 in the locked test, the hardest
at 0.55.

| | |
|---|---|
| Average difference on hard groups | **−0.2 percentage points** |
| 95% confidence interval | **−1.3 to +0.7 points** |
| Groups favouring model / k-mer / tie | 4 / 3 / 6 |

**The earlier correlation did not replicate.** No regrouping of the data was used to try to
rescue it, and the largest single effect on any hard group (*Bradyrhizobium*, −5.0 points)
favours k-mer counting.

One unresolved observation: taking the *best* of the four models for each group gives +0.7
points (interval −0.1 to +1.4). But taking the maximum of four noisy measurements is
upward-biased by construction, so that is the explanation to assume. Telling real
model-complementarity apart from selection bias would need a separate study with the
model-choice rule fixed in advance. **This is a hypothesis, not a finding.**

### 4.4 The conclusion survives imperfect sequences — and read length matters more than the model

All results above use clean reference sequences. A fourth experiment degraded them in six
realistic ways — sequencing errors at two rates, 10 bases lost from each end to simulate
imperfect primer binding, the natural variation between gene copies inside one bacterium, and
a combined worst case. Both methods received identical damaged sequences.

No condition changed the answer. But the degradation itself was revealing:

| damage | accuracy lost (both methods) |
|---|---|
| 0.1% sequencing error | ~0.3 points |
| 1% sequencing error | ~2–3 points |
| **10 bases trimmed each end** | **~4 points** |
| gene-copy variation | ~3 points |
| error + trimming combined | ~5 points |

**Losing read length costs more than a tenfold increase in sequencing error.** For anyone
choosing a 16S strategy in practice, primer placement and read coverage matter more than
which classification method is used — a bigger effect than anything separating the two
methods in this entire study.

### 4.5 Where 16S difficulty actually lives

Screening 53 candidate bacterial groups for difficulty produced a result independent of any
model, and arguably the most directly useful one here.

Groups chosen *specifically because* their 16S is reputed to be highly conserved —
*Rickettsia*, *Chlamydia*, *Xylella*, *Synechococcus*/*Prochlorococcus*, *Desulfovibrio*,
*Thermus* — turned out **perfectly separable**, at 100% accuracy.

Genuine difficulty was concentrated almost entirely in **subspecies comparisons and tight
species complexes**: *Mycobacterium abscessus* subspecies (55%), the *Xanthomonas campestris*
group (63%), the *Pseudomonas fluorescens* group (77%), the *Burkholderia cepacia* complex
(81%).

This is the third independent time in this project that the textbook claim "16S cannot
separate these organisms" has failed to hold up given enough isolates and the right region.
Within the groups screened here, **the lowest accuracies occurred in subspecies and
species-complex comparisons, and genera with a reputation for conserved 16S were separable**.
Taxonomic rank was not a prespecified variable and was not formally modelled, so this is a
description of these panels rather than a general rule. It is offered as a hypothesis worth
testing directly: that rank predicts 16S difficulty better than a genus's reputation does.

### 4.6 Layer depth barely matters — which sharpens the readout finding

Section 4.2 claims that how you read a model matters far more than which model you use. That
claim had an obvious weak point: every result used the **final layer** of each network, so
the conclusion could itself have been limited by that one choice.

This was tested directly on the two models whose internals can be read safely. Both were
first checked to confirm the extracted final layer exactly reproduces the standard output —
the check whose absence caused this analysis to be deferred earlier.

| model | layers tested | best layer | gain over the final layer |
|---|---|---|---|
| HyenaDNA | 5 of 10 | layer 8 | **+0.16 points** (interval −0.07 to +0.43) |
| NTv2 | 5 of 30 | **the final layer** | **0.00 points** |

And those figures are *upper bounds*: the best layer was picked by looking at the same data,
which biases it upward.

**Layer choice is a minor factor.** Accuracy varied by about 0.5 points across five depths,
against 21–27 points across readout *types*. **Readout type matters roughly fifty times more
than layer depth.** This strengthens rather than weakens section 4.2: what matters is whether
you train a classifier on the representation at all, not which layer you draw it from.

DNABERT-2 and DNABERT-S remain untested here, and that is stated as a limitation rather than
glossed. The attempt was not made because extracting their intermediate layers correctly
requires reimplementing part of the model, and a silently wrong representation would be worse
than an acknowledged gap.

## 5. What this study does and does not show

**Supported:**
- For separating a bacterial species from its close relatives using one 16S section,
  pretrained DNA language models perform equivalently to k-mer counting, within about half a
  percentage point.
- This holds across **58 independent bacterial groups** in three experiments, spanning easy,
  intermediate and genuinely hard problems, and survives realistic sequence degradation.
- DNABERT-S, trained specifically to distinguish species, did not outperform DNABERT-2, which
  was not. Its training used long genome fragments; that ability does not transfer
  automatically to short, highly similar 16S sections.
- HyenaDNA performed slightly but reliably worse than k-mer counting (−0.9 points) — clear
  statistically, below the threshold of practical importance.

**Not supported:**
- *That these models contain no useful information about 16S sequences.* They do: all four
  reach 95–96%, against an empirical ceiling of 97.0% for this reference population. That
  ceiling is reference-specific, not a universal maximum: it is computed on the same isolate
  unit, the same species-balanced weighting, and the same regions as the reported accuracy,
  and it would change with different reference sampling or class weights.
- *That DNA language models are not useful for genomics.* One task, one family of reading
  methods, four models.
- *That k-mer counting is better.* The finding is equivalence, not superiority.
- *That the models help on hard problems.* Explicitly tested in experiment 3 and not found.

## 6. Practical implications

**For anyone building a tool to recommend 16S regions:** there is no accuracy justification
for including a DNA language model. k-mer counting matches it and is far cheaper — 36 to 280
seconds and a 447 MB model file, versus 0.19 seconds and no model.

**For anyone planning 16S sequencing:** read length and primer placement matter more than the
choice of classification method, and whether your organisms are separable depends more on
whether you are distinguishing species or subspecies than on the genus's reputation.

**For anyone benchmarking DNA language models:** include a random-projection control, test
more than one way of reading the model, and check whether any result survives removing its
single largest contributing group. In this project all three checks changed a conclusion.

**For the specific project this grew out of:** the planned next stages were a recommendation
tool and then fine-tuning a model. The evidence supports building the tool with classical
methods. It does not support fine-tuning: the precondition for that work — classical methods
succeeding where the pretrained models fail — was looked for directly and not observed.

## 7. Honest limitations

Two limitations in the earlier version of this report have since been addressed; they are
marked below.

- ~~Half the test groups were too easy.~~ **Addressed** by experiment 3, which selected 13
  groups specifically for difficulty and reached the same conclusion.
- ~~Clean reference sequences only.~~ **Addressed** by experiment 4, which degraded sequences
  six ways and reached the same conclusion. The error model used (uniform random base changes)
  is still a crude stand-in for real sequencing artefacts.
- ~~Only the final layer of each model was tested.~~ **Partly addressed** (see 4.6). Tested
  for HyenaDNA and NTv2 across five depths each: the final layer was as good as any
  intermediate one. Still untested for DNABERT-2 and DNABERT-S, which store intermediate
  layers in a compressed internal format that cannot be extracted without reimplementing part
  of the model.
- **Four models.** Newer and larger ones exist and were not tested. Given that model choice
  matters about four times less than reading method, this is a lower priority than it looks.
- **One task.** Closed-group identification, where the right answer is always among the
  listed candidates.
- **The very hardest comparisons are excluded by design.** Groups need at least 10 target
  isolates and 2 relatives with 5 each. Organisms too rare to meet that — often the most
  clinically awkward to identify — cannot be assessed this way.

## 8. Summary

Pretrained DNA language models do not currently offer an accuracy advantage over k-mer
counting for identifying bacteria from 16S sections. This was tested with the analysis plan
fixed in advance, confirmed on an independent set of bacteria, confirmed again on groups
chosen specifically because they are difficult, and confirmed again on deliberately degraded
sequences. 58 bacterial groups in total.

The more useful findings are methodological and practical: a random projection of k-mer
counts matches every learned representation, while permuted labels collapse every method to
chance; how a model is read matters about four times more than which model is chosen, on a
like-for-like comparison; read length matters more than either; and the organisms that are
genuinely hard to tell apart by 16S are not the ones with that reputation.

Three times during this work a confident-looking conclusion turned out to be an artifact —
once from a poor reading method, once from a single small group driving an average, and once
from a correlation that did not replicate. Each is documented here rather than removed,
because the checks that caught them are the part most worth reusing.
