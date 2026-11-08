# Do pretrained DNA language models help identify bacteria from 16S sequences?

**Final report — 8 November 2026**

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

The test covered **53 groups of bacteria in total**, each consisting of one target species
and three to six of its closest relatives — for example *Staphylococcus aureus* against six
related staphylococci, or *Brucella melitensis* against three other *Brucella* species. The
groups spanned six major bacterial phyla. Sequences came from the NCBI public database:
19,289 records in total.

**Two separate sets of bacteria were used.** The first set (31 groups) was used to make every
decision about how to run the analysis. The second set (22 groups, from entirely different
bacterial families) was used only once, at the end, to check those decisions. The analysis
plan for the second set was written down and locked before the data was downloaded — this is
verifiable in the version history, which shows the plan fixed 94 minutes before the first
sequence arrived.

## 3. The main result

On the locked test set of 20 analysable groups:

| | |
|---|---|
| Average difference (DNA language model minus k-mer counting) | **−0.2 percentage points** |
| 95% confidence interval | **−0.6 to +0.2 percentage points** |
| Groups where the two methods differed by more than 1 point | **3 of 20** |

**The two approaches perform the same.** The confidence interval is narrow enough to rule out
any advantage for the language models larger than about 0.2 percentage points, and any
disadvantage larger than about 0.6 points. The threshold agreed in advance, below which a
difference would not change any practical decision, was 1 percentage point — so the entire
interval sits comfortably inside it.

This is not a case of "we couldn't find a difference because the study was too small". The
study was designed in advance to be precise enough to detect a 1-point difference, and it
achieved better precision than planned. The finding is that there is nothing meaningful there
to find.

The result also **replicated**. The first set of 31 groups gave +0.1 points; the second set of
22 unrelated groups gave −0.2 points. Different bacteria, same conclusion.

### Accuracy of each method

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

## 4. Three findings that matter more than the headline

### 4.1 A random method scores at the top of the table

The top row above is a control: k-mer counts deliberately scrambled through a random
projection, destroying their biological meaning. It scored 96.2% — nominally higher than every
language model.

This is not a flaw in the models. It shows that **at this level of difficulty, most of the
measured accuracy comes from the flexibility of the classifier, not from the biological
content of the representation**. Any future claim that one representation is better than
another on this kind of task should be checked against a random control first. Without it, two
methods can appear to differ when neither is contributing anything specific.

### 4.2 How you read the model matters far more than which model you use

A neural network's output can be turned into a usable representation in several ways. Three
were tested, from simplest to most involved. The differences were striking:

- Across these three reading methods, accuracy varied by **21 to 27 percentage points**.
- Across the four different models, using the best reading method, accuracy varied by
  **0.7 percentage points**.

**The choice of how to read the model matters roughly thirty times more than the choice of
model.** This has a direct practical consequence: any study reporting that a DNA language
model "failed" at a task, based on a single simple reading method, may be measuring its own
method rather than the model.

This project made exactly that error earlier, and the record of it is kept in this repository.
An earlier analysis concluded that DNABERT-2 performed no better than chance (50%) on
*Staphylococcus aureus*. Re-examined with a better reading method and more data, the same model
on the same bacteria scored 96.4%. The model had the information all along; the earlier method
could not reach it.

### 4.3 The models may do relatively better on hard problems — but this is unconfirmed

Half the bacterial groups tested turned out to be easy: k-mer counting was already over 99%
accurate, leaving no room for any method to show an advantage.

Looking only at the hard groups, there is a suggestive pattern. The harder the problem, the
better the language models did relative to k-mer counting (correlation 0.60, p = 0.006).

**This must be treated cautiously.** It was noticed after the main result, not predicted in
advance, and rests on only five difficult groups. On those five groups the difference was
+0.3 points with a confidence interval from −0.5 to +1.2 — consistent with a real small
benefit, and equally consistent with none.

It is the single most promising lead for future work, and it is the reason this report does
not conclude that DNA language models are useless here.

## 5. What this study does and does not show

**Supported:**
- For separating a bacterial species from its close relatives using one 16S section, pretrained
  DNA language models perform equivalently to k-mer counting, within half a percentage point.
- This holds across 45 independent bacterial groups in two separate experiments.
- DNABERT-S, which was specifically trained to distinguish species, did not outperform
  DNABERT-2, which was not. Its training used long genome fragments, and that ability does not
  automatically carry over to short, highly similar 16S sections.
- HyenaDNA performed slightly but reliably worse than k-mer counting (−0.9 points). The
  difference is statistically clear but below the threshold of practical importance.

**Not supported:**
- *That these models contain no useful information about 16S sequences.* They do: all four
  reach 95–96% accuracy, close to the maximum achievable for these data (97%, set by the fact
  that some species share identical sequences in the section being read).
- *That DNA language models are not useful for genomics.* This study tested one task, one
  family of reading methods, four models, and clean reference sequences.
- *That k-mer counting is better.* The finding is equivalence, not superiority.

## 6. Practical implications

**For anyone building a tool to recommend 16S regions:** there is currently no accuracy
justification for including a DNA language model. k-mer counting matches it and is far
cheaper. In this study, producing the language model representations took 36 to 280 seconds
per model and required a 447 MB model file; the k-mer equivalent took 0.19 seconds and
required no model at all.

**For anyone benchmarking DNA language models:** include a random-projection control, and test
more than one way of reading the model. On this evidence, either omission can change the
conclusion.

**For the specific project this grew out of:** the planned next stages were a recommendation
tool and then fine-tuning a model. The evidence supports building the tool with classical
methods. It does not support fine-tuning, because the precondition for that work — classical
methods succeeding where the pretrained models fail — was not observed.

## 7. Honest limitations

- **Half the test groups were too easy** to distinguish between methods. A future study should
  deliberately select difficult groups. This is the main design weakness.
- **Only the final layer of each model was tested.** Intermediate layers might work better. The
  attempt to test this was abandoned rather than botched: the DNABERT-2 implementation stores
  intermediate layers in a compressed internal format, and extracting them correctly would have
  required reimplementing part of the model, risking a silently wrong result. Four different
  ways of summarising the final layer were tested instead, and none helped.
- **Clean reference sequences only.** Real sequencing introduces errors, primer mismatches, and
  complications from bacteria carrying several slightly different copies of the 16S gene. None
  of that was simulated.
- **Four models.** Newer and larger models exist and were not tested.
- **One task.** Closed-group identification, where the correct answer is always among the
  listed candidates.

## 8. Summary

Pretrained DNA language models do not currently offer an accuracy advantage over k-mer counting
for identifying bacteria from 16S sections. This was tested carefully, with the analysis plan
fixed in advance, and confirmed on an independent set of bacteria.

The more useful findings are methodological: a random control scores as highly as any real
method on easy problems, and how a model is read matters about thirty times more than which
model is chosen. Both are easy to overlook, and overlooking either can produce a confident
conclusion that is wrong — as happened once in this project, and is documented here rather
than hidden.
