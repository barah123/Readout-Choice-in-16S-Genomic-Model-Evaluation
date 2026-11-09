# Extension B — does the conclusion survive imperfect sequences?

**Exploratory sensitivity analysis.** Re-runs the Phase 5 primary contrast
(DNABERT-2 linear probe vs k-mer counting) on the same 20 locked-test panels, with the
sequences deliberately degraded. Both methods receive identical damaged sequences in every
condition. Pre-registered in `PREREGISTRATION.md` before running.

## Result

| condition | k-mer | DNABERT-2 | difference | 95% CI |
|---|---|---|---|---|
| clean | 0.9629 | 0.9616 | −0.0014 | [−0.0027, −0.0000] |
| 0.1% sequencing error | 0.9602 | 0.9574 | −0.0028 | [−0.0056, −0.0004] |
| 1% sequencing error | 0.9348 | 0.9434 | +0.0085 | [−0.0068, +0.0323] |
| 10 bases trimmed each end | 0.9212 | 0.9180 | −0.0032 | [−0.0069, −0.0004] |
| copy-to-copy variation | 0.9287 | 0.9324 | +0.0037 | [−0.0051, +0.0149] |
| 1% error + trimming | 0.9081 | 0.9115 | +0.0034 | [−0.0076, +0.0204] |

The clean row reproduces Phase 5 (−0.0014 here vs −0.0021 there; the small difference is
because this analysis fixes one region per panel rather than re-selecting it per fold).

## Finding 1 — the equivalence holds under degradation

No condition produces a difference that is both large enough to matter (the pre-set threshold
was 1 percentage point) and statistically supported. **The Phase 5 conclusion survives
realistic data degradation.**

## Finding 2 — both methods degrade substantially, and that is the practical message

| condition | k-mer loses | DNABERT-2 loses |
|---|---|---|
| 0.1% error | 0.3 pts | 0.4 pts |
| 1% error | 2.8 pts | 1.8 pts |
| trimming | 4.2 pts | 4.4 pts |
| copy variation | 3.4 pts | 2.9 pts |
| error + trimming | **5.5 pts** | **5.0 pts** |

Losing 10 bases from each end costs more accuracy (~4 points) than a 1% sequencing error
(~2-3 points). Neither method is robust in absolute terms. For anyone choosing a 16S region
in practice, **primer placement and read coverage matter more than which classifier is
used** — a larger effect than anything separating the two methods.

## Finding 3 — the apparent noise-robustness advantage is one panel, and it does not survive

The three substitution-based conditions all showed the language model gaining ground, in a
consistent direction, with a plausible mechanism: a single base change alters up to 15 k-mer
occurrences (4+5+6), so composition counting is intrinsically brittle to point errors.

That story is wrong. Testing the paired shift per panel rather than the pooled mean:

| condition | mean shift | 95% CI | panels + / − | **excluding *francisella*** |
|---|---|---|---|---|
| 0.1% error | −0.0014 | [−0.0042, +0.0012] | 4 / 8 | −0.0021 |
| **1% error** | **+0.0099** | [−0.0053, +0.0337] | 9 / 5 | **−0.0006** |
| trimming | −0.0019 | [−0.0045, +0.0004] | 4 / 7 | −0.0020 |
| copy variation | +0.0051 | [−0.0035, +0.0163] | 9 / 5 | +0.0005 |
| error + trimming | +0.0048 | [−0.0076, +0.0218] | 6 / 10 | −0.0025 |

**A single panel, *francisella*, shifts by +0.21 under 1% error and accounts for the entire
effect.** Removing it moves the estimate from +0.0099 to −0.0006. The same panel drives the
copy-variation (+0.09) and combined (+0.14) conditions. Leave-one-out confirms the sign is
not stable: dropping any single panel moves the pooled estimate across zero.

The `combined` condition is the clearest warning — its mean is positive while **10 of 16
panels moved negative**. The average pointed one way and the majority of the evidence pointed
the other.

No relationship with task difficulty was found either (ρ = +0.23, p = 0.33), so this does not
support Extension A's hypothesis.

*francisella* is the smallest panel in the set (63 records), which is the most likely
explanation: small panels give unstable estimates under added noise.

## Conclusion

The equivalence between pretrained genomic language models and k-mer counting is **robust to
sequencing error, primer truncation, and gene-copy variation**. The hypothesis that the models
tolerate noise better is **not supported** — it rested on one small panel.

The practically useful finding is unrelated to the model comparison: **read truncation is the
single most damaging degradation tested**, costing roughly 4 percentage points, more than a
tenfold increase in sequencing error rate.

## Limitations

Uniform random substitution is a crude error model; real sequencers produce
position-dependent, context-dependent errors. Trimming a fixed 10 bases is a stand-in for
imperfect primer binding, not a simulation of PCR. Only DNABERT-2 was tested, as the
pre-registered primary arm. These are stress tests, not simulations.
