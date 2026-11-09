| Model | R0 geometry (pair AUC) | R1 nearest reference (% bal. acc.) | R2 linear probe (% bal. acc.) | R2 − R1 (pts) | 95% CI |
|---|---|---|---|---|---|
| DNABERT-2 | 0.741 | 91.8 | 94.6 | +2.83 | [+1.37, +4.57] |
| DNABERT-S | 0.722 | 91.3 | 94.2 | +2.92 | [+1.46, +4.64] |
| HyenaDNA | 0.670 | 90.9 | 93.9 | +3.02 | [+1.44, +4.88] |
| NTv2 | 0.724 | 91.2 | 94.0 | +2.79 | [+1.29, +4.52] |
| k-mer counting | — | — | 94.7 | — | — |

Development set, 25 groups. **R0 is a pair-level AUC and R1/R2 are isolate-level balanced accuracy; R0 must not be subtracted from R1 or R2.** The same-endpoint readout effect is R1 to R2: +2.8 to +3.0 pts, against a 0.73 pt spread across the four models at R2 (~4x). R0 is reported because geometric separation and supervised performance can disagree, not as a comparable accuracy.
