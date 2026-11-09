| Set | Method | n groups | Accuracy of chosen region (%) | Regret (pts) | Picked the best region |
|---|---|---|---|---|---|
| locked test | k-mer counting | 20 | 96.1 | 0.65 | 37% |
| locked test | DNABERT-2 | 20 | 96.1 | 0.57 | 44% |
| locked test | DNABERT-S | 20 | 96.1 | 0.67 | 35% |
| locked test | HyenaDNA | 20 | 95.1 | 1.61 | 29% |
| locked test | NTv2 | 20 | 95.5 | 1.20 | 34% |
| hard panels | k-mer counting | 13 | 81.7 | 1.95 | 60% |
| hard panels | DNABERT-2 | 13 | 81.5 | 2.02 | 60% |
| hard panels | DNABERT-S | 13 | 81.6 | 1.89 | 62% |
| hard panels | HyenaDNA | 13 | 81.7 | 1.37 | 66% |
| hard panels | NTv2 | 13 | 82.1 | 1.61 | 62% |

Each method ranks regions using training isolates only, commits to its top choice, and is scored on held-out isolates. Regret is the accuracy lost relative to the best region identifiable in hindsight.
