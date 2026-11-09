# Extension C — are intermediate layers better than the final layer?

**Pre-registered before running.** Closes the strongest internal tension in this work: the
headline claim is that *how you read a model* dominates, yet every result used the final
layer only.

**Answer: no. The final layer is as good as any intermediate layer, for both models tested.**

---

## 1. Correctness gate

The reason this analysis was deferred in Phase 6 was the risk of extracting a silently wrong
representation. Both models were therefore checked before any result was computed, with the
check enforced in code (the run aborts if it fails):

> the extracted final layer must be numerically identical to the standard forward output

| model | hidden states exposed | gate |
|---|---|---|
| HyenaDNA | 10 | **passed** (`allclose`, atol 1e-6) |
| NTv2 | 30 | **passed** (`allclose`, atol 1e-6) |

DNABERT-2 and DNABERT-S still fail this standard — they run their encoder on unpadded
sequences internally, so intermediate layers come back shaped `(total_tokens, hidden)` rather
than `(batch, seq, hidden)`, and reconstructing the padded form means reproducing their
unpad/repad index logic. They remain untested, and that remains a stated limitation.

## 2. Result

13 screened-hard panels (the only set with headroom), 5 depths per model, everything else
frozen. Mean balanced accuracy across panels:

**HyenaDNA** (k-mer baseline 0.8280)

| layer | accuracy | vs k-mer |
|---|---|---|
| 2 | 0.8286 | +0.0007 |
| 4 | 0.8293 | +0.0013 |
| 7 | 0.8311 | +0.0032 |
| 8 | 0.8340 | +0.0060 |
| **9 (final, used in the frozen analysis)** | **0.8324** | **+0.0044** |

Best depth 8, gain over the final layer **+0.0016, 95% CI [−0.0007, +0.0043]**.

**NTv2** (k-mer baseline 0.8280)

| layer | accuracy | vs k-mer |
|---|---|---|
| 7 | 0.8222 | −0.0057 |
| 14 | 0.8211 | −0.0068 |
| 22 | 0.8263 | −0.0017 |
| 26 | 0.8233 | −0.0046 |
| **29 (final, used in the frozen analysis)** | **0.8264** | **−0.0015** |

Best depth **is** the final layer; gain **+0.0000**.

## 3. Reading this correctly

The quoted "gain over final layer" is an **upper bound, not an estimate**. The best depth was
chosen by looking at these same panels, so taking the maximum over five noisy measurements is
upward-biased by construction — the same selection effect flagged for the best-of-4 model
comparison in Extension A. Even accepting that inflated figure, the gain is +0.0016 for
HyenaDNA with an interval spanning zero, and exactly zero for NTv2.

For NTv2 the deepest layer was genuinely best, so there is no selection concern at all there.

## 4. What this settles, and what it does not

**Settled, for these two models:** the frozen choice of the final layer was not a mistake and
did not understate them. Outcome 1 of the three pre-registered possibilities. The main
equivalence result does not need re-running.

**Also settled:** layer choice is a *minor* factor. The spread across five depths is 0.005
(HyenaDNA) and 0.005 (NTv2). For comparison, the spread across readout *types* measured
earlier — raw geometry vs nearest-reference vs linear probe — was 0.21 to 0.27. **Readout
type matters roughly fifty times more than layer depth.** That sharpens rather than
undermines the original claim: the thing that matters is whether you train a probe at all,
not which layer you draw from.

**Not settled:** DNABERT-2 and DNABERT-S, which carry the primary arm of the main analysis.
The honest position is that the layer question has been answered for half the models, using
the two where it could be answered safely, and that the answer was negative in both cases.
Extending it to the MosaicBERT pair requires implementing and validating their unpad/repad
logic — worthwhile, but it was not worth risking a wrong number to avoid an honest
limitation.

## 5. Incidental observation

On `mycobacterium_abscessus_subsp` — the hardest panel in the project, with only 3 unique
sequences across all isolates — every model at every layer returns exactly 0.577, identical
to k-mer and identical to the ambiguity ceiling. When three sequences must separate two
subspecies, no representation can help: the information is absent from the data, not from
the model. This is the clearest single illustration of the ceiling concept used throughout
this project.
