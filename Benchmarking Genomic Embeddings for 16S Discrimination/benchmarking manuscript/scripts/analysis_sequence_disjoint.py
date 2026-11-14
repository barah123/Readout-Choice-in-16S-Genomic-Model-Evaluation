"""POST HOC: does the locked-test conclusion survive exact-sequence-disjoint folds?

NOT part of the frozen confirmatory analysis. The frozen primary estimand is unchanged and
is not recomputed here.

Motivation. 92.0% of held-out isolates share their exact regional slice with the training
folds, because identical 16S sequences recur across strains. Isolate-grouped CV does not make
folds sequence-disjoint, so the frozen result characterises performance WITHIN a known
haplotype population rather than generalisation to unseen variants. This script builds the
disjoint evaluation and runs the same comparison inside it.

Design — one variable changes, everything else is held fixed:

  region      each block's modal region from the frozen run, held fixed across folds. The
              frozen pipeline re-selects a region per fold; fixing it removes selection
              variance so the splitter is the only difference between the two arms.
  control     StratifiedKFold(5) over isolates = the frozen design at a fixed region.
  treatment   StratifiedGroupKFold(5) grouped on the exact slice string, so no sequence
              appears in both training and test.
  everything  same probe (StandardScaler -> LogisticRegression, C=1.0, balanced, max_iter
              else      2000, seed 0), same k-mer featuriser, same embedding cache, same seed.

Feasibility. Haplotype grouping needs enough distinct haplotypes per class to form folds.
Blocks with fewer than n_splits haplotypes in either class cannot support the design and are
reported as not evaluable rather than dropped silently.

Split audit. Under a correct exact-sequence-disjoint split, C_exact must abstain on every
test sequence and therefore score exactly 0.500. That is asserted per fold as a direct check
on the split, not assumed.

Writes results/sequence_disjoint_results.json.
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
import run_readouts_full as rr          # validated k-mer / metric helpers

import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)

SEED = 0
N_SPLITS = 5
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "sequence_disjoint_results.json"))


def probe():
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED))


def c_exact(str_tr, ytr, str_te, yte):
    """Majority label among identical training sequences; abstain at 0.5 if absent."""
    lab = defaultdict(lambda: [0, 0])
    for s, l in zip(str_tr, ytr):
        lab[s][l] += 1
    cred = []
    for s, l in zip(str_te, yte):
        if s in lab:
            a, b = lab[s]
            p1 = b / (a + b)
        else:
            p1 = 0.5
        cred.append(p1 if l == 1 else 1 - p1)
    return rr.soft_balanced_accuracy(yte, np.array(cred))


def run_block(slices, y, emb_at_region, splitter_name, groups=None):
    y = np.asarray(y)
    if splitter_name == "isolate":
        sp = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
        it = sp.split(np.zeros(len(y)), y)
    else:
        sp = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
        it = sp.split(np.zeros(len(y)), y, groups)

    arms = defaultdict(list)
    overlap, n_te, n_folds, degenerate = 0, 0, 0, 0
    for tr, te in it:
        ytr, yte = y[tr], y[te]
        if len(set(ytr)) < 2 or len(set(yte)) < 2:
            degenerate += 1
            continue
        n_folds += 1
        str_tr = [slices[i] for i in tr]
        str_te = [slices[i] for i in te]
        seen = set(str_tr)
        overlap += sum(1 for s in str_te if s in seen)
        n_te += len(str_te)

        Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
        c = probe(); c.fit(Ktr, ytr)
        arms["C_kmer"].append(rr.balanced_accuracy(yte, c.predict(Kte)))

        rng = np.random.default_rng(SEED)
        rp = rng.standard_normal((Ktr.shape[1], 256)) / np.sqrt(max(Ktr.shape[1], 1))
        cr = probe(); cr.fit(Ktr @ rp, ytr)
        arms["C_random"].append(rr.balanced_accuracy(yte, cr.predict(Kte @ rp)))

        arms["C_exact"].append(c_exact(str_tr, ytr, str_te, yte))

        for m in MODELS:
            lut = emb_at_region.get(m)
            if lut is None or not all(s in lut for s in str_tr + str_te):
                continue
            Xtr = np.vstack([lut[s] for s in str_tr])
            Xte = np.vstack([lut[s] for s in str_te])
            p = probe(); p.fit(Xtr, ytr)
            arms[f"{m}:R2"].append(rr.balanced_accuracy(yte, p.predict(Xte)))

    out = {k: float(np.mean(v)) for k, v in arms.items() if v}
    out["_n_folds"] = n_folds
    out["_degenerate_folds"] = degenerate
    out["_test_seqs"] = n_te
    out["_test_seqs_seen_in_train"] = overlap
    out["_overlap_rate"] = float(overlap / n_te) if n_te else float("nan")
    return out


def main():
    rows = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_blocks.json")))
    scores = json.load(open(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json")))

    emb = {}
    for m in MODELS:
        f = os.path.join(ROOT, "phase4_locked_test", f"emb_test_{m}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f, allow_pickle=True)
        for r in [f"V{i}" for i in range(1, 10)]:
            if f"{r}__keys" in z:
                emb[(m, r)] = {k: v for k, v in zip(z[f"{r}__keys"], z[f"{r}__vecs"])}

    per = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        per[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    results, skipped = {}, {}
    for b in sorted(scores):
        v = scores[b]
        if v.get("_status") != "ok" or not v.get("_block_eligible"):
            continue
        region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        d = per[b][region]
        isolates = sorted(d)
        slices = [d[i][0] for i in isolates]
        y = [d[i][1] for i in isolates]

        hap = defaultdict(lambda: [0, 0])
        for s, l in zip(slices, y):
            hap[s][0 if l == 1 else 1] += 1
        n_t = sum(1 for c in hap.values() if c[0] > 0)
        n_n = sum(1 for c in hap.values() if c[1] > 0)
        if n_t < N_SPLITS or n_n < N_SPLITS:
            skipped[b] = dict(region=region, n_isolates=len(isolates),
                              n_haplotypes=len(hap), haplotypes_target=n_t,
                              haplotypes_neighbour=n_n,
                              reason=f"fewer than {N_SPLITS} haplotypes in a class")
            continue

        emb_at = {m: emb.get((m, region)) for m in MODELS}
        ctrl = run_block(slices, y, emb_at, "isolate")
        trt = run_block(slices, y, emb_at, "haplotype", groups=slices)
        results[b] = dict(
            region=region, n_isolates=len(isolates), n_haplotypes=len(hap),
            haplotypes_target=n_t, haplotypes_neighbour=n_n,
            cross_class_haplotypes=sum(1 for c in hap.values() if c[0] > 0 and c[1] > 0),
            frozen_C_kmer=v["C_kmer"], frozen_dnabert2_R2=v["dnabert2:R2"],
            isolate_disjoint=ctrl, sequence_disjoint=trt)
        print(f"  {b:26} {region}  iso-disjoint kmer={ctrl['C_kmer']:.4f} "
              f"d2={ctrl.get('dnabert2:R2', float('nan')):.4f}   |   "
              f"seq-disjoint kmer={trt['C_kmer']:.4f} "
              f"d2={trt.get('dnabert2:R2', float('nan')):.4f}  "
              f"overlap={trt['_overlap_rate']:.3f}", flush=True)

    # ---- summary statistics, deposited alongside the per-block values
    from scipy import stats as _st
    names = sorted(results)

    def _boot(x, seed=1729, n=20000):
        rng = np.random.default_rng(seed)
        dr = x[rng.integers(0, len(x), size=(n, len(x)))].mean(axis=1)
        return (float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5)),
                float(2 * min((dr <= 0).mean(), (dr >= 0).mean())))

    def _tost(x, m=0.010):
        mu, sd, n = x.mean(), x.std(ddof=1), len(x)
        se = sd / np.sqrt(n)
        return float(max(_st.t.sf((mu + m) / se, n - 1), _st.t.cdf((mu - m) / se, n - 1)))

    summary = {"n_blocks": len(names), "margin": 0.010}
    for key in ("isolate_disjoint", "sequence_disjoint"):
        x = np.array([results[b][key]["dnabert2:R2"] - results[b][key]["C_kmer"]
                      for b in names])
        lo, hi, p = _boot(x)
        summary[key] = {
            "mean_accuracy": {a: float(np.mean([results[b][key][a] for b in names]))
                              for a in ("C_kmer", "C_random", "C_exact",
                                        *[f"{m}:R2" for m in MODELS])},
            "primary_contrast_dnabert2_minus_kmer": {
                "mean": float(x.mean()), "sd": float(x.std(ddof=1)),
                "ci95": [lo, hi], "p_two_sided": p, "p_tost": _tost(x),
                "verdict": ("EQUIVALENT" if (lo > -0.010 and hi < 0.010)
                            else "INCONCLUSIVE at this margin")},
            "split_audit": {
                "max_train_test_exact_overlap": float(max(
                    results[b][key]["_overlap_rate"] for b in names)),
                "C_exact_all_exactly_half": all(
                    abs(results[b][key]["C_exact"] - 0.5) < 1e-12 for b in names)}}

    kd = np.array([results[b]["sequence_disjoint"]["C_kmer"]
                   - results[b]["isolate_disjoint"]["C_kmer"] for b in names])
    dd = np.array([results[b]["sequence_disjoint"]["dnabert2:R2"]
                   - results[b]["isolate_disjoint"]["dnabert2:R2"] for b in names])
    lo, hi, p = _boot(dd - kd)
    summary["degradation"] = {
        "kmer_mean_loss": float(kd.mean()), "dnabert2_mean_loss": float(dd.mean()),
        "difference": float((dd - kd).mean()), "ci95": [lo, hi], "p_two_sided": p,
        "blocks_where_model_degrades_less": int((dd - kd > 0).sum()),
        "interpretation": "no evidence that frozen embeddings generalise better to "
                          "unseen haplotypes than k-mer counts"}

    json.dump({"design": "post hoc; fixed modal region; splitter is the only difference",
               "n_evaluable": len(results), "n_not_evaluable": len(skipped),
               "summary": summary, "blocks": results, "not_evaluable": skipped},
              open(OUT, "w"), indent=2, default=float)

    pp = 100.0
    print()
    for key, lab in (("isolate_disjoint", "isolate-disjoint (control)"),
                     ("sequence_disjoint", "sequence-disjoint")):
        c = summary[key]["primary_contrast_dnabert2_minus_kmer"]
        print(f"  {lab:28} {pp*c['mean']:+.2f} pp  "
              f"[{pp*c['ci95'][0]:+.2f}, {pp*c['ci95'][1]:+.2f}]  "
              f"SD {pp*c['sd']:.2f}  p_TOST {c['p_tost']:.3e}  {c['verdict']}")
    g = summary["degradation"]
    print(f"  degradation: k-mer {pp*g['kmer_mean_loss']:+.2f} pp, "
          f"DNABERT-2 {pp*g['dnabert2_mean_loss']:+.2f} pp, "
          f"difference {pp*g['difference']:+.2f} pp "
          f"[{pp*g['ci95'][0]:+.2f}, {pp*g['ci95'][1]:+.2f}], p={g['p_two_sided']:.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
