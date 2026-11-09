"""Extension B: does the Phase 5 equivalence survive realistic, imperfect sequences?

Sensitivity analysis of the Phase 5 primary contrast (DNABERT-2:R2 vs C_kmer) on the same
20 locked-test panels, with the extracted sequences deliberately degraded.

Both arms receive IDENTICAL perturbed sequences in every condition, so no arm is advantaged.
Perturbation seeds are fixed. EXPLORATORY: this re-uses already-unsealed data to test the
robustness of an existing conclusion, not to establish a new one.

Conditions (declared in PREREGISTRATION.md):
  clean      - baseline, reproduces Phase 5
  err_0.1%   - high-quality Illumina substitution rate
  err_1%     - poor-quality run
  trim_10    - imperfect primer binding, 10 bases lost each end
  copy_var   - intragenomic copy variation: 1-3 substitutions per isolate
  combined   - err_1% + trim_10
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from transformers import AutoModel, AutoTokenizer

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
import run_readouts_full as rr

SEED = 0
MODEL_DIR = "../../glm_prototype_bundle (1)/dnabert2_local"
BASES = np.array(list("ACGT"))


def perturb(seq, condition, rng):
    s = list(seq)
    if condition in ("err_0.1%", "err_1%", "combined"):
        rate = 0.001 if condition == "err_0.1%" else 0.01
        n = rng.binomial(len(s), rate)
        for pos in rng.choice(len(s), size=n, replace=False) if n else []:
            cur = s[pos]
            alt = [b for b in "ACGT" if b != cur]
            s[pos] = alt[rng.integers(0, 3)]
    if condition == "copy_var":
        n = int(rng.integers(1, 4))
        n = min(n, len(s))
        for pos in rng.choice(len(s), size=n, replace=False):
            cur = s[pos]
            alt = [b for b in "ACGT" if b != cur]
            s[pos] = alt[rng.integers(0, 3)]
    out = "".join(s)
    if condition in ("trim_10", "combined"):
        out = out[10:-10] if len(out) > 25 else out
    return out


def probe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=2000, class_weight="balanced",
                                            random_state=SEED))


def main():
    rows = json.load(open("../phase4_locked_test/test_blocks.json"))
    scores = json.load(open("../phase4_locked_test/locked_test_scores.json"))
    blocks = {k: v for k, v in scores.items()
              if v.get("_status") == "ok" and v.get("_block_eligible")}

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    tok = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
    model = AutoModel.from_pretrained(MODEL_DIR, trust_remote_code=True)
    model.eval()
    for _n, mod in model.named_modules():
        if hasattr(mod, "rebuild_alibi_tensor") and getattr(mod, "alibi", "x") is None:
            mod.rebuild_alibi_tensor(size=mod._current_alibi_size)

    def embed(seqs, batch_size=16):
        outs = []
        for i in range(0, len(seqs), batch_size):
            b = seqs[i:i + batch_size]
            enc = tok(b, return_tensors="pt", padding=True, truncation=True, max_length=256)
            with torch.no_grad():
                h = model(**enc)[0]
            am = enc["attention_mask"].bool()
            m = am.clone()
            L = am.sum(1)
            m[:, 0] = False
            m[torch.arange(m.size(0)), (L - 1).clamp(min=0)] = False
            empty = m.sum(1) == 0
            if empty.any():
                m[empty] = am[empty]
            mf = m.unsqueeze(-1).float()
            outs.append(((h * mf).sum(1) / mf.sum(1).clamp(min=1)).float().numpy())
        return np.vstack(outs)

    conditions = ["clean", "err_0.1%", "err_1%", "trim_10", "copy_var", "combined"]
    results = defaultdict(dict)

    for cond in conditions:
        deltas, kmers, d2s = [], [], []
        for bid, v in blocks.items():
            region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
            d = by_block[bid][region]
            isolates = sorted(d)
            y = np.array([d[i][1] for i in isolates])
            if len(np.unique(y)) < 2 or min(np.bincount(y)) < 5:
                continue
            rng = np.random.default_rng(abs(hash((bid, cond))) % (2**32))
            slices = [perturb(d[i][0], cond, rng) for i in isolates]

            uniq = sorted(set(slices))
            E = embed(uniq)
            lut = {s: E[j] for j, s in enumerate(uniq)}

            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
            ks, ds = [], []
            for tr, te in skf.split(np.zeros(len(y)), y):
                ytr, yte = y[tr], y[te]
                if len(set(ytr)) < 2 or len(set(yte)) < 2:
                    continue
                str_tr = [slices[i] for i in tr]
                str_te = [slices[i] for i in te]
                Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
                c = probe(); c.fit(Ktr, ytr)
                ks.append(rr.balanced_accuracy(yte, c.predict(Kte)))
                Xtr = np.vstack([lut[s] for s in str_tr])
                Xte = np.vstack([lut[s] for s in str_te])
                p = probe(); p.fit(Xtr, ytr)
                ds.append(rr.balanced_accuracy(yte, p.predict(Xte)))
            if not ks:
                continue
            km, dm = float(np.mean(ks)), float(np.mean(ds))
            kmers.append(km); d2s.append(dm); deltas.append(dm - km)
            results[cond][bid] = dict(C_kmer=km, dnabert2=dm, delta=dm - km)

        deltas = np.array(deltas)
        rngb = np.random.default_rng(1729)
        dr = deltas[rngb.integers(0, len(deltas), size=(2000, len(deltas)))].mean(axis=1)
        lo, hi = np.percentile(dr, [2.5, 97.5])
        results[cond]["_summary"] = dict(
            n=len(deltas), kmer_mean=float(np.mean(kmers)),
            dnabert2_mean=float(np.mean(d2s)), delta_mean=float(deltas.mean()),
            ci=[float(lo), float(hi)])
        print(f"{cond:11} n={len(deltas):3}  C_kmer={np.mean(kmers):.4f}  "
              f"DNABERT-2={np.mean(d2s):.4f}  delta={deltas.mean():+.4f} "
              f"[{lo:+.4f}, {hi:+.4f}]", flush=True)

    json.dump(results, open("robustness_results.json", "w"), indent=2, default=float)

    print()
    print("=" * 92)
    print("How much each arm degrades relative to its own clean baseline")
    print("=" * 92)
    base = results["clean"]["_summary"]
    print(f"{'condition':12}{'C_kmer drop':>14}{'DNABERT-2 drop':>17}{'delta shift':>14}")
    for cond in conditions:
        s = results[cond]["_summary"]
        print(f"{cond:12}{s['kmer_mean']-base['kmer_mean']:>+14.4f}"
              f"{s['dnabert2_mean']-base['dnabert2_mean']:>+17.4f}"
              f"{s['delta_mean']-base['delta_mean']:>+14.4f}")
    print("\nwrote robustness_results.json")


if __name__ == "__main__":
    main()
