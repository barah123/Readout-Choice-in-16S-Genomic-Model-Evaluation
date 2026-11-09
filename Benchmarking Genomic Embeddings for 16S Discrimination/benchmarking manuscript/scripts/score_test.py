"""Phase 4 step 4: produce locked-test predictions under the FROZEN plan.

Implements FREEZE.md §3 exactly:
  - outer 5-fold stratified split over ISOLATES, defined ONCE per block and applied
    consistently across regions (so region selection and evaluation use the same split);
  - within each outer fold, the region is chosen by INNER grouped CV on TRAINING ISOLATES
    ONLY, scored with the classical arm (C_kmer) — no query labels touch the choice;
  - the chosen region is then applied IDENTICALLY to every arm;
  - arms: DNABERT-2:R2 (primary model arm), DNABERT-S/HyenaDNA/NTv2:R2 (secondary),
    C_kmer (primary comparator), C_exact / C_align / C_random (reference).

Writes per-block arm scores plus the full per-region exploratory cell matrix.
This script COMPUTES predictions; it does not compute the primary contrast or any
equivalence test — those are Phase 5's prespecified summaries.
"""
import json
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
import run_readouts_full as rr   # reuse validated kmer/identity/metric helpers

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 0            # frozen
N_SPLITS = 5        # frozen
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
REGIONS = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"]


def probe():
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED))


def kmer_score_inner(slices, y, n_splits=3):
    """Inner CV score for region selection: C_kmer only, training isolates only."""
    y = np.asarray(y)
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < n_splits:
        return float("nan")
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    out = []
    for tr, te in skf.split(np.zeros(len(y)), y):
        Ktr, Kte = rr.kmer_fit_transform([slices[i] for i in tr], [slices[i] for i in te])
        c = probe(); c.fit(Ktr, y[tr])
        out.append(rr.balanced_accuracy(y[te], c.predict(Kte)))
    return float(np.mean(out)) if out else float("nan")


def score_block(block_id, per_region, emb, meta):
    """per_region: {region: {isolate: (slice, label)}}"""
    # union of isolates with a consistent label
    labels = {}
    for region, d in per_region.items():
        for iso, (_, lab) in d.items():
            labels[iso] = lab
    isolates = sorted(labels)
    y_all = np.array([labels[i] for i in isolates])
    res = {"_block": block_id, **meta}
    if len(np.unique(y_all)) < 2 or min(np.bincount(y_all)) < N_SPLITS:
        res["_status"] = "insufficient_for_outer_split"
        return res

    eligible_regions = [r for r in REGIONS
                        if meta.get("_eligible_regions") and r in meta["_eligible_regions"]]
    if not eligible_regions:
        res["_status"] = "no_eligible_region"
        return res

    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    arm_scores = defaultdict(list)
    chosen_regions = []

    for tr_idx, te_idx in skf.split(np.zeros(len(isolates)), y_all):
        tr_iso = {isolates[i] for i in tr_idx}
        te_iso = {isolates[i] for i in te_idx}

        # ---- region selection: training isolates only, classical arm only
        best_r, best_s = None, -np.inf
        for r in eligible_regions:
            d = per_region[r]
            sl = [d[i][0] for i in isolates if i in tr_iso and i in d]
            yy = [d[i][1] for i in isolates if i in tr_iso and i in d]
            if len(set(yy)) < 2:
                continue
            s = kmer_score_inner(sl, yy)
            if s == s and s > best_s:
                best_s, best_r = s, r
        if best_r is None:
            continue
        chosen_regions.append(best_r)

        # ---- evaluate every arm at the chosen region, same split
        d = per_region[best_r]
        tr = [i for i in isolates if i in tr_iso and i in d]
        te = [i for i in isolates if i in te_iso and i in d]
        ytr = np.array([d[i][1] for i in tr])
        yte = np.array([d[i][1] for i in te])
        if len(set(ytr)) < 2 or len(set(yte)) < 2:
            continue
        str_tr = [d[i][0] for i in tr]
        str_te = [d[i][0] for i in te]

        # C_kmer (primary comparator)
        Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
        c = probe(); c.fit(Ktr, ytr)
        arm_scores["C_kmer"].append(rr.balanced_accuracy(yte, c.predict(Kte)))

        # C_random (flexibility reference)
        rng = np.random.default_rng(SEED)
        rp = rng.standard_normal((Ktr.shape[1], 256)) / np.sqrt(max(Ktr.shape[1], 1))
        cr = probe(); cr.fit(Ktr @ rp, ytr)
        arm_scores["C_random"].append(rr.balanced_accuracy(yte, cr.predict(Kte @ rp)))

        # C_exact
        lab = defaultdict(lambda: [0, 0])
        for s, l in zip(str_tr, ytr):
            lab[s][l] += 1
        cred = []
        for s, l in zip(str_te, yte):
            if s in lab:
                a, b = lab[s]; p1 = b / (a + b)
            else:
                p1 = 0.5
            cred.append(p1 if l == 1 else 1 - p1)
        arm_scores["C_exact"].append(rr.soft_balanced_accuracy(yte, np.array(cred)))

        # C_align
        uniq = sorted(set(str_tr) | set(str_te))
        pos = {s: k for k, s in enumerate(uniq)}
        D = rr.identity_dist(uniq)
        sub = D[np.ix_([pos[s] for s in str_tr], [pos[s] for s in str_te])]
        cp = rr.nn_credit(sub, ytr)
        arm_scores["C_align"].append(
            rr.soft_balanced_accuracy(yte, np.where(yte == 1, cp, 1 - cp)))

        # model R2 probes
        for m in MODELS:
            lut = emb.get((m, best_r))
            if lut is None:
                continue
            try:
                Xtr = np.vstack([lut[s] for s in str_tr])
                Xte = np.vstack([lut[s] for s in str_te])
            except KeyError:
                continue
            p = probe(); p.fit(Xtr, ytr)
            arm_scores[f"{m}:R2"].append(rr.balanced_accuracy(yte, p.predict(Xte)))

    for k, v in arm_scores.items():
        res[k] = float(np.mean(v)) if v else float("nan")
    res["_chosen_regions"] = chosen_regions
    res["_n_outer_folds"] = len(chosen_regions)
    res["_status"] = "ok" if chosen_regions else "no_usable_fold"
    return res


def main():
    rows = json.load(open("test_blocks.json"))
    summ = json.load(open("test_block_summary.json"))

    emb_files = {m: f"emb_test_{m}.npz" for m in MODELS}
    emb = {}
    for m, f in emb_files.items():
        if not os.path.exists(f):
            print(f"[skip] {f} missing", file=sys.stderr)
            continue
        z = np.load(f, allow_pickle=True)
        for r in REGIONS:
            if f"{r}__keys" in z:
                emb[(m, r)] = {k: v for k, v in zip(z[f"{r}__keys"], z[f"{r}__vecs"])}
    if not emb:
        sys.exit("no test embedding caches found")

    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])

    jobs = []
    for bid, per_region in by_block.items():
        info = summ.get(bid, {})
        meta = {"_family": info.get("family"), "_rank": info.get("rank"),
                "_expected_class": info.get("expected_class"),
                "_eligible_regions": info.get("eligible_regions", []),
                "_block_eligible": info.get("block_eligible", False),
                "_ceilings": {r: v.get("ceiling_bal_isolate")
                              for r, v in info.get("regions", {}).items()}}
        jobs.append((bid, per_region, emb, meta))

    print(f"scoring {len(jobs)} locked-test blocks under the frozen plan", flush=True)
    out = Parallel(n_jobs=4, verbose=10)(delayed(score_block)(*j) for j in jobs)
    json.dump({r["_block"]: r for r in out}, open("locked_test_scores.json", "w"),
              indent=2, default=float)
    ok = sum(1 for r in out if r.get("_status") == "ok" and r.get("_block_eligible"))
    print(f"\nwrote locked_test_scores.json  ({ok} eligible blocks scored)")
    print("Primary contrast and equivalence test are Phase 5, not computed here.")


if __name__ == "__main__":
    main()
