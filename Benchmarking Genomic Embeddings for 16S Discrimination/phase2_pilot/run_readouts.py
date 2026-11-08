"""Phase 2 step 5: the readout hierarchy (protocol §10) and classical controls (§11),
evaluated under isolate-grouped CV (§12.1) with matched ceilings (§8).

Readouts implemented
  R0  intrinsic geometry : within- vs between-class cosine separation AUC (no supervision)
  R1  reference lookup   : tie-aware 1-NN on frozen embeddings, isolate-grouped CV
  R2  linear probe       : regularized logistic regression, isolate-grouped CV
Classical controls
  C_exact   : exact slice lookup with explicit ambiguity (majority label of matching slice)
  C_align   : tie-aware 1-NN on pairwise sequence identity
  C_kmer    : canonical k-mer composition + the SAME linear probe used for R2
  C_random  : random projection of matched dim + same probe (flexibility reference, §11)

Design decisions (declared before results):
- CV unit is the isolate group, via StratifiedGroupKFold - no isolate's records ever straddle
  a fold. This is strictly stronger than Phase 0-3a's leave-one-haplotype-out.
- k-mer vocabulary and any scaler are fit on TRAINING FOLDS ONLY, fixing the transductive
  leak found in Phase 1's audit (audit_report.md finding 7).
- Primary metric: species-balanced accuracy, isolate-weighted (each isolate contributes its
  modal slice once), matching the ceiling_bal_isolate computed in build_blocks.py.
- Chance is 0.5 by construction for balanced accuracy; ceilings come from block_summary.json.
"""
import json
import sys
import warnings
from collections import defaultdict

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=ConvergenceWarning)

SEED = 0
N_SPLITS = 5
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]


# ---------------------------------------------------------------- metrics
def balanced_accuracy(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    accs = []
    for cls in (0, 1):
        m = y_true == cls
        if m.sum():
            accs.append(float((y_pred[m] == cls).mean()))
    return float(np.mean(accs)) if accs else float("nan")


def soft_balanced_accuracy(y_true, credit):
    """credit[i] in [0,1] = fraction of tied nearest neighbours matching the true label."""
    y_true, credit = np.asarray(y_true), np.asarray(credit)
    accs = []
    for cls in (0, 1):
        m = y_true == cls
        if m.sum():
            accs.append(float(credit[m].mean()))
    return float(np.mean(accs)) if accs else float("nan")


# ---------------------------------------------------------------- features
def kmer_fit_transform(train_seqs, test_seqs, ks=(4, 5, 6)):
    """Vocabulary fit on TRAIN ONLY (Phase 1 audit finding 7)."""
    vocab = sorted({(k, s[i:i + k]) for k in ks for s in train_seqs
                    for i in range(len(s) - k + 1) if "N" not in s[i:i + k]})
    index = {kv: j for j, kv in enumerate(vocab)}

    def encode(seqs):
        M = np.zeros((len(seqs), len(vocab)))
        for r, s in enumerate(seqs):
            for k in ks:
                for i in range(len(s) - k + 1):
                    j = index.get((k, s[i:i + k]))
                    if j is not None:
                        M[r, j] += 1
            for k in ks:
                cols = [j for (kk, _), j in index.items() if kk == k]
                if cols:
                    tot = M[r, cols].sum()
                    if tot:
                        M[r, cols] /= tot
        return M

    return encode(train_seqs), encode(test_seqs)


def identity_matrix(seqs):
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                         match_score=2, mismatch_score=-1)
    uniq = sorted(set(seqs))
    idx = {s: i for i, s in enumerate(uniq)}
    n = len(uniq)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = uniq[i], uniq[j]
            aln = al.align(a, b)[0]
            matches = 0
            for tb, qb in zip(aln.aligned[0], aln.aligned[1]):
                matches += sum(1 for x, y in zip(a[int(tb[0]):int(tb[1])],
                                                 b[int(qb[0]):int(qb[1])]) if x == y)
            d = 1 - 2 * matches / (len(a) + len(b))
            D[i, j] = D[j, i] = d
    return D, idx


# ---------------------------------------------------------------- readouts
def r0_geometry_auc(X, y, rng):
    """Within- vs between-class cosine separation, pair-level AUC. No supervision.
    Pairs are dependent; this is a descriptive geometry statistic only (protocol §13)."""
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    S = Xn @ Xn.T
    n = len(y)
    iu = np.triu_indices(n, 1)
    same = (y[:, None] == y[None, :])[iu]
    sim = S[iu]
    pos, neg = sim[same], sim[~same]
    if not len(pos) or not len(neg):
        return float("nan")
    # subsample for tractability, fixed seed
    if len(pos) > 20000:
        pos = rng.choice(pos, 20000, replace=False)
    if len(neg) > 20000:
        neg = rng.choice(neg, 20000, replace=False)
    # AUC = P(within-class similarity > between-class similarity)
    pos_sorted = np.sort(pos)
    ranks = np.searchsorted(pos_sorted, neg, side="left")
    auc = 1 - ranks.mean() / len(pos_sorted)
    return float(auc)


def nn_credit(D_train_test, y_train):
    """Tie-aware 1-NN: fraction of tied nearest training neighbours with each label."""
    credit_pos = []
    for col in range(D_train_test.shape[1]):
        d = D_train_test[:, col]
        mn = d.min()
        tied = y_train[d <= mn + 1e-9]
        credit_pos.append(float((tied == 1).mean()))
    return np.array(credit_pos)


def evaluate_region(slices, y, groups, emb_lookup, rng):
    """Returns {method: balanced_accuracy} for one panel x region."""
    y = np.asarray(y)
    groups = np.asarray(groups)
    res = {}
    n_classes = len(np.unique(y))
    if n_classes < 2:
        return res

    n_splits = min(N_SPLITS, np.bincount(y).min(), len(np.unique(groups)))
    if n_splits < 2:
        return res
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    folds = list(sgkf.split(np.zeros(len(y)), y, groups))

    # ---- geometry (R0): no CV needed, it is unsupervised and descriptive
    for mname, M in emb_lookup.items():
        X = np.vstack([M[s] for s in slices])
        res[f"{mname}:R0_geometry_auc"] = r0_geometry_auc(X, y, rng)

    # ---- distance-based R1 / C_align need a precomputed distance over unique slices
    uniq = sorted(set(slices))
    pos_of = {s: i for i, s in enumerate(uniq)}
    slice_pos = np.array([pos_of[s] for s in slices])

    D_ident, _ = identity_matrix(uniq)
    D_emb = {}
    for mname, M in emb_lookup.items():
        U = np.vstack([M[s] for s in uniq])
        Un = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-12)
        D_emb[mname] = 1 - Un @ Un.T

    acc = defaultdict(list)
    for tr, te in folds:
        ytr, yte = y[tr], y[te]
        if len(np.unique(ytr)) < 2:
            continue
        str_tr = [slices[i] for i in tr]
        str_te = [slices[i] for i in te]

        # --- C_exact: majority label of an identical training slice; abstain -> 0.5 credit
        lab_of = defaultdict(lambda: [0, 0])
        for s, lab in zip(str_tr, ytr):
            lab_of[s][lab] += 1
        cr = []
        for s, lab in zip(str_te, yte):
            if s in lab_of:
                a, b = lab_of[s]
                p1 = b / (a + b)
            else:
                p1 = 0.5
            cr.append(p1 if lab == 1 else 1 - p1)
        acc["C_exact"].append(soft_balanced_accuracy(yte, np.array(cr)))

        # --- C_align: tie-aware 1-NN on identity distance
        sub = D_ident[np.ix_(slice_pos[tr], slice_pos[te])]
        cp = nn_credit(sub, ytr)
        acc["C_align"].append(soft_balanced_accuracy(yte, np.where(yte == 1, cp, 1 - cp)))

        # --- C_kmer: train-only vocab + linear probe
        Ktr, Kte = kmer_fit_transform(str_tr, str_te)
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000, class_weight="balanced",
                                               random_state=SEED))
        clf.fit(Ktr, ytr)
        acc["C_kmer"].append(balanced_accuracy(yte, clf.predict(Kte)))

        # --- C_random: matched-dim random projection of the k-mer space
        rp = rng.standard_normal((Ktr.shape[1], 256)) / np.sqrt(Ktr.shape[1])
        clf_r = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, class_weight="balanced",
                                                 random_state=SEED))
        clf_r.fit(Ktr @ rp, ytr)
        acc["C_random"].append(balanced_accuracy(yte, clf_r.predict(Kte @ rp)))

        # --- per-model R1 and R2
        for mname, M in emb_lookup.items():
            subm = D_emb[mname][np.ix_(slice_pos[tr], slice_pos[te])]
            cpm = nn_credit(subm, ytr)
            acc[f"{mname}:R1_nn"].append(
                soft_balanced_accuracy(yte, np.where(yte == 1, cpm, 1 - cpm)))

            Xtr = np.vstack([M[s] for s in str_tr])
            Xte = np.vstack([M[s] for s in str_te])
            p = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, class_weight="balanced",
                                                 random_state=SEED))
            p.fit(Xtr, ytr)
            acc[f"{mname}:R2_linear"].append(balanced_accuracy(yte, p.predict(Xte)))

    for k, v in acc.items():
        res[k] = float(np.mean(v)) if v else float("nan")
    res["_n_folds"] = len(folds)
    return res


def main():
    rows = json.load(open("blocks.json"))
    summary = json.load(open("block_summary.json"))
    rng = np.random.default_rng(SEED)

    emb = {}
    for m in MODELS:
        try:
            z = np.load(f"emb_{m}.npz", allow_pickle=True)
        except FileNotFoundError:
            print(f"[skip] emb_{m}.npz not found", file=sys.stderr)
            continue
        emb[m] = z
    if not emb:
        sys.exit("no embedding caches found")

    by_pr = defaultdict(list)
    for r in rows:
        by_pr[(r["panel"], r["region"])].append(r)

    out = {}
    for (panel, region), rs in sorted(by_pr.items()):
        lookup = {}
        ok = True
        for m, z in emb.items():
            kk, vv = f"{region}__keys", f"{region}__vecs"
            if kk not in z:
                ok = False
                break
            lookup[m] = {k: v for k, v in zip(z[kk], z[vv])}
        if not ok:
            continue

        slices = [r["slice"] for r in rs]
        y = [r["is_target"] for r in rs]
        groups = [r["isolate_group"] for r in rs]
        res = evaluate_region(slices, y, groups, lookup, rng)
        if not res:
            continue
        res["_ceiling_bal"] = summary[panel][region]["ceiling_bal_isolate"]
        res["_n_target_isolates"] = summary[panel][region]["target_isolates"]
        out[f"{panel}|{region}"] = res
        print(f"done {panel} {region}", flush=True)

    json.dump(out, open("readout_results.json", "w"), indent=2)
    print("wrote readout_results.json")


if __name__ == "__main__":
    main()
