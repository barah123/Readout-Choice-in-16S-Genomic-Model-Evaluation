"""Readout hierarchy + classical controls across all blocks, parallelised over block x region.

Identical methodology to run_readouts.py (validated in the 2-panel pilot); the changes are
scale (31 blocks), parallelism, and that only regions passing the pre-declared eligibility
gate are scored as confirmatory - ineligible ones are still computed but tagged, never
silently dropped (protocol §19: keep biological-no-signal, insufficient-data, technical-
failure and excluded-by-protocol as distinct codes).
"""
import json
import sys
import warnings
from collections import defaultdict

import numpy as np
from Bio.Align import PairwiseAligner
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 0
N_SPLITS = 5
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]


def balanced_accuracy(y, p):
    y, p = np.asarray(y), np.asarray(p)
    a = [float((p[y == c] == c).mean()) for c in (0, 1) if (y == c).sum()]
    return float(np.mean(a)) if a else float("nan")


def soft_balanced_accuracy(y, credit):
    y, credit = np.asarray(y), np.asarray(credit)
    a = [float(credit[y == c].mean()) for c in (0, 1) if (y == c).sum()]
    return float(np.mean(a)) if a else float("nan")


def kmer_fit_transform(train_seqs, test_seqs, ks=(4, 5, 6)):
    vocab = sorted({(k, s[i:i + k]) for k in ks for s in train_seqs
                    for i in range(len(s) - k + 1) if "N" not in s[i:i + k]})
    index = {kv: j for j, kv in enumerate(vocab)}
    cols_by_k = defaultdict(list)
    for (k, _), j in index.items():
        cols_by_k[k].append(j)

    def enc(seqs):
        M = np.zeros((len(seqs), len(vocab)))
        for r, s in enumerate(seqs):
            for k in ks:
                for i in range(len(s) - k + 1):
                    j = index.get((k, s[i:i + k]))
                    if j is not None:
                        M[r, j] += 1
            for k in ks:
                c = cols_by_k.get(k)
                if c:
                    t = M[r, c].sum()
                    if t:
                        M[r, c] /= t
        return M
    return enc(train_seqs), enc(test_seqs)


def identity_dist(uniq):
    al = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                         match_score=2, mismatch_score=-1)
    n = len(uniq)
    D = np.zeros((n, n))
    for i in range(n):
        a = uniq[i]
        for j in range(i + 1, n):
            b = uniq[j]
            aln = al.align(a, b)[0]
            m = 0
            for tb, qb in zip(aln.aligned[0], aln.aligned[1]):
                m += sum(1 for x, y in zip(a[int(tb[0]):int(tb[1])],
                                           b[int(qb[0]):int(qb[1])]) if x == y)
            d = 1 - 2 * m / (len(a) + len(b))
            D[i, j] = D[j, i] = d
    return D


def geometry_auc(X, y, rng):
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    S = Xn @ Xn.T
    iu = np.triu_indices(len(y), 1)
    same = (y[:, None] == y[None, :])[iu]
    sim = S[iu]
    pos, neg = sim[same], sim[~same]
    if not len(pos) or not len(neg):
        return float("nan")
    if len(pos) > 20000:
        pos = rng.choice(pos, 20000, replace=False)
    if len(neg) > 20000:
        neg = rng.choice(neg, 20000, replace=False)
    ps = np.sort(pos)
    return float(1 - np.searchsorted(ps, neg, side="left").mean() / len(ps))


def nn_credit(D_tr_te, y_tr):
    out = []
    for c in range(D_tr_te.shape[1]):
        d = D_tr_te[:, c]
        tied = y_tr[d <= d.min() + 1e-9]
        out.append(float((tied == 1).mean()))
    return np.array(out)


def eval_cell(block_id, region, slices, y, groups, emb, meta):
    y = np.asarray(y); groups = np.asarray(groups)
    rng = np.random.default_rng(SEED)
    res = {"_block": block_id, "_region": region, **meta}
    if len(np.unique(y)) < 2:
        res["_status"] = "single_class"
        return res
    ns = min(N_SPLITS, int(np.bincount(y).min()), len(np.unique(groups)))
    if ns < 2:
        res["_status"] = "insufficient_groups"
        return res

    uniq = sorted(set(slices))
    pos_of = {s: i for i, s in enumerate(uniq)}
    sp = np.array([pos_of[s] for s in slices])
    try:
        D_ident = identity_dist(uniq)
    except Exception as e:
        res["_status"] = f"identity_failed:{e}"
        return res

    D_emb, Xfull = {}, {}
    for m, lut in emb.items():
        try:
            U = np.vstack([lut[s] for s in uniq])
        except KeyError:
            continue
        Un = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-12)
        D_emb[m] = 1 - Un @ Un.T
        Xfull[m] = U[sp]
        res[f"{m}:R0_geometry_auc"] = geometry_auc(U[sp], y, rng)

    folds = list(StratifiedGroupKFold(n_splits=ns, shuffle=True,
                                      random_state=SEED).split(np.zeros(len(y)), y, groups))
    acc = defaultdict(list)
    for tr, te in folds:
        ytr, yte = y[tr], y[te]
        if len(np.unique(ytr)) < 2:
            continue
        str_tr = [slices[i] for i in tr]
        str_te = [slices[i] for i in te]

        lab = defaultdict(lambda: [0, 0])
        for s, l in zip(str_tr, ytr):
            lab[s][l] += 1
        cr = []
        for s, l in zip(str_te, yte):
            if s in lab:
                a, b = lab[s]
                p1 = b / (a + b)
            else:
                p1 = 0.5
            cr.append(p1 if l == 1 else 1 - p1)
        acc["C_exact"].append(soft_balanced_accuracy(yte, np.array(cr)))

        cp = nn_credit(D_ident[np.ix_(sp[tr], sp[te])], ytr)
        acc["C_align"].append(soft_balanced_accuracy(yte, np.where(yte == 1, cp, 1 - cp)))

        Ktr, Kte = kmer_fit_transform(str_tr, str_te)
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000, class_weight="balanced",
                                               random_state=SEED))
        clf.fit(Ktr, ytr)
        acc["C_kmer"].append(balanced_accuracy(yte, clf.predict(Kte)))

        rp = rng.standard_normal((Ktr.shape[1], 256)) / np.sqrt(max(Ktr.shape[1], 1))
        clf_r = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, class_weight="balanced",
                                                 random_state=SEED))
        clf_r.fit(Ktr @ rp, ytr)
        acc["C_random"].append(balanced_accuracy(yte, clf_r.predict(Kte @ rp)))

        for m in D_emb:
            cpm = nn_credit(D_emb[m][np.ix_(sp[tr], sp[te])], ytr)
            acc[f"{m}:R1_nn"].append(
                soft_balanced_accuracy(yte, np.where(yte == 1, cpm, 1 - cpm)))
            p = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, class_weight="balanced",
                                                 random_state=SEED))
            p.fit(Xfull[m][tr], ytr)
            acc[f"{m}:R2_linear"].append(balanced_accuracy(yte, p.predict(Xfull[m][te])))

    for k, v in acc.items():
        res[k] = float(np.mean(v)) if v else float("nan")
    res["_n_folds"] = len(folds)
    res["_status"] = "ok"
    return res


def main():
    rows = json.load(open("blocks_full.json"))
    summary = json.load(open("block_summary_full.json"))

    emb_store = {}
    for m in MODELS:
        try:
            emb_store[m] = np.load(f"emb_full_{m}.npz", allow_pickle=True)
        except FileNotFoundError:
            print(f"[skip] emb_full_{m}.npz missing", file=sys.stderr)
    if not emb_store:
        sys.exit("no embedding caches")

    by = defaultdict(list)
    for r in rows:
        by[(r["block_id"], r["region"])].append(r)

    lut_cache = {}
    for m, z in emb_store.items():
        for region in {k.split("__")[0] for k in z.files if "__keys" in k}:
            lut_cache[(m, region)] = {k: v for k, v in
                                      zip(z[f"{region}__keys"], z[f"{region}__vecs"])}

    jobs = []
    for (bid, region), rs in sorted(by.items()):
        # Restrict each cell's lookup to just the slices that cell uses. Passing the whole
        # per-region table would re-pickle several MB per model to every worker, per cell.
        need = {r["slice"] for r in rs}
        emb = {}
        for m in emb_store:
            full = lut_cache.get((m, region))
            if full is None:
                continue
            sub = {s: full[s] for s in need if s in full}
            if sub:
                emb[m] = sub
        if not emb:
            continue
        binfo = summary.get(bid, {}).get("regions", {}).get(region, {})
        meta = {
            "_ceiling_bal": binfo.get("ceiling_bal_isolate"),
            "_target_isolates": binfo.get("target_isolates"),
            "_eligible": binfo.get("eligible", False),
            "_family": summary.get(bid, {}).get("family"),
            "_rank": summary.get(bid, {}).get("rank"),
            "_expected_class": summary.get(bid, {}).get("expected_class"),
        }
        jobs.append((bid, region, [r["slice"] for r in rs], [r["is_target"] for r in rs],
                     [r["isolate_group"] for r in rs], emb, meta))

    print(f"{len(jobs)} block x region cells to evaluate", flush=True)
    out = Parallel(n_jobs=6, verbose=10)(delayed(eval_cell)(*j) for j in jobs)
    json.dump({f"{r['_block']}|{r['_region']}": r for r in out},
              open("readout_results_full.json", "w"), indent=2)
    print("wrote readout_results_full.json")


if __name__ == "__main__":
    main()
