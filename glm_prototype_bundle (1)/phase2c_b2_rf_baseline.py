"""Phase 2c: the real B2 baseline (k-mer + random forest, plan §4.2), plus a
leakage-controlled re-evaluation of the region-level metrics from Phases 2/2b.

Why this script exists
----------------------
MASTER_PLAN §4.1 called for replacing the cheap raw-distance k-mer nearest-neighbour
proxy used in Phase 2b with the random-forest classifier the plan actually specifies.
While setting that up, the per-region haplotype structure turned out to invalidate the
region-level LOO-NN metric used in BOTH Phase 2 and Phase 2b:

  * The anthracis class collapses to 1-3 unique slice haplotypes per region (V1 and V8:
    all 9 anthracis sequences are one identical string).
  * Every region contains haplotypes shared ACROSS the class boundary, so the label is
    not a function of the slice. This caps achievable accuracy below 1 for any method.
  * 81-89% of samples therefore have TIED nearest neighbours at distance 0, and the
    reported LOO-NN accuracy is determined by argmin/argmax index order.
  * In V1, V4 and V7 the Phase 2 embedding accuracy EXCEEDS the Bayes-optimal accuracy,
    which is only possible via that tie-breaking artifact.

So this script does three things, in order:
  1. ORACLE CEILING  - per region, the best accuracy/balanced-accuracy any classifier
     could achieve given cross-class-identical slices. Reported as a first-class result.
  2. ERRATUM         - recompute LOO-NN under explicit tie policies to quantify how
     unstable the Phase 2 / 2b region rankings are.
  3. REAL B2 + HEAD-TO-HEAD - random forest on k-mer features vs random forest on
     DNABERT-2 embeddings, under IDENTICAL leave-one-haplotype-out CV so that
     string-identical slices never straddle a train/test split. Plus a label-permutation
     null, because with this few haplotypes a point estimate alone is uninterpretable.
     B3 (per-site Shannon entropy) and B4 (length / variable-site density) are added as
     formal region-level rankers, as §4.1 item 1 asks.

Slicing/QA logic is copied verbatim from embed_and_score.py (the canonical version) per
the CLAUDE.md convention, so every row stays joinable by (region, id).
"""
import json
import itertools
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio.Align import PairwiseAligner
from sklearn.ensemble import RandomForestClassifier
from joblib import Parallel, delayed

RNG_SEED = 0
rng_global = np.random.default_rng(RNG_SEED)

records = json.load(open("cereus_group_16s.json"))
ec_ref = json.load(open("ecoli_16s.txt"))["seq"]
aligner = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                          match_score=2, mismatch_score=-1)
REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
           "V8": (1243, 1294), "V9": (1435, 1465)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}
CEREUS_GROUP = ("B_cereus", "B_thuringiensis", "B_weihenstephanensis")


def slice_regions(seq):
    aln = aligner.align(ec_ref, seq)[0]
    ref_to_query = {}
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1 = int(t_blk[0]), int(t_blk[1]); q0, q1 = int(q_blk[0]), int(q_blk[1])
        for i in range(t0, t1):
            ref_to_query[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in REGIONS.items():
        qs = [ref_to_query[p] for p in range(s - 1, e) if p in ref_to_query]
        sl = seq[min(qs):max(qs) + 1] if qs else ""
        out[name] = sl if sl and 0.75 * EXPECTED[name] <= len(sl) <= 1.30 * EXPECTED[name] else None
    return out


def identity(a, b):
    """2*matches/(len_a+len_b). Matching positions, NOT aligned positions - see the
    Phase 2 erratum in prototype_report.md §2."""
    if not a or not b:
        return None
    aln = aligner.align(a, b)[0]
    matches = 0
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        seg_t = a[int(t_blk[0]):int(t_blk[1])]
        seg_q = b[int(q_blk[0]):int(q_blk[1])]
        matches += sum(1 for x, y in zip(seg_t, seg_q) if x == y)
    return 2 * matches / (len(a) + len(b))


# ---- same QA as Phases 0-2b ----
seqs, seen = [], set()
for grp, hdr, seq in records:
    if grp == "E_coli_ref" or seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    seqs.append((grp, hdr.split()[0], slice_regions(seq)))
print(f"{len(seqs)} sequences after QA (expect 43, matching Phases 0-2b)", flush=True)

# ---- MANDATORY identical-sequence QA check (CLAUDE.md convention) ----
_a = seqs[0][2]["V4"]
_b = seqs[1][2]["V4"]
assert abs(identity(_a, _a) - 1.0) < 1e-12, "identity(x,x) != 1.0 - scoring function is broken"
if _a != _b:
    assert identity(_a, _b) < 1.0, "identity() returned 1.0 for non-identical slices"
print("identity() identical-sequence QA check: PASS", flush=True)

groups = np.array([s[0] for s in seqs])
ids = np.array([s[1] for s in seqs])

# ---- cached DNABERT-2 embeddings ----
# Uses the Phase 2d CORRECTED cache. The original dnabert2_embeddings.npz leaks the class
# label through batch composition (phase2d_report.md §2) and must not be used for scoring.
z = np.load("dnabert2_embeddings_fixed.npz", allow_pickle=True)
assert list(z["ids"]) == list(ids), "embedding cache order does not match QA'd sequence order"
print("embedding cache alignment check: PASS", flush=True)


def kmer_matrix(slices, ks=(4, 5, 6)):
    """Dense k-mer frequency matrix over the vocabulary present in `slices`.
    Frequencies normalised within each k separately, then concatenated."""
    vocab = sorted({(k, s[i:i + k])
                    for k in ks for s in slices
                    for i in range(len(s) - k + 1) if "N" not in s[i:i + k]})
    index = {kv: j for j, kv in enumerate(vocab)}
    M = np.zeros((len(slices), len(vocab)))
    for r, s in enumerate(slices):
        for k in ks:
            n = len(s) - k + 1
            if n <= 0:
                continue
            for i in range(n):
                km = s[i:i + k]
                if "N" in km:
                    continue
                M[r, index[(k, km)]] += 1
    # normalise within each k block separately, then concatenate
    for k in ks:
        cols = [j for j, (kk, _) in enumerate(vocab) if kk == k]
        if not cols:
            continue
        block = M[:, cols]
        tot = block.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1.0
        M[:, cols] = block / tot
    return M, vocab


def per_site_entropy(slices):
    """B3: mean per-column Shannon entropy over a gapless left-anchored stack.
    Slices within a region differ in length by a few nt; columns are taken up to the
    shortest slice so no padding character inflates entropy."""
    if not slices:
        return float("nan")
    L = min(len(s) for s in slices)
    if L == 0:
        return float("nan")
    ents = []
    for c in range(L):
        col = [s[c] for s in slices]
        n = len(col)
        counts = {}
        for b in col:
            counts[b] = counts.get(b, 0) + 1
        ps = np.array([v / n for v in counts.values()])
        ents.append(float(-(ps * np.log2(ps)).sum()))
    return float(np.mean(ents))


def variable_site_density(slices):
    """B4 component: fraction of columns (up to shortest slice) that are polymorphic."""
    if not slices:
        return float("nan")
    L = min(len(s) for s in slices)
    if L == 0:
        return float("nan")
    var = sum(1 for c in range(L) if len({s[c] for s in slices}) > 1)
    return var / L


def balanced_accuracy(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    accs = []
    for cls in (0, 1):
        m = y_true == cls
        if m.sum() == 0:
            continue
        accs.append((y_pred[m] == cls).mean())
    return float(np.mean(accs)) if accs else float("nan")


def loho_cv(X, y, hap, seed=0, n_estimators=500):
    """Leave-one-haplotype-out CV. All sequences sharing a slice haplotype are held out
    together, so string-identical slices can never straddle a train/test split.

    Returns (balanced_accuracy, n_predicted) or (nan, 0) if not estimable.
    A haplotype whose label is ambiguous (appears in both classes) is kept - excluding it
    would inflate the score by removing exactly the hard cases.
    """
    hap = np.asarray(hap)
    uniq = np.unique(hap)
    # need >=2 haplotypes in each class among the TRAINING remainder for every fold
    preds = np.full(len(y), -1)
    for h in uniq:
        test = hap == h
        train = ~test
        if len(np.unique(y[train])) < 2:
            return float("nan"), 0, "train fold lost a class"
        clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed,
                                     class_weight="balanced", n_jobs=1)
        clf.fit(X[train], y[train])
        preds[test] = clf.predict(X[test])
    if (preds < 0).any():
        return float("nan"), 0, "unpredicted rows"
    return balanced_accuracy(y, preds), int(len(y)), "ok"


def nn_expected_acc(D, labels):
    """Tie-aware LOO nearest-neighbour: expected balanced accuracy averaged over uniform
    random tie-breaking, rather than whichever neighbour argmin happens to return."""
    labels = np.asarray(labels)
    n = len(labels)
    per = np.zeros(n)
    for i in range(n):
        d = D[i].copy()
        d[i] = np.inf
        mn = d.min()
        tied = np.where(d <= mn + 1e-9)[0]
        per[i] = (labels[tied] == labels[i]).mean()
    accs = []
    for cls in (0, 1):
        m = labels == cls
        if m.sum():
            accs.append(per[m].mean())
    return float(np.mean(accs))


# ======================================================================
results = {}
print("\n" + "=" * 100)
print("1. ORACLE CEILING - best achievable given cross-class-identical slices")
print("=" * 100)
print(f"{'reg':5} {'N':>3} {'n_an':>5} {'n_cg':>5} {'hap_an':>7} {'hap_cg':>7} {'hap_ambig':>10} "
      f"{'oracle_acc':>11} {'oracle_bal':>11}")
for name in REGIONS:
    an = [i for i in range(len(seqs)) if groups[i] == "B_anthracis" and seqs[i][2][name]]
    cg = [i for i in range(len(seqs)) if groups[i] in CEREUS_GROUP and seqs[i][2][name]]
    if not an or not cg:
        continue
    sl = {i: seqs[i][2][name] for i in an + cg}
    hap = defaultdict(lambda: [0, 0])
    for i in an:
        hap[sl[i]][0] += 1
    for i in cg:
        hap[sl[i]][1] += 1
    N = len(an) + len(cg)
    # Ceiling for PLAIN accuracy: per haplotype predict the majority label.
    err = sum(min(a, c) for a, c in hap.values())
    oracle_acc = 1 - err / N
    # Ceiling for BALANCED accuracy: a different rule maximises it. Each anthracis sequence
    # is worth 1/(2*n_an) and each cereus-group sequence 1/(2*n_cg), so for an ambiguous
    # haplotype the balanced-optimal call can be the *minority* label. Using the majority
    # rule here (as an earlier draft did) understates the ceiling.
    oracle_bal = sum(max(a / (2 * len(an)), c / (2 * len(cg))) for a, c in hap.values())
    n_hap_an = len({sl[i] for i in an})
    n_hap_cg = len({sl[i] for i in cg})
    n_ambig = sum(1 for a, c in hap.values() if a > 0 and c > 0)
    results[name] = dict(n=N, n_an=len(an), n_cg=len(cg), hap_an=n_hap_an, hap_cg=n_hap_cg,
                         hap_ambig=n_ambig, oracle_acc=oracle_acc, oracle_bal_acc=oracle_bal)
    print(f"{name:5} {N:3d} {len(an):5d} {len(cg):5d} {n_hap_an:7d} {n_hap_cg:7d} {n_ambig:10d} "
          f"{oracle_acc:11.3f} {oracle_bal:11.3f}")

# ======================================================================
print("\n" + "=" * 100)
print("2. ERRATUM - Phase 2/2b LOO-NN vs its own ceiling, and vs tie policy")
print("=" * 100)
phase2 = json.load(open("phase2_summary.json"))
phase2b = json.load(open("phase2b_kmer_region_summary.json"))
print(f"{'reg':5} {'rep_embed':>10} {'rep_kmer':>9} {'oracle_acc':>11} {'%tied':>7} "
      f"{'embed_tieaware_bal':>19} {'kmer_tieaware_bal':>18}")
for name in REGIONS:
    if name not in results:
        continue
    an = [i for i in range(len(seqs)) if groups[i] == "B_anthracis" and seqs[i][2][name]]
    cg = [i for i in range(len(seqs)) if groups[i] in CEREUS_GROUP and seqs[i][2][name]]
    idx = an + cg
    labels = np.array([0] * len(an) + [1] * len(cg))
    slices = [seqs[i][2][name] for i in idx]

    E = z[name][idx]
    En = E / np.linalg.norm(E, axis=1, keepdims=True)
    D_emb = 1 - (En @ En.T)
    M, _ = kmer_matrix(slices)
    D_km = np.sqrt(((M[:, None, :] - M[None, :, :]) ** 2).sum(-1))

    n = len(idx)
    ntied = 0
    for i in range(n):
        d = D_emb[i].copy(); d[i] = np.inf
        if (d <= d.min() + 1e-9).sum() > 1:
            ntied += 1
    emb_ta = nn_expected_acc(D_emb, labels)
    km_ta = nn_expected_acc(D_km, labels)
    rep_e = phase2[name]["loo_acc"]
    flag = "  <-- reported EXCEEDS oracle" if rep_e > results[name]["oracle_acc"] + 1e-9 else ""
    results[name].update(pct_tied=100 * ntied / n, embed_nn_tieaware_bal=emb_ta,
                         kmer_nn_tieaware_bal=km_ta,
                         reported_embed_loo=rep_e, reported_kmer_loo=phase2b[name]["loo_acc"],
                         reported_exceeds_oracle=bool(rep_e > results[name]["oracle_acc"] + 1e-9))
    print(f"{name:5} {rep_e:10.3f} {phase2b[name]['loo_acc']:9.3f} {results[name]['oracle_acc']:11.3f} "
          f"{100*ntied/n:6.0f}% {emb_ta:19.3f} {km_ta:18.3f}{flag}")

# ======================================================================
print("\n" + "=" * 100)
print("3. REAL B2 (k-mer + random forest) vs RF on DNABERT-2 embeddings")
print("   Leave-one-haplotype-out CV; balanced accuracy; 200-draw label-permutation null")
print("=" * 100)
N_PERM = 200
print(f"{'reg':5} {'hap_an':>7} {'RF_kmer_bal':>12} {'RF_embed_bal':>13} {'perm_mean':>10} "
      f"{'perm_p95':>9} {'p_kmer':>8} {'p_embed':>8}  note")
for name in REGIONS:
    if name not in results:
        continue
    an = [i for i in range(len(seqs)) if groups[i] == "B_anthracis" and seqs[i][2][name]]
    cg = [i for i in range(len(seqs)) if groups[i] in CEREUS_GROUP and seqs[i][2][name]]
    idx = an + cg
    y = np.array([0] * len(an) + [1] * len(cg))
    slices = [seqs[i][2][name] for i in idx]
    hap_codes = np.array([hash(s) for s in slices])

    # not estimable if either class has <2 haplotypes: holding out the only haplotype of a
    # class removes that class from training entirely.
    hap_an = len({slices[i] for i in range(len(an))})
    hap_cg = len({slices[i] for i in range(len(an), len(idx))})
    if hap_an < 2 or hap_cg < 2:
        results[name].update(rf_kmer_bal=None, rf_embed_bal=None,
                             estimable=False, reason=f"only {hap_an} anthracis haplotype(s)")
        print(f"{name:5} {hap_an:7d} {'n/e':>12} {'n/e':>13} {'':>10} {'':>9} {'':>8} {'':>8}  "
              f"NOT ESTIMABLE ({hap_an} anthracis haplotype(s))")
        continue

    Xk, _ = kmer_matrix(slices)
    Xe = z[name][idx]
    bal_k, nk, sk = loho_cv(Xk, y, hap_codes, seed=RNG_SEED)
    bal_e, ne, se = loho_cv(Xe, y, hap_codes, seed=RNG_SEED)

    # permutation null: shuffle labels at the HAPLOTYPE level, preserving haplotype blocks
    uniq_h = np.unique(hap_codes)

    def one_perm(draw, X):
        # Free sequence-level label permutation: H0 = label independent of slice sequence.
        # Labels are NOT collapsed within a haplotype, because doing so would erase the
        # cross-class ambiguity (one string carrying both labels) that dominates this
        # panel - the null must face the same ambiguity the observed data does.
        rng = np.random.default_rng(RNG_SEED + draw)
        yp = rng.permutation(y)
        if len(np.unique(yp)) < 2:
            return np.nan
        b, _, _ = loho_cv(X, yp, hap_codes, seed=RNG_SEED, n_estimators=100)
        return b

    perm_k = np.array([s for s in Parallel(n_jobs=8)(
        delayed(one_perm)(d, Xk) for d in range(N_PERM)) if not np.isnan(s)])
    perm_e = np.array([s for s in Parallel(n_jobs=8)(
        delayed(one_perm)(d, Xe) for d in range(N_PERM)) if not np.isnan(s)])
    perm_scores = perm_k
    # (B+1)/(M+1) permutation p-value - never reports p=0, which 200 draws cannot support
    pm = float(perm_k.mean()) if len(perm_k) else float("nan")
    p95 = float(np.percentile(perm_k, 95)) if len(perm_k) else float("nan")
    p_k = float((np.sum(perm_k >= bal_k) + 1) / (len(perm_k) + 1)) if len(perm_k) else float("nan")
    p_e = float((np.sum(perm_e >= bal_e) + 1) / (len(perm_e) + 1)) if len(perm_e) else float("nan")

    results[name].update(rf_kmer_bal=bal_k, rf_embed_bal=bal_e, estimable=True,
                         perm_mean_kmer=pm, perm_p95_kmer=p95,
                         perm_mean_embed=float(perm_e.mean()) if len(perm_e) else float("nan"),
                         p_kmer=p_k, p_embed=p_e, n_perm=int(len(perm_k)))
    note = "" if (p_k < 0.05 or p_e < 0.05) else "neither beats permutation null"
    print(f"{name:5} {hap_an:7d} {bal_k:12.3f} {bal_e:13.3f} {pm:10.3f} {p95:9.3f} "
          f"{p_k:8.3f} {p_e:8.3f}  {note}")

# ======================================================================
print("\n" + "=" * 100)
print("4. B3 (per-site Shannon entropy) and B4 (length / variable-site density)")
print("=" * 100)
print(f"{'reg':5} {'B3_entropy_all':>15} {'B3_entropy_an':>14} {'B4_mean_len':>12} {'B4_var_site_density':>20}")
for name in REGIONS:
    if name not in results:
        continue
    an = [i for i in range(len(seqs)) if groups[i] == "B_anthracis" and seqs[i][2][name]]
    cg = [i for i in range(len(seqs)) if groups[i] in CEREUS_GROUP and seqs[i][2][name]]
    all_sl = [seqs[i][2][name] for i in an + cg]
    an_sl = [seqs[i][2][name] for i in an]
    e_all = per_site_entropy(all_sl)
    e_an = per_site_entropy(an_sl)
    mlen = float(np.mean([len(s) for s in all_sl]))
    vsd = variable_site_density(all_sl)
    results[name].update(b3_entropy_all=e_all, b3_entropy_anthracis=e_an,
                         b4_mean_len=mlen, b4_var_site_density=vsd)
    print(f"{name:5} {e_all:15.4f} {e_an:14.4f} {mlen:12.1f} {vsd:20.4f}")

# ======================================================================
print("\n" + "=" * 100)
print("5. Benjamini-Hochberg across all (region x arm) cells, per plan §6")
print("=" * 100)
tests = [(f"{n}:{arm}", results[n][f"p_{arm}"])
         for n in REGIONS if results.get(n, {}).get("estimable")
         for arm in ("kmer", "embed")]
tests.sort(key=lambda t: t[1])
m = len(tests)
alpha = 0.05
survivors = set()
kmax = 0
for i, (lab, p) in enumerate(tests, start=1):
    if p <= i / m * alpha:
        kmax = i
for i, (lab, p) in enumerate(tests, start=1):
    if i <= kmax:
        survivors.add(lab)
qvals = {}
prev = 1.0
for i in range(m, 0, -1):
    lab, p = tests[i - 1]
    prev = min(prev, p * m / i)
    qvals[lab] = prev
print(f"{m} tests, alpha={alpha}, permutation p floor = {1/(N_PERM+1):.4f}")
print(f"{'cell':12} {'bal_acc':>8} {'p':>7} {'q(BH)':>7} {'survives':>9} {'vs chance 0.5':>14} {'vs oracle':>10}")
for lab, p in tests:
    reg, arm = lab.split(":")
    bal = results[reg][f"rf_{arm}_bal"]
    orc = results[reg]["oracle_bal_acc"]
    print(f"{lab:12} {bal:8.3f} {p:7.3f} {qvals[lab]:7.3f} {str(lab in survivors):>9} "
          f"{bal-0.5:+14.3f} {bal-orc:+10.3f}")
results["_bh"] = dict(n_tests=m, alpha=alpha, survivors=sorted(survivors),
                      qvalues=qvals, perm_p_floor=1 / (N_PERM + 1))
if survivors:
    print(f"\nBH survivors: {sorted(survivors)}")
else:
    print("\nBH survivors: NONE - no (region x arm) cell shows separability above its "
          "permutation null after correction")

with open("phase2c_b2_rf_summary.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nsaved phase2c_b2_rf_summary.json")
