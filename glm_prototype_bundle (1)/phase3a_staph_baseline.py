"""Phase 3a: full baseline (B1-B4 + gLM) on the S. aureus vs CoNS panel.

MASTER_PLAN.md §4 item 2's recommended next step: move off the B. anthracis/cereus
expected-null target to one where real partial separation is published (Chakravorty et al.
2007: V1 best for S. aureus vs coagulase-negative staphylococci). Run it with the
methodology Phase 2c/2d had to retrofit onto the anthracis panel, built in from the start:
leave-one-haplotype-out CV, a sequence-level label-permutation null per arm, BH correction,
oracle-ceiling reporting, and a DNABERT-2 embedding that is a pure function of the sequence
by construction (dedupe-then-broadcast, per phase2d_reembed_fixed.py).

Panel: target = Staphylococcus aureus; near-neighbor class = pooled coagulase-negative
staphylococci (S. epidermidis, S. haemolyticus, S. saprophyticus, S. hominis,
S. lugdunensis, S. capitis), fetched in fetch_staph.py -> staph_group_16s.json (787 raw
records -> 195 unique near-full-length sequences after QA). The haplotype pre-check
(haplotype_gate_staph.py) found V1 has ZERO cross-class haplotype ambiguity and oracle
accuracy = 1.000 - unlike the anthracis panel, this target can show a real winner.

Four arms, all scored under IDENTICAL LOHO-CV so that string-identical slices never
straddle a train/test split:
  B1  - identity: tie-aware 1-NN on pairwise sequence identity (the plan's simplest baseline)
  B2  - k-mer (4/5/6-mer) + random forest
  gLM - DNABERT-2 embedding (pure function of sequence) + random forest
B3 (per-site Shannon entropy) and B4 (length / variable-site density) are diagnostic
rankers, not scored classifiers, matching their role in Phase 2c.

Slicing/identity/k-mer/entropy/RF-CV functions are copied verbatim from
phase2c_b2_rf_baseline.py (the canonical corrected version) per the CLAUDE.md convention.
"""
import json
import numpy as np
import torch
from collections import defaultdict
from Bio.Align import PairwiseAligner
from sklearn.ensemble import RandomForestClassifier
from joblib import Parallel, delayed
from transformers import AutoModel, AutoTokenizer

RNG_SEED = 0
MODEL_DIR = "dnabert2_local"
N_PERM = 200

staph_records = json.load(open("staph_group_16s.json"))
ec_ref = json.load(open("ecoli_16s.txt"))["seq"]
aligner = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                          match_score=2, mismatch_score=-1)
REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173),
           "V8": (1243, 1294), "V9": (1435, 1465)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}
TARGET = "S_aureus"


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


# ---- QA (same filter as the cereus-group panel: dedup, 1450-1600nt full-length) ----
seqs, seen = [], set()
for grp, hdr, seq in staph_records:
    if seq in seen or not (1450 <= len(seq) <= 1600):
        continue
    seen.add(seq)
    seqs.append((grp, hdr.split()[0], slice_regions(seq)))
print(f"{len(seqs)} unique near-full-length sequences after QA (of {len(staph_records)} raw)", flush=True)
groups = np.array([s[0] for s in seqs])
ids = np.array([s[1] for s in seqs])

# ---- MANDATORY identical-sequence QA check (CLAUDE.md convention) ----
_a = seqs[0][2]["V4"] or seqs[0][2]["V2"]
assert abs(identity(_a, _a) - 1.0) < 1e-12, "identity(x,x) != 1.0 - scoring function is broken"
print("identity() identical-sequence QA check: PASS", flush=True)


def kmer_matrix(slices, ks=(4, 5, 6)):
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
    hap = np.asarray(hap)
    uniq = np.unique(hap)
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


def loho_nn_cv(D, y, hap):
    """Leave-one-haplotype-out tie-aware 1-NN on a precomputed distance matrix D. For each
    held-out haplotype group, the nearest neighbour is found only among the REMAINING pool
    (so a haplotype's own members never serve as their own match) - the identity-baseline
    analogue of loho_cv for a trained classifier. Ties get fractional credit (expected
    accuracy under uniform random tie-break), matching nn_expected_acc's convention."""
    hap = np.asarray(hap)
    y = np.asarray(y)
    n = len(y)
    per = np.full(n, np.nan)
    for h in np.unique(hap):
        test = np.where(hap == h)[0]
        train = np.where(hap != h)[0]
        if len(train) == 0:
            continue
        for i in test:
            d = D[i, train]
            mn = d.min()
            tied = train[d <= mn + 1e-9]
            per[i] = (y[tied] == y[i]).mean()
    mask = ~np.isnan(per)
    accs = []
    for cls in (0, 1):
        m = mask & (y == cls)
        if m.sum():
            accs.append(per[m].mean())
    return float(np.mean(accs)) if accs else float("nan")


# ======================================================================
# Region-indexed index sets, reused across every section below
# ======================================================================
region_idx = {}
for name in REGIONS:
    tgt = [i for i in range(len(seqs)) if groups[i] == TARGET and seqs[i][2][name]]
    cns = [i for i in range(len(seqs)) if groups[i] != TARGET and seqs[i][2][name]]
    region_idx[name] = (tgt, cns)

# ======================================================================
print("\n" + "=" * 100)
print("1. ORACLE CEILING - best achievable given cross-class-identical slices")
print("=" * 100)
results = {}
print(f"{'reg':5} {'N':>4} {'n_tgt':>5} {'n_cns':>5} {'hap_tgt':>7} {'hap_cns':>7} {'hap_ambig':>10} "
      f"{'oracle_acc':>11} {'oracle_bal':>11}")
for name in REGIONS:
    tgt, cns = region_idx[name]
    if not tgt or not cns:
        continue
    sl = {i: seqs[i][2][name] for i in tgt + cns}
    hap = defaultdict(lambda: [0, 0])
    for i in tgt:
        hap[sl[i]][0] += 1
    for i in cns:
        hap[sl[i]][1] += 1
    N = len(tgt) + len(cns)
    err = sum(min(a, c) for a, c in hap.values())
    oracle_acc = 1 - err / N
    oracle_bal = sum(max(a / (2 * len(tgt)), c / (2 * len(cns))) for a, c in hap.values())
    n_hap_tgt = len({sl[i] for i in tgt})
    n_hap_cns = len({sl[i] for i in cns})
    n_ambig = sum(1 for a, c in hap.values() if a > 0 and c > 0)
    results[name] = dict(n=N, n_tgt=len(tgt), n_cns=len(cns), hap_tgt=n_hap_tgt, hap_cns=n_hap_cns,
                         hap_ambig=n_ambig, oracle_acc=oracle_acc, oracle_bal_acc=oracle_bal)
    print(f"{name:5} {N:4d} {len(tgt):5d} {len(cns):5d} {n_hap_tgt:7d} {n_hap_cns:7d} {n_ambig:10d} "
          f"{oracle_acc:11.3f} {oracle_bal:11.3f}")

# ======================================================================
print("\n" + "=" * 100)
print("2. B1 - identity baseline: tie-aware 1-NN on pairwise sequence identity, LOHO-CV")
print("=" * 100)
print(f"{'reg':5} {'B1_bal':>8} {'perm_mean':>10} {'p_B1':>8}  note")
for name in REGIONS:
    if name not in results:
        continue
    tgt, cns = region_idx[name]
    idx = tgt + cns
    y = np.array([0] * len(tgt) + [1] * len(cns))
    slices = [seqs[i][2][name] for i in idx]
    hap_codes = np.array([hash(s) for s in slices])

    n = len(idx)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = 1 - identity(slices[i], slices[j])
            D[i, j] = D[j, i] = d

    hap_tgt = len({slices[i] for i in range(len(tgt))})
    hap_cns = len({slices[i] for i in range(len(tgt), n)})
    if hap_tgt < 2 or hap_cns < 2:
        results[name].update(b1_bal=None, estimable_b1=False)
        print(f"{name:5} {'n/e':>8} {'':>10} {'':>8}  NOT ESTIMABLE ({hap_tgt} target haplotype(s))")
        continue

    bal_b1 = loho_nn_cv(D, y, hap_codes)

    def one_perm_b1(draw):
        rng = np.random.default_rng(RNG_SEED + draw)
        yp = rng.permutation(y)
        if len(np.unique(yp)) < 2:
            return np.nan
        return loho_nn_cv(D, yp, hap_codes)

    perm_b1 = np.array([s for s in Parallel(n_jobs=8)(
        delayed(one_perm_b1)(d) for d in range(N_PERM)) if not np.isnan(s)])
    pm = float(perm_b1.mean()) if len(perm_b1) else float("nan")
    p_b1 = float((np.sum(perm_b1 >= bal_b1) + 1) / (len(perm_b1) + 1)) if len(perm_b1) else float("nan")
    results[name].update(b1_bal=bal_b1, estimable_b1=True, perm_mean_b1=pm, p_b1=p_b1, n_perm_b1=int(len(perm_b1)))
    note = "" if p_b1 < 0.05 else "does not beat permutation null"
    print(f"{name:5} {bal_b1:8.3f} {pm:10.3f} {p_b1:8.3f}  {note}")

# ======================================================================
print("\n" + "=" * 100)
print("3. B2 (k-mer + random forest), LOHO-CV, 200-draw label-permutation null")
print("=" * 100)
print(f"{'reg':5} {'hap_tgt':>7} {'RF_kmer_bal':>12} {'perm_mean':>10} {'p_kmer':>8}  note")
for name in REGIONS:
    if name not in results:
        continue
    tgt, cns = region_idx[name]
    idx = tgt + cns
    y = np.array([0] * len(tgt) + [1] * len(cns))
    slices = [seqs[i][2][name] for i in idx]
    hap_codes = np.array([hash(s) for s in slices])
    hap_tgt = len({slices[i] for i in range(len(tgt))})
    hap_cns = len({slices[i] for i in range(len(tgt), len(idx))})

    if hap_tgt < 2 or hap_cns < 2:
        results[name].update(rf_kmer_bal=None, estimable=False, reason=f"only {hap_tgt} target haplotype(s)")
        print(f"{name:5} {hap_tgt:7d} {'n/e':>12} {'':>10} {'':>8}  NOT ESTIMABLE ({hap_tgt} target haplotype(s))")
        continue

    Xk, _ = kmer_matrix(slices)
    bal_k, nk, sk = loho_cv(Xk, y, hap_codes, seed=RNG_SEED)

    def one_perm(draw, X):
        rng = np.random.default_rng(RNG_SEED + draw)
        yp = rng.permutation(y)
        if len(np.unique(yp)) < 2:
            return np.nan
        b, _, _ = loho_cv(X, yp, hap_codes, seed=RNG_SEED, n_estimators=100)
        return b

    perm_k = np.array([s for s in Parallel(n_jobs=8)(
        delayed(one_perm)(d, Xk) for d in range(N_PERM)) if not np.isnan(s)])
    pm = float(perm_k.mean()) if len(perm_k) else float("nan")
    p_k = float((np.sum(perm_k >= bal_k) + 1) / (len(perm_k) + 1)) if len(perm_k) else float("nan")
    results[name].update(rf_kmer_bal=bal_k, estimable=True, perm_mean_kmer=pm, p_kmer=p_k, n_perm=int(len(perm_k)))
    note = "" if p_k < 0.05 else "does not beat permutation null"
    print(f"{name:5} {hap_tgt:7d} {bal_k:12.3f} {pm:10.3f} {p_k:8.3f}  {note}")

# ======================================================================
print("\n" + "=" * 100)
print("4. B3 (per-site Shannon entropy) and B4 (length / variable-site density)")
print("=" * 100)
print(f"{'reg':5} {'B3_entropy_all':>15} {'B3_entropy_tgt':>15} {'B4_mean_len':>12} {'B4_var_site_density':>20}")
for name in REGIONS:
    if name not in results:
        continue
    tgt, cns = region_idx[name]
    all_sl = [seqs[i][2][name] for i in tgt + cns]
    tgt_sl = [seqs[i][2][name] for i in tgt]
    e_all = per_site_entropy(all_sl)
    e_tgt = per_site_entropy(tgt_sl)
    mlen = float(np.mean([len(s) for s in all_sl]))
    vsd = variable_site_density(all_sl)
    results[name].update(b3_entropy_all=e_all, b3_entropy_target=e_tgt, b4_mean_len=mlen, b4_var_site_density=vsd)
    print(f"{name:5} {e_all:15.4f} {e_tgt:15.4f} {mlen:12.1f} {vsd:20.4f}")

# ======================================================================
print("\n" + "=" * 100)
print("5. gLM arm - DNABERT-2 embedding (pure function of sequence) + random forest")
print("=" * 100)
tok = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModel.from_pretrained(MODEL_DIR, trust_remote_code=True)
model.eval()
for _n, module in model.named_modules():
    if hasattr(module, "rebuild_alibi_tensor") and getattr(module, "alibi", "x") is None:
        module.rebuild_alibi_tensor(size=module._current_alibi_size)


def embed_batch(seqlist, batch_size=16):
    outs = []
    for i in range(0, len(seqlist), batch_size):
        batch = seqlist[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=256)
        with torch.no_grad():
            out = model(**enc)[0]
        am = enc["attention_mask"].bool()
        mask = am.clone()
        lengths = am.sum(1)
        mask[:, 0] = False
        rows = torch.arange(mask.size(0))
        mask[rows, lengths - 1] = False
        m = mask.unsqueeze(-1).float()
        outs.append(((out * m).sum(1) / m.sum(1).clamp(min=1)).numpy())
    return np.vstack(outs)


embeddings = {}
for name in REGIONS:
    present = [seqs[i][2][name] if seqs[i][2][name] else "" for i in range(len(seqs))]
    ok = [i for i, s in enumerate(present) if s]
    E = np.full((len(seqs), 768), np.nan)
    if ok:
        uniq = sorted({present[i] for i in ok})
        U = embed_batch(uniq)
        lookup = {s: U[j] for j, s in enumerate(uniq)}
        for i in ok:
            E[i] = lookup[present[i]]
    embeddings[name] = E
    print(f"embedded {name}: {len(ok)}/{len(seqs)} slices ({len({present[i] for i in ok})} unique)", flush=True)

# embedding purity check (CLAUDE.md convention): identical slices -> identical embeddings
worst = 0.0
for name in REGIONS:
    E = embeddings[name]
    byslice = defaultdict(list)
    for i in range(len(seqs)):
        if seqs[i][2][name]:
            byslice[seqs[i][2][name]].append(i)
    for s, idxs in byslice.items():
        if len(idxs) > 1:
            sub = E[idxs]
            worst = max(worst, float(np.abs(sub - sub[0]).max()))
assert worst == 0.0, f"embedding is batch-dependent: max within-haplotype deviation {worst:.3e}"
print(f"embedding purity check (identical slices -> identical embeddings): PASS (deviation = {worst:.3e})")

np.savez_compressed("staph_dnabert2_embeddings.npz", **embeddings,
                    ids=ids, groups=groups)
print("saved staph_dnabert2_embeddings.npz")

print(f"\n{'reg':5} {'RF_embed_bal':>13} {'perm_mean':>10} {'p_embed':>8}  note")
for name in REGIONS:
    if name not in results or not results[name].get("estimable"):
        print(f"{name:5} {'n/e':>13} {'':>10} {'':>8}  NOT ESTIMABLE (see B2 section above)")
        continue
    tgt, cns = region_idx[name]
    idx = tgt + cns
    y = np.array([0] * len(tgt) + [1] * len(cns))
    slices = [seqs[i][2][name] for i in idx]
    hap_codes = np.array([hash(s) for s in slices])
    Xe = embeddings[name][idx]

    bal_e, ne, se = loho_cv(Xe, y, hap_codes, seed=RNG_SEED)

    def one_perm_e(draw):
        rng = np.random.default_rng(RNG_SEED + draw)
        yp = rng.permutation(y)
        if len(np.unique(yp)) < 2:
            return np.nan
        b, _, _ = loho_cv(Xe, yp, hap_codes, seed=RNG_SEED, n_estimators=100)
        return b

    perm_e = np.array([s for s in Parallel(n_jobs=8)(
        delayed(one_perm_e)(d) for d in range(N_PERM)) if not np.isnan(s)])
    pm = float(perm_e.mean()) if len(perm_e) else float("nan")
    p_e = float((np.sum(perm_e >= bal_e) + 1) / (len(perm_e) + 1)) if len(perm_e) else float("nan")
    results[name].update(rf_embed_bal=bal_e, perm_mean_embed=pm, p_embed=p_e, n_perm_embed=int(len(perm_e)))
    note = "" if p_e < 0.05 else "does not beat permutation null"
    print(f"{name:5} {bal_e:13.3f} {pm:10.3f} {p_e:8.3f}  {note}")

# ======================================================================
print("\n" + "=" * 100)
print("6. Benjamini-Hochberg across all (region x arm) cells, per plan §6")
print("=" * 100)
tests = []
for n in REGIONS:
    r = results.get(n, {})
    if r.get("estimable_b1"):
        tests.append((f"{n}:B1_identity", r["p_b1"]))
    if r.get("estimable"):
        tests.append((f"{n}:B2_kmer", r["p_kmer"]))
        tests.append((f"{n}:gLM_embed", r["p_embed"]))
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

bal_lookup = {}
for n in REGIONS:
    r = results.get(n, {})
    if "b1_bal" in r and r["b1_bal"] is not None:
        bal_lookup[f"{n}:B1_identity"] = r["b1_bal"]
    if r.get("rf_kmer_bal") is not None:
        bal_lookup[f"{n}:B2_kmer"] = r["rf_kmer_bal"]
    if r.get("rf_embed_bal") is not None:
        bal_lookup[f"{n}:gLM_embed"] = r["rf_embed_bal"]

print(f"{m} tests, alpha={alpha}, permutation p floor = {1/(N_PERM+1):.4f}")
print(f"{'cell':15} {'bal_acc':>8} {'p':>7} {'q(BH)':>7} {'survives':>9} {'vs chance 0.5':>14} {'vs oracle':>10}")
for lab, p in tests:
    reg = lab.split(":")[0]
    bal = bal_lookup[lab]
    orc = results[reg]["oracle_bal_acc"]
    print(f"{lab:15} {bal:8.3f} {p:7.3f} {qvals[lab]:7.3f} {str(lab in survivors):>9} "
          f"{bal-0.5:+14.3f} {bal-orc:+10.3f}")
results["_bh"] = dict(n_tests=m, alpha=alpha, survivors=sorted(survivors), qvalues=qvals,
                      perm_p_floor=1 / (N_PERM + 1))
if survivors:
    print(f"\nBH survivors: {sorted(survivors)}")
else:
    print("\nBH survivors: NONE")

with open("phase3a_staph_summary.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nsaved phase3a_staph_summary.json")
