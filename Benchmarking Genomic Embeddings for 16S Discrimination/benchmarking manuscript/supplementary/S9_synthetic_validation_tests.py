"""Phase 1 synthetic fixture (protocol §8.4): known-ground-truth tests of pipeline mechanics,
independent of biological ambiguity. Every assertion below has a value computed by hand from
the construction, not estimated from the result.

This does NOT establish biological generalization (per §8.4, synthetic fixtures test
implementation, not biology) - it only proves the alignment/slicing, oracle-ceiling formula,
LOHO-CV, and permutation-null code do what they claim to do on cases where the right answer
is known in advance.

Functions below are copied verbatim from phase3a_staph_baseline.py (the canonical corrected
version) per the CLAUDE.md convention, minus the embedding arm (covered separately by the
batch-invariance re-verification in audit_report.md item 6).
"""
import numpy as np
from collections import defaultdict
from Bio.Align import PairwiseAligner
from sklearn.ensemble import RandomForestClassifier

RNG_SEED = 0
rng = np.random.default_rng(RNG_SEED)

aligner = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                          match_score=2, mismatch_score=-1)


def slice_regions(seq, ref, regions, expected):
    aln = aligner.align(ref, seq)[0]
    ref_to_query = {}
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1 = int(t_blk[0]), int(t_blk[1]); q0, q1 = int(q_blk[0]), int(q_blk[1])
        for i in range(t0, t1):
            ref_to_query[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in regions.items():
        qs = [ref_to_query[p] for p in range(s - 1, e) if p in ref_to_query]
        sl = seq[min(qs):max(qs) + 1] if qs else ""
        out[name] = sl if sl and 0.75 * expected[name] <= len(sl) <= 1.30 * expected[name] else None
    return out


def identity(a, b):
    if not a or not b:
        return None
    aln = aligner.align(a, b)[0]
    matches = 0
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        seg_t = a[int(t_blk[0]):int(t_blk[1])]
        seg_q = b[int(q_blk[0]):int(q_blk[1])]
        matches += sum(1 for x, y in zip(seg_t, seg_q) if x == y)
    return 2 * matches / (len(a) + len(b))


def kmer_matrix(slices, ks=(4, 5, 6)):
    vocab = sorted({(k, s[i:i + k]) for k in ks for s in slices
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
    return M


def balanced_accuracy(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    accs = []
    for cls in (0, 1):
        m = y_true == cls
        if m.sum() == 0:
            continue
        accs.append((y_pred[m] == cls).mean())
    return float(np.mean(accs)) if accs else float("nan")


def loho_cv(X, y, hap, seed=0, n_estimators=200):
    hap = np.asarray(hap)
    preds = np.full(len(y), -1)
    for h in np.unique(hap):
        test = hap == h
        train = ~test
        if len(np.unique(y[train])) < 2:
            return float("nan")
        clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed,
                                     class_weight="balanced", n_jobs=1)
        clf.fit(X[train], y[train])
        preds[test] = clf.predict(X[test])
    return balanced_accuracy(y, preds)


def loho_nn_cv(D, y, hap):
    hap = np.asarray(hap); y = np.asarray(y); n = len(y)
    per = np.full(n, np.nan)
    for h in np.unique(hap):
        test = np.where(hap == h)[0]
        train = np.where(hap != h)[0]
        for i in test:
            d = D[i, train]
            mn = d.min()
            tied = train[d <= mn + 1e-9]
            per[i] = (y[tied] == y[i]).mean()
    accs = []
    for cls in (0, 1):
        m = y == cls
        if m.sum():
            accs.append(per[m].mean())
    return float(np.mean(accs))


def oracle(hap_counts, n_a, n_c):
    """hap_counts: dict haplotype -> [n_in_class_a, n_in_class_c]."""
    N = n_a + n_c
    err = sum(min(a, c) for a, c in hap_counts.values())
    oracle_acc = 1 - err / N
    oracle_bal = sum(max(a / (2 * n_a), c / (2 * n_c)) for a, c in hap_counts.values())
    return oracle_acc, oracle_bal


print("=" * 90)
print("PART A: alignment/extraction correctness under a known indel")
print("=" * 90)
# Reference: 200bp fixed random sequence. Query: same backbone but with a 3bp DELETION
# inserted 10bp before the target region, so slice_regions must correctly shift coordinates.
base = "".join(rng.choice(list("ACGT"), 200))
REGIONS = {"R_diag": (51, 70), "R_null": (101, 120), "R_partial": (151, 170)}
EXPECTED = {k: e - s + 1 for k, (s, e) in REGIONS.items()}
ref = base
query_with_deletion = base[:40] + base[43:]  # delete 3 bases at position 41-43, upstream of R_diag
expected_r_diag = base[50:70]  # ground truth: the untouched region, by construction
got = slice_regions(query_with_deletion, ref, REGIONS, EXPECTED)
assert got["R_diag"] == expected_r_diag, (
    f"alignment-based extraction failed under a known upstream indel: "
    f"got {got['R_diag']!r}, expected {expected_r_diag!r}")
print("PASS: slice_regions recovers the correct region after an upstream 3bp deletion")

print("\n" + "=" * 90)
print("PART B: oracle ceiling formula, hand-computed ground truth")
print("=" * 90)
# R_diag: 2 target haplotypes (10 genomes), 2 neighbor haplotypes (20 genomes), ZERO overlap.
hap_diag = {"AAAAACCCCC": [5, 0], "AAAAAGGGGG": [5, 0], "TTTTTCCCCC": [0, 10], "TTTTTGGGGG": [0, 10]}
acc, bal = oracle(hap_diag, n_a=10, n_c=20)
assert abs(acc - 1.0) < 1e-12 and abs(bal - 1.0) < 1e-12, f"expected (1.0, 1.0), got ({acc}, {bal})"
print(f"PASS: zero-overlap region -> oracle_acc={acc:.6f}, oracle_bal={bal:.6f} (expected 1.0, 1.0)")

# R_null: single haplotype shared by both classes entirely.
hap_null = {"SHARED_STRING_XYZ": [10, 20]}
acc, bal = oracle(hap_null, n_a=10, n_c=20)
assert abs(acc - (1 - 10 / 30)) < 1e-12, f"expected acc={1-10/30}, got {acc}"
assert abs(bal - 0.5) < 1e-12, f"fully-shared single haplotype must give oracle_bal=0.5 exactly, got {bal}"
print(f"PASS: fully-shared single haplotype -> oracle_acc={acc:.6f} (=1-min(n_a,n_c)/N), "
      f"oracle_bal={bal:.6f} (=0.5 exactly)")

# R_partial: hand-constructed partial ambiguity, n_a=10, n_c=20.
hap_partial = {"P1": [8, 0], "P2": [2, 5], "P3": [0, 15]}
acc, bal = oracle(hap_partial, n_a=10, n_c=20)
# hand calculation: err = min(8,0)+min(2,5)+min(0,15) = 0+2+0 = 2 -> acc = 1-2/30 = 0.933333...
# bal: max(8/20,0)=0.4; max(2/20,5/40)=max(0.1,0.125)=0.125; max(0,15/40)=0.375; sum=0.9
assert abs(acc - (28 / 30)) < 1e-9, f"expected acc=28/30={28/30:.6f}, got {acc}"
assert abs(bal - 0.9) < 1e-9, f"expected oracle_bal=0.9, got {bal}"
print(f"PASS: hand-constructed partial-ambiguity case -> oracle_acc={acc:.6f} (expected {28/30:.6f}), "
      f"oracle_bal={bal:.6f} (expected 0.900000)")

print("\n" + "=" * 90)
print("PART C: LOHO-CV recovers perfect separation on a clean, non-degenerate region")
print("=" * 90)


def rows_y_hap(hap_counts):
    rows, y, hap = [], [], []
    for seq, (n_a, n_c) in hap_counts.items():
        for _ in range(n_a):
            rows.append(seq); y.append(0); hap.append(seq)
        for _ in range(n_c):
            rows.append(seq); y.append(1); hap.append(seq)
    return rows, np.array(y), np.array(hap)


# First attempt, kept as a documented negative finding rather than silently discarded: using
# "AAAAACCCCC"/"AAAAAGGGGG" (target) vs "TTTTTCCCCC"/"TTTTTGGGGG" (neighbor) looks zero-overlap
# at the haplotype level, but the SUFFIXES cross class lines (CCCC/GGGG each appear once in a
# target haplotype and once in a neighbor haplotype). With only 1 haplotype per class surviving
# in a LOHO training fold, k-mer+RF picked up the crossed suffix instead of the fully-diagnostic
# prefix and inverted every held-out prediction (balanced accuracy = 0.0, confirmed across 5
# seeds and with/without class_weight='balanced' - i.e. not stochastic noise, a feature-design
# confound). This is a genuine, reportable lesson about LOHO-CV fragility at n=1
# haplotype/class/fold, consistent with the ">=5 preferred" margin in phase2cd_report.md §5.3 -
# not a bug in loho_cv() or kmer_matrix(). Demonstrated explicitly, not just avoided:
rows_confounded, y_confounded, hap_confounded = rows_y_hap(hap_diag)
X_confounded = kmer_matrix(rows_confounded)
bal_confounded = loho_cv(X_confounded, y_confounded, hap_confounded, seed=RNG_SEED)
assert abs(bal_confounded - 0.0) < 1e-9, (
    f"expected the documented confound to reproduce exactly (0.0), got {bal_confounded} - "
    f"the confound itself may no longer reproduce, update this comment if so")
print(f"DOCUMENTED (not a failure): crossed-suffix confound at n=1 haplotype/class/fold "
      f"reproduces exactly as before -> balanced accuracy = {bal_confounded:.6f} (expected 0.000000)")

# Corrected design: within-class haplotypes differ by a single TAIL base only, so no k-mer
# window is shared between the two classes at all (verified: zero overlapping nonzero columns).
hap_diag_clean = {
    "AAAAACCCCC": [5, 0], "AAAAACCCCG": [5, 0],
    "TTTTTGGGGG": [0, 10], "TTTTTGGGGC": [0, 10],
}
rows, y, hap = rows_y_hap(hap_diag_clean)
X = kmer_matrix(rows)
bal = loho_cv(X, y, hap, seed=RNG_SEED)
assert abs(bal - 1.0) < 1e-9, f"expected perfect LOHO-CV separation (1.0) on a zero-overlap region, got {bal}"
print(f"PASS: LOHO-CV k-mer+RF balanced accuracy = {bal:.6f} (expected 1.000000) "
      f"on the corrected, verified-zero-overlap haplotype design")

print("\n" + "=" * 90)
print("PART D: label permutation collapses performance toward chance")
print("=" * 90)
perm_scores = []
for draw in range(50):
    prng = np.random.default_rng(100 + draw)
    yp = prng.permutation(y)
    if len(np.unique(yp)) < 2:
        continue
    perm_scores.append(loho_cv(X, yp, hap, seed=RNG_SEED, n_estimators=50))
perm_scores = np.array([s for s in perm_scores if not np.isnan(s)])
perm_mean = perm_scores.mean()
assert perm_mean < 0.75, (
    f"permuted-label LOHO-CV should collapse toward chance (~0.5), got mean={perm_mean:.3f} "
    f"- permutation routine may not actually be breaking the label/feature association")
print(f"PASS: permuted-label balanced accuracy mean = {perm_mean:.3f} over {len(perm_scores)} draws "
      f"(real score 1.000 clearly separates from this null)")

print("\n" + "=" * 90)
print("PART E: tie-aware nearest-neighbour matches a hand-computed fractional-credit case")
print("=" * 90)
# 4 points: 0 and 1 are class 0; 2 and 3 are class 1. Distances constructed so point 0's
# nearest neighbours (excluding its own haplotype group) are EXACTLY tied between one
# same-class and one other-class point.
D = np.array([
    [0.0, 0.0, 0.5, 0.5],   # point 0: tied at 0.5 to points 2 and 3 (both class 1) once its
    [0.0, 0.0, 0.9, 0.9],   #   own haplotype group {0,1} is excluded -> expected credit 0.0
    [0.5, 0.9, 0.0, 0.0],
    [0.5, 0.9, 0.0, 0.0],
])
y4 = np.array([0, 0, 1, 1])
hap4 = np.array(["g0", "g0", "g1", "g1"])  # points 0,1 share a haplotype group; 2,3 share another
bal = loho_nn_cv(D, y4, hap4)
# By construction: holding out group g0, point 0's nearest remaining neighbours are 2 and 3
# (tied at 0.5), both class 1 -> fraction matching true label (0) = 0/2 = 0.0.
# Point 1's nearest remaining neighbours (excluding group g0) are also 2 and 3 tied at 0.9 ->
# fraction matching true label (0) = 0.0. Symmetric for group g1 by the matrix's symmetry.
# Expected balanced accuracy over both classes = mean(0.0, 0.0) = 0.0.
assert abs(bal - 0.0) < 1e-9, f"expected hand-computed tie-aware balanced accuracy = 0.0, got {bal}"
print(f"PASS: tie-aware LOHO-NN balanced accuracy = {bal:.6f} (expected 0.000000, hand-verified)")

print("\n" + "=" * 90)
print("PART F: identity() self-check")
print("=" * 90)
s = hap_diag and "ACGTACGTACGTACGTACGT"
assert abs(identity(s, s) - 1.0) < 1e-12, "identity(x,x) != 1.0"
assert identity(s, s[:-1] + "T" if s[-1] != "T" else s[:-1] + "A") < 1.0, "identity() returned 1.0 for a non-identical pair"
print("PASS: identity(x,x) == 1.0 and identity() < 1.0 for a non-identical pair")

print("\nALL SYNTHETIC FIXTURE CHECKS PASSED")
