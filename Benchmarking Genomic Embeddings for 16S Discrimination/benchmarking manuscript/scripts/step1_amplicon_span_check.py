"""Extended Step 1: does the gLM picture hold at AMPLICON length?

gLM evaluated strict hypervariable spans of 31-107 nt. regioneer's unit is the assay: a
primer-defined amplicon of 253-490 nt, primer-trimmed. This re-slices the SAME cached
full-length records at amplicon coordinates and recomputes, per span:

  - distinct observed haplotypes per class (the eligibility unit proposed for section 7.2)
  - target-neighbour collisions (the fraction of target isolates sharing an observed
    sequence with any declared neighbour)
  - the empirical ambiguity ceiling
  - a resolvability flag, from which the WITHIN-PANEL ICC across candidate assays is
    estimated -- the parameter the block-count calculation is most sensitive to
  - k-mer (B4) balanced accuracy under the same probe and isolate-grouped folds, so that
    length transfer is tested rather than inferred from ceilings alone

Amplicon coordinates are E. coli J01695 positions of the INSERT between standard primers,
i.e. after primer removal, matching protocol 1 section 5's "observed sequence".

Writes results/step1_amplicon_spans.json.
"""
import json, os, sys, warnings
from collections import Counter, defaultdict

import numpy as np
from Bio.Align import PairwiseAligner
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "phase2_pilot"))
import run_readouts_full as rr
warnings.filterwarnings("ignore", category=ConvergenceWarning)

_EC = json.load(open(os.path.join(ROOT, "..", "glm_prototype_bundle (1)", "ecoli_16s.txt")))["seq"]
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "step1_amplicon_spans.json"))

# strict hypervariable regions (gLM's unit)
REGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
           "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294),
           "V9": (1435, 1465)}
# primer-trimmed amplicon inserts (regioneer's unit)
AMPLICONS = {"A_V4_515F_806R":   (534, 786),
             "A_V3V4_341F_806R": (358, 786),
             "A_V1V2_27F_338R":  (28, 321),
             "A_V1V3_27F_534R":  (28, 515),
             "A_V6V8_939F_1378R": (958, 1358)}
SPANS = {**{k: v for k, v in REGIONS.items()}, **AMPLICONS}
EXPECTED = {k: e - s + 1 for k, (s, e) in SPANS.items()}
SEED, N_SPLITS = 0, 5


def slice_one(seq):
    al = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                         match_score=2, mismatch_score=-1)
    aln = al.align(_EC, seq)[0]
    r2q = {}
    for t_blk, q_blk in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1, q0 = int(t_blk[0]), int(t_blk[1]), int(q_blk[0])
        for i in range(t0, t1):
            r2q[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in SPANS.items():
        qs = [r2q[p] for p in range(s - 1, e) if p in r2q]
        sl = seq[min(qs):max(qs) + 1] if qs else ""
        out[name] = sl if sl and 0.75 * EXPECTED[name] <= len(sl) <= 1.30 * EXPECTED[name] else None
    return out


def ceiling_bal(slices, y):
    y = np.asarray(y)
    n_t, n_n = int((y == 1).sum()), int((y == 0).sum())
    if not n_t or not n_n: return float("nan")
    c = defaultdict(lambda: [0, 0])
    for s, lab in zip(slices, y): c[s][0 if lab == 1 else 1] += 1
    return sum(max(a / (2 * n_t), b / (2 * n_n)) for a, b in c.values())


def probe():
    return make_pipeline(StandardScaler(), LogisticRegression(
        max_iter=2000, class_weight="balanced", random_state=SEED))


def kmer_acc(slices, y):
    y = np.asarray(y)
    if len(set(y)) < 2 or min(np.bincount(y)) < N_SPLITS: return float("nan")
    out = []
    for tr, te in StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED).split(np.zeros(len(y)), y):
        if len(set(y[tr])) < 2 or len(set(y[te])) < 2: continue
        Ktr, Kte = rr.kmer_fit_transform([slices[i] for i in tr], [slices[i] for i in te])
        c = probe(); c.fit(Ktr, y[tr])
        out.append(rr.balanced_accuracy(y[te], c.predict(Kte)))
    return float(np.mean(out)) if out else float("nan")


def main():
    summ = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_block_summary.json")))
    elig = sorted(b for b, v in summ.items() if v.get("block_eligible"))
    raw = os.path.join(ROOT, "phase4_locked_test", "raw_species_test")

    recs = defaultdict(list)
    for f in os.listdir(raw):
        if not f.endswith(".json"): continue
        for r in json.load(open(os.path.join(raw, f))):
            if r.get("block_id") in elig: recs[r["block_id"]].append(r)
    print(f"blocks {len(recs)}, records {sum(len(v) for v in recs.values())}", flush=True)

    results = {}
    for b in elig:
        rs = recs.get(b, [])
        if not rs: continue
        sl = Parallel(n_jobs=8, batch_size=16)(delayed(slice_one)(r["sequence"]) for r in rs)
        per_span = {}
        for span in SPANS:
            by_iso = defaultdict(list)
            for r, s in zip(rs, sl):
                if s[span]: by_iso[r["isolate_group"]].append((s[span], r["is_target"]))
            if not by_iso: continue
            iso = sorted(by_iso)
            seqs = [Counter(x[0] for x in by_iso[i]).most_common(1)[0][0] for i in iso]
            y = [int(by_iso[i][0][1]) for i in iso]
            if len(set(y)) < 2: continue
            hap = defaultdict(lambda: [0, 0])
            for s_, l in zip(seqs, y): hap[s_][0 if l else 1] += 1
            nbr = {s_ for s_, l in zip(seqs, y) if not l}
            tgt = [s_ for s_, l in zip(seqs, y) if l]
            cb = ceiling_bal(seqs, y)
            per_span[span] = dict(
                expected_len=EXPECTED[span], n_isolates=len(iso),
                median_len=int(np.median([len(s_) for s_ in seqs])),
                n_haplotypes=len(hap),
                hap_target=sum(1 for c in hap.values() if c[0] > 0),
                hap_neighbour=sum(1 for c in hap.values() if c[1] > 0),
                collision_fraction=sum(1 for s_ in tgt if s_ in nbr) / max(len(tgt), 1),
                ceiling_bal=cb, resolvable_099=int(cb >= 0.99),
                kmer_acc=kmer_acc(seqs, y))
        results[b] = per_span
        a = per_span.get("A_V3V4_341F_806R", {}); v4 = per_span.get("V4", {})
        print(f"  {b:26} V4 strict: hap={v4.get('n_haplotypes','-'):>4} "
              f"coll={v4.get('collision_fraction',float('nan')):.2f} "
              f"ceil={v4.get('ceiling_bal',float('nan')):.3f}   |  V3V4 amplicon: "
              f"hap={a.get('n_haplotypes','-'):>4} coll={a.get('collision_fraction',float('nan')):.2f} "
              f"ceil={a.get('ceiling_bal',float('nan')):.3f}", flush=True)

    json.dump({"spans": {k: {"coords": v, "expected_len": EXPECTED[k]} for k, v in SPANS.items()},
               "blocks": results}, open(OUT, "w"), indent=2, default=float)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
