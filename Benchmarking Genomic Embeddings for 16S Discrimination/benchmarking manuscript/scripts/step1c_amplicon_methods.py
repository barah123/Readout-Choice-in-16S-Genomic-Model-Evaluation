"""Step 1c: does the gLM k-mer/encoder equivalence TRANSFER to amplicon length?

Step 1 (section 2B of the amendment) showed headroom grows with span length -- 1.58 pp across
strict regions, 2.19 pp across amplicons. Ceilings bound what ANY method can do; they do not
say whether methods differ. This runs the comparison itself at amplicon length, like for like
with gLM: same probe, same isolate-grouped folds, same pooling, conditional endpoint (no PCR).

Caches the re-sliced sequences so the alignment pass is paid once.
Writes results/step1c_amplicon_methods.json.
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
import embed_all_models as eam
warnings.filterwarnings("ignore", category=ConvergenceWarning)

_EC = json.load(open(os.path.join(ROOT, "..", "glm_prototype_bundle (1)", "ecoli_16s.txt")))["seq"]
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
CACHE = os.path.join(RES, "step1c_slices.json")
OUT = os.path.join(RES, "step1c_amplicon_methods.json")

SPANS = {"V2": (137, 242), "V4": (576, 682),
         "A_V4_515F_806R": (534, 786), "A_V1V2_27F_338R": (28, 321),
         "A_V6V8_939F_1378R": (958, 1358), "A_V3V4_341F_806R": (358, 786),
         "A_V1V3_27F_534R": (28, 515)}
EXPECTED = {k: e - s + 1 for k, (s, e) in SPANS.items()}
SEED, N_SPLITS = 0, 5


def slice_one(seq):
    al = PairwiseAligner(mode="global", open_gap_score=-5, extend_gap_score=-2,
                         match_score=2, mismatch_score=-1)
    aln = al.align(_EC, seq)[0]
    r2q = {}
    for t, q in zip(aln.aligned[0], aln.aligned[1]):
        t0, t1, q0 = int(t[0]), int(t[1]), int(q[0])
        for i in range(t0, t1):
            r2q[i] = q0 + (i - t0)
    out = {}
    for name, (s, e) in SPANS.items():
        qs = [r2q[p] for p in range(s - 1, e) if p in r2q]
        sl = seq[min(qs):max(qs) + 1] if qs else ""
        out[name] = sl if sl and 0.75 * EXPECTED[name] <= len(sl) <= 1.30 * EXPECTED[name] else None
    return out


def probe():
    return make_pipeline(StandardScaler(), LogisticRegression(
        max_iter=2000, class_weight="balanced", random_state=SEED))


def build_cache():
    summ = json.load(open(os.path.join(ROOT, "phase4_locked_test", "test_block_summary.json")))
    elig = sorted(b for b, v in summ.items() if v.get("block_eligible"))
    raw = os.path.join(ROOT, "phase4_locked_test", "raw_species_test")
    recs = defaultdict(list)
    for f in os.listdir(raw):
        if f.endswith(".json"):
            for r in json.load(open(os.path.join(raw, f))):
                if r.get("block_id") in elig:
                    recs[r["block_id"]].append(r)
    data = {}
    for b in elig:
        rs = recs.get(b, [])
        if not rs: continue
        sl = Parallel(n_jobs=8, batch_size=16)(delayed(slice_one)(r["sequence"]) for r in rs)
        per = {}
        for span in SPANS:
            by = defaultdict(list)
            for r, s in zip(rs, sl):
                if s[span]: by[r["isolate_group"]].append((s[span], int(r["is_target"])))
            if not by: continue
            iso = sorted(by)
            per[span] = {"iso": iso,
                         "seq": [Counter(x[0] for x in by[i]).most_common(1)[0][0] for i in iso],
                         "y": [by[i][0][1] for i in iso]}
        data[b] = per
        print(f"  sliced {b}", flush=True)
    json.dump(data, open(CACHE, "w"))
    return data


def main():
    data = json.load(open(CACHE)) if os.path.exists(CACHE) else build_cache()
    uniq = sorted({s for pb in data.values() for sp in pb.values() for s in sp["seq"]})
    print(f"\nunique slices to embed: {len(uniq)}", flush=True)
    tok, model, spec = eam.load("dnabert2")
    E = eam.embed_unique(uniq, tok, model, spec, batch_size=16)
    lut = {s: E[i] for i, s in enumerate(uniq)}
    print(f"embedded {E.shape}", flush=True)

    res = defaultdict(dict)
    for b, pb in data.items():
        for span, v in pb.items():
            y = np.array(v["y"]); seqs = v["seq"]
            if len(set(y)) < 2 or min(np.bincount(y)) < N_SPLITS: continue
            km, d2 = [], []
            for tr, te in StratifiedKFold(N_SPLITS, shuffle=True,
                                          random_state=SEED).split(np.zeros(len(y)), y):
                if len(set(y[tr])) < 2 or len(set(y[te])) < 2: continue
                str_tr = [seqs[i] for i in tr]; str_te = [seqs[i] for i in te]
                Ktr, Kte = rr.kmer_fit_transform(str_tr, str_te)
                c = probe(); c.fit(Ktr, y[tr])
                km.append(rr.balanced_accuracy(y[te], c.predict(Kte)))
                p = probe(); p.fit(np.vstack([lut[s] for s in str_tr]), y[tr])
                d2.append(rr.balanced_accuracy(y[te], p.predict(np.vstack([lut[s] for s in str_te]))))
            if km:
                res[span][b] = {"C_kmer": float(np.mean(km)), "dnabert2": float(np.mean(d2))}

    def boot(x, n=20000, seed=1729):
        rng = np.random.default_rng(seed)
        dr = x[rng.integers(0, len(x), size=(n, len(x)))].mean(axis=1)
        return float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))

    summary = {}
    print(f"\n{'span':22}{'nt':>5}{'n':>4}{'k-mer':>9}{'DNABERT-2':>11}{'delta pp':>10}"
          f"{'95% CI pp':>20}  verdict")
    for span in SPANS:
        bs = sorted(res[span])
        if len(bs) < 5: continue
        k = np.array([res[span][b]["C_kmer"] for b in bs])
        d = np.array([res[span][b]["dnabert2"] for b in bs])
        x = d - k; lo, hi = boot(x)
        v = "EQUIVALENT" if (lo > -0.01 and hi < 0.01) else "INCONCLUSIVE at +-1pp"
        summary[span] = {"n": len(bs), "kmer": float(k.mean()), "dnabert2": float(d.mean()),
                         "delta": float(x.mean()), "ci95": [lo, hi], "verdict": v,
                         "expected_len": EXPECTED[span]}
        print(f"{span:22}{EXPECTED[span]:>5}{len(bs):>4}{100*k.mean():>9.2f}{100*d.mean():>11.2f}"
              f"{100*x.mean():>+10.2f}  [{100*lo:+7.2f},{100*hi:+7.2f}]  {v}")
    json.dump({"summary": summary, "per_block": {k: v for k, v in res.items()}},
              open(OUT, "w"), indent=2, default=float)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
