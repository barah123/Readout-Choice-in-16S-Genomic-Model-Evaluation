"""Extension C: is an intermediate layer better than the final layer?

Pre-registered in PREREGISTRATION.md. Tests HyenaDNA and NTv2 (both verified to expose
intermediate layers correctly) on the 13 screened-hard panels, at 5 depths each, with
everything else frozen: same masked-mean pooling, same linear probe, same isolate-grouped
5-fold CV, same seed. Only the layer changes.

A correctness gate aborts the run unless the extracted final layer reproduces the standard
forward output - the check whose absence caused this analysis to be deferred for DNABERT-2.

NTv2 lives in `glm-ntv2`, which is pinned to transformers 4.34.1 and has no sklearn. Rather
than add packages to a verified environment, embedding and scoring are separate steps, the
same split used everywhere else in this project.

Run:  conda run -n glm-benchmark python extension_C_layers.py hyenadna
      conda run -n glm-ntv2      python extension_C_layers.py ntv2 embed
      conda run -n glm-benchmark python extension_C_layers.py ntv2 score
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import torch

sys.path.insert(0, os.path.abspath("../phase2_pilot"))

SEED = 0
N_DEPTHS = 5

SPECS = {
    "hyenadna": dict(path="LongSafari/hyenadna-medium-160k-seqlen-hf", max_length=2048),
    "ntv2": dict(path="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
                 max_length=512),
}


def load(name):
    spec = SPECS[name]
    if name == "ntv2":
        from transformers import AutoModelForMaskedLM, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(spec["path"], trust_remote_code=True)
        full = AutoModelForMaskedLM.from_pretrained(spec["path"], trust_remote_code=True)
        full.eval()
        return tok, full.esm, spec
    from transformers import AutoModel, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(spec["path"], trust_remote_code=True)
    mod = AutoModel.from_pretrained(spec["path"], trust_remote_code=True)
    mod.eval()
    return tok, mod, spec


def correctness_gate(tok, enc_model, spec):
    """ABORT unless the extracted final layer equals the standard forward output."""
    enc = tok(["ACGT" * 30], return_tensors="pt", truncation=True,
              max_length=spec["max_length"])
    with torch.no_grad():
        base = enc_model(**enc)[0]
        out = enc_model(**enc, output_hidden_states=True)
    hs = out.hidden_states
    ok = torch.allclose(hs[-1], base, atol=1e-6)
    if not ok:
        raise RuntimeError("CORRECTNESS GATE FAILED: extracted final layer does not match "
                           "the standard forward output. Refusing to produce numbers from a "
                           "representation that may be silently wrong.")
    return len(hs)


def embed_layers(seqs, tok, enc_model, spec, depths, batch_size=16):
    acc = defaultdict(list)
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i + batch_size]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True,
                  max_length=spec["max_length"])
        with torch.no_grad():
            hs = enc_model(**enc, output_hidden_states=True).hidden_states
        am = (enc["attention_mask"].bool() if "attention_mask" in enc
              else torch.ones(hs[-1].shape[:2], dtype=torch.bool))
        mask = am.clone()
        L = am.sum(1)
        mask[:, 0] = False
        mask[torch.arange(mask.size(0)), (L - 1).clamp(min=0)] = False
        empty = mask.sum(1) == 0
        if empty.any():
            mask[empty] = am[empty]
        m = mask.unsqueeze(-1).float()
        for d in depths:
            h = hs[d]
            acc[d].append(((h * m).sum(1) / m.sum(1).clamp(min=1)).float().numpy())
    return {d: np.vstack(v) for d, v in acc.items()}


def panel_inputs():
    """The (block, region, isolates, labels, slices) the frozen analysis used."""
    rows = json.load(open("hard_blocks.json"))
    summ = json.load(open("hard_block_summary.json"))
    scores = json.load(open("hard_scores.json"))
    eligible = [b for b, v in summ.items() if v.get("block_eligible")]
    by_block = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by_block[r["block_id"]][r["region"]][r["isolate_group"]] = (r["slice"], r["is_target"])
    out = []
    for bid in sorted(eligible):
        v = scores.get(bid, {})
        if v.get("_status") != "ok":
            continue
        region = Counter(v["_chosen_regions"]).most_common(1)[0][0]
        d = by_block[bid][region]
        isolates = sorted(d)
        y = np.array([d[i][1] for i in isolates])
        if len(np.unique(y)) < 2 or min(np.bincount(y)) < 5:
            continue
        out.append((bid, region, [d[i][0] for i in isolates], y))
    return out


def do_embed(name):
    """Embed every hard panel at the chosen depths and cache to disk."""
    tok, enc_model, spec = load(name)
    n_hs = correctness_gate(tok, enc_model, spec)
    last = n_hs - 1
    depths = sorted({int(round(last * f)) for f in (0.25, 0.5, 0.75, 0.9, 1.0)})
    print(f"{name}: {n_hs} hidden-state entries; depths {depths} (last={last})", flush=True)
    store = {"__depths__": np.array(depths), "__last__": np.array([last])}
    for bid, region, slices, _ in panel_inputs():
        uniq = sorted(set(slices))
        mats = embed_layers(uniq, tok, enc_model, spec, depths)
        store[f"{bid}__keys"] = np.array(uniq, dtype=object)
        for dep, M in mats.items():
            store[f"{bid}__L{dep}"] = M
        print(f"  {bid:34} {len(uniq)} unique slices", flush=True)
    np.savez_compressed(f"extC_emb_{name}.npz", **store)
    print(f"saved extC_emb_{name}.npz")


def main(name, cached=False):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    import run_readouts_full as rr

    def probe():
        return make_pipeline(StandardScaler(),
                             LogisticRegression(max_iter=2000, class_weight="balanced",
                                                random_state=SEED))

    if cached:
        z = np.load(f"extC_emb_{name}.npz", allow_pickle=True)
        depths = [int(d) for d in z["__depths__"]]
        last = int(z["__last__"][0])
        tok = enc_model = spec = None
        print(f"{name}: using cached embeddings, depths {depths} (last={last}, "
              f"frozen analysis used the final layer)\n", flush=True)
    else:
        tok, enc_model, spec = load(name)
        n_hs = correctness_gate(tok, enc_model, spec)
        last = n_hs - 1
        depths = sorted({int(round(last * f)) for f in (0.25, 0.5, 0.75, 0.9, 1.0)})
        z = None
        print(f"{name}: {n_hs} hidden-state entries; testing depths {depths} "
              f"(last={last}, frozen analysis used the final layer)\n", flush=True)

    out = {}
    for bid, region, slices, y in panel_inputs():
        uniq = sorted(set(slices))
        if cached:
            keys = list(z[f"{bid}__keys"])
            luts = {dep: {k: v for k, v in zip(keys, z[f"{bid}__L{dep}"])}
                    for dep in depths}
        else:
            mats = embed_layers(uniq, tok, enc_model, spec, depths)
            luts = {dep: {s: M[j] for j, s in enumerate(uniq)} for dep, M in mats.items()}

        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
        per = defaultdict(list)
        kmer = []
        for tr, te in skf.split(np.zeros(len(y)), y):
            ytr, yte = y[tr], y[te]
            if len(set(ytr)) < 2 or len(set(yte)) < 2:
                continue
            s_tr = [slices[i] for i in tr]
            s_te = [slices[i] for i in te]
            Ktr, Kte = rr.kmer_fit_transform(s_tr, s_te)
            c = probe(); c.fit(Ktr, ytr)
            kmer.append(rr.balanced_accuracy(yte, c.predict(Kte)))
            for dep, lut in luts.items():
                Xtr = np.vstack([lut[s] for s in s_tr])
                Xte = np.vstack([lut[s] for s in s_te])
                p = probe(); p.fit(Xtr, ytr)
                per[dep].append(rr.balanced_accuracy(yte, p.predict(Xte)))
        out[bid] = {"region": region, "C_kmer": float(np.mean(kmer)),
                    **{f"layer_{dep}": float(np.mean(per[dep])) for dep in depths}}
        best = max(depths, key=lambda dd: np.mean(per[dd]))
        print(f"{bid:34} kmer={np.mean(kmer):.3f}  "
              + "  ".join(f"L{dep}={np.mean(per[dep]):.3f}" for dep in depths)
              + f"   best=L{best}", flush=True)

    json.dump({"model": name, "depths": depths, "last": last, "per_block": out},
              open(f"extension_C_{name}.json", "w"), indent=2, default=float)

    print()
    print("=" * 96)
    print(f"{name}: mean over {len(out)} hard panels")
    print("=" * 96)
    means = {dep: np.mean([out[b][f"layer_{dep}"] for b in out]) for dep in depths}
    km = np.mean([out[b]["C_kmer"] for b in out])
    for dep in depths:
        tag = "  <- FINAL (used in the frozen analysis)" if dep == last else ""
        print(f"  layer {dep:3}: {means[dep]:.4f}   vs k-mer {km:.4f}  "
              f"delta {means[dep]-km:+.4f}{tag}")
    best_dep = max(depths, key=lambda d: means[d])
    gain = means[best_dep] - means[last]
    print(f"\n  best depth = {best_dep}; gain over final layer = {gain:+.4f}")

    d = np.array([out[b][f"layer_{best_dep}"] - out[b][f"layer_{last}"] for b in out])
    rng = np.random.default_rng(1729)
    dr = d[rng.integers(0, len(d), size=(2000, len(d)))].mean(axis=1)
    lo, hi = np.percentile(dr, [2.5, 97.5])
    print(f"  95% CI on that gain = [{lo:+.4f}, {hi:+.4f}]")
    print("  NOTE: best depth was chosen using these same panels, so this gain is "
          "upward-biased\n        by selection and is an upper bound, not an estimate.")


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "hyenadna"
    mode = sys.argv[2] if len(sys.argv) > 2 else "run"
    if mode == "embed":
        do_embed(name)
    elif mode == "score":
        main(name, cached=True)
    else:
        main(name, cached=False)
