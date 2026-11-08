"""Embed every unique slice in the full 31-block set, reusing the validated embedding path
from embed_all_models.py (unique-slice dedup -> pure function of sequence by construction,
label-independent batch ordering, masked mean pooling excluding each row's own specials).

Run:  conda run -n glm-benchmark python embed_full.py dnabert2 dnabert_s hyenadna
      conda run -n glm-ntv2      python embed_full.py ntv2
"""
import json
import sys
import time

import numpy as np

from embed_all_models import REGIONS, embed_unique, load


def main(models):
    rows = json.load(open("blocks_full.json"))
    by_region = {}
    for r in rows:
        by_region.setdefault(r["region"], set()).add(r["slice"])
    print({k: len(v) for k, v in sorted(by_region.items())}, flush=True)

    for name in models:
        tok, model, spec = load(name)
        store, meta = {}, {}
        t0 = time.time()
        for region in REGIONS:
            uniq = sorted(by_region.get(region, []))
            if not uniq:
                continue
            E = embed_unique(uniq, tok, model, spec)
            assert E.shape[0] == len(uniq)
            store[f"{region}__keys"] = np.array(uniq, dtype=object)
            store[f"{region}__vecs"] = E
            meta[region] = dict(n_unique=len(uniq), dim=int(E.shape[1]))
            print(f"{name} {region}: {len(uniq)} unique -> dim {E.shape[1]} "
                  f"({time.time()-t0:.0f}s elapsed)", flush=True)
        np.savez_compressed(f"emb_full_{name}.npz", **store)
        meta["_elapsed_s"] = round(time.time() - t0, 1)
        meta["_model_path"] = spec["path"]
        json.dump(meta, open(f"emb_full_{name}_meta.json", "w"), indent=2)
        print(f"saved emb_full_{name}.npz ({meta['_elapsed_s']}s)\n", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["dnabert2"])
