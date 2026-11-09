"""Phase 4 step 3: embed locked-test slices with the four frozen checkpoints.

Reuses embed_all_models.load / embed_unique unchanged — identical pooling, identical
unique-slice dedup, identical label-independent batch ordering as development.

Run:  conda run -n glm-benchmark python embed_test.py dnabert2 dnabert_s hyenadna
      conda run -n glm-ntv2      python embed_test.py ntv2
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
from embed_all_models import REGIONS, embed_unique, load


def main(models):
    rows = json.load(open("test_blocks.json"))
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
            assert np.isfinite(E).all(), f"{name} {region}: non-finite embeddings"
            store[f"{region}__keys"] = np.array(uniq, dtype=object)
            store[f"{region}__vecs"] = E
            meta[region] = dict(n_unique=len(uniq), dim=int(E.shape[1]))
            print(f"{name} {region}: {len(uniq)} unique -> dim {E.shape[1]} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        np.savez_compressed(f"emb_test_{name}.npz", **store)
        meta["_elapsed_s"] = round(time.time() - t0, 1)
        meta["_model_path"] = spec["path"]
        json.dump(meta, open(f"emb_test_{name}_meta.json", "w"), indent=2)
        print(f"saved emb_test_{name}.npz ({meta['_elapsed_s']}s)\n", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["dnabert2"])
