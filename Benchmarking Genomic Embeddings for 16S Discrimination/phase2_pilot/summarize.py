"""Phase 2 step 6: render the pilot results, with protocol §22's decision table applied
per cell rather than a single global verdict.

Deliberately NOT done here: confirmatory p-values, BH/Holm correction, equivalence testing.
Protocol §18 reserves those for Phase 5 on locked data; this is a development-phase pilot,
and computing them now then reusing the same panels later would invalidate them.
"""
import json
from collections import defaultdict

import numpy as np

res = json.load(open("readout_results.json"))
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
CLASSICAL = ["C_exact", "C_align", "C_kmer", "C_random"]


def fmt(v):
    return "  n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.3f}"


panels = sorted({k.split("|")[0] for k in res})

for panel in panels:
    keys = sorted([k for k in res if k.startswith(panel + "|")],
                  key=lambda k: int(k.split("|V")[1]))
    print("=" * 118)
    print(f"PANEL: {panel}")
    print("=" * 118)

    # ---- classical controls vs ceiling
    print("\nCLASSICAL CONTROLS (isolate-grouped CV, species-balanced accuracy)")
    print(f"{'region':8}{'ceiling':>9}{'n_tgt_iso':>11}" + "".join(f"{c:>11}" for c in CLASSICAL))
    for k in keys:
        r = res[k]
        reg = k.split("|")[1]
        print(f"{reg:8}{fmt(r['_ceiling_bal']):>9}{r['_n_target_isolates']:>11}"
              + "".join(f"{fmt(r.get(c)):>11}" for c in CLASSICAL))

    # ---- per readout level
    for level, label in [("R0_geometry_auc", "R0  intrinsic geometry (pair AUC, unsupervised)"),
                         ("R1_nn", "R1  nearest-reference lookup (frozen embeddings)"),
                         ("R2_linear", "R2  linear probe on frozen embeddings")]:
        print(f"\n{label}")
        print(f"{'region':8}{'ceiling':>9}" + "".join(f"{m:>12}" for m in MODELS)
              + f"{'C_kmer':>11}")
        for k in keys:
            r = res[k]
            reg = k.split("|")[1]
            print(f"{reg:8}{fmt(r['_ceiling_bal']):>9}"
                  + "".join(f"{fmt(r.get(f'{m}:{level}')):>12}" for m in MODELS)
                  + f"{fmt(r.get('C_kmer')):>11}")
    print()

# ------------------------------------------------------------------ headline contrast
print("=" * 118)
print("HEADLINE CONTRAST: best frozen-model readout vs best classical, per panel x region")
print("=" * 118)
print(f"{'panel':16}{'region':8}{'ceiling':>9}{'best_classical':>24}{'best_model_readout':>34}{'delta':>9}")
rows = []
for k in sorted(res):
    r = res[k]
    panel, reg = k.split("|")
    cls = {c: r.get(c) for c in CLASSICAL if c != "C_random" and r.get(c) == r.get(c)}
    mdl = {f"{m}:{lv}": r.get(f"{m}:{lv}")
           for m in MODELS for lv in ("R1_nn", "R2_linear")
           if r.get(f"{m}:{lv}") == r.get(f"{m}:{lv}")}
    if not cls or not mdl:
        continue
    bc = max(cls, key=cls.get)
    bm = max(mdl, key=mdl.get)
    delta = mdl[bm] - cls[bc]
    rows.append((panel, reg, delta))
    print(f"{panel:16}{reg:8}{fmt(r['_ceiling_bal']):>9}"
          f"{bc + ' ' + fmt(cls[bc]):>24}{bm + ' ' + fmt(mdl[bm]):>34}{delta:>+9.3f}")

print()
for panel in panels:
    d = [x[2] for x in rows if x[0] == panel]
    if d:
        wins = sum(1 for v in d if v > 0.01)
        losses = sum(1 for v in d if v < -0.01)
        print(f"{panel}: model-best beats classical-best in {wins}/{len(d)} regions, "
              f"loses in {losses}, within +/-0.01 in {len(d)-wins-losses}. "
              f"mean delta = {np.mean(d):+.3f}")

print("\nNote: development-phase pilot. No confirmatory p-values or multiplicity correction "
      "here by design (protocol §18 reserves those for Phase 5 on locked data).")
