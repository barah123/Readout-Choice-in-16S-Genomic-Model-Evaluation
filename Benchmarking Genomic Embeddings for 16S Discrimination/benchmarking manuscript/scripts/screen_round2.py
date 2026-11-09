"""Extension A, screening round 2.

WHY THIS EXISTS (amendment logged in PREREGISTRATION.md before any model touched these
panels): round 1 screened 30 candidates with a deliberately cheap sample (retmax=100).
11 returned INSUFFICIENT DATA - a data-collection shortfall, not a finding about difficulty.
Round 2 re-screens those 11 at retmax=300, and adds 12 further candidates of the type round 1
showed to be reliably hard: subspecies-level and tight species-complex panels.

SELECTION CRITERIA ARE UNCHANGED (ceiling <= 0.97 AND k-mer <= 0.95). Screening remains
classical-only; no genomic language model is involved anywhere in this file.
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
import build_blocks_full as bb
import fetch_blocks

from screen_hard_panels import CANDIDATES, CEILING_MAX, KMER_MAX, kmer_cv

fetch_blocks.OUTDIR = "raw_screen2"
fetch_blocks.RETMAX = 300          # deeper sample than round 1

RETRY = ["ehrlichia_anaplasma", "coxiella_rickettsiella", "treponema", "mycoplasma",
         "geobacter", "nitrosomonas", "methylococcus", "arcobacter",
         "akkermansia_verruco", "tannerella", "francisella_tularensis_subsp"]

NEW = [
    dict(block_id="lactobacillus_casei_group", family="Lactobacillaceae*",
         target="Lacticaseibacillus casei",
         neighbours=["Lacticaseibacillus paracasei", "Lacticaseibacillus rhamnosus",
                     "Lacticaseibacillus zeae"]),
    dict(block_id="lactococcus_lactis_subsp", family="Streptococcaceae*",
         target="Lactococcus lactis subsp. lactis",
         neighbours=["Lactococcus lactis subsp. cremoris",
                     "Lactococcus lactis subsp. hordniae"]),
    dict(block_id="bacillus_subtilis_subsp", family="Bacillaceae*",
         target="Bacillus subtilis subsp. subtilis",
         neighbours=["Bacillus subtilis subsp. spizizenii",
                     "Bacillus subtilis subsp. inaquosorum",
                     "Bacillus amyloliquefaciens"]),
    dict(block_id="bifidobacterium_animalis_subsp", family="Bifidobacteriaceae*",
         target="Bifidobacterium animalis subsp. lactis",
         neighbours=["Bifidobacterium animalis subsp. animalis",
                     "Bifidobacterium pseudolongum"]),
    dict(block_id="salmonella_subspecies", family="Enterobacteriaceae*",
         target="Salmonella enterica subsp. enterica",
         neighbours=["Salmonella enterica subsp. salamae",
                     "Salmonella enterica subsp. arizonae",
                     "Salmonella enterica subsp. diarizonae",
                     "Salmonella enterica subsp. houtenae"]),
    dict(block_id="streptococcus_dysgalactiae_subsp", family="Streptococcaceae*",
         target="Streptococcus dysgalactiae subsp. equisimilis",
         neighbours=["Streptococcus dysgalactiae subsp. dysgalactiae",
                     "Streptococcus canis", "Streptococcus pyogenes"]),
    dict(block_id="campylobacter_jejuni_subsp", family="Campylobacteraceae*",
         target="Campylobacter jejuni subsp. jejuni",
         neighbours=["Campylobacter jejuni subsp. doylei", "Campylobacter coli"]),
    dict(block_id="clostridium_botulinum_groups", family="Clostridiaceae*",
         target="Clostridium botulinum",
         neighbours=["Clostridium sporogenes", "Clostridium novyi",
                     "Clostridium haemolyticum"]),
    dict(block_id="burkholderia_cepacia_complex", family="Burkholderiaceae*",
         target="Burkholderia cenocepacia",
         neighbours=["Burkholderia multivorans", "Burkholderia vietnamiensis",
                     "Burkholderia ambifaria", "Burkholderia stabilis"]),
    dict(block_id="yersinia_enterocolitica_subsp", family="Yersiniaceae*",
         target="Yersinia enterocolitica subsp. enterocolitica",
         neighbours=["Yersinia enterocolitica subsp. palearctica",
                     "Yersinia kristensenii", "Yersinia intermedia"]),
    dict(block_id="xanthomonas_campestris_pv", family="Xanthomonadaceae*",
         target="Xanthomonas campestris",
         neighbours=["Xanthomonas axonopodis", "Xanthomonas vesicatoria",
                     "Xanthomonas arboricola"]),
    dict(block_id="pseudomonas_fluorescens_group", family="Pseudomonadaceae*",
         target="Pseudomonas fluorescens",
         neighbours=["Pseudomonas putida", "Pseudomonas chlororaphis",
                     "Pseudomonas protegens", "Pseudomonas synxantha"]),
]


def screen_one(cand):
    bid = cand["block_id"]
    recs = []
    for sp, is_t in [(cand["target"], True)] + [(n, False) for n in cand["neighbours"]]:
        try:
            r, _ = fetch_blocks.fetch_species(bid, sp, is_t)
            recs.extend(r)
        except Exception as e:
            print(f"  [fetch fail] {sp}: {e}", flush=True)
    if not recs:
        return dict(status="no_records", family=cand["family"])

    sliced = Parallel(n_jobs=8, batch_size=16)(
        delayed(bb._slice_one)(r["sequence"]) for r in recs)
    for r, sl in zip(recs, sliced):
        r["slices"] = sl

    best = None
    for region in bb.REGIONS:
        present = [r for r in recs if r["slices"][region]]
        if len(present) < 20:
            continue
        by_iso = defaultdict(list)
        for r in present:
            by_iso[r["isolate_group"]].append(r)
        iso = []
        for k, rs in by_iso.items():
            modal = Counter(x["slices"][region] for x in rs).most_common(1)[0][0]
            iso.append((rs[0]["is_target"], modal))
        n_t = sum(1 for t, _ in iso if t)
        n_n = len(iso) - n_t
        if n_t < 10 or n_n < 10:
            continue
        cnt = defaultdict(lambda: [0, 0])
        for t, sl in iso:
            cnt[sl][0 if t else 1] += 1
        _, ceil = bb.ceilings(cnt, n_t, n_n)
        km = kmer_cv([sl for _, sl in iso], [int(t) for t, _ in iso])
        if km != km:
            continue
        if best is None or km > best["kmer"]:
            best = dict(region=region, ceiling=ceil, kmer=km,
                        n_target=n_t, n_neighbour=n_n)

    if best is None:
        return dict(status="insufficient_data", family=cand["family"])
    hard = bool(best["ceiling"] <= CEILING_MAX and best["kmer"] <= KMER_MAX)
    return dict(status="screened", family=cand["family"], hard=hard,
                target=cand["target"], neighbours=cand["neighbours"], **best)


def main():
    os.makedirs("raw_screen2", exist_ok=True)
    by_id = {c["block_id"]: c for c in CANDIDATES}
    todo = [by_id[b] for b in RETRY if b in by_id] + NEW
    print(f"round 2: {len(todo)} candidates "
          f"({len([b for b in RETRY if b in by_id])} re-screened deeper, {len(NEW)} new)\n",
          flush=True)

    results = {}
    for cand in todo:
        res = screen_one(cand)
        results[cand["block_id"]] = res
        if res["status"] == "screened":
            print(f"{cand['block_id']:34} reg={res['region']:3} "
                  f"ceiling={res['ceiling']:.3f} kmer={res['kmer']:.3f} "
                  f"n_t={res['n_target']:4} {'<== HARD' if res['hard'] else ''}", flush=True)
        else:
            print(f"{cand['block_id']:34} {res['status'].upper()}", flush=True)

    json.dump(results, open("screen_results_round2.json", "w"), indent=2, default=float)
    r1 = json.load(open("screen_results.json"))
    hard1 = [k for k, v in r1.items() if v.get("hard")]
    hard2 = [k for k, v in results.items() if v.get("hard")]
    print(f"\nround 1 hard: {len(hard1)}  round 2 hard: {len(hard2)}  "
          f"TOTAL: {len(hard1) + len(hard2)}")
    print("target >=12 for planned precision; <8 => report underpowered")
    for h in hard1 + hard2:
        print(f"  {h}")


if __name__ == "__main__":
    main()
