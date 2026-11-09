"""Extension A step 1: screen candidate panels for DIFFICULTY, using classical methods only.

Selection criteria (declared in PREREGISTRATION.md before running):
  a panel is HARD if, at its best region,
    (1) observed-sequence ambiguity ceiling <= 0.97, AND
    (2) k-mer + linear probe balanced accuracy <= 0.95.

No genomic language model is used anywhere in this file. Selection therefore cannot favour
the models.

Candidates are genera/complexes chosen for documented 16S conservation, plus several
subspecies-level panels (inherently hard by construction). Families already used in Phases
2-6 are avoided where the biology allows; where a hard group sits in a used family that
overlap is recorded in the `family` field and reported.

Screening uses a cheap sample (retmax=100). Panels passing the screen get a full fetch later.
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.abspath("../phase2_pilot"))
import build_blocks_full as bb
import fetch_blocks
import run_readouts_full as rr

fetch_blocks.OUTDIR = "raw_screen"
fetch_blocks.RETMAX = 100          # cheap screening sample
SEED = 0
CEILING_MAX = 0.97                 # pre-declared
KMER_MAX = 0.95                    # pre-declared

CANDIDATES = [
    # --- obligate intracellular / highly conserved 16S ---
    dict(block_id="rickettsia", family="Rickettsiaceae", target="Rickettsia rickettsii",
         neighbours=["Rickettsia conorii", "Rickettsia prowazekii", "Rickettsia typhi",
                     "Rickettsia felis"]),
    dict(block_id="ehrlichia_anaplasma", family="Anaplasmataceae",
         target="Ehrlichia chaffeensis",
         neighbours=["Ehrlichia canis", "Anaplasma phagocytophilum",
                     "Anaplasma marginale", "Wolbachia pipientis"]),
    dict(block_id="chlamydia", family="Chlamydiaceae", target="Chlamydia trachomatis",
         neighbours=["Chlamydia muridarum", "Chlamydia suis", "Chlamydia psittaci",
                     "Chlamydia pneumoniae"]),
    dict(block_id="coxiella_rickettsiella", family="Coxiellaceae",
         target="Coxiella burnetii",
         neighbours=["Rickettsiella grylli", "Diplorickettsia massiliensis"]),
    # --- spirochaetes ---
    dict(block_id="treponema", family="Spirochaetaceae", target="Treponema pallidum",
         neighbours=["Treponema denticola", "Treponema phagedenis",
                     "Treponema socranskii"]),
    dict(block_id="brachyspira", family="Brachyspiraceae",
         target="Brachyspira hyodysenteriae",
         neighbours=["Brachyspira pilosicoli", "Brachyspira innocens",
                     "Brachyspira murdochii"]),
    # --- mollicutes ---
    dict(block_id="mycoplasma", family="Mycoplasmataceae", target="Mycoplasma pneumoniae",
         neighbours=["Mycoplasma genitalium", "Mycoplasma hominis",
                     "Ureaplasma urealyticum", "Mycoplasma fermentans"]),
    dict(block_id="spiroplasma", family="Spiroplasmataceae", target="Spiroplasma citri",
         neighbours=["Spiroplasma melliferum", "Spiroplasma kunkelii",
                     "Spiroplasma apis"]),
    # --- plant-pathogen complexes, close species ---
    dict(block_id="pectobacterium_dickeya", family="Pectobacteriaceae",
         target="Pectobacterium carotovorum",
         neighbours=["Pectobacterium atrosepticum", "Dickeya dadantii",
                     "Dickeya solani", "Pectobacterium brasiliense"]),
    dict(block_id="erwinia_pantoea", family="Erwiniaceae", target="Erwinia amylovora",
         neighbours=["Erwinia pyrifoliae", "Pantoea agglomerans",
                     "Pantoea ananatis", "Erwinia tasmaniensis"]),
    dict(block_id="xylella", family="Xanthomonadaceae*", target="Xylella fastidiosa",
         neighbours=["Xanthomonas oryzae", "Xanthomonas citri",
                     "Xylella taiwanensis"]),
    # --- conserved environmental groups ---
    dict(block_id="bradyrhizobium", family="Nitrobacteraceae",
         target="Bradyrhizobium japonicum",
         neighbours=["Bradyrhizobium elkanii", "Bradyrhizobium diazoefficiens",
                     "Nitrobacter winogradskyi"]),
    dict(block_id="synechococcus_prochlorococcus", family="Synechococcaceae",
         target="Prochlorococcus marinus",
         neighbours=["Synechococcus elongatus", "Synechococcus sp.",
                     "Cyanobium gracile"]),
    dict(block_id="nostoc_anabaena", family="Nostocaceae", target="Nostoc punctiforme",
         neighbours=["Anabaena variabilis", "Nostoc commune",
                     "Aphanizomenon flos-aquae"]),
    dict(block_id="desulfovibrio", family="Desulfovibrionaceae",
         target="Desulfovibrio vulgaris",
         neighbours=["Desulfovibrio desulfuricans", "Desulfovibrio piger",
                     "Desulfovibrio gigas"]),
    dict(block_id="geobacter", family="Geobacteraceae", target="Geobacter sulfurreducens",
         neighbours=["Geobacter metallireducens", "Geobacter uraniireducens",
                     "Pelobacter carbinolicus"]),
    dict(block_id="nitrosomonas", family="Nitrosomonadaceae",
         target="Nitrosomonas europaea",
         neighbours=["Nitrosomonas eutropha", "Nitrosospira multiformis",
                     "Nitrosomonas communis"]),
    dict(block_id="rhodobacter_paracoccus", family="Rhodobacteraceae",
         target="Paracoccus denitrificans",
         neighbours=["Rhodobacter sphaeroides", "Rhodobacter capsulatus",
                     "Paracoccus yeei"]),
    dict(block_id="halomonas", family="Halomonadaceae", target="Halomonas elongata",
         neighbours=["Halomonas halodenitrificans", "Chromohalobacter salexigens",
                     "Halomonas meridiana"]),
    dict(block_id="alteromonas", family="Alteromonadaceae",
         target="Alteromonas macleodii",
         neighbours=["Pseudoalteromonas haloplanktis", "Pseudoalteromonas atlantica",
                     "Marinobacter hydrocarbonoclasticus"]),
    dict(block_id="methylococcus", family="Methylococcaceae",
         target="Methylococcus capsulatus",
         neighbours=["Methylomonas methanica", "Methylomicrobium album",
                     "Methylobacter tundripaludum"]),
    dict(block_id="thermus_deinococcus", family="Thermaceae",
         target="Thermus thermophilus",
         neighbours=["Thermus aquaticus", "Deinococcus radiodurans",
                     "Meiothermus ruber"]),
    dict(block_id="arcobacter", family="Arcobacteraceae", target="Arcobacter butzleri",
         neighbours=["Arcobacter cryaerophilus", "Arcobacter skirrowii",
                     "Arcobacter nitrofigilis"]),
    dict(block_id="akkermansia_verruco", family="Akkermansiaceae",
         target="Akkermansia muciniphila",
         neighbours=["Verrucomicrobium spinosum", "Akkermansia glycaniphila"]),
    dict(block_id="eggerthella_atopobium", family="Eggerthellaceae",
         target="Eggerthella lenta",
         neighbours=["Atopobium vaginae", "Olsenella uli", "Gordonibacter pamelaeae"]),
    dict(block_id="tannerella", family="Tannerellaceae", target="Tannerella forsythia",
         neighbours=["Parabacteroides distasonis", "Parabacteroides merdae"]),
    dict(block_id="myxococcus", family="Myxococcaceae", target="Myxococcus xanthus",
         neighbours=["Corallococcus coralloides", "Stigmatella aurantiaca",
                     "Myxococcus fulvus"]),
    # --- subspecies-level panels: hard by construction ---
    dict(block_id="mycobacterium_abscessus_subsp", family="Mycobacteriaceae*",
         target="Mycobacterium abscessus subsp. abscessus",
         neighbours=["Mycobacterium abscessus subsp. massiliense",
                     "Mycobacterium abscessus subsp. bolletii"]),
    dict(block_id="francisella_tularensis_subsp", family="Francisellaceae*",
         target="Francisella tularensis subsp. tularensis",
         neighbours=["Francisella tularensis subsp. holarctica",
                     "Francisella tularensis subsp. mediasiatica"]),
    dict(block_id="lactobacillus_delbrueckii_subsp", family="Lactobacillaceae*",
         target="Lactobacillus delbrueckii subsp. bulgaricus",
         neighbours=["Lactobacillus delbrueckii subsp. lactis",
                     "Lactobacillus delbrueckii subsp. delbrueckii"]),
]


def probe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=2000, class_weight="balanced",
                                            random_state=SEED))


def kmer_cv(slices, y, n_splits=5):
    y = np.asarray(y)
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < 3:
        return float("nan")
    ns = min(n_splits, int(min(np.bincount(y))))
    skf = StratifiedKFold(n_splits=ns, shuffle=True, random_state=SEED)
    out = []
    for tr, te in skf.split(np.zeros(len(y)), y):
        Ktr, Kte = rr.kmer_fit_transform([slices[i] for i in tr], [slices[i] for i in te])
        c = probe(); c.fit(Ktr, y[tr])
        out.append(rr.balanced_accuracy(y[te], c.predict(Kte)))
    return float(np.mean(out)) if out else float("nan")


def main():
    os.makedirs("raw_screen", exist_ok=True)
    results = {}

    for cand in CANDIDATES:
        bid = cand["block_id"]
        recs = []
        for sp, is_t in [(cand["target"], True)] + [(n, False) for n in cand["neighbours"]]:
            try:
                r, _ = fetch_blocks.fetch_species(bid, sp, is_t)
                recs.extend(r)
            except Exception as e:
                print(f"  [fetch fail] {sp}: {e}", flush=True)
        if not recs:
            results[bid] = dict(status="no_records", **{k: cand[k] for k in ("family",)})
            print(f"{bid:34} NO RECORDS", flush=True)
            continue

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
            # "best region" = the one the classical arm does best at
            if best is None or km > best["kmer"]:
                best = dict(region=region, ceiling=ceil, kmer=km,
                            n_target=n_t, n_neighbour=n_n)

        if best is None:
            results[bid] = dict(status="insufficient_data", family=cand["family"])
            print(f"{bid:34} INSUFFICIENT DATA", flush=True)
            continue

        hard = bool(best["ceiling"] <= CEILING_MAX and best["kmer"] <= KMER_MAX)
        results[bid] = dict(status="screened", family=cand["family"], hard=hard,
                            target=cand["target"], neighbours=cand["neighbours"], **best)
        print(f"{bid:34} reg={best['region']:3} ceiling={best['ceiling']:.3f} "
              f"kmer={best['kmer']:.3f} n_t={best['n_target']:4} "
              f"{'<== HARD' if hard else ''}", flush=True)

    json.dump(results, open("screen_results.json", "w"), indent=2, default=float)
    hard = [k for k, v in results.items() if v.get("hard")]
    print(f"\nselected HARD panels: {len(hard)}/{len(CANDIDATES)}")
    print(f"target >=12 for planned precision; <8 => underpowered")
    for h in hard:
        print(f"  {h}")


if __name__ == "__main__":
    main()
