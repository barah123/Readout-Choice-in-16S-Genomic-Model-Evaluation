"""Phase 2 step 2: quantify what isolate-aware grouping changes vs. the Phase 0-3a
string-dedup approach. Directly tests Phase 1 audit findings 1 and 4.

Key question: how many records collapse to the same isolate but NOT the same sequence string?
Those are exactly the cases the old exact-string dedup treated as independent replicates and
could let straddle a train/test split.
"""
import json
from collections import Counter, defaultdict

for path, name in [("staph_panel_meta.json", "S. aureus / CoNS"),
                   ("cereus_panel_meta.json", "B. cereus group")]:
    recs = json.load(open(path))
    print("=" * 94)
    print(f"{name}  ({path})")
    print("=" * 94)
    print(f"raw records (1350-1600nt, ACGTN only): {len(recs)}")

    with_strain = sum(1 for r in recs if r["isolate_inferred"])
    print(f"  with a usable strain field: {with_strain} "
          f"({100*with_strain/len(recs):.0f}%)  |  without: {len(recs)-with_strain}")

    # how the two dedup policies differ
    uniq_seq = len({r["sequence"] for r in recs})
    uniq_iso = len({r["isolate_group"] for r in recs})
    print(f"  unique sequence strings : {uniq_seq}")
    print(f"  unique isolate groups   : {uniq_iso}")

    # the dangerous case: one isolate, multiple DIFFERENT sequence strings
    by_iso = defaultdict(set)
    for r in recs:
        by_iso[r["isolate_group"]].add(r["sequence"])
    multi_seq_iso = {k: v for k, v in by_iso.items() if len(v) > 1}
    n_records_affected = sum(
        sum(1 for r in recs if r["isolate_group"] == k) for k in multi_seq_iso)
    print(f"  isolates carrying >1 distinct sequence string: {len(multi_seq_iso)}"
          f"  (covering {n_records_affected} records)")
    if multi_seq_iso:
        ex = sorted(multi_seq_iso.items(), key=lambda kv: -len(kv[1]))[:5]
        for k, v in ex:
            print(f"      {k}: {len(v)} distinct strings")
    print("  ^ these are the records old string-dedup treated as independent but are not")

    # converse: one sequence string shared across multiple isolates (biologically real,
    # must be preserved, NOT merged)
    by_seq = defaultdict(set)
    for r in recs:
        by_seq[r["sequence"]].add(r["isolate_group"])
    shared = {k: v for k, v in by_seq.items() if len(v) > 1}
    print(f"  sequence strings shared by >1 isolate: {len(shared)} "
          f"(biologically real; preserved, not merged)")

    # per-species isolate counts - the real support for splitting
    per_sp = defaultdict(set)
    for r in recs:
        per_sp[r["species_group"]].add(r["isolate_group"])
    print("\n  independent isolates per species (the unit that governs split feasibility):")
    for sp, isos in sorted(per_sp.items(), key=lambda kv: -len(kv[1])):
        n_rec = sum(1 for r in recs if r["species_group"] == sp)
        print(f"    {sp:24} isolates={len(isos):4}  records={n_rec:4}")
    print()
