"""Audit the isolate-grouping key: coverage, how much the grouping binds, and how much of
it depends on strain-string formatting that was deliberately left un-normalised.

Produces every isolate-grouping number quoted in Methods ("Isolate grouping"), so each is
reproducible from the fetched records rather than retyped. Run from anywhere.

The aggressive-renormalisation arm is a sensitivity check only: it upper-cases each strain
string and deletes every non-alphanumeric character, then re-derives the groups. The
difference bounds the effect of the formatting component of normalisation. It does not bound
catalogue equivalence (one isolate deposited under an ATCC number and its DSM equivalent),
which cannot be resolved from the records themselves.
"""
import glob
import json
import os
import re
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

SETS = {
    "development": "phase2_pilot/raw_species/*.json",
    "locked_test": "phase4_locked_test/raw_species_test/*.json",
    "hard_panels": "phase7_extensions/raw_species_hard/*.json",
}


def aggressive(s):
    return re.sub(r"[^A-Z0-9]", "", s.upper())


def main():
    out, tot = {}, Counter()
    for name, pat in SETS.items():
        recs = []
        for f in sorted(glob.glob(os.path.join(ROOT, pat))):
            recs.extend(json.load(open(f)))
        if not recs:
            print(f"[skip] {name}: no records under {pat}")
            continue

        groups = Counter(r["isolate_group"] for r in recs)
        multi = {k: v for k, v in groups.items() if v > 1}
        agg = {f'{r["species_query"]}|{aggressive(r["strain"])}' if r.get("strain")
               else r["isolate_group"] for r in recs}

        acc_sp = defaultdict(set)
        for r in recs:
            acc_sp[r["accession"]].add((r["block_id"], r["species_query"]))
        cross = [a for a, v in acc_sp.items() if len({sp for _, sp in v}) > 1]
        same_block = [a for a in cross if len({b for b, _ in acc_sp[a]}) == 1]

        d = dict(
            records=len(recs),
            accessions=len(acc_sp),
            records_with_strain=sum(1 for r in recs if r.get("strain")),
            isolate_groups=len(groups),
            groups_after_aggressive_norm=len(agg),
            groups_with_multiple_records=len(multi),
            records_in_multi_record_groups=sum(multi.values()),
            largest_group=max(groups.values()),
            cross_species_accessions=len(cross),
            cross_species_within_one_block=len(same_block),
        )
        d["strain_coverage_pct"] = round(100 * d["records_with_strain"] / d["records"], 2)
        d["groups_merged_by_aggressive_norm"] = (d["isolate_groups"]
                                                 - d["groups_after_aggressive_norm"])
        out[name] = d
        for k in ("records", "records_with_strain", "isolate_groups",
                  "groups_after_aggressive_norm", "groups_with_multiple_records",
                  "records_in_multi_record_groups", "cross_species_accessions",
                  "cross_species_within_one_block"):
            tot[k] += d[k]
        print(f"{name:12} " + "  ".join(f"{k}={v}" for k, v in d.items()))

    tot = dict(tot)
    tot["strain_coverage_pct"] = round(100 * tot["records_with_strain"] / tot["records"], 2)
    tot["groups_merged_by_aggressive_norm"] = (tot["isolate_groups"]
                                               - tot["groups_after_aggressive_norm"])
    tot["pct_groups_merged_by_aggressive_norm"] = round(
        100 * tot["groups_merged_by_aggressive_norm"] / tot["isolate_groups"], 2)
    out["total_analysed"] = tot
    print("\ntotal        " + "  ".join(f"{k}={v}" for k, v in tot.items()))

    dest = os.path.join(HERE, "..", "results", "isolate_grouping_audit.json")
    json.dump(out, open(os.path.abspath(dest), "w"), indent=2)
    print(f"\nwrote {os.path.abspath(dest)}")


if __name__ == "__main__":
    main()
