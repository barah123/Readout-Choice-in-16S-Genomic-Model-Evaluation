"""AQ06 ambiguous-base audit: how much of the empirical ambiguity ceiling is an artefact of N?

N is the only non-ACGT character present in the retrieved records. A slice carrying an N is a
distinct string, so it forms its own haplotype; splitting a haplotype can only DECREASE
sum(min(a, b)) and INCREASE sum(max(a/2n_pos, b/2n_neg)), i.e. it can only make the ceiling
look more generous. This script measures that, per confirmatory block, by recomputing the
ceiling at the region each block actually used with the N-bearing isolate units removed.
"""
import json
import os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
TEST = os.path.join(ROOT, "phase4_locked_test")


def ceilings(counts, n_pos, n_neg):
    total = n_pos + n_neg
    if not total or not n_pos or not n_neg:
        return float("nan"), float("nan")
    err = sum(min(a, b) for a, b in counts.values())
    bal = sum(max(a / (2 * n_pos), b / (2 * n_neg)) for a, b in counts.values())
    return 1 - err / total, bal


def cell_ceiling(units):
    """units: list of (slice, label)"""
    counts = defaultdict(lambda: [0, 0])
    n_pos = n_neg = 0
    for s, lab in units:
        counts[s][0 if lab else 1] += 1   # index 0 = target, matching build_blocks_full
        n_pos += bool(lab)
        n_neg += not lab
    return ceilings(counts, n_pos, n_neg)


def main():
    rows = json.load(open(os.path.join(TEST, "test_blocks.json")))
    sc = json.load(open(os.path.join(TEST, "locked_test_scores.json")))
    elig = sorted(b for b, v in sc.items()
                  if v.get("_status") == "ok" and v.get("_block_eligible"))

    by_cell = defaultdict(list)
    for r in rows:
        by_cell[(r["block_id"], r["region"])].append((r["slice"], r["is_target"]))

    out, deltas = {}, []
    for b in elig:
        region = Counter(sc[b]["_chosen_regions"]).most_common(1)[0][0]
        units = by_cell[(b, region)]
        keep = [u for u in units if "N" not in u[0]]
        a_all, b_all = cell_ceiling(units)
        a_kp, b_kp = cell_ceiling(keep)
        out[b] = {"region": region, "n_units": len(units),
                  "n_units_with_N": len(units) - len(keep),
                  "ceiling_bal_as_computed": b_all, "ceiling_bal_N_units_removed": b_kp,
                  "delta_pp": (b_all - b_kp) * 100}
        deltas.append(out[b]["delta_pp"])
        print(f"{b:26} {region}  units={len(units):4} withN={len(units)-len(keep):3}  "
              f"ceiling {b_all*100:6.2f} -> {b_kp*100:6.2f}  delta {(b_all-b_kp)*100:+.3f} pp")

    summary = {"n_blocks": len(elig),
               "max_abs_delta_pp": max(abs(d) for d in deltas),
               "mean_delta_pp": sum(deltas) / len(deltas),
               "blocks_with_any_N_unit": sum(1 for v in out.values() if v["n_units_with_N"])}
    print("\n" + json.dumps(summary, indent=2))
    dest = os.path.join(HERE, "..", "results", "ambiguous_base_audit.json")
    json.dump({"per_block": out, "summary": summary}, open(dest, "w"), indent=2)
    print(f"wrote {os.path.abspath(dest)}")


if __name__ == "__main__":
    main()
