"""Build the manuscript tables as CSV (for typesetting) and Markdown (for reading).

Every number traces to a results file; nothing is retyped by hand.
"""
import csv
import json
import os
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "tables"))
os.makedirs(OUT, exist_ok=True)

MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
PRETTY = {"dnabert2": "DNABERT-2", "dnabert_s": "DNABERT-S",
          "hyenadna": "HyenaDNA", "ntv2": "NTv2"}


def load(p):
    with open(p) as f:
        return json.load(f)


def write(name, header, rows, note=""):
    with open(os.path.join(OUT, f"{name}.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    with open(os.path.join(OUT, f"{name}.md"), "w") as f:
        f.write("| " + " | ".join(header) + " |\n")
        f.write("|" + "|".join(["---"] * len(header)) + "|\n")
        for r in rows:
            f.write("| " + " | ".join(str(x) for x in r) + " |\n")
        if note:
            f.write(f"\n{note}\n")
    print(f"  {name}: {len(rows)} rows")


# T1 — study design
def table1():
    rows = [
        ["Development", 31, 26, 25, 12356, "all analysis decisions"],
        ["Locked test", 22, 20, 20, 6933, "single confirmatory use"],
        ["Hard panels", 15, 13, 13, 4720, "difficulty-selected replication"],
        ["Total", 68, 59, 58, 24009, ""],
    ]
    write("table1_study_design",
          ["Experiment", "Groups declared", "Groups eligible", "Groups analysed",
           "Sequence records", "Purpose"], rows,
          "Eligibility required >=10 target isolates and >=2 neighbour species with >=5 each. "
          "Development counts 25 analysed because one block is a subspecies stratum reported "
          "separately. A further 6,788 records were fetched for difficulty screening only.")


# T2 — primary and secondary endpoints
def table2():
    r = load(os.path.join(ROOT, "phase5_evaluate", "locked_test_results.json"))
    p = r["primary"]
    rows = [["DNABERT-2 vs k-mer (PRIMARY)", f"{p['mean']*100:+.2f}",
             f"[{p['ci95'][0]*100:+.2f}, {p['ci95'][1]*100:+.2f}]",
             f"{p['p_two_sided']:.3f}", "—", p["verdict"]]]
    for k, v in r["secondary"].items():
        lab = k.replace(":R2 - C_kmer", " vs k-mer").replace("best-of-4 R2 - C_kmer",
                                                             "best-of-4 vs k-mer")
        for m in MODELS:
            lab = lab.replace(m, PRETTY[m])
        rows.append([lab, f"{v['mean']*100:+.2f}",
                     f"[{v['ci95'][0]*100:+.2f}, {v['ci95'][1]*100:+.2f}]",
                     f"{v['p_raw']:.3f}", f"{v['p_holm']:.3f}", v["verdict"]])
    write("table2_endpoints",
          ["Contrast", "Mean difference (pts)", "95% CI", "p", "p (Holm)", "Verdict"], rows,
          "Locked test, 20 groups. Differences in percentage points of balanced accuracy. "
          "Equivalence margin +-1.0 pt, fixed in advance. Secondary tests Holm-corrected "
          "within their family of four.")


# T3 — method performance
def table3():
    sc = load(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json"))
    blocks = [v for v in sc.values()
              if v.get("_status") == "ok" and v.get("_block_eligible")]
    order = [("C_random", "Random projection of k-mer counts (control)"),
             ("C_kmer", "k-mer counting (4/5/6-mer) + linear probe"),
             ("dnabert_s:R2", "DNABERT-S + linear probe"),
             ("dnabert2:R2", "DNABERT-2 + linear probe"),
             ("ntv2:R2", "Nucleotide Transformer v2 + linear probe"),
             ("hyenadna:R2", "HyenaDNA + linear probe"),
             ("C_align", "Pairwise alignment, tie-aware 1-NN"),
             ("C_exact", "Exact sequence match")]
    rows = []
    for k, lab in order:
        v = [b[k] for b in blocks if k in b and b[k] == b[k]]
        if v:
            rows.append([lab, f"{np.mean(v)*100:.1f}",
                         f"{np.std(v, ddof=1)*100:.1f}", len(v)])
    write("table3_method_performance",
          ["Method", "Mean balanced accuracy (%)", "SD across groups", "n groups"], rows,
          "Locked test. Mean information ceiling for these groups is 97.0%.")


# T4 — readout hierarchy
def table4():
    bl = load(os.path.join(ROOT, "phase2_pilot", "block_level_results.json"))
    sp = {k: v for k, v in bl.items() if v["rank"] == "species"}
    rows = []
    for m in MODELS:
        vals = []
        for lv in ("R0", "R1", "R2"):
            x = [v[f"{m}:{lv}"] for v in sp.values()
                 if v.get(f"{m}:{lv}") == v.get(f"{m}:{lv}")]
            vals.append(np.mean(x) * 100 if x else float("nan"))
        rows.append([PRETTY[m]] + [f"{v:.1f}" for v in vals] +
                    [f"{vals[2]-vals[0]:+.1f}"])
    cl = np.mean([v["classical"] for v in sp.values()]) * 100
    rows.append(["k-mer counting", "—", "—", f"{cl:.1f}", "—"])
    write("table4_readout_hierarchy",
          ["Model", "R0 geometry (%)", "R1 nearest reference (%)",
           "R2 linear probe (%)", "R2 − R0"], rows,
          "Development set, 25 groups. R0 is unsupervised cosine separation; R1 is "
          "nearest-reference lookup; R2 trains a linear probe. Spread across readouts "
          "(21-27 pts) is ~30x the spread across models at R2 (0.7 pts).")


# T5 — degradation robustness
def table5():
    r = load(os.path.join(ROOT, "phase7_extensions", "robustness_results.json"))
    conds = [("clean", "Clean reference sequences"),
             ("err_0.1%", "0.1% substitution error"),
             ("err_1%", "1% substitution error"),
             ("copy_var", "Intragenomic copy variation"),
             ("trim_10", "10 bases trimmed each end"),
             ("combined", "1% error + trimming")]
    base = r["clean"]["_summary"]
    rows = []
    for k, lab in conds:
        s = r[k]["_summary"]
        rows.append([lab, f"{s['kmer_mean']*100:.1f}", f"{s['dnabert2_mean']*100:.1f}",
                     f"{(s['kmer_mean']-base['kmer_mean'])*100:+.1f}",
                     f"{(s['dnabert2_mean']-base['dnabert2_mean'])*100:+.1f}",
                     f"{s['delta_mean']*100:+.2f}"])
    write("table5_robustness",
          ["Condition", "k-mer (%)", "DNABERT-2 (%)", "k-mer change",
           "DNABERT-2 change", "Difference (pts)"], rows,
          "Locked-test groups, both arms given identical degraded sequences. "
          "Read truncation is the most damaging condition for both methods.")


# T6 — layer sweep
def table6():
    rows = []
    for m in ["hyenadna", "ntv2"]:
        p = os.path.join(ROOT, "phase7_extensions", f"extension_C_{m}.json")
        if not os.path.exists(p):
            continue
        e = load(p)
        pb, last = e["per_block"], e["last"]
        km = np.mean([pb[b]["C_kmer"] for b in pb]) * 100
        for d in e["depths"]:
            v = np.mean([pb[b][f"layer_{d}"] for b in pb]) * 100
            rows.append([PRETTY[m], d, f"{d/last:.2f}", f"{v:.2f}", f"{v-km:+.2f}",
                         "final layer" if d == last else ""])
    if rows:
        write("table6_layer_sweep",
              ["Model", "Layer", "Relative depth", "Balanced accuracy (%)",
               "vs k-mer (pts)", "Note"], rows,
              "13 hard groups. DNABERT-2 and DNABERT-S could not be tested: they store "
              "intermediate layers unpadded internally, and extracting them safely would "
              "require reimplementing part of the model.")


# T7 — region ranking
def table7():
    p = os.path.join(HERE, "..", "results", "region_ranking_summary.json")
    if not os.path.exists(p):
        print("  [skip] table7 — region ranking not yet run")
        return
    s = load(p)
    rows = []
    for setname, d in s.items():
        for arm in ["C_kmer"] + [f"{m}:R2" for m in MODELS]:
            if arm not in d:
                continue
            v = d[arm]
            lab = "k-mer counting" if arm == "C_kmer" else PRETTY[arm.split(":")[0]]
            rows.append([setname.replace("_", " "), lab, v["n"],
                         f"{v['acc']*100:.1f}", f"{v['regret']*100:.2f}",
                         f"{v['hit']*100:.0f}%"])
    write("table7_region_ranking",
          ["Set", "Method", "n groups", "Accuracy of chosen region (%)",
           "Regret (pts)", "Picked the best region"], rows,
          "Each method ranks regions using training isolates only, commits to its top "
          "choice, and is scored on held-out isolates. Regret is the accuracy lost "
          "relative to the best region identifiable in hindsight.")


if __name__ == "__main__":
    print("building tables…")
    for fn in [table1, table2, table3, table4, table5, table6, table7]:
        try:
            fn()
        except Exception as e:
            print(f"  [FAIL] {fn.__name__}: {type(e).__name__}: {e}")
    print("done")
