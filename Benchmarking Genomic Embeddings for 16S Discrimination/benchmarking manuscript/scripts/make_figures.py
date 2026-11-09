"""Publication-ready figures for the benchmark manuscript.

Writes PNG (300 dpi) and PDF (vector) for each figure into ../figures/.
House style: no chart junk, no gridlines competing with data, colourblind-safe palette,
consistent fonts, every panel readable at single-column width.
"""
import json
import os
import sys
from collections import Counter, defaultdict

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------- house style
mpl.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.4,
})

# colourblind-safe (Okabe-Ito)
C = {
    "kmer":    "#0072B2",   # blue
    "model":   "#D55E00",   # vermillion
    "random":  "#999999",   # grey
    "align":   "#56B4E9",   # sky
    "exact":   "#CC79A7",   # pink
    "ceiling": "#000000",
    "dnabert2": "#D55E00", "dnabert_s": "#E69F00",
    "hyenadna": "#009E73", "ntv2": "#CC79A7",
}
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
PRETTY = {"dnabert2": "DNABERT-2", "dnabert_s": "DNABERT-S",
          "hyenadna": "HyenaDNA", "ntv2": "NTv2"}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


def load(p):
    with open(p) as f:
        return json.load(f)


# ================================================================ Figure 1
def fig1_main_result():
    """Forest plot: per-group difference on the locked test, plus the three replications."""
    r = load(os.path.join(ROOT, "phase5_evaluate", "locked_test_results.json"))
    pb = r["per_block"]
    names = sorted(pb, key=lambda b: pb[b]["delta"])
    d = np.array([pb[b]["delta"] * 100 for b in names])

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(7.2, 4.6), gridspec_kw={"width_ratios": [2.6, 1]})

    y = np.arange(len(names))
    cols = [C["model"] if v > 0 else C["kmer"] for v in d]
    ax.barh(y, d, color=cols, height=0.68)
    ax.axvline(0, color="k", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ") for n in names])
    ax.set_xlabel("DNABERT-2 − k-mer  (percentage points)")
    ax.set_title("a   Per-group difference, locked test (n = 20)", loc="left",
                 fontweight="bold")
    ax.axvspan(-1, 1, color="k", alpha=0.055, zorder=0)
    ax.text(0.02, 0.055, "shaded band: ±1 pt, the threshold fixed in advance\n"
            "for a difference large enough to matter",
            transform=ax.transAxes, fontsize=7, ha="left", va="bottom", color="#444")
    ax.set_ylim(-1.3, len(names) - 0.3)

    # replication panel
    reps = [
        ("Development\n(25 groups)", 0.14, None, None),
        ("Locked test\n(20 groups)", r["primary"]["mean"] * 100,
         r["primary"]["ci95"][0] * 100, r["primary"]["ci95"][1] * 100),
    ]
    try:
        a = load(os.path.join(ROOT, "phase7_extensions", "extension_A_results.json"))
        reps.append(("Hard panels\n(13 groups)", a["mean"] * 100,
                     a["ci"][0] * 100, a["ci"][1] * 100))
    except FileNotFoundError:
        pass

    yy = np.arange(len(reps))[::-1]
    for i, (lab, m, lo, hi) in zip(yy, reps):
        if lo is not None:
            ax2.plot([lo, hi], [i, i], color="k", lw=1.5, solid_capstyle="round")
            ax2.plot([lo, lo], [i - .09, i + .09], color="k", lw=1.2)
            ax2.plot([hi, hi], [i - .09, i + .09], color="k", lw=1.2)
        ax2.plot(m, i, "o", ms=7, color=C["model"], zorder=5,
                 markeredgecolor="white", markeredgewidth=1)
    ax2.axvline(0, color="k", lw=0.9)
    ax2.axvspan(-1, 1, color="k", alpha=0.055, zorder=0)
    ax2.set_yticks(yy)
    ax2.set_yticklabels([x[0] for x in reps])
    ax2.set_ylim(-0.6, len(reps) - 0.4)
    ax2.set_xlabel("mean difference (pts)")
    ax2.set_title("b   Replication", loc="left", fontweight="bold")

    fig.tight_layout(w_pad=2.2)
    save(fig, "fig1_main_result")


# ================================================================ Figure 2
def fig2_readout_dominates():
    """Readout vs model vs layer — with metrics kept strictly separate.

    R0 is a PAIR-LEVEL AUC. R1 and R2 are ISOLATE-LEVEL BALANCED ACCURACY. An earlier
    version of this figure plotted all three on one accuracy axis and derived a ~30x
    ratio from the R0-R2 gap. That subtraction is invalid across different metrics.
    R0 is now shown on its own axis and explicitly labelled; every quantitative
    comparison uses R1 vs R2, which share an endpoint.
    """
    bl = load(os.path.join(ROOT, "phase2_pilot", "block_level_results.json"))
    sp = {k: v for k, v in bl.items() if v["rank"] == "species"}

    def mean_of(key):
        v = [b[key] for b in sp.values() if b.get(key) == b.get(key)]
        return float(np.mean(v)) if v else float("nan")

    def paired_ci(a_key, b_key, seed=1729):
        d = np.array([b[b_key] - b[a_key] for b in sp.values()
                      if b.get(a_key) == b.get(a_key) and b.get(b_key) == b.get(b_key)])
        rng = np.random.default_rng(seed)
        dr = d[rng.integers(0, len(d), size=(4000, len(d)))].mean(axis=1)
        return d.mean(), np.percentile(dr, 2.5), np.percentile(dr, 97.5), len(d)

    fig, axes = plt.subplots(1, 4, figsize=(9.6, 3.2),
                             gridspec_kw={"width_ratios": [0.8, 1, 1, 0.95]})

    # (a) R0 — DIFFERENT METRIC, own axis
    ax = axes[0]
    vals = [mean_of(f"{m}:R0") for m in MODELS]
    ax.bar(range(4), vals, color=[C[m] for m in MODELS], width=0.62)
    ax.set_xticks(range(4))
    ax.set_xticklabels([PRETTY[m] for m in MODELS], rotation=38, ha="right", fontsize=7)
    ax.set_ylim(0.5, 0.85)
    ax.set_ylabel("pair AUC")
    ax.set_title("a   R0 geometry\n(pair AUC — different\nmetric, not comparable)",
                 loc="left", fontweight="bold", fontsize=8.6)

    # (b) R1 -> R2, SAME endpoint
    ax = axes[1]
    for m in MODELS:
        ax.plot([0, 1], [mean_of(f"{m}:R1"), mean_of(f"{m}:R2")], "o-",
                color=C[m], label=PRETTY[m], ms=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["R1\nnearest ref.", "R2\nlinear probe"])
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylim(0.90, 0.96)
    ax.set_ylabel("balanced accuracy")
    mu, lo, hi, n = paired_ci("dnabert2:R1", "dnabert2:R2")
    ax.set_title(f"b   Readout, same endpoint\n+{mu*100:.1f} pts "
                 f"[{lo*100:+.1f}, {hi*100:+.1f}]", loc="left", fontweight="bold",
                 fontsize=8.6)
    ax.legend(loc="lower right", fontsize=7)

    # (c) model choice at R2
    ax = axes[2]
    vals = [mean_of(f"{m}:R2") for m in MODELS]
    ax.bar(range(4), vals, color=[C[m] for m in MODELS], width=0.62)
    ax.set_xticks(range(4))
    ax.set_xticklabels([PRETTY[m] for m in MODELS], rotation=38, ha="right", fontsize=7)
    ax.set_ylim(0.90, 0.96)
    ax.set_ylabel("balanced accuracy")
    ax.set_title(f"c   Model choice at R2\nspread {(max(vals)-min(vals))*100:.2f} pts",
                 loc="left", fontweight="bold", fontsize=8.6)

    # (d) same axis, balanced-accuracy points only
    ax = axes[3]
    readout = np.mean([mean_of(f"{m}:R2") - mean_of(f"{m}:R1") for m in MODELS]) * 100
    model = (max(vals) - min(vals)) * 100
    layer = np.nan
    lv = []
    for m in ["hyenadna", "ntv2"]:
        pth = os.path.join(ROOT, "phase7_extensions", f"extension_C_{m}.json")
        if os.path.exists(pth):
            e = load(pth); pb = e["per_block"]
            vv = [np.mean([pb[b][f"layer_{d}"] for b in pb]) for d in e["depths"]]
            lv.append((max(vv) - min(vv)) * 100)
    if lv:
        layer = float(np.mean(lv))

    labs = ["readout\n(R1 to R2)", "model\nchoice", "layer depth\n(hard panels)"]
    bars = [readout, model, layer]
    cols = ["#444444", C["kmer"], "#BBBBBB"]
    yy = np.arange(3)[::-1]
    ax.barh(yy, bars, color=cols, height=0.58)
    for y_, v in zip(yy, bars):
        if v == v:
            ax.text(v + 0.09, y_, f"{v:.2f}", va="center", fontsize=8, fontweight="bold")
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=7.6)
    ax.set_xlabel("spread (balanced-accuracy pts)")
    ax.set_xlim(0, max(v for v in bars if v == v) * 1.38)
    ax.set_title(f"d   Same metric, same axis\nreadout ~{readout/model:.0f}x model",
                 loc="left", fontweight="bold", fontsize=8.6)
    ax.text(0.02, -0.52, "grey: measured on the hard panels,\nnot the development set",
            transform=ax.transAxes, fontsize=6.6, color="#555", va="top")

    fig.tight_layout(w_pad=1.5)
    save(fig, "fig2_readout_dominates")


# ================================================================ Figure 3
def fig3_methods_and_random():
    """All methods incl. the random-projection control, with the ceiling marked."""
    r = load(os.path.join(ROOT, "phase5_evaluate", "locked_test_results.json"))
    order = [("C_random", "Random projection\nof k-mer counts", C["random"]),
             ("C_kmer", "k-mer counting", C["kmer"]),
             ("dnabert_s:R2", "DNABERT-S", C["dnabert_s"]),
             ("dnabert2:R2", "DNABERT-2", C["dnabert2"]),
             ("ntv2:R2", "NTv2", C["ntv2"]),
             ("hyenadna:R2", "HyenaDNA", C["hyenadna"]),
             ("C_align", "Alignment 1-NN", C["align"]),
             ("C_exact", "Exact match", C["exact"])]
    scores = load(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json"))
    blocks = [v for v in scores.values()
              if v.get("_status") == "ok" and v.get("_block_eligible")]
    means, labels, colors = [], [], []
    for key, lab, col in order:
        vals = [b[key] for b in blocks if key in b and b[key] == b[key]]
        if vals:
            means.append(np.mean(vals) * 100); labels.append(lab); colors.append(col)
    ceil = np.mean([v["ceiling"] for v in r["per_block"].values()]) * 100

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    y = np.arange(len(means))[::-1]
    ax.barh(y, means, color=colors, height=0.68)
    ax.axvline(ceil, color="k", ls="--", lw=1.1)
    ax.text(ceil - 0.3, y[0] + 0.75, f"ceiling {ceil:.1f}%", fontsize=7.5,
            ha="right", va="center")
    for yi, v in zip(y, means):
        ax.text(v - 0.4, yi, f"{v:.1f}", va="center", ha="right",
                color="white", fontsize=7.5, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlim(88, 99)
    ax.set_xlabel("balanced accuracy (%), locked test")
    ax.set_title("A random projection of k-mer counts matches\nevery learned representation",
                 loc="left", fontweight="bold", fontsize=9.5)

    # label-permutation control: the test that shows the task is actually learned
    pp = os.path.join(RES, "label_permutation_control.json")
    if os.path.exists(pp):
        s = load(pp)["summary"]
        perm = np.mean([v["permuted"] for v in s.values()]) * 100
        ax.axvline(perm, color="#B00020", ls="-", lw=1.4, zorder=4)
        ax.text(perm + 0.25, y[-1] - 0.2,
                f"permuted labels\n{perm:.1f}% (chance)",
                fontsize=7, color="#B00020", va="center")
        ax.set_xlim(min(48, perm - 2), 99)
    fig.tight_layout()
    save(fig, "fig3_methods_and_random")


# ================================================================ Figure 4
def fig4_difficulty():
    """Difficulty: saturation on the locked test, and the hard-panel test."""
    scores = load(os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json"))
    lt = {k: v for k, v in scores.items()
          if v.get("_status") == "ok" and v.get("_block_eligible")}
    err = np.array([(1 - v["C_kmer"]) * 100 for v in lt.values()])
    dlt = np.array([(v["dnabert2:R2"] - v["C_kmer"]) * 100 for v in lt.values()])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2))

    ax = axes[0]
    ax.scatter(err, dlt, s=34, color=C["model"], edgecolor="white", linewidth=0.7, zorder=3)
    ax.axhline(0, color="k", lw=0.9)
    ax.axvline(1, color=C["random"], ls=":", lw=1.1)
    ax.text(1.25, ax.get_ylim()[1] * 0.82, "10 of 20 groups\nsaturated", fontsize=7.3,
            color="#444", va="top")
    ax.set_xlabel("k-mer error (pts) — how much room there is")
    ax.set_ylabel("DNABERT-2 − k-mer (pts)")
    ax.set_title("a   Locked test: half the groups\n     leave no room to differ",
                 loc="left", fontweight="bold")

    ax = axes[1]
    try:
        a = load(os.path.join(ROOT, "phase7_extensions", "extension_A_results.json"))
        pp = a["per_panel"]
        names = sorted(pp, key=lambda b: pp[b]["delta"])
        d = np.array([pp[b]["delta"] * 100 for b in names])
        yy = np.arange(len(names))
        ax.barh(yy, d, color=[C["model"] if v > 0 else C["kmer"] for v in d], height=0.66)
        ax.axvline(0, color="k", lw=0.9)
        ax.set_yticks(yy)
        short = {"mycobacterium_abscessus_subsp": "M. abscessus subsp.",
                 "clostridium_botulinum_groups": "C. botulinum groups",
                 "burkholderia_cepacia_complex": "B. cepacia complex",
                 "pseudomonas_fluorescens_group": "P. fluorescens group",
                 "lactobacillus_delbrueckii_subsp": "L. delbrueckii subsp.",
                 "streptococcus_dysgalactiae_subsp": "S. dysgalactiae subsp.",
                 "lactococcus_lactis_subsp": "L. lactis subsp.",
                 "bacillus_subtilis_subsp": "B. subtilis subsp.",
                 "lactobacillus_casei_group": "L. casei group",
                 "xanthomonas_campestris_pv": "X. campestris group",
                 "pectobacterium_dickeya": "Pectobacterium / Dickeya",
                 "nostoc_anabaena": "Nostoc / Anabaena",
                 "bradyrhizobium": "Bradyrhizobium"}
        ax.set_yticklabels([short.get(n, n.replace("_", " ")) for n in names], fontsize=7.2)
        ax.set_xlabel("DNABERT-2 − k-mer (pts)")
        ax.set_title(f"b   Hard groups only (n = {a['n']}):\n     mean "
                     f"{a['mean']*100:+.1f} pts, CI "
                     f"[{a['ci'][0]*100:+.1f}, {a['ci'][1]*100:+.1f}]",
                     loc="left", fontweight="bold")
    except FileNotFoundError:
        pass
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig4_difficulty")


# ================================================================ Figure 5
def fig5_robustness():
    """Degradation: both arms fall together; truncation hurts most."""
    r = load(os.path.join(ROOT, "phase7_extensions", "robustness_results.json"))
    conds = ["clean", "err_0.1%", "err_1%", "copy_var", "trim_10", "combined"]
    nice = {"clean": "clean", "err_0.1%": "0.1% error", "err_1%": "1% error",
            "copy_var": "copy variation", "trim_10": "10 bp trimmed\neach end",
            "combined": "1% error +\ntrimming"}
    km = [r[c]["_summary"]["kmer_mean"] * 100 for c in conds]
    d2 = [r[c]["_summary"]["dnabert2_mean"] * 100 for c in conds]

    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    x = np.arange(len(conds))
    w = 0.38
    ax.bar(x - w / 2, km, w, label="k-mer counting", color=C["kmer"])
    ax.bar(x + w / 2, d2, w, label="DNABERT-2", color=C["model"])
    ax.set_xticks(x); ax.set_xticklabels([nice[c] for c in conds], fontsize=7.4)
    ax.set_ylim(86, 98)
    ax.set_ylabel("balanced accuracy (%)")
    ax.legend(loc="upper right")
    ax.set_title("Both methods degrade together; losing read length costs most",
                 loc="left", fontweight="bold", fontsize=9.5)
    fig.tight_layout()
    save(fig, "fig5_robustness")


# ================================================================ Figure 6
def fig6_region_ranking():
    """Task T2: can a method pick the right region in advance?"""
    p = os.path.join(RES, "region_ranking_summary.json")
    if not os.path.exists(p):
        print("  [skip] fig6 — region ranking not yet run")
        return
    s = load(p)
    sets = [k for k in ("locked_test", "hard_panels") if k in s]
    fig, axes = plt.subplots(1, len(sets), figsize=(3.6 * len(sets), 3.3), squeeze=False)
    for ax, setname in zip(axes[0], sets):
        d = s[setname]
        arms = [a for a in ["C_kmer"] + [f"{m}:R2" for m in MODELS] if a in d]
        reg = [d[a]["regret"] * 100 for a in arms]
        labs = ["k-mer" if a == "C_kmer" else PRETTY[a.split(":")[0]] for a in arms]
        cols = [C["kmer"] if a == "C_kmer" else C[a.split(":")[0]] for a in arms]
        y = np.arange(len(arms))[::-1]
        ax.barh(y, reg, color=cols, height=0.66)
        ax.set_yticks(y); ax.set_yticklabels(labs)
        ax.set_xlabel("regret (pts) — lower is better")
        ax.set_title(f"{setname.replace('_', ' ')}  (n = {d[arms[0]]['n']})",
                     loc="left", fontweight="bold")
    fig.suptitle("Region ranking: accuracy lost by choosing the wrong region",
                 x=0.01, ha="left", fontweight="bold", fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93], w_pad=2.0)
    save(fig, "fig6_region_ranking")


# ================================================================ Figure 7
def fig7_ceiling():
    """Where 16S difficulty lives: achieved accuracy against the information ceiling."""
    pts = []
    for setname, pth, sumpth in [
        ("locked test", os.path.join(ROOT, "phase4_locked_test", "locked_test_scores.json"),
         None),
        ("hard panels", os.path.join(ROOT, "phase7_extensions", "hard_scores.json"), None)]:
        if not os.path.exists(pth):
            continue
        sc = load(pth)
        for b, v in sc.items():
            if v.get("_status") != "ok" or not v.get("_block_eligible"):
                continue
            ch = Counter(v["_chosen_regions"]).most_common(1)[0][0]
            ceil = (v.get("_ceilings") or {}).get(ch)
            if ceil is None or ceil != ceil:
                continue
            pts.append((setname, b, ceil * 100, v["C_kmer"] * 100,
                        v["dnabert2:R2"] * 100))

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    for setname, mk in [("locked test", "o"), ("hard panels", "s")]:
        sub = [p for p in pts if p[0] == setname]
        if not sub:
            continue
        ax.scatter([p[2] for p in sub], [p[3] for p in sub], marker=mk, s=38,
                   color=C["kmer"], edgecolor="white", linewidth=0.7,
                   label=f"k-mer, {setname}", zorder=3, alpha=0.9)
    lim = [50, 102]
    ax.plot(lim, lim, "k--", lw=1, zorder=1)
    ax.text(99, 96, "ceiling\n(no method can\nexceed this)", fontsize=7, ha="right",
            color="#444")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("information ceiling (%)")
    ax.set_ylabel("achieved balanced accuracy (%)")
    ax.legend(loc="lower right")
    ax.set_title("Methods sit close to the information ceiling",
                 loc="left", fontweight="bold", fontsize=9.5)
    fig.tight_layout()
    save(fig, "fig7_ceiling")


if __name__ == "__main__":
    print("building figures…")
    for fn in [fig1_main_result, fig2_readout_dominates, fig3_methods_and_random,
               fig4_difficulty, fig5_robustness, fig6_region_ranking, fig7_ceiling]:
        try:
            fn()
        except Exception as e:
            print(f"  [FAIL] {fn.__name__}: {type(e).__name__}: {e}")
    print("done")
