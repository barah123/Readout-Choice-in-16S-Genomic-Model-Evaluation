"""Publication-ready figures for the benchmark manuscript.

Writes PNG (300 dpi) and PDF (vector) for each figure into ../figures/.
House style: no chart junk, no gridlines competing with data, colourblind-safe palette,
consistent fonts, every panel readable at single-column width.
"""
import csv
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
    "axes.labelweight": "bold",
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
    "exact":   "#4D4D4D",   # dark grey — was #CC79A7, which collided with ntv2 below
    "ceiling": "#000000",
    "dnabert2": "#D55E00", "dnabert_s": "#E69F00",
    "hyenadna": "#009E73", "ntv2": "#CC79A7",
}
MODELS = ["dnabert2", "dnabert_s", "hyenadna", "ntv2"]
PRETTY = {"dnabert2": "DNABERT-2", "dnabert_s": "DNABERT-S",
          "hyenadna": "HyenaDNA", "ntv2": "NTv2"}


PLOS_NUM = {"fig1_main_result": 1, "fig2_readout_dominates": 2,
            "fig3_methods_and_random": 3, "fig4_difficulty": 4,
            "fig5_robustness": 5, "fig6_region_ranking": 6, "fig7_ceiling": 7}


def save(fig, name):
    """Working copies as PNG/PDF, plus a PLOS submission TIFF.

    PLOS Computational Biology requires TIFF or EPS, 300-600 dpi, width between 2.63 and
    7.5 inches, height at most 8.75 inches, and LZW compression for TIFF. Files are named
    Fig1.tiff ... Fig7.tiff as the journal expects.
    """
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    n = PLOS_NUM.get(name)
    if n:
        tif_dir = os.path.join(FIG, "plos_tiff")
        os.makedirs(tif_dir, exist_ok=True)
        path = os.path.join(tif_dir, f"Fig{n}.tiff")
        # Render to RGBA then flatten onto white. Saving straight to TIFF leaves an alpha
        # channel, which PLOS production does not expect and which can render
        # unpredictably in print.
        import io
        from PIL import Image
        buf = io.BytesIO()
        fig.savefig(buf, dpi=300, format="png", facecolor="white")
        buf.seek(0)
        im = Image.open(buf)
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA")
            bg = Image.new("RGB", im.size, (255, 255, 255))
            bg.paste(im, mask=im.split()[-1])
            im = bg
        else:
            im = im.convert("RGB")
        im.save(path, format="TIFF", compression="tiff_lzw", dpi=(300, 300))
        buf.close()
        # Measure the SAVED raster, not fig.get_size_inches(): savefig runs with
        # bbox_inches="tight", so the file is larger than the declared figure size and a
        # figure declared inside the limit can land outside it.
        w, h = im.size[0] / 300.0, im.size[1] / 300.0
        ok_w = 2.63 <= w <= 7.5
        ok_h = h <= 8.75
        mb = os.path.getsize(path) / 1e6
        flag = "" if (ok_w and ok_h and mb < 10) else "   <-- CHECK"
        print(f"  Fig{n}.tiff  {w:.2f} x {h:.2f} in  {mb:.1f} MB"
              f"  width_ok={ok_w} height_ok={ok_h}{flag}")
    plt.close(fig)


def load(p):
    with open(p) as f:
        return json.load(f)


def _boot_ci(d, seed=1729, n_boot=4000):
    rng = np.random.default_rng(seed)
    draws = d[rng.integers(0, len(d), size=(n_boot, len(d)))].mean(axis=1)
    return float(d.mean()), float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def development_primary_contrast():
    """The MATCHED development contrast: DNABERT-2:R2 minus C_kmer, paired over blocks.

    An earlier version of Figure 1b hardcoded 0.14, which is the stored `model - classical`
    value in block_level_results.json. That quantity is a DOUBLE ORACLE: the best of four
    models per block (DNABERT-2 in 20/25, DNABERT-S 3, NTv2 1, HyenaDNA 1) AND the best of
    two classical comparators per block (C_kmer 17/25, C_align 8/25). Plotting it beside the
    locked-test and hard-panel rows put a twice-selected number on a shared axis with two
    pre-named contrasts, and it was the only positive row.

    This recomputes the same contrast the other two rows use, at each block's selected region.
    """
    bl = load(os.path.join(ROOT, "phase2_pilot", "block_level_results.json"))
    rr = load(os.path.join(ROOT, "phase2_pilot", "readout_results_full.json"))
    sp = {k: v for k, v in bl.items() if v["rank"] == "species"}
    d = []
    for b, v in sp.items():
        cell = rr.get(f"{b}|{v['region']}")
        if cell and "C_kmer" in cell and "dnabert2:R2_linear" in cell:
            d.append(cell["dnabert2:R2_linear"] - cell["C_kmer"])
    return np.array(d)


# ================================================================ Figure 1
def fig1_main_result():
    """Forest plot: per-group difference on the locked test, plus the three replications."""
    r = load(os.path.join(ROOT, "phase5_evaluate", "locked_test_results.json"))
    pb = r["per_block"]
    names = sorted(pb, key=lambda b: pb[b]["delta"])
    d = np.array([pb[b]["delta"] * 100 for b in names])

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(6.9, 4.6), gridspec_kw={"width_ratios": [2.6, 1]})

    y = np.arange(len(names))
    cols = [C["model"] if v > 0 else C["kmer"] for v in d]
    ax.barh(y, d, color=cols, height=0.68)
    ax.axvline(0, color="k", lw=0.9)
    # Seven of the 20 blocks differ by exactly zero and would otherwise render as no bar at
    # all, indistinguishable from a missing block. Mark them as points on the zero line.
    # (comamonadaceae at -0.002 and proteus_morganellaceae at +0.008 pp are near-zero but
    # not exact, and keep their hairline bars.)
    zero = [i for i, v in enumerate(d) if abs(v) < 1e-12]
    if zero:
        ax.plot([0] * len(zero), [y[i] for i in zero], "o", color="white",
                markeredgecolor="#333333", markeredgewidth=1.1, ms=4.5, zorder=5,
                ls="none", label="difference exactly zero")
        ax.legend(loc="upper left", fontsize=7.5, handletextpad=0.3, borderpad=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ") for n in names])
    ax.set_xlabel("DNABERT-2 − k-mer  (percentage points)")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=12)
    ax.axvspan(-1, 1, color="k", alpha=0.055, zorder=0)
    ax.set_ylim(-1.3, len(names) - 0.3)

    # replication panel — every row is the SAME contrast, DNABERT-2:R2 minus C_kmer
    dev = development_primary_contrast()
    dm, dlo, dhi = _boot_ci(dev)
    reps = [
        ("Development\n(25 groups,\nexploratory)", dm * 100, dlo * 100, dhi * 100),
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
    ax2.set_xlabel("mean difference  (percentage points)")
    ax2.set_title("b", loc="left", fontweight="bold", fontsize=12)

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

    fig, axgrid = plt.subplots(2, 2, figsize=(7.2, 6.2))
    axes = [axgrid[0, 0], axgrid[0, 1], axgrid[1, 0], axgrid[1, 1]]

    # (a) R0 — DIFFERENT METRIC, own axis
    ax = axes[0]
    vals = [mean_of(f"{m}:R0") for m in MODELS]
    ax.bar(range(4), vals, color=[C[m] for m in MODELS], width=0.62)
    ax.set_xticks(range(4))
    ax.set_xticklabels([PRETTY[m] for m in MODELS], rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0.5, 0.85)
    ax.set_ylabel("pair AUC")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=12)

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
    # Not drawn on the plot (the annotation was removed), but NOT dead: this is the source
    # of the paired effect quoted in the Fig 2 legend — +2.83 pts, 95% CI [+1.37, +4.57].
    mu, lo, hi, n = paired_ci("dnabert2:R1", "dnabert2:R2")
    ax.set_title("b", loc="left", fontweight="bold", fontsize=12)
    ax.legend(loc="lower right", fontsize=8)

    # (c) model choice at R2
    ax = axes[2]
    vals = [mean_of(f"{m}:R2") for m in MODELS]
    ax.bar(range(4), vals, color=[C[m] for m in MODELS], width=0.62)
    ax.set_xticks(range(4))
    ax.set_xticklabels([PRETTY[m] for m in MODELS], rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0.90, 0.96)
    ax.set_ylabel("balanced accuracy")
    ax.set_title("c", loc="left", fontweight="bold", fontsize=12)

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
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=8)
    # Every bar is a difference between two balanced accuracies: R2 minus R1 for the readout
    # bar, best minus worst model for the model bar, best minus worst depth for the layer bar.
    ax.set_xlabel("difference in balanced accuracy\n(percentage points)")
    ax.set_xlim(0, max(v for v in bars if v == v) * 1.38)
    ax.set_title("d", loc="left", fontweight="bold", fontsize=12)

    fig.tight_layout(w_pad=2.0, h_pad=2.4)
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
    # Fold-matched, matching Figure 7 and the aggregation of the accuracies it bounds. The
    # whole-block ceiling this previously averaged does not bound a mean of per-fold
    # accuracies; both happen to round to the same displayed value here.
    cpath = os.path.join(RES, "fold_matched_ceiling_audit.csv")
    if os.path.exists(cpath):
        with open(cpath) as f:
            crows = [x for x in csv.DictReader(f) if x["set"] == "locked test"]
        ceil = float(np.mean([float(x["ceiling_fold_matched"]) for x in crows])) * 100
        ceil_label = f"fold-matched ambiguity ceiling, {ceil:.1f}%"
    else:
        ceil = np.mean([v["ceiling"] for v in r["per_block"].values()]) * 100
        ceil_label = f"ambiguity ceiling, {ceil:.1f}%"

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    y = np.arange(len(means))[::-1]
    ax.barh(y, means, color=colors, height=0.68)
    # Both reference lines go in a legend rather than on-plot text: the permutation line
    # previously carried no label at all, and a legend labels both consistently.
    ax.axvline(ceil, color="k", ls="--", lw=1.1, zorder=4, label=ceil_label)
    for yi, v in zip(y, means):
        ax.text(v - 0.4, yi, f"{v:.1f}", va="center", ha="right",
                color="white", fontsize=7.5, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlim(88, 99)
    ax.set_xlabel("balanced accuracy (%), locked test")

    # label-permutation control: the test that shows the task is actually learned
    pp = os.path.join(RES, "label_permutation_control.json")
    if os.path.exists(pp):
        perm_sum = load(pp)["summary"]
        perm = np.mean([v["permuted"] for v in perm_sum.values()]) * 100
        ax.axvline(perm, color="#B00020", ls="-", lw=1.4, zorder=4,
                   label=f"label-permuted floor, {perm:.1f}%")
        ax.set_xlim(min(48, perm - 2), 99)
    # above the axes: every bar spans the full width, so an in-axes legend would sit on data
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1, fontsize=7.5,
              framealpha=0, handlelength=2.2, labelspacing=0.35)
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
    ax.set_xlabel("k-mer error (percentage points) — how much room there is")
    ax.set_ylabel("DNABERT-2 − k-mer  (percentage points)")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=12)

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
        ax.set_xlabel("DNABERT-2 − k-mer  (percentage points)")
        ax.set_title("b", loc="left", fontweight="bold", fontsize=12)
        # the pooled estimate, with the equivalence margin it fails to stay inside
        ax.axvspan(-1, 1, color="k", alpha=0.055, zorder=0)
        m, lo, hi = a["mean"] * 100, a["ci"][0] * 100, a["ci"][1] * 100
        yp = -1.15
        ax.plot([lo, hi], [yp, yp], color="k", lw=1.6, solid_capstyle="round", zorder=5)
        for e in (lo, hi):
            ax.plot([e, e], [yp - .22, yp + .22], color="k", lw=1.2, zorder=5)
        ax.plot(m, yp, "o", ms=6.5, color=C["model"], markeredgecolor="white",
                markeredgewidth=1, zorder=6)
        ax.set_ylim(yp - 0.9, len(names) - 0.35)
        ax.axhline(-0.5, color="#BBBBBB", lw=0.7)
        ax.set_yticks(list(yy) + [yp])
        ax.set_yticklabels([short.get(n, n.replace("_", " ")) for n in names]
                           + [f"pooled  {m:+.2f}\n[{lo:+.2f}, {hi:+.2f}]"], fontsize=7.2)
        for t in ax.get_yticklabels()[-1:]:
            t.set_fontweight("bold")
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

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(5.6, 5.2),
                                  gridspec_kw={"height_ratios": [1.45, 1]}, sharex=True)
    x = np.arange(len(conds))
    w = 0.38
    ax.bar(x - w / 2, km, w, label="k-mer counting", color=C["kmer"])
    ax.bar(x + w / 2, d2, w, label="DNABERT-2", color=C["model"])
    ax.set_ylim(86, 98)
    ax.set_ylabel("balanced accuracy (%)")
    ax.legend(loc="upper right")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=12)

    # (b) the paired block-level difference with its 95% bootstrap interval. The per-arm
    # bars above carry no uncertainty because the analysis is paired; the interval belongs
    # on the difference, which is what the frozen plan estimates.
    mu = np.array([r[c]["_summary"]["delta_mean"] * 100 for c in conds])
    lo = np.array([r[c]["_summary"]["ci"][0] * 100 for c in conds])
    hi = np.array([r[c]["_summary"]["ci"][1] * 100 for c in conds])
    ax2.axhspan(-1, 1, color="k", alpha=0.055, zorder=0)
    ax2.axhline(0, color="k", lw=0.9, zorder=1)
    for xi, m_, l_, h_ in zip(x, mu, lo, hi):
        ax2.plot([xi, xi], [l_, h_], color="k", lw=1.5, solid_capstyle="round", zorder=3)
        ax2.plot([xi - .09, xi + .09], [l_, l_], color="k", lw=1.1, zorder=3)
        ax2.plot([xi - .09, xi + .09], [h_, h_], color="k", lw=1.1, zorder=3)
    ax2.plot(x, mu, "o", ms=6, color=C["model"], markeredgecolor="white",
             markeredgewidth=1, zorder=5, ls="none")
    ax2.set_xticks(x); ax2.set_xticklabels([nice[c] for c in conds], fontsize=7.4)
    ax2.set_ylabel("DNABERT-2 − k-mer\n(percentage points)")
    ax2.set_title("b", loc="left", fontweight="bold", fontsize=12)

    fig.tight_layout(h_pad=1.4)
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
    nice = {"locked_test": "locked test", "hard_panels": "hard panels"}
    fig, axes = plt.subplots(2, len(sets), figsize=(3.6 * len(sets), 5.4), squeeze=False)
    letters = iter("abcd")
    # Row 1: accuracy actually achieved at the selected region. Row 2: regret against that
    # method's own best region. Accuracy was previously omitted, so the figure showed the
    # shortfall without the level it is a shortfall from.
    for row, key, xlab in [(0, "acc", "balanced accuracy (%) at the selected region"),
                           (1, "regret", "regret (percentage points) — lower is better")]:
        for col, setname in enumerate(sets):
            ax = axes[row][col]
            d = s[setname]
            arms = [a for a in ["C_kmer"] + [f"{m}:R2" for m in MODELS] if a in d]
            val = [d[a][key] * 100 for a in arms]
            labs = ["k-mer" if a == "C_kmer" else PRETTY[a.split(":")[0]] for a in arms]
            cols = [C["kmer"] if a == "C_kmer" else C[a.split(":")[0]] for a in arms]
            y = np.arange(len(arms))[::-1]
            if key == "acc":
                # Dot plot, not bars: these accuracies sit near 96% and differ by fractions
                # of a point, so any readable axis is truncated — and a truncated baseline
                # on a BAR chart misrepresents the ratios. Dots carry no such implication.
                for yi, v_, c_ in zip(y, val, cols):
                    ax.plot([min(val) - 10, v_], [yi, yi], color="#DDDDDD", lw=1.2, zorder=1)
                    ax.plot(v_, yi, "o", ms=7, color=c_, markeredgecolor="white",
                            markeredgewidth=1, zorder=3)
            else:
                ax.barh(y, val, color=cols, height=0.66)
            ax.set_yticks(y); ax.set_yticklabels(labs)
            ax.set_xlabel(xlab, fontsize=8)
            if key == "acc":
                lo = min(val)
                ax.set_xlim(lo - 1.2, min(100, max(val) + 0.8))
                n = d[arms[0]]["n"]
                ax.set_title(f"{next(letters)}", loc="left", fontweight="bold", fontsize=12)
                ax.text(1.0, 1.06, f"{nice[setname]} ({n} blocks)", transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=8.5, fontweight="bold")
            else:
                ax.set_title(f"{next(letters)}", loc="left", fontweight="bold", fontsize=12)
    fig.tight_layout(w_pad=2.0, h_pad=2.6)
    save(fig, "fig6_region_ranking")


# ================================================================ Figure 7
def fig7_ceiling():
    """Achieved accuracy against the ambiguity ceiling, on a matched aggregation.

    The reported accuracy is the mean over five outer folds of per-fold balanced accuracy.
    The ceiling published alongside it was evaluated once on the WHOLE block, so the two
    were not aggregated over the same thing and a fold mean can exceed a whole-population
    bound: 3 of 66 cells did, and the figure nevertheless annotated the diagonal "no method
    can exceed this". The ceiling plotted here is fold-matched — the same formula applied to
    each held-out fold at that fold's chosen region, then averaged over folds — which bounds
    all 66 cells and differs from the whole-block value by +0.001 pp on average
    (scripts/audit_fold_matched_ceiling.py).

    Both arms are now drawn. Earlier the DNABERT-2 value was computed and silently discarded,
    so the figure showed only the k-mer baseline.
    """
    path = os.path.join(RES, "fold_matched_ceiling_audit.csv")
    if not os.path.exists(path):
        print("  [skip] fig7 — run scripts/audit_fold_matched_ceiling.py first")
        return
    with open(path) as f:
        rows = list(csv.DictReader(f))

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    lim = [50, 102]
    ax.plot(lim, lim, "k--", lw=1, zorder=1)

    for setname, mk in [("locked test", "o"), ("hard panels", "s")]:
        sub = [r for r in rows if r["set"] == setname]
        if not sub:
            continue
        x = [float(r["ceiling_fold_matched"]) * 100 for r in sub]
        ax.scatter(x, [float(r["C_kmer"]) * 100 for r in sub], marker=mk, s=34,
                   color=C["kmer"], edgecolor="white", linewidth=0.6, zorder=3, alpha=0.9,
                   label=f"k-mer, {setname}")
        ax.scatter(x, [float(r["dnabert2_R2"]) * 100 for r in sub], marker=mk, s=34,
                   color=C["model"], edgecolor="white", linewidth=0.6, zorder=3, alpha=0.9,
                   label=f"DNABERT-2, {setname}")

    # The diagonal is labelled away from the data, which previously overprinted the
    # top-right cluster and was barely legible.
    # The region above the diagonal is empty by construction, so the label goes there and
    # cannot overprint either the data or the legend.
    ax.annotate("fold-matched ambiguity ceiling\n(attainable maximum)",
                xy=(66, 66), xytext=(53, 82), fontsize=7.5, color="#444", ha="left",
                arrowprops=dict(arrowstyle="-", lw=0.7, color="#888",
                                shrinkA=2, shrinkB=3))
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("empirical ambiguity ceiling (%), fold-matched")
    ax.set_ylabel("achieved balanced accuracy (%)")
    ax.legend(loc="lower right", fontsize=7.5)

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
