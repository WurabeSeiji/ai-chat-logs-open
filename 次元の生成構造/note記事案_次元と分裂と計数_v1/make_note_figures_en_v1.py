#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""English versions of the five note-article figures (saved as *_en.png)."""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle

import make_note_figures_v1 as ja

BASE = os.path.dirname(os.path.abspath(__file__))
EXP = ja.EXP
FIG = ja.FIG

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

C_BLUE, C_GREEN, C_RED, C_ORANGE, C_GRAY = (
    ja.C_BLUE, ja.C_GREEN, ja.C_RED, ja.C_ORANGE, ja.C_GRAY)


def fig1_en():
    ns = np.arange(3, 13)
    m = ns * (ns - 1) // 2
    rank = 2 * np.minimum(ns, m // 2)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ns, m, "o-", color=C_BLUE, lw=2,
            label="number of relation waves (explodes as N squared)")
    ax.plot(ns, rank, "s-", color=C_GREEN, lw=2,
            label="rotational degrees of freedom (grows only linearly)")
    ax.axhline(3, color="k", ls="--", lw=1.5,
               label="uniquely readable spatial directions (fixed at three)")
    ax.set_xlabel("number of bodies N")
    ax.set_ylabel("count")
    ax.set_title("Relations explode, rotations grow slowly,\nreadable directions stop at three")
    ax.set_xticks(ns)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    ax.annotate("66 relation waves at N=12", xy=(12, 66), xytext=(8.6, 60),
                arrowprops=dict(arrowstyle="->", color=C_BLUE), color=C_BLUE)
    ax.annotate("still three directions", xy=(11, 3), xytext=(9.8, 12),
                arrowprops=dict(arrowstyle="->", color="k"))
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_relations_explode_directions_three_en.png"), dpi=160)
    plt.close(fig)


def fig2_en():
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Overflowing relation waves organize into a few rotation planes")

    rng = np.random.default_rng(7)
    t = np.linspace(0, 2 * np.pi, 120)
    for i in range(15):
        y0 = 0.9 + i * 0.55
        ph = rng.uniform(0, 2 * np.pi)
        ax.plot(0.5 + 2.2 * t / t[-1], y0 + 0.18 * np.sin(3 * t + ph),
                color=C_BLUE, lw=1.0, alpha=0.8)
    ax.text(1.6, 9.5, "15 relation waves (6 bodies)", ha="center", fontsize=11, color=C_BLUE)

    ax.add_patch(FancyArrowPatch((3.2, 5.0), (5.0, 5.0),
                                 arrowstyle="-|>", mutation_scale=25, color="k"))
    ax.text(4.1, 5.35, "organize", ha="center", fontsize=11)

    labels = ["Rotation plane 1", "Rotation plane 2", "Rotation plane 3"]
    for k, lab in enumerate(labels):
        cy = 7.6 - k * 2.3
        e = Ellipse((6.9, cy), 3.0, 1.5, facecolor="#e8f0fe",
                    edgecolor=C_GREEN, lw=2)
        ax.add_patch(e)
        th = np.linspace(0, 2 * np.pi, 60)
        ax.plot(6.9 + 1.05 * np.cos(th), cy + 0.5 * np.sin(th),
                color=C_GREEN, lw=1, ls=":")
        ax.add_patch(FancyArrowPatch((6.9 + 1.05, cy), (6.9 + 1.02, cy + 0.09),
                                     arrowstyle="-|>", mutation_scale=14, color=C_GREEN))
        ax.text(6.9, cy, lab, ha="center", va="center", fontsize=11)
    ax.text(6.9, 1.0, "Each plane has the same form as the two-wave relation (AB)\n= one particle-like unit",
            ha="center", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_planes_as_registers_en.png"), dpi=160)
    plt.close(fig)


def fig3_en():
    path = os.path.join(EXP, "spontaneous_splitting_result_v1", "dormant_fraction_curves_v1.csv")
    taus3, f3 = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            taus3.append(float(row["tau"]))
            f3.append(float(row["f_median_delta_1e-05"]))
    taus3 = np.array(taus3)
    f3 = np.array(f3)

    curves = ja._load_largeN_curves()

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    ax.semilogy(taus3, f3, color=C_BLUE, lw=1.6,
                label="3 bodies (3 relation waves, seed $10^{-10}$ in energy ratio)")
    big_colors = {300: C_ORANGE, 1000: C_RED}
    for n_big, taus_b, f_b in curves:
        m_big = n_big * (n_big - 1) // 2
        ax.semilogy(taus_b, f_b, color=big_colors.get(n_big, C_GREEN), lw=1.8,
                    label=f"{n_big} bodies ({m_big:,} relation waves, seed $10^{{-30}}$)")
    ax.set_xlabel("time (steps)")
    ax.set_ylabel("fraction transferred to new waves (log scale)")
    ax.set_title("A nearly zero component rises geometrically (measured data)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=9.5)
    ax.set_ylim(1e-31, 2.0)
    if curves:
        n_big, taus_b, f_b = curves[-1]
        i_mid = int(np.argmin(np.abs(f_b - 1e-15)))
        ax.annotate("climbs 20+ orders of magnitude\nat a constant rate",
                    xy=(taus_b[i_mid], f_b[i_mid]), xytext=(taus_b[i_mid] + 1100, 1e-19),
                    arrowprops=dict(arrowstyle="->", color="k"), fontsize=10.5)
        ax.annotate("starts at $10^{-30}$ — a dormant\ncomponent far below measurement",
                    xy=(taus_b[0], f_b[0]), xytext=(400, 1e-24),
                    arrowprops=dict(arrowstyle="->", color="k"), fontsize=10.5)
    ax.annotate("stabilizes at a constant fraction\n= the new wave takes hold",
                xy=(taus3[-80], f3[-80]), xytext=(700, 4e-5),
                arrowprops=dict(arrowstyle="->", color="k"), fontsize=10.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_dormant_growth_en.png"), dpi=160)
    plt.close(fig)


def fig4_en():
    rows = ["zero-closure stratum\n(conserved value = 0)", "stratum with large\nconserved value"]
    ns = [3, 4, 5, 6]
    classes = [
        ["bounded\nmixing", "bounded\nmixing", "expansion", "expansion"],
        ["recurrence", "expansion", "expansion", "expansion"],
    ]
    color = {"expansion": "#f4b6b6", "bounded\nmixing": "#d9d9d9", "recurrence": "#b6c8f4"}

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.set_xlim(0, len(ns))
    ax.set_ylim(0, len(rows))
    for i, row in enumerate(classes):
        for j, cls in enumerate(row):
            y = len(rows) - 1 - i
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor=color[cls],
                                   edgecolor="white", lw=3))
            ax.text(j + 0.5, y + 0.5, cls, ha="center", va="center", fontsize=11)
    ax.set_xticks([j + 0.5 for j in range(len(ns))])
    ax.set_xticklabels([f"{n} bodies" for n in ns], fontsize=11)
    ax.set_yticks([len(rows) - 1 - i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels(rows, fontsize=10)
    ax.set_title("Three fates of splitting — decided by the number of bodies\nand the conserved quantity")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.text(0.5, 0.02,
             "expansion = a new wave takes hold and spreads / bounded mixing = spreads but never takes hold\n"
             "recurrence = keeps returning deep toward the original wave",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(os.path.join(FIG, "fig4_three_fates_en.png"), dpi=160)
    plt.close(fig)


def fig5_en():
    import json
    with open(os.path.join(EXP, "counting_ceiling_result_v1", "summary_v1.json")) as fh:
        d = json.load(fh)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.axis("off")
    ax1.set_title("Existence layer: no ceiling on the number of waves")
    levels = [(8.6, "1 term", 1), (6.8, "2 terms", 2), (5.0, "4 terms", 4), (3.2, "8 terms", 8)]
    t = np.linspace(0, 2 * np.pi, 100)
    for y0, lab, k in levels:
        ax1.text(0.4, y0, lab, fontsize=10.5, va="center")
        width = 6.6 / k
        for i in range(k):
            x0 = 2.6 + i * width
            ax1.plot(x0 + width * 0.85 * t / t[-1],
                     y0 + 0.35 * np.sin(k * t / 2 + i), color=C_BLUE, lw=1.4)
    for ya, yb in [(8.2, 7.4), (6.4, 5.6), (4.6, 3.8)]:
        ax1.add_patch(FancyArrowPatch((1.2, ya), (1.2, yb),
                                      arrowstyle="-|>", mutation_scale=14, color=C_GRAY))
    ax1.text(5.5, 1.9, "... can be rewritten into any number of terms", fontsize=10.5, ha="center")
    ax1.text(5.5, 0.7, "all legitimate readings of one and the same closure",
             fontsize=10, ha="center", color=C_GRAY)

    u3 = d["u3"]
    ns = sorted(int(k) for k in u3)
    filled = {n: {fam: u3[str(n)][fam] / u3[str(n)]["ceiling"]
                  for fam in ("thermal", "frozen", "parent")} for n in ns}
    fams = [("thermal", "thermalized (mixed waves)", C_BLUE),
            ("frozen", "frozen state", C_ORANGE),
            ("parent", "single wave (parent)", C_GREEN)]
    x = np.arange(len(ns))
    w = 0.26
    for k, (fam, lab, col) in enumerate(fams):
        vals = [filled[n].get(fam, 0.0) for n in ns]
        ax2.bar(x + (k - 1) * w, vals, w, color=col, label=lab)
    ax2.axhline(1.0, color=C_RED, ls="--", lw=1.8)
    ax2.text(len(ns) - 0.55, 1.03, "ceiling allowed by the convention", color=C_RED,
             fontsize=10, ha="right")
    ax2.set_ylim(0, 1.15)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{n} bodies" for n in ns])
    ax2.set_ylabel("readout count / ceiling")
    ax2.set_title("Readout layer: stops at 10-40% of the ceiling", fontsize=11.5)
    ax2.legend(loc="center left", fontsize=9)
    ax2.grid(axis="y", alpha=0.3)

    fig.suptitle("Existence opens without limit; only the readout saturates", fontsize=14, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(FIG, "fig5_existence_vs_readout_en.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig1_en()
    fig2_en()
    fig3_en()
    fig4_en()
    fig5_en()
    print("English figures written to", FIG)
