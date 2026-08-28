"""記事図1（対数目盛のインフレーション図）の修正版再生成：原本 make_note_figs_en.py の Fig 1 節をそのまま抜き出し。usage: OUT ZIPDIR"""
# -*- coding: utf-8 -*-
"""English-labelled versions of the note figures (same data as make_note_figs.py)."""
import os, sys, csv, zipfile, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]; ZIPDIR = sys.argv[2]
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.titlesize": 15, "axes.titleweight": "bold",
                     "axes.labelsize": 13, "mathtext.fontset": "dejavusans", "axes.spines.top": False, "axes.spines.right": False})
RED, BLUE, GREEN, GRAY, NAVY, ORANGE = "#d7263d", "#1f5fd8", "#2e8b57", "#777777", "#12245e", "#e8871e"

def read(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try: out[k] = np.array([float(v) for v in vals])
        except ValueError: out[k] = vals
    return out

# ---- Fig 1: inflation, log scale, N=5 and N=16
z = zipfile.ZipFile(os.path.join(ZIPDIR, "complex_simplex_decompactification_N5_N16_20260826.zip"))
def load(n):
    rows = list(csv.DictReader(io.StringIO(z.read(n).decode())))
    return {k: np.array([float(r[k]) for r in rows]) for k in ("step", "H_perp", "H_parallel", "H_total")}
d5, d16 = load("results/N5_geometry_summary.csv"), load("results/N16_geometry_summary.csv")
fig, ax = plt.subplots(figsize=(10, 6.2), dpi=170)
ax.semilogy(d16["step"], d16["H_perp"], color=RED, lw=2.4, label=r"N = 16: amplitude$^2$ in the new directions, $H_\perp$")
ax.semilogy(d5["step"], d5["H_perp"], color=ORANGE, lw=2.0, label=r"N = 5: amplitude$^2$ in the new directions, $H_\perp$")
ax.semilogy(d16["step"], d16["H_parallel"], color=BLUE, lw=1.8, alpha=0.9, label=r"N = 16: amplitude$^2$ in the original plane, $H_\parallel$")
ax.set_ylim(1e-34, 5); ax.set_xlim(0, 800)
h0, hf = d16["H_perp"][0], d16["H_perp"][-1]
ax.axhline(h0, color=GRAY, lw=0.8, ls=":"); ax.axhline(hf, color=GRAY, lw=0.8, ls=":")
ax.annotate("", xy=(520, hf), xytext=(520, h0), arrowprops=dict(arrowstyle="<->", color=NAVY, lw=2.2))
ax.text(540, 1e-15, "about 31 orders\nof magnitude\n($10^{31}$ times)", color=NAVY, fontsize=15, weight="bold", va="center")
ax.text(200, hf * 0.12, f"N = 16 final {hf:.2f}  (N = 5: 0.4)", fontsize=11, color=RED, va="top")
ax.text(795, h0 * 0.25, "start $10^{-32}$ (N = 16) / $10^{-31}$ (N = 5)", ha="right", fontsize=11, color=RED)
ax.text(60, 1e-8, "straight line on a log scale\n= exponential growth at a constant rate", fontsize=10, color=RED, rotation=66, ha="center", va="center")
ax.set_xlabel("step"); ax.set_ylabel("squared amplitude (total = 1, log scale)")
ax.set_title("Inflation: the amplitude in the new directions climbs 31 orders of magnitude")
ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(0.99, 0.08), fontsize=10); fig.tight_layout()
fig.savefig(os.path.join(OUT, "note_fig0_inflation_log_en.png")); plt.close(fig)

