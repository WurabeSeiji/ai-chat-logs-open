#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第4章導入用の図:
  (a) A system carries a ledger of phase pairs
  (b) Readouts split into phase differences and phase rates
  (c) The origin is inherited by genealogy, not chosen from outside
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch

C_TEXT = "#222222"
C_GRID = "#d0d0d0"
C_LEDGER = "#f4f7fb"
C_A = "#1f77b4"
C_B = "#ff7f0e"
C_RATE = "#2ca02c"
C_DIFF = "#9467bd"
C_WARN = "#d62728"


def add_arrow(ax, start, end, color="#333333", rad=0.0, lw=1.4, ms=12):
    arr = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=color,
        connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)
    return arr


def draw_dial(ax, center, label, angle, color=C_A, radius=0.17):
    cx, cy = center
    ax.add_patch(Circle((cx, cy), radius, fc="white", ec=color, lw=1.2))
    ax.plot([cx, cx + radius * 0.78 * np.cos(angle)],
            [cy, cy + radius * 0.78 * np.sin(angle)],
            color=color, lw=1.8)
    ax.text(cx, cy - radius - 0.08, label, ha="center", va="top",
            fontsize=9, color=C_TEXT)


def draw_ledger(ax, xy, title, color=C_A, phase_shift=0.0):
    x, y = xy
    box = Rectangle((x - 0.78, y - 0.52), 1.56, 1.04,
                    fc=C_LEDGER, ec=color, lw=1.3)
    ax.add_patch(box)
    ax.text(x, y + 0.40, title, ha="center", va="center",
            fontsize=10, color=color, weight="bold")
    labels = ["x", "y", "z", "t", "R", "Q"]
    positions = [(-0.50, 0.12), (0, 0.12), (0.50, 0.12),
                 (-0.50, -0.27), (0, -0.27), (0.50, -0.27)]
    for i, (lab, (dx, dy)) in enumerate(zip(labels, positions)):
        draw_dial(ax, (x + dx, y + dy),
                  lab, phase_shift + 0.7 * i + 0.35, color=color, radius=0.12)


fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2))

# ---------- (a) ledger ----------
ax = axes[0]
ax.set_title("(a) Internal ledger", fontsize=12)
draw_ledger(ax, (0, 0), "system", C_A, phase_shift=0.2)
ax.text(0, -0.92, "one conjugate pair per axis", ha="center", fontsize=10)
ax.text(0, -1.17, "absolute phases are not read", ha="center",
        fontsize=10, color=C_WARN)
ax.set_xlim(-1.25, 1.25)
ax.set_ylim(-1.35, 1.08)
ax.axis("off")

# ---------- (b) readout split ----------
ax = axes[1]
ax.set_title("(b) What can be read", fontsize=12)
draw_ledger(ax, (-0.78, 0.32), "A", C_A, phase_shift=0.1)
draw_ledger(ax, (0.88, 0.32), "B", C_B, phase_shift=0.75)
add_arrow(ax, (-0.12, 0.05), (0.22, 0.05), C_DIFF, lw=1.6, ms=13)
ax.text(0.05, -0.17, r"$\Delta\theta$", ha="center", fontsize=12, color=C_DIFF)
ax.text(0.05, -0.46, "position / time\n(relational values)",
        ha="center", fontsize=9, color=C_DIFF)
add_arrow(ax, (-0.78, -0.75), (-0.22, -0.75), C_RATE, lw=1.6, ms=13)
ax.text(-0.50, -0.62, r"$d\theta/ds$", ha="center", fontsize=12, color=C_RATE)
ax.text(-0.50, -1.05, "mass / energy /\nangular momentum / charge",
        ha="center", fontsize=9, color=C_RATE)
ax.set_xlim(-1.75, 1.85)
ax.set_ylim(-1.28, 1.15)
ax.axis("off")

# ---------- (c) genealogy ----------
ax = axes[2]
ax.set_title("(c) Origin by genealogy", fontsize=12)
nodes = {
    "root": (0, 0.86),
    "A": (-0.65, 0.20),
    "B": (0.65, 0.20),
    "A1": (-1.00, -0.54),
    "A2": (-0.30, -0.54),
    "B1": (0.65, -0.54),
}
for parent, child in [("root", "A"), ("root", "B"),
                      ("A", "A1"), ("A", "A2"), ("B", "B1")]:
    add_arrow(ax, nodes[parent], nodes[child], "#777777", lw=1.1, ms=10)
for name, (x, y) in nodes.items():
    color = C_A if name.startswith("A") else C_B if name.startswith("B") else C_RATE
    ax.add_patch(Circle((x, y), 0.16, fc="white", ec=color, lw=1.4))
    ax.text(x, y, name, ha="center", va="center", fontsize=9, color=C_TEXT)
ax.text(0, 1.14, "first: no difference", ha="center", fontsize=10, color=C_RATE)
ax.text(0, -0.93, "splits inherit the ledger;\ndifferences accumulate",
        ha="center", fontsize=10)
ax.text(0, -1.25, "origin is not selected from outside", ha="center",
        fontsize=9, color=C_WARN)
ax.set_xlim(-1.35, 1.35)
ax.set_ylim(-1.42, 1.30)
ax.axis("off")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch4_ledger.{ext}", dpi=200, bbox_inches="tight")
print("done: fig_text_ch4_ledger (png+svg)")
