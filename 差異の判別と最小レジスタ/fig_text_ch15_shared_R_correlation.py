#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chapter 15 figure:
  Two externally separated partial systems can still share an R-ledger entry.
  The figure emphasizes shared readout constraints, not a signal traveling
  between A and B.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Arc


C_TEXT = "#222222"
C_MUTED = "#6d7480"
C_R = "#1f77b4"
C_A = "#d62728"
C_B = "#2ca02c"
C_SHARE = "#6f42c1"
C_FACE = "#f7f9fc"
C_NOTE = "#8a4b00"


def arrow(ax, start, end, color="#444444", lw=1.5, ms=12,
          rad=0.0, linestyle="-", alpha=1.0):
    arr = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=ms, lw=lw,
        color=color, linestyle=linestyle, alpha=alpha,
        connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)
    return arr


def label_box(ax, xy, w, h, title, lines, ec, fc=C_FACE):
    x, y = xy
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        fc=fc, ec=ec, lw=1.4
    )
    ax.add_patch(box)
    ax.text(x, y + h * 0.30, title, ha="center", va="center",
            fontsize=11.2, color=ec, weight="bold")
    for i, line in enumerate(lines):
        ax.text(x, y + h * 0.05 - i * 0.24, line,
                ha="center", va="center", fontsize=8.4, color=C_TEXT)
    return box


fig, ax = plt.subplots(figsize=(12.2, 6.2))
ax.set_xlim(0, 12.2)
ax.set_ylim(0, 6.2)
ax.set_aspect("equal")
ax.axis("off")

ax.text(6.1, 5.86, "Shared R ledger and separated readouts",
        ha="center", fontsize=15, color=C_TEXT, weight="bold")
ax.text(6.1, 5.55,
        "A and B are far apart on the external arc, but their records are constrained by the same shared R entry",
        ha="center", fontsize=9.4, color=C_MUTED)

# Shared R circle.
O = np.array([6.1, 3.05])
R = 1.95
ax.add_patch(Circle(O, R, fill=False, ec=C_R, lw=2.4))
ax.add_patch(Circle(O, 0.055, fc=C_R, ec=C_R))
ax.text(O[0], O[1] + 0.18, "shared R ledger", ha="center",
        fontsize=11.0, color=C_R, weight="bold")
ax.text(O[0], O[1] - 0.10, "same radius / conserved slot",
        ha="center", fontsize=8.6, color=C_TEXT)

# Two separated partial systems on the same R circle.
theta_A = np.deg2rad(142)
theta_B = np.deg2rad(28)
A = O + R * np.array([np.cos(theta_A), np.sin(theta_A)])
B = O + R * np.array([np.cos(theta_B), np.sin(theta_B)])

ax.plot([O[0], A[0]], [O[1], A[1]], color=C_R, lw=1.5, alpha=0.65)
ax.plot([O[0], B[0]], [O[1], B[1]], color=C_R, lw=1.5, alpha=0.65)
ax.add_patch(Circle(A, 0.12, fc=C_A, ec="white", lw=1.2, zorder=5))
ax.add_patch(Circle(B, 0.12, fc=C_B, ec="white", lw=1.2, zorder=5))
ax.text(A[0] - 0.35, A[1] + 0.18, "partial system A",
        ha="right", fontsize=9.5, color=C_A, weight="bold")
ax.text(B[0] + 0.35, B[1] + 0.18, "partial system B",
        ha="left", fontsize=9.5, color=C_B, weight="bold")

# Emphasize separated external positions.
ax.add_patch(Arc(O, 2 * R + 0.18, 2 * R + 0.18,
                 theta1=28, theta2=142, color=C_MUTED, lw=1.2,
                 linestyle="--"))
ax.text(O[0], O[1] + R + 0.34, "externally separated positions",
        ha="center", fontsize=8.8, color=C_MUTED)

# Readout boxes.
boxA = label_box(ax, (1.75, 1.35), 2.55, 1.22, "readout A",
                 ["local record", "R_A read here"], C_A)
boxB = label_box(ax, (10.45, 1.35), 2.55, 1.22, "readout B",
                 ["local record", "R_B read here"], C_B)
arrow(ax, (A[0] - 0.05, A[1] - 0.10), (2.65, 1.95),
      color=C_A, lw=1.5, ms=12, rad=0.16)
arrow(ax, (B[0] + 0.05, B[1] - 0.10), (9.55, 1.95),
      color=C_B, lw=1.5, ms=12, rad=-0.16)

# Shared ledger entries as the source of correlation.
ledger = FancyBboxPatch(
    (3.05, 0.25), 6.10, 1.24,
    boxstyle="round,pad=0.05,rounding_size=0.09",
    fc="#fbf8ff", ec=C_SHARE, lw=1.5
)
ax.add_patch(ledger)
ax.text(6.10, 1.22, "shared entries that can leave correlation",
        ha="center", va="center", fontsize=10.2,
        color=C_SHARE, weight="bold")
entries = [
    "R conservation\ncommon radius",
    "phase relation",
    "internal sign structure",
    "readout condition",
]
for i, entry in enumerate(entries):
    x = 3.70 + i * 1.50
    ax.text(x, 0.78, entry, ha="center", va="center",
            fontsize=7.2, color=C_TEXT, linespacing=1.05)
    ax.plot([x - 0.48, x + 0.48], [0.54, 0.54],
            color="#d8c6f1", lw=1.0)

arrow(ax, (3.10, 1.10), (3.95, 0.92), color=C_SHARE,
      lw=1.2, ms=10, rad=-0.16, linestyle="--", alpha=0.9)
arrow(ax, (9.10, 1.10), (8.23, 0.92), color=C_SHARE,
      lw=1.2, ms=10, rad=0.16, linestyle="--", alpha=0.9)
ax.text(6.1, 0.06,
        "Correlation is read after comparing local records; the diagram does not represent a signal from A to B.",
        ha="center", fontsize=8.6, color=C_NOTE)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch15_shared_R_correlation.{ext}",
                dpi=200, bbox_inches="tight")
print("done: fig_text_ch15_shared_R_correlation (png+svg)")
