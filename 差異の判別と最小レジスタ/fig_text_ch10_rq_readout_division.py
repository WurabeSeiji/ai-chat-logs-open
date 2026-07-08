#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chapter 10 figure:
  R-side readout changes common scale / norm.
  Q-side readout changes signed internal phase.
"""
import numpy as np
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Arc


C_TEXT = "#222222"
C_MUTED = "#777777"
C_LEDGER = "#f5f7fb"
C_R = "#1f77b4"
C_Q_POS = "#d62728"
C_Q_NEG = "#2ca02c"
C_GRID = "#cfd6df"
C_AXIS = "#555555"


def arrow(ax, start, end, color=C_AXIS, lw=1.4, ms=12, rad=0.0):
    arr = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=ms, lw=lw,
        color=color, connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)
    return arr


def draw_small_dial(ax, center, angle, color, label):
    cx, cy = center
    r = 0.14
    ax.add_patch(Circle((cx, cy), r, fc="white", ec=color, lw=1.2))
    ax.plot([cx, cx + 0.78 * r * np.cos(angle)],
            [cy, cy + 0.78 * r * np.sin(angle)],
            color=color, lw=1.6)
    ax.text(cx, cy - 0.23, label, ha="center", va="top",
            fontsize=8.5, color=C_TEXT)


def draw_ledger(ax, x, y, w=2.45, h=1.35):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h,
                           fc=C_LEDGER, ec="#9aa8b7", lw=1.3))
    ax.text(x, y + 0.48, "same extended ledger", ha="center",
            va="center", fontsize=11, color=C_TEXT, weight="bold")
    labs = ["x", "y", "z", "t", "R", "Q"]
    xs = [-0.85, -0.50, -0.15, 0.20, 0.55, 0.90]
    for i, (lab, dx) in enumerate(zip(labs, xs)):
        col = C_R if lab == "R" else C_Q_POS if lab == "Q" else "#6a7480"
        draw_small_dial(ax, (x + dx, y - 0.12), 0.45 + i * 0.58, col, lab)
    ax.text(x, y - 0.55, "role names, not privileged labels",
            ha="center", va="center", fontsize=8.5, color=C_MUTED)


def draw_grid(ax, center, scale, color=C_R, alpha=1.0):
    cx, cy = center
    span = 1.3 * scale
    step = span / 4
    for i in range(5):
        xx = cx - span / 2 + i * step
        yy = cy - span / 2 + i * step
        ax.plot([xx, xx], [cy - span / 2, cy + span / 2],
                color=C_GRID, lw=1.0, alpha=alpha)
        ax.plot([cx - span / 2, cx + span / 2], [yy, yy],
                color=C_GRID, lw=1.0, alpha=alpha)
    ax.add_patch(Rectangle((cx - span / 2, cy - span / 2), span, span,
                           fc="none", ec=color, lw=1.6, alpha=alpha))
    for ang in np.linspace(0.15, 2 * np.pi + 0.15, 6, endpoint=False):
        arrow(ax, (cx, cy),
              (cx + 0.43 * span * np.cos(ang), cy + 0.43 * span * np.sin(ang)),
              color=color, lw=1.2, ms=9)


def draw_phase_dial(ax, center, angle, color, direction_label, rad):
    cx, cy = center
    r = 0.48
    ax.add_patch(Circle((cx, cy), r, fc="#fffdf8", ec="#a6a6a6", lw=1.2))
    ax.add_patch(Circle((cx, cy), 0.055, fc=color, ec=color, lw=1.0))
    ax.plot([cx, cx + 0.78 * r * np.cos(angle)],
            [cy, cy + 0.78 * r * np.sin(angle)],
            color=color, lw=2.4)
    ax.add_patch(Arc((cx, cy), 1.18 * r, 1.18 * r,
                     theta1=25 if rad > 0 else 205,
                     theta2=145 if rad > 0 else 325,
                     color=color, lw=1.8))
    if rad > 0:
        arrow(ax, (cx + 0.02, cy + 0.31), (cx - 0.17, cy + 0.33),
              color=color, lw=1.6, ms=10)
    else:
        arrow(ax, (cx - 0.02, cy - 0.31), (cx + 0.17, cy - 0.33),
              color=color, lw=1.6, ms=10)
    ax.text(cx, cy - 0.72, direction_label, ha="center",
            va="center", fontsize=9.5, color=color, weight="bold")


fig, ax = plt.subplots(figsize=(12.0, 6.0))
ax.set_xlim(0, 12)
ax.set_ylim(0, 6)
ax.axis("off")

ax.text(6, 5.68, "R/Q readout division: same ledger, different faces",
        ha="center", fontsize=15, color=C_TEXT, weight="bold")

# Center ledger
draw_ledger(ax, 6.0, 3.70)
arrow(ax, (4.72, 3.55), (3.05, 3.55), color=C_R, lw=1.7, ms=14)
arrow(ax, (7.28, 3.55), (8.95, 3.55), color=C_Q_POS, lw=1.7, ms=14)

# R panel
ax.text(2.2, 5.05, "R-side readout", ha="center",
        fontsize=12.5, color=C_R, weight="bold")
ax.text(2.2, 4.72, "norm / scale / curvature side", ha="center",
        fontsize=9.5, color=C_TEXT)
draw_grid(ax, (1.55, 2.85), 0.72, color=C_R, alpha=0.65)
draw_grid(ax, (2.85, 2.85), 1.08, color=C_R, alpha=1.0)
arrow(ax, (2.05, 1.63), (2.55, 1.63), color=C_R, lw=1.6, ms=12)
ax.text(2.30, 1.36, "common scale changes", ha="center",
        fontsize=10, color=C_TEXT)
ax.text(2.30, 1.08, "all components are read through one size ledger",
        ha="center", fontsize=8.5, color=C_MUTED)

# Q panel
ax.text(9.8, 5.05, "Q-side readout", ha="center",
        fontsize=12.5, color=C_Q_POS, weight="bold")
ax.text(9.8, 4.72, "signed internal phase side", ha="center",
        fontsize=9.5, color=C_TEXT)
draw_phase_dial(ax, (9.05, 2.95), np.deg2rad(45), C_Q_POS, "+ phase", 1)
draw_phase_dial(ax, (10.55, 2.95), np.deg2rad(135), C_Q_NEG, "- phase", -1)
ax.plot([8.45, 11.15], [1.62, 1.62], color="#a6a6a6", lw=1.0)
ax.text(9.80, 1.36, "size can stay fixed; signed phase differs",
        ha="center", fontsize=10, color=C_TEXT)
ax.text(9.80, 1.08, "direction on an internal circle becomes readable",
        ha="center", fontsize=8.5, color=C_MUTED)

# Bottom note
ax.text(6, 0.48,
        "This is a division of readout roles, not a derivation of GR/EM field equations.",
        ha="center", fontsize=9.5, color="#8a4b00")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch10_rq_readout_division.{ext}",
                dpi=200, bbox_inches="tight")
print("done: fig_text_ch10_rq_readout_division (png+svg)")
