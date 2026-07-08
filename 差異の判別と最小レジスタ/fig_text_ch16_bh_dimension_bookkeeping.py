#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chapter 16 figure:
  Dimension bookkeeping for a 4D spatial black-hole-type structure in 5D+
  spacetime and an ordinary 3D spatial black-hole readout in 3+1 spacetime.
  The panels are juxtaposed without asserting an identification map.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, FancyArrowPatch


C_TEXT = "#222222"
C_MUTED = "#6d7480"
C_LEFT = "#6f42c1"
C_RIGHT = "#1f77b4"
C_HORIZON = "#d62728"
C_FACE_L = "#fbf8ff"
C_FACE_R = "#f7fbff"
C_NOTE = "#8a4b00"


def arrow(ax, start, end, color="#444444", lw=1.2, ms=11,
          rad=0.0, linestyle="-", alpha=1.0):
    arr = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=ms, lw=lw,
        color=color, linestyle=linestyle, alpha=alpha,
        connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)
    return arr


def rounded_panel(ax, xy, w, h, title, subtitle, ec, fc):
    x, y = xy
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.08,rounding_size=0.12",
        fc=fc, ec=ec, lw=1.6
    )
    ax.add_patch(p)
    ax.text(x + w / 2, y + h - 0.38, title,
            ha="center", va="center", fontsize=12.2,
            color=ec, weight="bold")
    ax.text(x + w / 2, y + h - 0.72, subtitle,
            ha="center", va="center", fontsize=8.6,
            color=C_MUTED)
    return p


fig, ax = plt.subplots(figsize=(12.4, 6.3))
ax.set_xlim(0, 12.4)
ax.set_ylim(0, 6.3)
ax.set_aspect("equal")
ax.axis("off")

ax.text(6.2, 5.92, "Dimension bookkeeping for black-hole-type readouts",
        ha="center", fontsize=15, color=C_TEXT, weight="bold")
ax.text(6.2, 5.58,
        "The panels are placed side by side as readout grammars; no identification map is asserted here.",
        ha="center", fontsize=9.2, color=C_MUTED)

# Panels
rounded_panel(
    ax, (0.55, 0.72), 5.28, 4.55,
    "5D+ spacetime",
    "4D spatial black-hole-type structure",
    C_LEFT, C_FACE_L
)
rounded_panel(
    ax, (6.58, 0.72), 5.28, 4.55,
    "3+1 spacetime",
    "3D spatial black-hole readout",
    C_RIGHT, C_FACE_R
)

# Separator
ax.plot([6.2, 6.2], [0.82, 5.10], color="#d0d6de", lw=1.2, linestyle="--")
ax.text(6.2, 0.48, "juxtaposition", ha="center", fontsize=8.4, color=C_MUTED)

# Left panel: 4D spatial object suggested by a 3D sphere plus an extra w axis.
O1 = np.array([3.18, 3.10])
ax.add_patch(Circle(O1, 1.05, fc="#ece2ff", ec=C_LEFT, lw=2.0, alpha=0.94))
ax.add_patch(Ellipse(O1, 2.10, 0.52, fill=False, ec=C_LEFT, lw=1.1, alpha=0.82))
ax.add_patch(Ellipse(O1, 0.70, 2.10, fill=False, ec=C_LEFT, lw=1.1, alpha=0.70))
ax.add_patch(Ellipse(O1, 1.55, 1.55, angle=45, fill=False, ec="#9c77d4", lw=1.0, alpha=0.65))
ax.add_patch(Circle(O1, 1.30, fill=False, ec=C_HORIZON, lw=2.0, linestyle="-"))
ax.text(O1[0], O1[1] + 1.58, "boundary response: S^3",
        ha="center", fontsize=8.8, color=C_HORIZON, weight="bold")
ax.text(O1[0], O1[1] - 1.58, "capacity / entropy side: ~ R^3",
        ha="center", fontsize=8.8, color=C_TEXT)
arrow(ax, (1.62, 2.02), (2.34, 2.46), color=C_LEFT, lw=1.1, ms=9, rad=-0.18)
ax.text(1.40, 1.84, "extra spatial\nledger direction",
        ha="center", va="center", fontsize=7.4, color=C_LEFT, linespacing=1.0)
ax.text(4.86, 1.72, "Tangherlini-type\ncomparison",
        ha="center", va="center", fontsize=7.8, color=C_MUTED, linespacing=1.0)

# Right panel: ordinary 3D spatial black hole readout with a 2-sphere horizon.
O2 = np.array([9.22, 3.10])
ax.add_patch(Circle(O2, 1.05, fc="#e7f2ff", ec=C_RIGHT, lw=2.0, alpha=0.96))
ax.add_patch(Ellipse(O2, 2.10, 0.56, fill=False, ec=C_RIGHT, lw=1.1, alpha=0.84))
ax.add_patch(Ellipse(O2, 0.72, 2.10, fill=False, ec=C_RIGHT, lw=1.1, alpha=0.70))
ax.add_patch(Circle(O2, 1.30, fill=False, ec=C_HORIZON, lw=2.0))
ax.text(O2[0], O2[1] + 1.58, "horizon readout: S^2",
        ha="center", fontsize=8.8, color=C_HORIZON, weight="bold")
ax.text(O2[0], O2[1] - 1.58, "area side: ~ R^2",
        ha="center", fontsize=8.8, color=C_TEXT)
ax.text(7.62, 1.82, "ordinary external\nspatial readout",
        ha="center", va="center", fontsize=7.8, color=C_MUTED, linespacing=1.0)
arrow(ax, (7.84, 2.05), (8.42, 2.52), color=C_RIGHT, lw=1.1, ms=9, rad=0.18)

# Shared vocabulary below.
shared_y = 0.98
for x, text, color in [
    (2.04, "internal modes", C_LEFT),
    (3.18, "horizon / boundary", C_HORIZON),
    (4.38, "few external hairs", C_MUTED),
    (8.08, "internal matter", C_RIGHT),
    (9.22, "horizon / boundary", C_HORIZON),
    (10.48, "few external hairs", C_MUTED),
]:
    ax.text(x, shared_y, text, ha="center", fontsize=7.1, color=color)

ax.text(6.2, 0.16,
        "This is a dimension-bookkeeping figure: it compares readout slots and boundary scaling, not a derivation of either black hole.",
        ha="center", fontsize=8.2, color=C_NOTE)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch16_bh_dimension_bookkeeping.{ext}",
                dpi=200, bbox_inches="tight")
print("done: fig_text_ch16_bh_dimension_bookkeeping (png+svg)")
