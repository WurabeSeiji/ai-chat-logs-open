#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chapter 12 figure:
  A local frame is transported through separated readout faces.
  After a closed loop, the returned frame differs from the starting frame.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Arc, Circle


C_TEXT = "#222222"
C_MUTED = "#777777"
C_A = "#1f77b4"
C_B = "#ff7f0e"
C_C = "#2ca02c"
C_CONN = "#6f42c1"
C_RETURN = "#d62728"
C_FACE = "#f5f7fb"
C_NOTE = "#8a4b00"


def arrow(ax, start, end, color="#444444", lw=1.6, ms=13,
          rad=0.0, style="-|>", linestyle="-", alpha=1.0):
    arr = FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=ms, lw=lw,
        color=color, linestyle=linestyle, alpha=alpha,
        connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)
    return arr


def draw_face(ax, center, title, subtitle, color):
    x, y = center
    w, h = 2.55, 1.10
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        fc=C_FACE, ec=color, lw=1.5
    )
    ax.add_patch(box)
    ax.text(x, y + 0.26, title, ha="center", va="center",
            fontsize=11.5, color=color, weight="bold")
    ax.text(x, y - 0.10, subtitle, ha="center", va="center",
            fontsize=8.2, color=C_TEXT)
    for i in range(3):
        ax.plot([x - 0.82, x + 0.82], [y - 0.34 + i * 0.16, y - 0.34 + i * 0.16],
                color="#d7dee8", lw=0.8)


def draw_frame(ax, origin, angle, color, label, dashed=False, alpha=1.0):
    x, y = origin
    length = 0.62
    e1 = (x + length * np.cos(angle), y + length * np.sin(angle))
    e2 = (x + length * np.cos(angle + np.pi / 2),
          y + length * np.sin(angle + np.pi / 2))
    linestyle = "--" if dashed else "-"
    arrow(ax, (x, y), e1, color=color, lw=2.0, ms=12,
          linestyle=linestyle, alpha=alpha)
    arrow(ax, (x, y), e2, color=color, lw=2.0, ms=12,
          linestyle=linestyle, alpha=alpha)
    ax.add_patch(Circle((x, y), 0.035, fc=color, ec=color, alpha=alpha))
    if label:
        ax.text(x, y - 0.36, label, ha="center", va="top",
                fontsize=8.4, color=color)


fig, ax = plt.subplots(figsize=(12.0, 6.2))
ax.set_xlim(0, 12)
ax.set_ylim(0, 6.2)
ax.axis("off")

ax.text(6.0, 5.86, "Frame transport and holonomy",
        ha="center", fontsize=15, color=C_TEXT, weight="bold")
ax.text(6.0, 5.55,
        "a local reference frame is carried through readout faces; the return mismatch is recorded",
        ha="center", fontsize=9.5, color=C_MUTED)

# Readout faces on a closed path.
A = (2.25, 4.35)
B = (9.75, 4.20)
C = (6.05, 1.42)
P = (2.25, 2.48)
draw_face(ax, A, "readout face A", "starting local ledger", C_A)
draw_face(ax, B, "readout face B", "different local convention", C_B)
draw_face(ax, C, "readout face C", "another local convention", C_C)

# Connection loop.
arrow(ax, (3.65, 4.35), (8.35, 4.22), color=C_CONN, lw=2.0, ms=16, rad=-0.05)
arrow(ax, (9.20, 3.57), (6.70, 2.00), color=C_CONN, lw=2.0, ms=16, rad=-0.10)
arrow(ax, (5.35, 1.96), (2.82, 2.55), color=C_CONN, lw=2.0, ms=16, rad=-0.08)
ax.text(6.0, 4.70, "transport by connection", ha="center",
        fontsize=9.6, color=C_CONN, weight="bold")
ax.text(8.70, 2.72, "compare local frames", ha="center",
        fontsize=8.6, color=C_CONN)
ax.text(3.35, 2.18, "return to A", ha="center",
        fontsize=8.6, color=C_CONN)

# Frames carried around the loop.
draw_frame(ax, (2.25, 4.23), np.deg2rad(0), C_A, "start frame")
draw_frame(ax, (9.75, 4.08), np.deg2rad(28), C_B, "after A -> B")
draw_frame(ax, (6.05, 1.30), np.deg2rad(58), C_C, "after B -> C")

# Returned frame at A compared with original.
compare_box = FancyBboxPatch(
    (P[0] - 1.12, P[1] - 0.78), 2.24, 1.40,
    boxstyle="round,pad=0.04,rounding_size=0.08",
    fc="#fffdf8", ec="#c9b27b", lw=1.2, linestyle="--"
)
ax.add_patch(compare_box)
ax.text(P[0], P[1] + 0.50, "compare again at A", ha="center",
        va="center", fontsize=9.2, color=C_TEXT, weight="bold")
draw_frame(ax, P, np.deg2rad(0), C_A, "",
           dashed=True, alpha=0.75)
draw_frame(ax, P, np.deg2rad(76), C_RETURN, "")
ax.text(P[0] - 0.70, P[1] - 0.50, "original", ha="center",
        va="center", fontsize=8.0, color=C_A)
ax.text(P[0] + 0.70, P[1] - 0.50, "returned", ha="center",
        va="center", fontsize=8.0, color=C_RETURN)
ax.add_patch(Arc(P, 1.05, 1.05, theta1=6, theta2=76,
                 color=C_RETURN, lw=2.0))
arrow(ax, (P[0] + 0.37, P[1] + 0.27), (P[0] + 0.18, P[1] + 0.48),
      color=C_RETURN, lw=1.6, ms=10)
ax.text(3.50, 2.70, "holonomy", ha="left", va="center",
        fontsize=11.5, color=C_RETURN, weight="bold")
ax.text(3.50, 2.42, "residual mismatch\nleft by the loop", ha="left",
        va="center", fontsize=8.6, color=C_TEXT)

# Safe interpretation note.
ax.text(6.0, 0.42,
        "This records a mismatch made visible by connection + closed transport; it is not a derivation of Yang-Mills or gravity equations.",
        ha="center", fontsize=9.0, color=C_NOTE)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch12_holonomy_transport.{ext}",
                dpi=200, bbox_inches="tight")
print("done: fig_text_ch12_holonomy_transport (png+svg)")
