#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chapter 11 figure:
  Two ledgers meet at a shared boundary.
  A connection aligns local frames / phase origins.
  A vertex map writes an external record.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Polygon


C_TEXT = "#222222"
C_MUTED = "#777777"
C_A = "#1f77b4"
C_B = "#ff7f0e"
C_BOUNDARY = "#6f42c1"
C_CONN = "#2ca02c"
C_VERTEX = "#d62728"
C_RECORD = "#8a5a00"
C_LEDGER = "#f5f7fb"


def arrow(ax, start, end, color="#444444", lw=1.6, ms=13, rad=0.0):
    arr = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=ms, lw=lw,
        color=color, connectionstyle=f"arc3,rad={rad}"
    )
    ax.add_patch(arr)
    return arr


def draw_ledger(ax, x, y, title, color, shift=0.0):
    w, h = 2.25, 1.35
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h,
                           fc=C_LEDGER, ec=color, lw=1.5))
    ax.text(x, y + 0.45, title, ha="center", va="center",
            fontsize=12, color=color, weight="bold")
    labels = ["x", "t", "R", "Q"]
    for i, lab in enumerate(labels):
        yy = y + 0.15 - i * 0.23
        ax.add_patch(Rectangle((x - 0.82, yy - 0.07), 1.64, 0.14,
                               fc="white", ec="#c9d1db", lw=0.7))
        ax.text(x - 0.68, yy, lab, ha="center", va="center",
                fontsize=8.5, color=C_TEXT)
        ax.plot([x - 0.47, x + 0.64], [yy, yy],
                color="#d8dee7", lw=1.0)
        ax.add_patch(Circle((x - 0.12 + 0.17 * ((i + shift) % 3), yy),
                            0.035, fc=color, ec=color, lw=0.8))
    ax.text(x, y - 0.86, "private entries", ha="center",
            fontsize=8.5, color=C_MUTED)


def draw_boundary(ax):
    ax.add_patch(Rectangle((3.55, 2.67), 4.90, 1.15,
                           fc="#f4efff", ec=C_BOUNDARY, lw=1.7))
    ax.text(6.0, 3.52, "shared boundary", ha="center",
            fontsize=12.5, color=C_BOUNDARY, weight="bold")
    for x in [4.35, 5.15, 6.0, 6.85, 7.65]:
        ax.plot([x, x], [2.82, 3.24], color="#b99fea", lw=1.0)
    ax.text(6.0, 2.86, "differences / coincidences / exchange marks",
            ha="center", fontsize=9.2, color=C_TEXT)


def draw_connection(ax):
    ax.add_patch(Rectangle((4.30, 1.88), 3.40, 0.50,
                           fc="#eef9ef", ec=C_CONN, lw=1.5))
    ax.text(6.0, 2.16, "connection filter", ha="center",
            fontsize=11, color=C_CONN, weight="bold")
    ax.text(6.0, 1.96, "align scale / phase origin / local frame",
            ha="center", fontsize=8.4, color=C_TEXT)


def draw_vertex(ax):
    diamond = Polygon([[6.0, 1.55], [6.65, 1.08],
                       [6.0, 0.62], [5.35, 1.08]],
                      closed=True, fc="#fff1f1", ec=C_VERTEX, lw=1.6)
    ax.add_patch(diamond)
    ax.text(6.0, 1.12, "vertex\nmap", ha="center", va="center",
            fontsize=10.5, color=C_VERTEX, weight="bold")


def draw_record(ax):
    x, y, w, h = 9.70, 0.76, 2.25, 1.12
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h,
                           fc="#fff8e8", ec=C_RECORD, lw=1.5))
    ax.text(x, y + 0.34, "external record", ha="center",
            fontsize=11.5, color=C_RECORD, weight="bold")
    for i, lab in enumerate(["coincidence", "delta", "exchange"]):
        yy = y + 0.10 - i * 0.24
        ax.plot([x - 0.78, x + 0.78], [yy, yy], color="#e4c88d", lw=1.0)
        ax.text(x, yy + 0.04, lab, ha="center", va="bottom",
                fontsize=8.0, color=C_TEXT)


fig, ax = plt.subplots(figsize=(12.0, 6.2))
ax.set_xlim(0, 12)
ax.set_ylim(0, 6.2)
ax.axis("off")

ax.text(6.0, 5.86, "Interaction as shared ledger readout",
        ha="center", fontsize=15, color=C_TEXT, weight="bold")
ax.text(6.0, 5.55,
        "not a force drawn on background time; a record produced through shared boundary + connection",
        ha="center", fontsize=9.5, color=C_MUTED)

draw_ledger(ax, 2.55, 4.55, "ledger A", C_A, shift=0)
draw_ledger(ax, 9.45, 4.55, "ledger B", C_B, shift=1)
draw_boundary(ax)
draw_connection(ax)
draw_vertex(ax)
draw_record(ax)

arrow(ax, (3.42, 4.04), (4.25, 3.80), color=C_A, lw=1.8, ms=14, rad=0.05)
arrow(ax, (8.58, 4.04), (7.75, 3.80), color=C_B, lw=1.8, ms=14, rad=-0.05)
arrow(ax, (6.0, 2.67), (6.0, 2.40), color=C_BOUNDARY, lw=1.7, ms=13)
arrow(ax, (6.0, 1.88), (6.0, 1.56), color=C_CONN, lw=1.7, ms=13)
arrow(ax, (6.65, 1.08), (8.55, 0.78), color=C_RECORD, lw=1.8, ms=14)

ax.text(3.28, 3.28, "same readout face", ha="right",
        fontsize=8.8, color=C_BOUNDARY)
ax.text(8.72, 3.28, "same readout face", ha="left",
        fontsize=8.8, color=C_BOUNDARY)
ax.text(3.1, 0.34,
        "Concrete laws, constants, and amplitudes are outside this figure.",
        ha="center", fontsize=9.0, color="#8a4b00")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch11_interaction_readout_flow.{ext}",
                dpi=200, bbox_inches="tight")
print("done: fig_text_ch11_interaction_readout_flow (png+svg)")
