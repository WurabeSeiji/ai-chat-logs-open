#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第3章導入用の図:
  (a) セル対 (a,b) と D(a,b)=(b,-a) の配線
  (b) D の反復: 90°遅延を2回で NOT
  (c) 整数格子で厳密に閉じる離散位相の代表 (Z4/Z6)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

C_A = "#1f77b4"
C_B = "#ff7f0e"
C_NEG = "#d62728"
C_GRID = "#cccccc"
C_TEXT = "#222222"
C_GREEN = "#2ca02c"


def add_box(ax, xy, text, fc, ec="#333333"):
    x, y = xy
    box = Rectangle((x - 0.28, y - 0.18), 0.56, 0.36,
                    facecolor=fc, edgecolor=ec, lw=1.2)
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center", fontsize=12, color=C_TEXT)
    return box


def add_arrow(ax, start, end, color="#333333", rad=0.0, lw=1.5):
    arr = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                          lw=lw, color=color,
                          connectionstyle=f"arc3,rad={rad}")
    ax.add_patch(arr)
    return arr


fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))

# ---------- (a) wiring ----------
ax = axes[0]
ax.set_title("(a) Two integer cells and one inverter", fontsize=12)
add_box(ax, (-1.35, 0.55), r"$a$", "#e9f2ff", C_A)
add_box(ax, (-1.35, -0.55), r"$b$", "#fff3e6", C_B)
add_box(ax, (1.35, 0.55), r"$b$", "#fff3e6", C_B)
add_box(ax, (1.35, -0.55), r"$-a$", "#ffecec", C_NEG)

add_arrow(ax, (-1.05, -0.55), (1.05, 0.55), C_B, rad=0.0, lw=1.7)
add_arrow(ax, (-1.05, 0.55), (0.15, -0.55), C_A, rad=0.0, lw=1.7)
inv = Circle((0.42, -0.55), 0.16, facecolor="white", edgecolor=C_NEG, lw=1.4)
ax.add_patch(inv)
ax.text(0.42, -0.55, "-", ha="center", va="center", fontsize=14, color=C_NEG)
add_arrow(ax, (0.58, -0.55), (1.05, -0.55), C_NEG, rad=0.0, lw=1.7)

ax.text(0, 1.08, r"$D(a,b)=(b,-a)$", ha="center", fontsize=13)
ax.text(0, -1.18, "swap + sign flip; no multiplier", ha="center", fontsize=10)
ax.set_xlim(-2.05, 2.05)
ax.set_ylim(-1.45, 1.45)
ax.axis("off")

# ---------- (b) rotation / NOT ----------
ax = axes[1]
ax.set_title("(b) D twice gives NOT", fontsize=12)
ax.axhline(0, color="#777777", lw=0.8)
ax.axvline(0, color="#777777", lw=0.8)
ax.text(1.72, -0.12, r"$a$", fontsize=11)
ax.text(0.08, 1.72, r"$b$", fontsize=11)
theta = np.deg2rad(35)
r = 1.25
pts = []
labels = [r"$(a,b)$", r"$D$", r"$D^2=\mathrm{NOT}$", r"$D^3$"]
for k in range(4):
    ang = theta - k * np.pi / 2
    pts.append(np.array([r * np.cos(ang), r * np.sin(ang)]))

circle = plt.Circle((0, 0), r, fill=False, ls="--", lw=1.0, color="#aaaaaa")
ax.add_patch(circle)
for i, p in enumerate(pts):
    ax.plot(p[0], p[1], "o", color=[C_A, C_B, C_NEG, C_GREEN][i], ms=6)
    offset = np.array([0.12 if p[0] >= 0 else -0.52,
                       0.12 if p[1] >= 0 else -0.22])
    ax.text(p[0] + offset[0], p[1] + offset[1], labels[i], fontsize=10)
for i in range(4):
    start = pts[i]
    end = pts[(i + 1) % 4]
    add_arrow(ax, start * 0.93, end * 0.93, color="#555555", rad=-0.25, lw=1.1)
ax.text(0, -1.72, r"$D^4=\mathrm{id}$", ha="center", fontsize=11)
ax.set_xlim(-1.85, 1.85)
ax.set_ylim(-1.85, 1.85)
ax.set_aspect("equal")
ax.axis("off")

# ---------- (c) allowed lattices ----------
ax = axes[2]
ax.set_title("(c) Exact finite rotations on integer lattices", fontsize=12)

def mini_circle(ax, center, n, title, color):
    cx, cy = center
    rr = 0.53
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(cx + rr * np.cos(t), cy + rr * np.sin(t),
            color="#999999", lw=0.9, ls="--")
    for k in range(n):
        ang = 2 * np.pi * k / n
        x = cx + rr * np.cos(ang)
        y = cy + rr * np.sin(ang)
        ax.plot([cx, x], [cy, y], color=C_GRID, lw=0.8)
        ax.plot(x, y, "o", color=color, ms=5)
    ax.text(cx, cy - 0.82, title, ha="center", fontsize=10)

mini_circle(ax, (-0.85, 0.35), 4, r"$\mathbb{Z}_4$: 90 deg", C_A)
mini_circle(ax, (0.85, 0.35), 6, r"$\mathbb{Z}_6$: 60 deg", C_B)
ax.text(0, -0.95, r"allowed orders: $\{1,2,3,4,6\}$",
        ha="center", fontsize=11)
ax.text(0, -1.25, r"45 deg needs $1/\sqrt{2}$ -> approximation",
        ha="center", fontsize=9, color=C_NEG)
ax.set_xlim(-1.75, 1.75)
ax.set_ylim(-1.45, 1.35)
ax.axis("off")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_text_ch3_cellpair.{ext}", dpi=200, bbox_inches="tight")
print("done: fig_text_ch3_cellpair (png+svg)")
