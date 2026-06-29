#!/usr/bin/env python3
"""
Experimental-setup schematic for the single-source double-slit configuration.

Conditions:
    Source     : single point source, wavelength lambda_0 = 1, speed c = 1
    Source->slit plane distance  L = 10
    Slit separation              W = 5   (slits at y = +/- W/2 = +/- 2.5)
    Half-angle (axis to one slit): theta = arctan((W/2)/L) ~ 14.0 deg
    Slit->screen distance        D >> L  (far-field / Fraunhofer)

Outputs: fig_setup_double_slit.png, fig_setup_double_slit.svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "mathtext.fontset": "cm",
    "axes.linewidth": 0.0,
})

# --- geometry -------------------------------------------------------------
L = 10.0                       # source -> slit plane
W = 5.0                        # slit separation
sy = W / 2.0                   # slit offset = 2.5
theta = np.degrees(np.arctan2(sy, L))   # half-angle ~ 14.04 deg
barrier_x = L
screen_x = 17.0                # drawn position (not to scale)

fig, ax = plt.subplots(figsize=(10.5, 5.6))

# optical axis
ax.plot([-0.5, screen_x + 0.3], [0, 0], ls="--", color="0.5", lw=1.0, zorder=1)

# --- source ---------------------------------------------------------------
ax.plot(0, 0, "o", color="black", ms=10, zorder=6)
ax.annotate(r"Source $S$" "\n" r"$\lambda_0 = 1,\ c = 1$",
            (0, 0), textcoords="offset points", xytext=(-6, -42),
            ha="center", va="top", fontsize=11)

# --- rays: source -> slits ------------------------------------------------
ax.plot([0, L], [0, sy], color="#1f5fbf", lw=1.6, zorder=3)
ax.plot([0, L], [0, -sy], color="#1f5fbf", lw=1.6, zorder=3)

# --- half-angle arc at source --------------------------------------------
arc = Arc((0, 0), 5.4, 5.4, angle=0, theta1=0, theta2=theta,
          color="#c0392b", lw=1.6, zorder=4)
ax.add_patch(arc)
ax.annotate(rf"$\theta \approx {theta:.1f}^\circ$",
            (3.0, 0.62), color="#c0392b", fontsize=11)

# --- barrier with two slits ----------------------------------------------
gap = 0.42
for (y0, y1) in [(-4.6, -sy - gap), (-sy + gap, sy - gap), (sy + gap, 4.6)]:
    ax.plot([barrier_x, barrier_x], [y0, y1], color="black", lw=5,
            solid_capstyle="butt", zorder=4)
# slit markers
for s, lab, dy in [(sy, r"Slit $1$", 10), (-sy, r"Slit $2$", -20)]:
    ax.plot(barrier_x, s, "o", mfc="white", mec="black", ms=7, zorder=5)
    ax.annotate(lab, (barrier_x, s), textcoords="offset points",
                xytext=(9, dy), fontsize=11)

# --- rays: slits -> screen (far-field, ~parallel) -------------------------
for s in [sy, -sy]:
    ax.plot([barrier_x, 13.0], [s, s], color="#1f5fbf", lw=1.2, alpha=0.6, zorder=2)
    ax.plot([14.2, screen_x], [s, s], color="#1f5fbf", lw=1.2, alpha=0.6, zorder=2)

# not-to-scale break on the axis (//)
for xb in (13.45, 13.75):
    ax.plot([xb - 0.12, xb + 0.12], [-0.32, 0.32], color="0.4", lw=1.2, zorder=3)

# --- screen ---------------------------------------------------------------
ax.plot([screen_x, screen_x], [-4.8, 4.8], color="black", lw=3.5, zorder=4)
ax.annotate("Observation screen\n" r"(far-field, $D \gg L$)",
            (screen_x, -4.8), textcoords="offset points", xytext=(0, -10),
            ha="center", va="top", fontsize=11)

# --- dimension: L ---------------------------------------------------------
yL = -5.6
ax.annotate("", (0, yL), (L, yL),
            arrowprops=dict(arrowstyle="<->", color="black", lw=1.1))
ax.text(L / 2, yL - 0.35, r"$L = 10$", ha="center", va="top", fontsize=11)

# --- dimension: W ---------------------------------------------------------
xW = barrier_x + 1.0
ax.annotate("", (xW, sy), (xW, -sy),
            arrowprops=dict(arrowstyle="<->", color="black", lw=1.1))
ax.text(xW + 0.2, 0, r"$W = 5$", ha="left", va="center", fontsize=11)

# --- frame ----------------------------------------------------------------
ax.set_xlim(-2.2, screen_x + 1.6)
ax.set_ylim(-7.2, 5.6)
ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

cap = (r"Single-source double-slit configuration. A point source $S$ "
       r"($\lambda_0=1$, $c=1$) illuminates two slits separated by $W=5$ at "
       r"distance $L=10$; each slit subtends a half-angle "
       r"$\theta=\arctan[(W/2)/L]\approx 14.0^\circ$ from the optical axis "
       r"(full subtended angle $2\theta\approx 28.1^\circ$). The screen lies "
       r"in the far field, $D\gg L$ (horizontal scale broken, not to scale).")
fig.text(0.5, 0.015, cap, ha="center", va="bottom", fontsize=9.5, wrap=True)

fig.subplots_adjust(left=0.02, right=0.98, top=0.99, bottom=0.16)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_setup_double_slit.{ext}"),
                dpi=200)
print("saved fig_setup_double_slit.png / .svg")
