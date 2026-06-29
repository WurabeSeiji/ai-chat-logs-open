#!/usr/bin/env python3
"""
Enlarged setup (source -> double slit only; screen removed) with the source
position uncertainty drawn in:
  - a circle of radius lambda/2 (the +/- lambda/2 position uncertainty),
  - the cos^2 source-position distribution (Model B) along the diameter that
    passes through the centre PARALLEL to the screen (the transverse axis).

Geometry kept at L=10, W=5 (lambda0=1). A zoomed inset shows the small
source region clearly.

Outputs: fig_setup_source_uncertainty.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle

plt.rcParams.update({
    "font.family": "serif", "font.size": 12, "mathtext.fontset": "cm",
})

# --- geometry -------------------------------------------------------------
L = 10.0
W = 5.0
lam0 = 1.0
half = lam0 / 2.0                 # lambda/2 = 0.5  (the uncertainty radius)
sy = W / 2.0                      # slit offset = 2.5
theta = np.degrees(np.arctan2(sy, L))
barrier_x = L
A_cos = 1.4                        # drawing amplitude of the cos^2 bulge


def draw_source_uncertainty(ax, A=A_cos, lw_scale=1.0):
    """Draw the lambda/2 circle, the transverse diameter, and the cos^2 curve."""
    circ = Circle((0, 0), half, fc="#cfe3ff", ec="#1f5fbf",
                  lw=1.6 * lw_scale, alpha=0.85, zorder=3)
    ax.add_patch(circ)
    # transverse diameter (parallel to screen) through the centre
    ax.plot([0, 0], [-half, half], ls=(0, (4, 3)), color="#1f5fbf",
            lw=1.2 * lw_scale, zorder=4)
    # cos^2 source-position distribution, bulging to the left (-x)
    yy = np.linspace(-half, half, 400)
    vv = np.cos(np.pi * yy / lam0)**2          # =1 at y=0, =0 at y=+/-lambda/2
    ax.fill_betweenx(yy, -A * vv, 0, color="#f0a500", alpha=0.22, zorder=2)
    ax.plot(-A * vv, yy, color="#e08000", lw=2.3 * lw_scale, zorder=5)


fig, ax = plt.subplots(figsize=(11.5, 6.2))

# optical axis
ax.plot([-2.6, barrier_x + 1.2], [0, 0], ls="--", color="0.55", lw=1.0, zorder=1)

# rays: source -> slits
ax.plot([0, L], [0, sy], color="#1f5fbf", lw=1.5, zorder=2)
ax.plot([0, L], [0, -sy], color="#1f5fbf", lw=1.5, zorder=2)

# source uncertainty (true scale)
draw_source_uncertainty(ax)
ax.plot(0, 0, "o", color="black", ms=5, zorder=6)

# half-angle arc
arc = Arc((0, 0), 4.6, 4.6, angle=0, theta1=0, theta2=theta,
          color="#c0392b", lw=1.5, zorder=4)
ax.add_patch(arc)
ax.annotate(rf"$\theta\approx{theta:.1f}^\circ$", (2.6, 0.5),
            color="#c0392b", fontsize=11)

# barrier with two slits (screen removed)
gap = 0.42
for (y0, y1) in [(-4.3, -sy - gap), (-sy + gap, sy - gap), (sy + gap, 4.3)]:
    ax.plot([barrier_x, barrier_x], [y0, y1], color="black", lw=5,
            solid_capstyle="butt", zorder=4)
for s, lab, dy in [(sy, r"Slit $1$", 9), (-sy, r"Slit $2$", -18)]:
    ax.plot(barrier_x, s, "o", mfc="white", mec="black", ms=7, zorder=5)
    ax.annotate(lab, (barrier_x, s), textcoords="offset points",
                xytext=(8, dy), fontsize=11)

# labels for the source-uncertainty elements
ax.annotate(r"$\pm\lambda/2$ position" "\n" "uncertainty (radius $\\lambda/2$)",
            (0, half), textcoords="offset points", xytext=(14, 30),
            fontsize=10, color="#1f5fbf",
            arrowprops=dict(arrowstyle="->", color="#1f5fbf", lw=1.0))
ax.annotate(r"$\cos^2$ source distribution" "\n" "(Model B), span $\\pm\\lambda/2$",
            (-A_cos, 0), textcoords="offset points", xytext=(-12, -54),
            fontsize=10, color="#b36b00", ha="center",
            arrowprops=dict(arrowstyle="->", color="#b36b00", lw=1.0))

# dimension L
yL = -4.9
ax.annotate("", (0, yL), (L, yL), arrowprops=dict(arrowstyle="<->", color="black", lw=1.1))
ax.text(L / 2, yL - 0.3, r"$L = 10$", ha="center", va="top", fontsize=11)
# dimension W
xW = barrier_x + 0.7
ax.annotate("", (xW, sy), (xW, -sy), arrowprops=dict(arrowstyle="<->", color="black", lw=1.1))
ax.text(xW + 0.2, 0, r"$W = 5$", ha="left", va="center", fontsize=11)

# source label
ax.annotate(r"Source $S$ ($\lambda_0=1,\ c=1$)", (0, -half),
            textcoords="offset points", xytext=(40, -16), ha="center", fontsize=10)

# --- frame ----------------------------------------------------------------
ax.set_xlim(-2.7, barrier_x + 1.7)
ax.set_ylim(-5.3, 5.0)
ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)

cap = (r"Enlarged source-to-slit region ($L=10$, $W=5$; screen removed). The "
       r"source position carries a $\pm\lambda/2$ uncertainty (blue circle, "
       r"radius $\lambda/2$). The $\cos^2$ source-position distribution "
       r"(Model B; orange) is drawn along the diameter through the centre, "
       r"parallel to the screen, spanning $\pm\lambda/2$ (peak at the centre, "
       r"zero at $\pm\lambda/2$).")
fig.text(0.5, 0.012, cap, ha="center", va="bottom", fontsize=9, wrap=True)

fig.subplots_adjust(left=0.02, right=0.98, top=0.99, bottom=0.12)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_setup_source_uncertainty.{ext}"), dpi=200)
print(f"theta = {theta:.4f} deg, lambda/2 = {half}")
print("saved fig_setup_source_uncertainty.png / .svg")
