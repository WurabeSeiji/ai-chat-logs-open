#!/usr/bin/env python3
"""
Pushforward demonstration (observation mode (a)): random sampling of the source
position y ~ P(y) = cos^2(pi y / lambda), one shift value u(y) read per trial,
histogrammed over many trials.

For each sampled y the central-fringe shift is the EXACT
    u(y) = 360 * (r_1(y) - r_2(y)) / lambda0   [degrees],
    r_k = sqrt(L^2 + (y - y_slit,k)^2).
The histogram of u is the pushforward rho(u) = P(y(u)) |dy/du|. Because the map
u(y) is nearly linear (paraxial), the histogram keeps the cos^2 SHAPE of P; the
small departure from the ideal cos^2(pi u/180) is the non-paraxial nonlinearity
(~ 1/2 tan^2(theta) ~ 3%).

This is mode (a): per-trial shift read out and histogrammed. It is NOT the
intensity-overlay (mode (b)), which would wash to a reduced-visibility fringe.

Outputs: fig_shift_histogram.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.size": 12, "mathtext.fontset": "cm",
})

# ---- geometry ------------------------------------------------------------
L = 10.0
W = 5.0
lam0 = 1.0
half = lam0 / 2.0
ysl1, ysl2 = +W / 2.0, -W / 2.0

# ---- sample y ~ P(y) = cos^2(pi y / lambda) on [-lambda/2, lambda/2] -------
np.random.seed(0)
N = 600000
yc = np.random.uniform(-half, half, size=3 * N)
acc = np.random.uniform(0.0, 1.0, size=3 * N) < np.cos(np.pi * yc / lam0)**2
y = yc[acc][:N]

# ---- exact per-trial shift u(y) (degrees) --------------------------------
r1 = np.sqrt(L**2 + (y - ysl1)**2)
r2 = np.sqrt(L**2 + (y - ysl2)**2)
u = 360.0 * (r1 - r2) / lam0          # = 2 pi (r1 - r2)/lambda0 in degrees

# ---- ideal (linear pushforward) cos^2 shape ------------------------------
ug = np.linspace(-90.0, 90.0, 600)
ideal = np.cos(np.pi * ug / 180.0)**2

# ---- plot ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11.5, 5.4))

bins = np.linspace(-95, 95, 96)
ax.hist(u, bins=bins, density=True, color="#9ecae1", edgecolor="#4a90c2",
        lw=0.5, label=r"histogram of shift $u$ over trials  ($\rho(u)$, $N=6\times10^5$)")

# normalise ideal to the same area as the histogram (density)
# cos^2 over [-90,90] deg has integral 90 (deg); density-normalised peak = 1/90
ideal_density = ideal / 90.0
ax.plot(ug, ideal_density, color="#d11f2d", lw=2.4, ls="--",
        label=r"ideal linear pushforward $\propto\cos^2(\pi u/180^\circ)$")

ax.axvline(-90, color="0.7", lw=0.8, ls=":")
ax.axvline(+90, color="0.7", lw=0.8, ls=":")
ax.set_xlim(-100, 100)
ax.set_ylim(0, None)
ax.set_xlabel(r"central-fringe shift  $u$  (degrees of $\Phi_0$)")
ax.set_ylabel(r"probability density  $\rho(u)$")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)

# quantify the non-paraxial edge compression
u_edge = 360.0 * (np.sqrt(L**2 + (half - ysl1)**2)
                  - np.sqrt(L**2 + (half - ysl2)**2)) / lam0
info = (rf"$L={L:.0f},\ W={W:.0f},\ \lambda_0={lam0:.0f}$" "\n"
        rf"linear edge: $u=\mp90.0^\circ$, exact edge: $u=\mp{abs(u_edge):.1f}^\circ$" "\n"
        rf"non-paraxial compression $0.5\,\tan^2\theta\approx{0.5*np.tan(np.arctan(0.25))**2*100:.1f}\%$")
ax.text(0.013, 0.97, info, transform=ax.transAxes, va="top", ha="left",
        fontsize=9, family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7"))

cap = (r"Observation mode (a): the source position $y$ is drawn from "
       r"$P(y)=\cos^2(\pi y/\lambda)$ each trial, and the single central-fringe "
       r"shift $u(y)$ is read and histogrammed. The histogram (the pushforward "
       r"$\rho(u)$) keeps the $\cos^2$ shape of $P$ because the geometric map "
       r"$u(y)$ is nearly linear; the slight edge compression vs the ideal "
       r"$\cos^2$ (dashed) is of the non-paraxial scale $\frac{1}{2}\tan^2\theta\approx3\%$. "
       r"No source-position measurement or conditioning is used -- only the "
       r"fringe shift is read.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.7, wrap=True)

fig.subplots_adjust(left=0.08, right=0.985, top=0.96, bottom=0.21)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_shift_histogram.{ext}"), dpi=200)

print(f"samples used: {len(y)}")
print(f"linear edge u = +/-90.0 deg, exact edge u = +/-{abs(u_edge):.2f} deg")
print(f"compression = {(90-abs(u_edge))/90*100:.2f}%  vs  0.5 tan^2(theta) = "
      f"{0.5*np.tan(np.arctan(0.25))**2*100:.2f}%")
print("saved fig_shift_histogram.png / .svg")
