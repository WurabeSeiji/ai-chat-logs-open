#!/usr/bin/env python3
"""
Paper 2, Figure 2 (N=1). Random sampling -> central cos^2 envelope.
EXACTLY Paper 1's calculation (same as fig_decomposition_static). Nothing added.
Same 9 sampling positions, same Phi0 grid, same weights, same envelope.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

L = 10.0
W = 5.0
lam0 = 1.0
half = lam0 / 2.0
ysl1, ysl2 = +W / 2.0, -W / 2.0

x_src = np.array([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
y_src = x_src * half
w_int = np.cos(np.pi * x_src / 2.0)**2

Phi0_deg = np.linspace(-720.0, 720.0, 16001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W / lam0)

def unit_intensity(y_i):
    r1 = np.sqrt(L**2 + (y_i - ysl1)**2)
    r2 = np.sqrt(L**2 + (y_i - ysl2)**2)
    Phi1 = 2.0 * np.pi * (r1 - ysl1 * s) / lam0
    Phi2 = 2.0 * np.pi * (r2 - ysl2 * s) / lam0
    psi = np.exp(1j * Phi1) + np.exp(1j * Phi2)
    return (psi * np.conj(psi)).real / 4.0

fig, ax = plt.subplots(figsize=(13.0, 5.8))

env_x = np.linspace(-90.0, 90.0, 600)
env_y = np.cos(np.pi * env_x / 180.0)**2
ax.plot(env_x, env_y, color="#f0a500", lw=3.2, zorder=5,
        label=r"$\cos^2$ envelope")

peaks_x, peaks_y = [], []
for x_i, y_i, w_i in zip(x_src, y_src, w_int):
    curve = w_i * unit_intensity(y_i)
    if x_i == 0.0:
        ax.plot(Phi0_deg, curve, color="#1f5fbf", lw=2.3, zorder=4,
                label=r"$x=0$ (centred, peak $1$)")
    else:
        ax.plot(Phi0_deg, curve, color="#2ca02c", lw=1.0, alpha=0.85, zorder=2)
    m = np.abs(Phi0_deg) < 180.0
    idx = np.argmax(curve[m])
    peaks_x.append(Phi0_deg[m][idx]); peaks_y.append(curve[m][idx])

ax.plot(peaks_x, peaks_y, 'o', color="#d11f2d", ms=5.5, zorder=6, label="fringe peaks")
ax.plot([], [], color="#2ca02c", lw=1.2,
        label=r"each source $x_i$ (peak $\cos^2(\pi x_i/2)$)")

for k in range(-2, 3):
    ax.axvline(360 * k, color="0.92", lw=0.8, zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)
ax.set_xlim(-720, 720)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-720, 721, 180))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s/\lambda_0$  (degrees)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=9.0, framealpha=0.96)
cap = (r"Paper 2, Fig. 2 ($N=1$). Exactly Paper 1's calculation. Same 9 sampling "
       r"positions $x_i\in[-0.8,0.8]$ and weights $\cos^2(\pi x_i/2)$; the fringe "
       r"peaks trace the $\cos^2$ envelope at the centre.")
fig.text(0.5, 0.01, cap, ha="center", va="bottom", fontsize=9, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.13)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_burst_N1_envelope.{ext}"), dpi=200)
print("saved fig_burst_N1_envelope.png / .svg")
