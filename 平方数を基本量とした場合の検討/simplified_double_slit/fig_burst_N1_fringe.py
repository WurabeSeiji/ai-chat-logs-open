#!/usr/bin/env python3
"""
Paper 2, Figure 1 (N=1), centred source.
EXACTLY Paper 1's calculation. Nothing added.
    psi = exp(i Phi1) + exp(i Phi2),  I = |psi|^2 / 4 = psi * conj(psi) / 4
    Phi_k = 2 pi ( r_k - y_slit,k s ) / lambda0,  r_k = sqrt(L^2+(y-y_slit,k)^2)
Same parameters, same Phi0 grid.
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
Phi0_deg = np.linspace(-720.0, 720.0, 16001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W / lam0)

y = 0.0
r1 = np.sqrt(L**2 + (y - ysl1)**2)
r2 = np.sqrt(L**2 + (y - ysl2)**2)
Phi1 = 2.0 * np.pi * (r1 - ysl1 * s) / lam0
Phi2 = 2.0 * np.pi * (r2 - ysl2 * s) / lam0
psi = np.exp(1j * Phi1) + np.exp(1j * Phi2)
I = (psi * np.conj(psi)).real / 4.0

fig, ax = plt.subplots(figsize=(13.0, 5.2))
ax.plot(Phi0_deg, I, color="#1f5fbf", lw=2.3, label=r"$N=1$ (centred source)")
for k in range(-2, 3):
    ax.axvline(360 * k, color="0.92", lw=0.8, zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)
ax.set_xlim(-720, 720)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-720, 721, 180))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s/\lambda_0$  (degrees)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=10, framealpha=0.96)
cap = (r"Paper 2, Fig. 1 ($N=1$), centred source. Exactly Paper 1's calculation "
       r"($\psi=e^{i\Phi_1}+e^{i\Phi_2}$, $I=|\psi|^2/4$), same parameters and "
       r"same $\Phi_0$ axis. Nothing added.")
fig.text(0.5, 0.01, cap, ha="center", va="bottom", fontsize=9, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.15)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_burst_N1_fringe.{ext}"), dpi=200)
print("saved fig_burst_N1_fringe.png / .svg")
