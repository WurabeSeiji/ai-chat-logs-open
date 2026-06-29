#!/usr/bin/env python3
"""
Diagnose the cage error: compare lambda=1/m vs lambda=13/m, both vs the SAME
base-phase axis phi (lambda=1) in degrees, -720..720.

A wave of wavelength lambda_m, as a function of base phase phi (phi_rad = 2 pi x,
x in lambda=1 units):  cos(2 pi x / lambda_m) = cos((1/lambda_m) phi_rad).
  lambda = 1/m   -> cos(m phi)              (cage = 1   -> spikes every 180 deg)
  lambda = 13/m  -> cos((m/13) phi)         (cage = 13  -> single spike over +-720)

m = 1,3,5,...,25 (13 odd harmonics), equal amplitude.

Outputs: fig_oddharm_cage_compare.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

m_list = np.arange(1, 26, 2)                 # 1,3,...,25
phi_deg = np.linspace(-720.0, 720.0, 20001)
phr = np.radians(phi_deg)

sum_1 = np.sum([np.cos(m * phr) for m in m_list], axis=0)          # lambda = 1/m
sum_13 = np.sum([np.cos((m / 13.0) * phr) for m in m_list], axis=0)  # lambda = 13/m

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13.0, 7.6), sharex=True)

ax1.plot(phi_deg, sum_1, color="#d11f2d", lw=1.3)
ax1.axhline(0, color="0.85", lw=0.8)
ax1.set_ylabel(r"$\sum_m\cos(m\varphi)$")
ax1.set_title(r"$\lambda=1/m$ (cage $=1$): dense spike train, $\pm13$ every $180^\circ$ "
              r"-- this is what I mistakenly plotted", fontsize=11)
for k in range(-4, 5):
    ax1.axvline(180 * k, color="0.92", lw=0.7, zorder=0)

ax2.plot(phi_deg, sum_13, color="#1f5fbf", lw=1.5)
ax2.axhline(0, color="0.85", lw=0.8)
ax2.set_ylabel(r"$\sum_m\cos((m/13)\varphi)$")
ax2.set_xlabel(r"base phase $\varphi$  (degrees, $\lambda=1$)")
ax2.set_xlim(-720, 720)
ax2.set_xticks(np.arange(-720, 721, 180))
ax2.set_title(r"$\lambda=13/m$ (cage $=13$): SINGLE isolated central spike "
              r"(next spike at $\varphi=\pm2340^\circ$) -- the correct one", fontsize=11)
for k in range(-4, 5):
    ax2.axvline(180 * k, color="0.92", lw=0.7, zorder=0)

cap = (r"Cage diagnosis. Same 13 odd harmonics, same axis. $\lambda=1/m$ (top, "
       r"cage 1) gives a dense periodic spike train; $\lambda=13/m$ (bottom, cage "
       r"13) gives a single isolated central spike over $\pm720^\circ$. Using the "
       r"fundamental phase as the axis collapses cage 13 to 1 -- that was my error.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.8, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.13, hspace=0.2)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_cage_compare.{ext}"), dpi=200)

i0 = np.argmin(np.abs(phi_deg))
print(f"lambda=1/m  : sum at 0 = {sum_1[i0]:.2f}, at 180 = {sum_1[np.argmin(np.abs(phi_deg-180))]:.2f}")
print(f"lambda=13/m : sum at 0 = {sum_13[i0]:.2f}, at 180 = {sum_13[np.argmin(np.abs(phi_deg-180))]:.2f}, "
      f"at 360 = {sum_13[np.argmin(np.abs(phi_deg-360))]:.2f}")
print("saved fig_oddharm_cage_compare.png / .svg")
