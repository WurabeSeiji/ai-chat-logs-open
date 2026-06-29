#!/usr/bin/env python3
"""
Source waveform only (no slits), matching paper_odd_harmonic_localization:

    S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi),   phi in [-pi/2, pi/2]

Odd harmonic ORDERS n = 1,3,5,...,N (here N=25 -> 13 harmonics).
"13" = NUMBER of odd harmonics = (N+1)/2 = peak height S_N(0). NO cage rescaling.
Closed form: S_N = sin((N+1)phi)/(2 sin phi); zeros exactly at phi = +-pi/2.

Basic domain is the half-wave interval [-pi/2,pi/2] = [-90deg,90deg]:
within one such period S_N is a single localized central peak.

Outputs: fig_oddharm_source_waveform.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

N = 99                                    # highest odd harmonic order
n_list = np.arange(1, N + 1, 2)           # 1,3,5,...,99  -> 50 harmonics
phi_deg = np.linspace(-90.0, 90.0, 18001)  # half-wave interval [-pi/2, pi/2]
phr = np.radians(phi_deg)

waves = [np.cos(n * phr) for n in n_list]
total = np.sum(waves, axis=0)

# verify against closed form sin((N+1)phi)/(2 sin phi)
with np.errstate(divide="ignore", invalid="ignore"):
    closed = np.where(np.abs(np.sin(phr)) < 1e-12,
                      (N + 1) / 2.0,
                      np.sin((N + 1) * phr) / (2.0 * np.sin(phr)))
max_diff = np.max(np.abs(total - closed))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.0, 7.6), sharex=True)

# --- panel 1: the 13 cos(n phi) waves ------------------------------------
cmap = plt.cm.viridis(np.linspace(0, 0.92, len(n_list)))
for w, n, c in zip(waves, n_list, cmap):
    ax1.plot(phi_deg, w, lw=0.9, color=c, alpha=0.85)
ax1.set_ylabel(r"$\cos(n\varphi)$")
ax1.set_ylim(-1.15, 1.15)
ax1.set_title(rf"{len(n_list)} odd-harmonic cos waves  $\cos(n\varphi),\ n=1,3,\dots,{N}$  "
              r"(all $=1$ at $\varphi=0$, all $=0$ at $\varphi=\pm90^\circ$)", fontsize=11)
ax1.axvline(0, color="0.9", lw=0.7, zorder=0)

# --- panel 2: the OBSERVABLE squared amplitude, peak-normalized -----------
peak = total[np.argmin(np.abs(phi_deg))]            # S_N(0) = (N+1)/2
I_hat = (total / peak) ** 2                          # |S_N|^2 normalized to 1.0
ax2.plot(phi_deg, I_hat, color="#1f5fbf", lw=1.8)
ax2.axhline(0, color="0.85", lw=0.8)
ax2.axvline(0, color="0.9", lw=0.7, zorder=0)
ax2.set_ylabel(r"$|S_N(\varphi)|^2$  (normalized to peak 1.0)")
ax2.set_xlabel(r"phase $\varphi$  (degrees;  basic domain $[-90^\circ,90^\circ]=[-\pi/2,\pi/2]$)")
ax2.set_xlim(-90, 90)
ax2.set_ylim(-0.03, 1.08)
ax2.set_xticks(np.arange(-90, 91, 15))
ax2.set_title(r"composite squared amplitude $\widehat{I}_N=|S_N/S_N(0)|^2$ "
              r"(the paper's Fig. 1 observable): single localized central peak, "
              r"zeros at $\varphi=\pm90^\circ$", fontsize=11)

cap = (r"Source waveform only (no slits), per paper_odd_harmonic_localization. "
       rf"$N={N}$, odd orders $n=1,3,\dots,{N}$ ({len(n_list)} harmonics; count$=$peak height "
       rf"$={(N+1)//2}$). "
       r"Basic domain is the half-wave interval $[-\pi/2,\pi/2]$ -- within one period "
       r"$S_N$ is a single localized peak with side-lobes, exactly zero at the ends. "
       r"NO cage rescaling. (direct sum vs closed form max diff "
       rf"$={max_diff:.1e}$).")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.8, wrap=True)
fig.subplots_adjust(left=0.085, right=0.985, top=0.95, bottom=0.14, hspace=0.18)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_source_waveform.{ext}"), dpi=200)

print(f"N={N}, odd orders n={list(n_list)} ({len(n_list)} harmonics)")
print(f"S_N(0) = {total[np.argmin(np.abs(phi_deg))]:.4f}  (expected {(N+1)/2})")
print(f"S_N(+-90) = {total[0]:.4e}, {total[-1]:.4e}  (expected 0)")
print(f"direct-sum vs closed-form max diff = {max_diff:.2e}")
print("saved fig_oddharm_source_waveform.png / .svg")
