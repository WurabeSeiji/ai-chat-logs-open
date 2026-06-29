#!/usr/bin/env python3
"""
EXACT two-slit interference of the localized odd-harmonic wave, nothing dropped.
Base n=27 (lambda=1), cage=27, L=270, W=135; source centred at y=0.
All 50 odd harmonics n=1,3,...,99 (lambda_n = 27/n).

Exact source->slit distance (transverse slit offset kept):
    r_k = sqrt(L^2 + (y - y_slit,k)^2) = sqrt(270^2 + 67.5^2) = 278.31
Phase of harmonic n through slit k:
    Phi_k^(n) = (2 pi n / 27)(r_k - y_slit,k s)
psi(s) = sum_n [exp(i Phi1^(n)) + exp(i Phi2^(n))],   I = |psi|^2.
Reference phase Phi0 = 2 pi W s (lambda=1 = n=27); range +-720 deg.

No approximation: L, the transverse slit offset W/2, and the per-harmonic source
phase are all kept. This is the honest exact result.

Outputs: fig_oddharm_interference_L270.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

cage = 27.0
L, W = 270.0, 135.0
ysl1, ysl2 = +W / 2.0, -W / 2.0
y = 0.0
N = 99
n_list = np.arange(1, N + 1, 2)               # 50 odd harmonics

Phi0_deg = np.linspace(-720.0, 720.0, 60001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W)  # Phi0 = 2 pi W s, lambda=1

r1 = np.sqrt(L**2 + (y - ysl1)**2)
r2 = np.sqrt(L**2 + (y - ysl2)**2)

psi = np.zeros_like(s, dtype=complex)
for n in n_list:
    k = 2.0 * np.pi * n / cage
    psi += np.exp(1j * k * (r1 - ysl1 * s)) + np.exp(1j * k * (r2 - ysl2 * s))
I = (psi * np.conj(psi)).real
I_norm = I / I.max()

n_waves = 2 * len(n_list)
n_pairs = n_waves * (n_waves - 1) // 2
loc = (I_norm[1:-1] > I_norm[:-2]) & (I_norm[1:-1] > I_norm[2:]) & (I_norm[1:-1] >= 0.5)
n_peaks = int(loc.sum())
ic = I_norm[np.argmin(np.abs(Phi0_deg))]

fig, ax = plt.subplots(figsize=(13.5, 5.6))
ax.plot(Phi0_deg, I_norm, color="#d11f2d", lw=1.4,
        label=rf"{len(n_list)} odd harmonics $\times$ 2 slits $=$ {n_waves} waves "
              rf"(exact, {n_pairs} pairs)")
for kk in range(-2, 3):
    ax.axvline(360 * kk, color="0.92", lw=0.8, zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)
ax.set_xlim(-720, 720)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-720, 721, 180))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s$  (degrees; $\lambda=1=n{=}27$)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
cap = (r"EXACT two-slit interference of the localized odd-harmonic wave: base $n=27$, "
       r"cage$=27$, $L=270$, $W=135$, $\lambda_n=27/n$, 50 harmonics. "
       r"$\Phi_k^{(n)}=(2\pi n/27)(r_k-y_{{\rm slit},k}s)$ with $r_k=\sqrt{L^2+(W/2)^2}=278.31$; "
       r"$L$, the slit offset $W/2$, and the per-harmonic source phase are all kept "
       r"(nothing dropped). " + rf"$I(0)={ic:.3f}$, {n_peaks} peaks $\geq0.5$.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.6, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.20)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_interference_L270.{ext}"), dpi=200)

print(f"EXACT: base n=27, cage=27, L=270, W=135, r_k={r1:.4f}, r_k/cage={r1/cage:.4f}")
print(f"I(0)={ic:.4f}, peaks>=0.5: {n_peaks}")
print("saved fig_oddharm_interference_L270.png / .svg")
