#!/usr/bin/env python3
"""
Odd-harmonic localized-wave two-slit interference, centred source (y=0).
FULL calculation: 13 odd harmonics x 2 slits = 26 waves, summed then squared.
No shortcut, no localization assumed -- it must emerge from the full sum.

Source = equal-amplitude odd harmonics with wavelengths lambda_m = 13/m,
m = 1,3,5,...,25 (cage lambda0 = 13; m=13 gives the base lambda=1).
Each harmonic m, slit k:
    Phi_k^(m) = (2 pi / lambda_m)(r_k - y_slit,k s) = (2 pi m / 13)(r_k - y_slit,k s)
    r_k = sqrt(L^2 + (y - y_slit,k)^2),  y = 0
psi(s) = sum_m [ exp(i Phi1^(m)) + exp(i Phi2^(m)) ],  I = |psi|^2.

Same geometry/axis as Paper 1: L=10, W=5, slits +-2.5,
Phi0 = 2 pi W s / lambda (lambda=1), range +-720 deg.

Outputs: fig_oddharm_M13_centered.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

# ---- geometry (Paper 1) --------------------------------------------------
L = 10.0
W = 5.0
ysl1, ysl2 = +W / 2.0, -W / 2.0
y = 0.0                                   # source fixed at centre
cage = 13.0                               # lambda0 = 13  (M=13)
m_list = np.arange(1, 26, 2)              # 1,3,5,...,25  -> 13 odd harmonics

Phi0_deg = np.linspace(-720.0, 720.0, 16001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W / 1.0)   # Phi0 = 2 pi W s / lambda, lambda=1

r1 = np.sqrt(L**2 + (y - ysl1)**2)
r2 = np.sqrt(L**2 + (y - ysl2)**2)

# ---- full 26-wave coherent sum, then square ------------------------------
psi = np.zeros_like(s, dtype=complex)
for m in m_list:
    kappa = 2.0 * np.pi * m / cage
    psi += np.exp(1j * kappa * (r1 - ysl1 * s))
    psi += np.exp(1j * kappa * (r2 - ysl2 * s))
I = (psi * np.conj(psi)).real
I_norm = I / I.max()

n_waves = 2 * len(m_list)
n_pairs = n_waves * (n_waves - 1) // 2

# ---- plot ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.0, 5.4))
ax.plot(Phi0_deg, I_norm, color="#1f5fbf", lw=1.6,
        label=rf"{len(m_list)} odd harmonics $\times$ 2 slits $=$ {n_waves} waves "
              rf"(full interference, {n_pairs} pairs)")
for k in range(-2, 3):
    ax.axvline(360 * k, color="0.92", lw=0.8, zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)
ax.set_xlim(-720, 720)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-720, 721, 180))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s/\lambda$  (degrees, $\lambda=1$)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
cap = (r"Odd-harmonic localized wave, centred source ($y=0$), full calculation. "
       r"Source = equal-amplitude odd harmonics $\lambda_m=13/m$, $m=1,3,\dots,25$ "
       r"(cage $\lambda_0=13$). Each harmonic through both slits with the Paper-1 "
       r"phase $\Phi_k^{(m)}=(2\pi m/13)(r_k-y_{{\rm slit},k}s)$; the 26 waves are "
       r"summed and squared ($I=|\sum e^{i\Phi}|^2$, 325 interference pairs). "
       r"No localization assumed.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.7, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.18)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_M13_centered.{ext}"), dpi=200)

ic = I_norm[np.argmin(np.abs(Phi0_deg))]
ipk = Phi0_deg[np.argmax(I_norm)]
print(f"waves={n_waves}, pairs={n_pairs}, r1=r2={r1:.5f}")
print(f"I at Phi0=0: {ic:.4f};  global peak at Phi0={ipk:.2f} deg")
print("saved fig_oddharm_M13_centered.png / .svg")
