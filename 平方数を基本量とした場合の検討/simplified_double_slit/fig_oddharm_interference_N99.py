#!/usr/bin/env python3
"""
Two-slit far-field interference of the localized odd-harmonic wave (N=99).

Source = equal-amplitude odd harmonics, wavelength lambda_n = 13/n
(n=13 -> lambda=1 = reference; fundamental n=1 -> lambda=13). Slit separation W=5.

In the far field each harmonic n gives the two-slit amplitude
    2 cos( (2 pi n / 13)(W/2) s ) = 2 cos( n theta_scr ),   theta_scr = pi W s / 13,
so the coherent sum over the odd harmonics reproduces the localized wave itself:
    psi(s) = sum_n 2 cos(n theta_scr) = 2 S_N(theta_scr),
    I(s)   = |psi|^2 = 4 S_N(theta_scr)^2,   S_N = sin((N+1)theta)/(2 sin theta).
The screen interference IS the squared isolated wave: a single central localized
fringe. Reference phase Phi0 = 2 pi W s (lambda=1 = n=13), theta_scr = Phi0/26.

All 50 odd harmonics n=1,3,...,99 are summed explicitly (none omitted); the result
is verified against the closed form 4 S_N^2 to machine precision.

Outputs: fig_oddharm_interference_N99.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

W = 5.0
ysl1, ysl2 = +W / 2.0, -W / 2.0
cage = 13.0
N = 99
n_list = np.arange(1, N + 1, 2)               # 1,3,...,99 -> 50 harmonics

Phi0_deg = np.linspace(-720.0, 720.0, 60001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W)  # Phi0 = 2 pi W s, lambda=1

# ---- full coherent sum over all harmonics x both slits, then square ------
psi = np.zeros_like(s, dtype=complex)
for n in n_list:
    kappa = 2.0 * np.pi * n / cage            # 2 pi / lambda_n
    psi += np.exp(1j * kappa * (-ysl1 * s)) + np.exp(1j * kappa * (-ysl2 * s))
I = (psi * np.conj(psi)).real
I_norm = I / I.max()

# ---- analytic check: I = 4 S_N(theta)^2,  theta = pi W s / 13 = Phi0/26 ---
theta = np.pi * W * s / cage
with np.errstate(divide="ignore", invalid="ignore"):
    S = np.where(np.abs(np.sin(theta)) < 1e-12, (N + 1) / 2.0,
                 np.sin((N + 1) * theta) / (2.0 * np.sin(theta)))
I_analytic = (S / S.max())**2
max_diff = np.max(np.abs(I_norm - I_analytic))

n_waves = 2 * len(n_list)
n_pairs = n_waves * (n_waves - 1) // 2
loc = (I_norm[1:-1] > I_norm[:-2]) & (I_norm[1:-1] > I_norm[2:]) & (I_norm[1:-1] >= 0.5)
n_peaks = int(loc.sum())

# ---- plot ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.5, 5.6))
ax.plot(Phi0_deg, I_norm, color="#1f5fbf", lw=1.8, zorder=3,
        label=rf"{len(n_list)} odd harmonics $\times$ 2 slits $=$ {n_waves} waves "
              rf"(full sum, {n_pairs} pairs)")
ax.plot(Phi0_deg, I_analytic, color="#f0a500", lw=1.0, ls="--", zorder=4,
        label=r"analytic $4\,S_N(\theta_{\rm scr})^2$")
for k in range(-2, 3):
    ax.axvline(360 * k, color="0.92", lw=0.8, zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)
ax.set_xlim(-720, 720)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-720, 721, 180))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s$  "
              r"(degrees; $\lambda=1=$ the $n{=}13$ wave; $\theta_{\rm scr}=\Phi_0/26$)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
cap = (r"Two-slit far-field interference of the localized odd-harmonic wave, $N=99$ "
       r"(all 50 harmonics $n=1,3,\dots,99$, $\lambda_n=13/n$, $W=5$). Each harmonic "
       r"gives $2\cos(n\theta_{\rm scr})$, so the coherent sum is $2S_N(\theta_{\rm scr})$ "
       r"and the screen intensity $I=4S_N(\theta_{\rm scr})^2$: a single central localized "
       rf"fringe. Full sum vs closed form: max diff $={max_diff:.1e}$.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.6, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.20)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_interference_N99.{ext}"), dpi=200)

ic = I_norm[np.argmin(np.abs(Phi0_deg))]
print(f"harmonics: n=1..{N} ({len(n_list)}), waves={n_waves}, pairs={n_pairs}")
print(f"I at Phi0=0: {ic:.4f};  peaks >= 0.5: {n_peaks}")
print(f"full sum vs analytic 4 S_N^2 : max diff = {max_diff:.2e}")
print("saved fig_oddharm_interference_N99.png / .svg")
