#!/usr/bin/env python3
"""
Same as fig_oddharm_interference_N99 but with the horizontal axis widened to
Phi0 in [-9360, 9360] deg  (theta_scr = Phi0/26 in [-360 deg, 360 deg]).

|S_N|^2 has period 180 deg in theta_scr, so I = 4 S_N^2 repeats every
180 deg * 26 = 4680 deg in Phi0:
    repeat peaks at Phi0 = 0, +-4680, +-9360 deg  (theta_scr = 0, +-180, +-360)
    zeros        at Phi0 = +-2340, +-7020 deg     (theta_scr = +-90, +-270)
Over this wide window the single localized fringe is seen repeating as a periodic
train of isolated fringes.

Source = equal-amplitude odd harmonics lambda_n = 13/n, n=1,3,...,99; W=5.
psi = sum_n 2 cos(n theta_scr) = 2 S_N(theta_scr),  I = 4 S_N^2.

Outputs: fig_oddharm_interference_N99_wide.png / .svg
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

Phi0_deg = np.linspace(-9360.0, 9360.0, 200001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W)  # Phi0 = 2 pi W s, lambda=1

# ---- full coherent sum over all harmonics x both slits, then square ------
psi = np.zeros_like(s, dtype=complex)
for n in n_list:
    kappa = 2.0 * np.pi * n / cage
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

# ---- plot ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.5, 5.6))
ax.plot(Phi0_deg, I_norm, color="#1f5fbf", lw=1.1, zorder=3,
        label=rf"{len(n_list)} odd harmonics $\times$ 2 slits $=$ {n_waves} waves "
              rf"(full sum, {n_pairs} pairs)")
ax.plot(Phi0_deg, I_analytic, color="#f0a500", lw=0.9, ls="--", zorder=4,
        label=r"analytic $4\,S_N(\theta_{\rm scr})^2$")
for x0 in (-9360, -4680, 0, 4680, 9360):
    ax.axvline(x0, color="0.85", lw=0.9, zorder=0)
for x0 in (-7020, -2340, 2340, 7020):
    ax.axvline(x0, color="0.93", lw=0.7, ls=":", zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)
ax.set_xlim(-9360, 9360)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-9360, 9361, 2340))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s$  "
              r"(degrees; $\lambda=1=$ the $n{=}13$ wave; $\theta_{\rm scr}=\Phi_0/26$)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
cap = (r"Two-slit far-field interference of the localized odd-harmonic wave, $N=99$, "
       r"widened to $\Phi_0\in[-9360^\circ,9360^\circ]$ ($\theta_{\rm scr}\in[-360^\circ,360^\circ]$). "
       r"$I=4S_N^2$ has period $4680^\circ$ in $\Phi_0$: a periodic train of isolated "
       r"fringes (peaks at $\Phi_0=0,\pm4680^\circ,\pm9360^\circ$; zeros at "
       rf"$\pm2340^\circ,\pm7020^\circ$). Full sum vs closed form: max diff $={max_diff:.1e}$.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.6, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.20)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_interference_N99_wide.{ext}"), dpi=200)

ic = I_norm[np.argmin(np.abs(Phi0_deg))]
i4680 = I_norm[np.argmin(np.abs(Phi0_deg - 4680))]
i2340 = I_norm[np.argmin(np.abs(Phi0_deg - 2340))]
print(f"range +-9360 deg; I(0)={ic:.4f}, I(4680)={i4680:.4f}, I(2340)={i2340:.2e}")
print(f"full sum vs analytic 4 S_N^2 : max diff = {max_diff:.2e}")
print("saved fig_oddharm_interference_N99_wide.png / .svg")
