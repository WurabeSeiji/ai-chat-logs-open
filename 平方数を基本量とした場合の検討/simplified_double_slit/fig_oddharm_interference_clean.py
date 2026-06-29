#!/usr/bin/env python3
"""
Diagnose the 23-peak scramble. Two-slit interference of the localized odd-harmonic
wave (N=99), with vs without the source-to-slit propagation phase (2 pi n/13) r_k.

Top  : WITH source phase (finite L, slits between spikes) -> scrambled multi-band.
Bottom: WITHOUT source phase = standard far-field (Fraunhofer) two slit.
        Each harmonic n gives 2 cos(n theta_scr), theta_scr = pi W s / 13, so
            psi = sum_n 2 cos(n theta_scr) = 2 S_N(theta_scr),  I = 4 S_N^2.
        -> a SINGLE central localized fringe = the squared isolated wave.
The bottom panel is overlaid with the analytic 4 S_N^2 to confirm the identity.

lambda_n = 13/n, n=1,3,...,99 (50 harmonics); W=5, L=10; reference n=13 (lambda=1).
Phi0 = 2 pi W s, range +-720 deg;  theta_scr = Phi0/26.

Outputs: fig_oddharm_interference_clean.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

L, W = 10.0, 5.0
ysl1, ysl2 = +W / 2.0, -W / 2.0
y = 0.0
cage = 13.0
N = 99
n_list = np.arange(1, N + 1, 2)               # 50 harmonics

Phi0_deg = np.linspace(-720.0, 720.0, 60001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W)  # Phi0 = 2 pi W s (lambda=1)
r1 = np.sqrt(L**2 + (y - ysl1)**2)
r2 = np.sqrt(L**2 + (y - ysl2)**2)

psi_src = np.zeros_like(s, dtype=complex)     # WITH source phase
psi_far = np.zeros_like(s, dtype=complex)     # WITHOUT (far field)
for n in n_list:
    kappa = 2.0 * np.pi * n / cage
    psi_src += np.exp(1j * kappa * (r1 - ysl1 * s)) + np.exp(1j * kappa * (r2 - ysl2 * s))
    psi_far += np.exp(1j * kappa * (-ysl1 * s)) + np.exp(1j * kappa * (-ysl2 * s))
I_src = (psi_src * np.conj(psi_src)).real
I_far = (psi_far * np.conj(psi_far)).real
I_src_n = I_src / I_src.max()
I_far_n = I_far / I_far.max()

# analytic check: I_far = 4 S_N(theta)^2,  theta = pi W s / 13 = Phi0/26
theta = np.pi * W * s / cage
with np.errstate(divide="ignore", invalid="ignore"):
    S = np.where(np.abs(np.sin(theta)) < 1e-12, (N + 1) / 2.0,
                 np.sin((N + 1) * theta) / (2.0 * np.sin(theta)))
I_analytic_n = (S / S.max())**2
max_diff = np.max(np.abs(I_far_n - I_analytic_n))

def count_peaks(arr, thr=0.5):
    loc = (arr[1:-1] > arr[:-2]) & (arr[1:-1] > arr[2:]) & (arr[1:-1] >= thr)
    return np.where(loc)[0] + 1

p_src = count_peaks(I_src_n)
p_far = count_peaks(I_far_n)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13.5, 8.0), sharex=True)

ax1.plot(Phi0_deg, I_src_n, color="#d11f2d", lw=1.2)
ax1.plot(Phi0_deg[p_src], I_src_n[p_src], 'o', color="k", ms=3.5,
         label=rf"{len(p_src)} peaks $\geq0.5$")
ax1.set_ylabel(r"$|\psi|^2$ (norm.)")
ax1.set_ylim(-0.03, 1.08)
ax1.set_title(r"WITH source$\rightarrow$slit phase $(2\pi n/13)\,r_k$ (finite $L$, slits "
              r"between spikes): scrambled multi-band -- this was my error", fontsize=11)
ax1.legend(loc="upper right", fontsize=9)
for k in range(-2, 3):
    ax1.axvline(360 * k, color="0.92", lw=0.8, zorder=0)

ax2.plot(Phi0_deg, I_far_n, color="#1f5fbf", lw=1.6, zorder=3,
         label=rf"far-field two slit ({len(p_far)} peak$\geq0.5$)")
ax2.plot(Phi0_deg, I_analytic_n, color="#f0a500", lw=1.0, ls="--", zorder=4,
         label=r"analytic $4\,S_N(\theta)^2$ (norm.)")
ax2.set_ylabel(r"$|\psi|^2$ (norm.)")
ax2.set_xlabel(r"Two-slit reference phase $\Phi_0=2\pi W s$ (deg; $\lambda=1=n{=}13$; "
               r"$\theta_{\rm scr}=\Phi_0/26$)")
ax2.set_xlim(-720, 720)
ax2.set_ylim(-0.03, 1.08)
ax2.set_xticks(np.arange(-720, 721, 180))
ax2.set_title(r"WITHOUT source phase = far-field: $\psi=2S_N(\theta_{\rm scr})$, "
              r"$I=4S_N^2$ -- a SINGLE central localized fringe (matches intuition)",
              fontsize=11)
ax2.legend(loc="upper right", fontsize=9)
for k in range(-2, 3):
    ax2.axvline(360 * k, color="0.92", lw=0.8, zorder=0)

cap = (r"Diagnosis of the multi-band scramble. Same 50 odd harmonics, $W=5$. "
       r"Top: including the per-harmonic source$\rightarrow$slit phase scrambles the "
       r"pattern into many peaks. Bottom: the standard far-field two-slit interference "
       r"of the localized wave is exactly $4\,S_N(\theta_{\rm scr})^2$ "
       rf"(analytic overlay, max diff $={max_diff:.1e}$): one central localized fringe.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.6, wrap=True)
fig.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.13, hspace=0.22)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_oddharm_interference_clean.{ext}"), dpi=200)

print(f"WITH source phase : I(0)={I_src_n[np.argmin(np.abs(Phi0_deg))]:.4f}, peaks>=0.5: {len(p_src)}")
print(f"far-field (clean) : I(0)={I_far_n[np.argmin(np.abs(Phi0_deg))]:.4f}, peaks>=0.5: {len(p_far)}")
print(f"far-field vs analytic 4 S_N^2 : max diff = {max_diff:.2e}")
print("saved fig_oddharm_interference_clean.png / .svg")
