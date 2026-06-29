#!/usr/bin/env python3
"""
Paper 2, Fig 2: the localized odd-harmonic source wave S_N (isolated peak) for
N=17, on the half-wave interval [-pi/2, pi/2]. Shows both the signed amplitude
S_N and the normalized squared amplitude |S_N|^2 (the localized read-out量).

S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi) = sin((N+1)phi)/(2 sin phi).

Output: fig_paper2_localized_wave_N17.png / .svg
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

N = 17
n_list = np.arange(1, N + 1, 2)            # 1,3,...,17 -> 9 harmonics
phi_deg = np.linspace(-90.0, 90.0, 8001)
phr = np.radians(phi_deg)
S = np.sum([np.cos(n * phr) for n in n_list], axis=0)
I = (S / S.max())**2                        # normalized |S_N|^2

fig, ax = plt.subplots(figsize=(11.0, 4.8))
ax.plot(phi_deg, S / S.max(), color="#9bbcdd", lw=1.2,
        label=r"amplitude $S_N/S_N(0)$")
ax.plot(phi_deg, I, color="#1f5fbf", lw=2.2,
        label=r"squared $|S_N|^2$ (normalised)")
ax.axhline(0, color="0.85", lw=0.8); ax.axvline(0, color="0.92", lw=0.7)
ax.set_xlim(-90, 90); ax.set_ylim(-0.35, 1.05)
ax.set_xticks(np.arange(-90, 91, 30))
ax.set_xlabel(r"phase $\varphi$ (deg; half-wave interval $[-90^\circ,90^\circ]=[-\pi/2,\pi/2]$)")
ax.set_ylabel("amplitude / intensity")
ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
cap = (rf"Localized odd-harmonic source wave $S_N(\varphi)=\sum_{{m}}\cos((2m{{+}}1)\varphi)$, "
       rf"$N={N}$ ({len(n_list)} harmonics). Central peak $S_N(0)=(N{{+}}1)/2$, zeros at "
       r"$\pm90^\circ$; the squared amplitude is a single localized peak (width $\sim1/(N{+}1)$).")
fig.text(0.5, 0.01, cap, ha="center", va="bottom", fontsize=8.6, wrap=True)
fig.subplots_adjust(left=0.08, right=0.985, top=0.96, bottom=0.22)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_paper2_localized_wave_N17.{ext}"), dpi=200)
print(f"S_N(0)={S.max():.3f} (expected {(N+1)/2}); saved fig_paper2_localized_wave_N17.png/.svg")
