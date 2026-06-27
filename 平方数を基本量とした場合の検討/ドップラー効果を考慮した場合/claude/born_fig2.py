#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 2 : the localization kernel S_N and its universal collapse.
Generates born_localization_kernel.png / .svg  (2 panels):
 (left)  S_N(phi) normalized by its peak (N+1)/2, vs phase in % of half-wavelength,
         for N = 9, 99, 999 : the main peak narrows as ~ 1/(N+1).
 (right) the same curves vs the scaled variable u = (N+1) phi : they COLLAPSE onto
         the universal form sin(u)/u (sinc), the truncated-reproducing-kernel limit.
Consistent with the starting paper (odd-harmonic isolated-peak wave, Zenodo v0.4).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def S_N(x, N):
    x = np.asarray(x, float)
    s = np.sin(x)
    return np.where(np.abs(s) < 1e-12, (N+1)/2.0, np.sin((N+1)*x)/(2.0*s))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.0))
colors = {9: "#d62728", 99: "#1f77b4", 999: "#2ca02c"}

# (left) normalized kernel vs phase (% of half-wavelength): narrowing ~ 1/(N+1)
phi = np.linspace(-np.pi/2, np.pi/2, 400001)
xpct = phi/np.pi*100
for N in (9, 99, 999):
    Shat = S_N(phi, N)/((N+1)/2.0)
    axL.plot(xpct, Shat, color=colors[N], lw=1.4, label=f"S_N/peak,  N={N}  (width ~ 1/{N+1})")
axL.set_xlim(-12, 12)
axL.set_title("(left) localized kernel narrows as ~ 1/(N+1)", fontsize=10)
axL.set_xlabel("phase φ  (% of half-wavelength)")
axL.set_ylabel("S_N(φ) / peak,   peak = (N+1)/2")
axL.axhline(0, color="0.6", lw=0.7); axL.grid(alpha=0.25)
axL.legend(fontsize=8.5, loc="upper right")

# (right) vs scaled u = (N+1) phi : universal collapse onto sin(u)/u
for N in (9, 99, 999):
    u = (N+1)*phi
    Shat = S_N(phi, N)/((N+1)/2.0)
    m = np.abs(u) <= 4*np.pi
    axR.plot(u[m], Shat[m], color=colors[N], lw=1.3, alpha=0.85, label=f"N={N}")
uu = np.linspace(-4*np.pi, 4*np.pi, 2000)
sinc = np.where(np.abs(uu) < 1e-9, 1.0, np.sin(uu)/uu)
axR.plot(uu, sinc, color="k", lw=1.0, ls="--", label="sin(u)/u  (universal limit)")
axR.set_title("(right) universal collapse onto sin(u)/u,  u = (N+1)φ", fontsize=10)
axR.set_xlabel("scaled phase  u = (N+1) φ")
axR.set_ylabel("S_N / peak")
axR.axhline(0, color="0.6", lw=0.7); axR.grid(alpha=0.25)
axR.legend(fontsize=8.5, loc="upper right")

fig.suptitle("The localized kernel: equal-amplitude odd-harmonic isolated-peak wave "
             "S_N(φ) = sin((N+1)φ)/(2 sinφ)", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig("born_localization_kernel.png", dpi=150, bbox_inches="tight")
fig.savefig("born_localization_kernel.svg", bbox_inches="tight")
plt.close(fig)
print("saved born_localization_kernel.png/.svg")
