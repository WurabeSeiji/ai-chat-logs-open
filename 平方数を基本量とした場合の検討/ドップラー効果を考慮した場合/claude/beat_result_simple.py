#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The beat-conjecture result is simply the SQUARED AMPLITUDE |psi_base|^2.
The beat (nu vs nu0=pi/2) is only the mechanism that scans the observer phase
uniformly (irrational ratio 2/pi -> Weyl); it does NOT change the result.

  base wave cos(phi)        -> observation distribution = cos^2(phi).
  base wave S_N(phi)        -> observation distribution = |S_N(phi)|^2 (localized).

x-axis = phase, y-axis = squared amplitude.  That is all.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = os.path.dirname(os.path.abspath(__file__))
phi_deg = np.linspace(-180, 180, 4001)
phi = np.deg2rad(phi_deg)

# --- simplest base: a single cosine ---
amp = np.cos(phi)
I_cos = np.cos(phi) ** 2                       # = beat-conjecture result

# --- localized base: combined odd harmonics n=1,3,5,7,9 ---
S = sum(np.cos(n * phi) for n in [1, 3, 5, 7, 9])
I_S = (S / 5.0) ** 2                            # normalized |S_N|^2

fig, ax = plt.subplots(figsize=(10.5, 5.6))
ax.plot(phi_deg, amp, color="0.65", lw=1.2, ls="--",
        label=r"base amplitude  $\cos\varphi$")
ax.plot(phi_deg, I_cos, color="#1f77b4", lw=2.2,
        label=r"beat result $=|\psi_{\rm base}|^2=\cos^2\varphi$")
ax.plot(phi_deg, I_S, color="#d62728", lw=1.6, alpha=0.85,
        label=r"localized base: $|S_N(\varphi)|^2$, $N=9$ (for reference)")
ax.axhline(0, color="0.7", lw=0.7)
ax.set_xlim(-180, 180)
ax.set_xticks(range(-180, 181, 45))
ax.set_xlabel(r"phase  $\varphi$  (degrees)")
ax.set_ylabel(r"squared amplitude  $|\psi_{\rm base}|^2$")
ax.set_title(r"Beat-conjecture result = squared amplitude. "
             r"$\nu=1,\ \nu_0=\pi/2$ only makes the scan uniform; "
             r"the output is just $|\psi_{\rm base}|^2$.")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.25)
fig.tight_layout()
base = os.path.join(OUTDIR, "beat_result_simple")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)
print("ratio nu/nu0 = 2/pi =", 1.0/(np.pi/2), "(irrational -> uniform scan)")
print("Outputs ->", base + ".png / .svg")
