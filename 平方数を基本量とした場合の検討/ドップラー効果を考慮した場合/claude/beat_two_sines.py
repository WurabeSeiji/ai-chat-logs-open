#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two sine waves and their exact beat (interference).
  blue : nu = 1     -> lambda = 1        ,  sin(2 pi x / 1)   = sin(2 pi x)
  red  : nu0 = pi/2 -> lambda0 = 2/pi    ,  sin(2 pi x / (2/pi)) = sin(pi^2 x)
  green: exact superposition = blue + red   (the beat / interference waveform)
x-axis = position in wavelength units, 0..5 (0.5 pitch).  y-axis = amplitude.
c = 1, both start at the origin (x=0).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = os.path.dirname(os.path.abspath(__file__))

lam  = 1.0          # nu = 1
lam0 = 2.0/np.pi    # nu0 = pi/2  -> lambda0 = 2/pi

x = np.linspace(0.0, 5.0, 200001)
blue  = np.sin(2*np.pi * x / lam)      # = sin(2 pi x)
red   = np.sin(2*np.pi * x / lam0)     # = sin(pi^2 x)
green = blue + red                      # exact interference / beat

fig, ax = plt.subplots(figsize=(13.0, 5.4))
ax.plot(x, blue,  color="#1f77b4", lw=1.4, label=r"$\nu=1,\ \lambda=1$")
ax.plot(x, red,   color="#d62728", lw=1.4, label=r"$\nu_0=\pi/2,\ \lambda_0=2/\pi$")
ax.plot(x, green, color="#2ca02c", lw=1.8, alpha=0.9,
        label=r"beat (interference) $=\sin(2\pi x)+\sin(\pi^2 x)$")

ax.axhline(0, color="0.6", lw=0.7)
ax.set_xlim(0, 5)
ax.set_ylim(-2.1, 2.1)
ax.set_xticks(np.arange(0, 5.01, 0.5))
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_xlabel(r"position  $x$  (wavelength units, $\lambda=1$)")
ax.set_ylabel("amplitude")
ax.set_title(r"Two sines and their exact beat:  "
             r"$\nu=1$ ($\lambda=1$, blue), $\nu_0=\pi/2$ ($\lambda_0=2/\pi$, red), "
             r"sum (green)")
ax.legend(loc="upper right", fontsize=9, ncol=3)
ax.grid(alpha=0.3)
fig.tight_layout()

base = os.path.join(OUTDIR, "beat_two_sines")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)
print("lambda=%.4f  lambda0=2/pi=%.4f" % (lam, lam0))
print("Outputs ->", base + ".png / .svg")
