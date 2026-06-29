#!/usr/bin/env python3
"""
Source-position read-out (STATIC source, pure geometry, NO Doppler).

Each source position x_i over +/- lambda/2 is drawn separately:
  - its exact far-field two-slit intensity (identical shape, only SHIFTED by the
    exact source-side path difference), scaled so the peak height equals the
    weight cos^2(pi x_i / 2);
  - GREEN for x_i != 0, BLUE for x_i = 0 (the centred-source fringe);
  - YELLOW cos^2 envelope cos^2(pi Phi0/180) over +/-90 deg, through which the
    fringe peaks run.

No moving source, no Doppler: wavelength is lambda0 everywhere; the photon
propagates at c along the geometric paths r_k. The source position only shifts
the fringe (geometric path difference). The cos^2 envelope is the read-out
(transcription) of the input source-position distribution, not a derivation.

Outputs: fig_decomposition_static.png / .svg
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.size": 12, "mathtext.fontset": "cm",
})

# ---- geometry ------------------------------------------------------------
L = 10.0
W = 5.0
lam0 = 1.0
half = lam0 / 2.0
ysl1, ysl2 = +W / 2.0, -W / 2.0

x_src = np.array([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
y_src = x_src * half
w_int = np.cos(np.pi * x_src / 2.0)**2        # peak height = intensity weight

Phi0_deg = np.linspace(-720.0, 720.0, 16001)
s = np.radians(Phi0_deg) / (2.0 * np.pi * W / lam0)

def unit_intensity(y_i):
    """Exact unit-peak STATIC two-slit intensity for a point source at y_i.
    Pure geometry: wavelength lambda0 on both legs, no Doppler factor."""
    r1 = np.sqrt(L**2 + (y_i - ysl1)**2)
    r2 = np.sqrt(L**2 + (y_i - ysl2)**2)
    Phi1 = 2.0 * np.pi * (r1 - ysl1 * s) / lam0
    Phi2 = 2.0 * np.pi * (r2 - ysl2 * s) / lam0
    psi = np.exp(1j * Phi1) + np.exp(1j * Phi2)
    return (psi * np.conj(psi)).real / 4.0     # peak 1

# ---- figure --------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.0, 5.8))

env_x = np.linspace(-90.0, 90.0, 600)
env_y = np.cos(np.pi * env_x / 180.0)**2
ax.plot(env_x, env_y, color="#f0a500", lw=3.2, zorder=5,
        label=r"$\cos^2$ envelope = read-out of source distribution")

print("  x_i     peak Phi0 (exact)   peak height   cos^2 at that Phi0   dev")
peaks_x, peaks_y = [], []
for x_i, y_i, w_i in zip(x_src, y_src, w_int):
    curve = w_i * unit_intensity(y_i)
    if x_i == 0.0:
        ax.plot(Phi0_deg, curve, color="#1f5fbf", lw=2.3, zorder=4,
                label=r"$x=0$  (centred source, peak $1$)")
    else:
        ax.plot(Phi0_deg, curve, color="#2ca02c", lw=1.0, alpha=0.85, zorder=2)
    m = np.abs(Phi0_deg) < 180.0
    idx = np.argmax(curve[m])
    px, py = Phi0_deg[m][idx], curve[m][idx]
    peaks_x.append(px); peaks_y.append(py)
    env_at = np.cos(np.pi * px / 180.0)**2
    print(f"{x_i:+.1f}   {px:+8.2f}            {py:.4f}        {env_at:.4f}"
          f"            {py-env_at:+.4f}")

ax.plot(peaks_x, peaks_y, 'o', color="#d11f2d", ms=5.5, zorder=6,
        label="fringe peaks (exact geometry)")
ax.plot([], [], color="#2ca02c", lw=1.2,
        label=r"each source $x_i$ separate (green, peak $\cos^2(\pi x_i/2)$)")

for k in range(-2, 3):
    ax.axvline(360 * k, color="0.92", lw=0.8, zorder=0)
ax.axhline(0, color="0.92", lw=0.8, zorder=0)

ax.set_xlim(-720, 720)
ax.set_ylim(-0.03, 1.08)
ax.set_xticks(np.arange(-720, 721, 180))
ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s/\lambda_0$  (degrees)")
ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
ax.legend(loc="upper right", fontsize=9.0, framealpha=0.96)

cap = (r"Static source, pure geometry (no Doppler). Each source point $x_i$ over "
       r"$\pm\lambda/2$ gives its own exact far-field two-slit intensity -- "
       r"identical in shape, only shifted by the geometric source-side path "
       r"difference -- scaled to peak height $\cos^2(\pi x_i/2)$ (green; centre "
       r"$x=0$ in blue). The fringe peaks trace the $\cos^2$ envelope (yellow), "
       r"which is the read-out (transcription) of the input source-position "
       r"distribution. Small departures are the exact (non-paraxial) path "
       r"geometry.")
fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.7, wrap=True)

fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.17)

outdir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg"):
    fig.savefig(os.path.join(outdir, f"fig_decomposition_static.{ext}"), dpi=200)
print("\nsaved fig_decomposition_static.png / .svg")
