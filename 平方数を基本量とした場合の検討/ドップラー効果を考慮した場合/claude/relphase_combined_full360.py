#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Coherent two-source interference of the COMBINED odd-harmonic wave
 (all five odd modes n = 1,3,5,7,9 summed into one wave), full range +/-360 deg.
 Reproduces the relative-phase-contrast reference figure (Fig.2):
   - two sources at a COMMON phase offset alpha = +/-15 deg (phase difference 30 deg)
   - psi_alpha(phi) = sum_n [cos(n phi - alpha) + cos(n phi + alpha)]
                    = 2 cos(alpha) * sum_n cos(n phi)            (n = 1,3,5,7,9)
   - intensity  I(phi) = psi_alpha(phi)^2 ,  normalized to max = 1
 The squared amplitude has period 180 deg, so the main peaks sit at
   phi = -360, -180, 0, +180, +360 deg ; the zeros at +/-90, +/-270 deg bound the
 central closed-system core |phi| <= 90 deg.  Strict left-right symmetry about 0.

 NOTE: this is the phase-only baseline (NO nu-nu0 beat, NO Doppler, NO tilt yet).
 c=1 / tilt theta=30 / Doppler enter in the NEXT step (beat-induced centre jitter).

 Author: N. Kihara (peer-review dialogue note), 2026.  Combined wave, beat OFF.
================================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------- parameters --------------------------------- #
MODES     = [1, 3, 5, 7, 9]        # the five odd modes, combined into ONE wave
ALPHA_DEG = 15.0                   # common phase offset of each source (+/-alpha)
PHI_MAX   = 360.0                  # full range +/- 360 deg
NPHI      = 400001
OUTDIR    = os.path.dirname(os.path.abspath(__file__))

alpha = np.deg2rad(ALPHA_DEG)

# ------------------------------- the wave ----------------------------------- #
def psi_two_source(phi_rad):
    """Two-source coherent sum: sum_n [cos(n phi - a) + cos(n phi + a)]."""
    out = np.zeros_like(phi_rad)
    for n in MODES:
        out += np.cos(n * phi_rad - alpha) + np.cos(n * phi_rad + alpha)
    return out

def psi_factored(phi_rad):
    """Equivalent closed factorization: 2 cos(alpha) * sum_n cos(n phi)."""
    s = np.zeros_like(phi_rad)
    for n in MODES:
        s += np.cos(n * phi_rad)
    return 2.0 * np.cos(alpha) * s

# ------------------------------- compute ------------------------------------ #
phi_deg = np.linspace(-PHI_MAX, PHI_MAX, NPHI)
phi = np.deg2rad(phi_deg)

psi = psi_two_source(phi)
psi_chk = psi_factored(phi)
fact_err = np.max(np.abs(psi - psi_chk))      # verify the 2cos(alpha) factorization

I = psi ** 2
I_norm = I / I.max()                           # normalize to max = 1

# locate the main peaks (should be at -360,-180,0,180,360)
peak_phis = [-360, -180, 0, 180, 360]
peak_vals = [I_norm[np.argmin(np.abs(phi_deg - p))] for p in peak_phis]
# zeros bounding the central core
core_zeros = [-90, 90]
zero_vals = [I_norm[np.argmin(np.abs(phi_deg - z))] for z in core_zeros]

# ------------------------------- figure ------------------------------------- #
fig, ax = plt.subplots(figsize=(13.2, 6.6))

# central closed-system core shading
ax.axvspan(-90, 90, color="#2ca02c", alpha=0.12,
           label=r"Central closed-system region ($|\varphi|\leq 90^\circ$)")

ax.plot(phi_deg, I_norm, color="#1f77b4", lw=1.1,
        label=r"Normalized Intensity $I(\varphi)$")

# mark the five main peaks
ax.plot(peak_phis, peak_vals, "v", color="#d62728", ms=8, zorder=6)
for p in peak_phis:
    ax.axvline(p, color="#d62728", lw=0.6, ls=":", alpha=0.5)

ax.set_xlim(-PHI_MAX, PHI_MAX)
ax.set_ylim(-0.02, 1.06)
ax.set_xticks(range(-360, 361, 90))
ax.set_xlabel(r"Phase  $\varphi$  (degrees)")
ax.set_ylabel(r"Normalized Intensity  $I(\varphi)\,/\,\max I$")
ax.grid(alpha=0.25)
ax.legend(loc="upper right", fontsize=9, framealpha=0.95)

setup = (
    "Setup (formal, phase-only model):\n"
    "- 5 odd modes n=1,3,5,7,9 combined into ONE wave (coherent)\n"
    "- two sources at common offset  alpha = +/-15 deg  (phase diff 30 deg)\n"
    r"- psi(phi) = sum_n [cos(n phi - a) + cos(n phi + a)] = 2 cos(15) sum_n cos(n phi)"
    "\n- intensity  I(phi) = psi(phi)^2 ,  normalized max = 1\n"
    "- I has period 180 deg -> peaks at -360,-180,0,180,360 deg\n"
    "- zeros at +/-90, +/-270 deg bound the central |phi|<=90 core\n"
    "- NO nu-nu0 beat yet (baseline; Doppler/tilt are the next step)"
)
ax.text(0.012, 0.975, setup, transform=ax.transAxes, va="top", ha="left",
        fontsize=7.6, family="monospace",
        bbox=dict(boxstyle="round", fc="#fffbe6", ec="0.7", alpha=0.95))

ax.text(0.0, 1.005, "Central core peak\n(single maximum, symmetric about 0)",
        transform=ax.get_xaxis_transform(), ha="center", va="bottom",
        fontsize=8, color="#2ca02c")

ax.text(0.985, 0.04,
        "Strict left-right symmetry about phi=0\n(each main peak is a single maximum)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

fig.suptitle("Coherent Two-Source Interference  (combined odd modes n=1,3,5,7,9)\n"
             "Sources at +/-15 deg phase offset (phase difference 30 deg)  |  "
             "Full range +/-360 deg", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])

base = os.path.join(OUTDIR, "relphase_combined_full360")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)

# ------------------------------- report ------------------------------------- #
print("Combined odd modes n =", MODES, " alpha = %.0f deg" % ALPHA_DEG)
print("factorization check  max|psi_sum - 2cos(a)*sum cos| = %.2e" % fact_err)
print("main-peak check (should all be 1.0):")
for p, v in zip(peak_phis, peak_vals):
    print("   phi = %+5d deg  ->  I_norm = %.6f" % (p, v))
print("core-boundary zeros (should be ~0):")
for z, v in zip(core_zeros, zero_vals):
    print("   phi = %+5d deg  ->  I_norm = %.3e" % (z, v))
print("Outputs ->", base + ".png / .svg")
