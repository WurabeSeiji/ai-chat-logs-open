#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Combined odd-harmonic two-source interference, full +/-360 deg, with the
 c=1 / tilt theta=30 deg / relativistic-Doppler line OVERLAID on the phase-only
 baseline.  NO nu-nu0 beat (no centre jitter): the Doppler here is the STATIC
 frequency shift of a source/observer relative velocity beta at the fixed tilt.
================================================================================

Baseline (phase-only, blue):
    psi(phi)   = sum_n [cos(n phi - a) + cos(n phi + a)] ,  n = 1,3,5,7,9, a=15deg
    I0(phi)    = psi(phi)^2  (normalized max 1)
    -> period 180 deg, peaks at -360,-180,0,180,360.

Doppler line (c=1, tilt theta=30 deg, beta, red):
    A relative velocity beta at angle theta from the screen normal gives the
    relativistic Doppler factor (at this fixed tilt)
        D = 1 / ( gamma (1 - beta cos theta) ) ,  gamma = 1/sqrt(1-beta^2).
    The observed frequency of every mode is scaled n -> n D, so the OBSERVED
    pattern is the baseline with phi -> D phi:
        psi_D(phi) = sum_n [cos(n D phi - a) + cos(n D phi + a)]
        ID(phi)    = psi_D(phi)^2  (normalized max 1).
    For theta=30 deg the source is approaching (cos30>0) => blueshift D>1 =>
    the whole interference pattern is COMPRESSED in phase: the main peaks move
    from +/-180,+/-360 to +/-180/D, +/-360/D.  (Still left-right SYMMETRIC,
    because at a single fixed tilt D is one constant; the chirp/asymmetry of the
    rigorous treatment appears only when the screen angle varies across the range.)

    NO beat: D is fixed (no nu-nu0 fluctuation), so the central peak stays at
    phi=0 and nothing jitters.

Author: N. Kihara (peer-review dialogue note), 2026.  Doppler overlay, beat OFF.
================================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------- parameters --------------------------------- #
MODES     = [1, 3, 5, 7, 9]
ALPHA_DEG = 15.0
THETA_DEG = 30.0                   # optical-axis tilt
BETA      = 0.6                    # source/observer relative speed (c=1), default
PHI_MAX   = 360.0
NPHI      = 600001
OUTDIR    = os.path.dirname(os.path.abspath(__file__))

alpha = np.deg2rad(ALPHA_DEG)
theta = np.deg2rad(THETA_DEG)
gamma = 1.0 / np.sqrt(1.0 - BETA**2)
D = 1.0 / (gamma * (1.0 - BETA * np.cos(theta)))     # Doppler factor at the tilt

# ------------------------------- waves -------------------------------------- #
def psi(phi_rad, scale=1.0):
    out = np.zeros_like(phi_rad)
    for n in MODES:
        out += np.cos(n * scale * phi_rad - alpha) + np.cos(n * scale * phi_rad + alpha)
    return out

phi_deg = np.linspace(-PHI_MAX, PHI_MAX, NPHI)
phi = np.deg2rad(phi_deg)

I0 = psi(phi, 1.0) ** 2
I0 /= I0.max()
ID = psi(phi, D) ** 2
ID /= ID.max()

base_peaks = [-360, -180, 0, 180, 360]
dop_peaks = [p / D for p in base_peaks]

# ------------------------------- figure ------------------------------------- #
fig, ax = plt.subplots(figsize=(13.4, 6.8))

ax.axvspan(-90, 90, color="#2ca02c", alpha=0.10,
           label=r"Central closed-system region ($|\varphi|\leq 90^\circ$)")

ax.plot(phi_deg, I0, color="#1f77b4", lw=1.1,
        label=r"baseline (phase-only, no Doppler)  $I_0(\varphi)$")
ax.plot(phi_deg, ID, color="#d62728", lw=1.1, alpha=0.9,
        label=r"c=1, tilt $\theta=30^\circ$, Doppler ($\beta=%.1f$, $D=%.3f$)  $I_D(\varphi)$"
              % (BETA, D))

ax.plot(base_peaks, [1]*len(base_peaks), "v", color="#1f77b4", ms=7, zorder=6)
ax.plot(dop_peaks, [1]*len(dop_peaks), "v", color="#d62728", ms=7, zorder=6)

ax.set_xlim(-PHI_MAX, PHI_MAX)
ax.set_ylim(-0.02, 1.07)
ax.set_xticks(range(-360, 361, 90))
ax.set_xlabel(r"Phase  $\varphi$  (degrees)")
ax.set_ylabel(r"Normalized Intensity  $I(\varphi)\,/\,\max I$")
ax.grid(alpha=0.25)
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)

setup = (
    "Combined 5 odd modes n=1,3,5,7,9 ; two sources alpha=+/-15 deg.\n"
    "BLUE  baseline : I0 = [sum_n 2cos(15) cos(n phi)]^2,  peaks +/-180,+/-360.\n"
    "RED   Doppler  : c=1, tilt theta=30 deg, beta=%.1f (gamma=%.3f).\n"
    "  Doppler factor D = 1/[gamma(1-beta cos theta)] = %.3f  (blueshift, approach).\n"
    "  every mode n -> nD  =>  pattern compressed: phi -> D phi.\n"
    "  peaks move to +/-180/D=%.1f, +/-360/D=%.1f deg.\n"
    "  NO nu-nu0 beat: D fixed, centre stays at phi=0 (no jitter).\n"
    "  Single fixed tilt => still left-right symmetric (chirp/asymmetry needs\n"
    "  the screen angle to vary across the range = the rigorous next step)."
    % (BETA, gamma, D, 180/D, 360/D)
)
ax.text(0.012, 0.975, setup, transform=ax.transAxes, va="top", ha="left",
        fontsize=7.3, family="monospace",
        bbox=dict(boxstyle="round", fc="#fffbe6", ec="0.7", alpha=0.95))

fig.suptitle("Coherent two-source interference (combined n=1,3,5,7,9): "
             "phase-only baseline vs c=1 / tilt 30 deg / Doppler (no beat)",
             fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])

base = os.path.join(OUTDIR, "relphase_combined_doppler_overlay")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)

# ------------------------------- report ------------------------------------- #
print("MODES =", MODES, " alpha=%.0f deg  theta=%.0f deg  beta=%.2f" %
      (ALPHA_DEG, THETA_DEG, BETA))
print("gamma = %.4f   Doppler factor D = %.4f   (lambda_obs = 1/D = %.4f)" %
      (gamma, D, 1.0/D))
print("baseline peaks  :", base_peaks)
print("Doppler  peaks  :", [round(p, 1) for p in dop_peaks])
print("Outputs ->", base + ".png / .svg")
