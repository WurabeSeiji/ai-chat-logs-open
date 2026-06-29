#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 EXACT (direction-dependent) relativistic Doppler on the combined odd-harmonic
 two-source interference, overlaid on the phase-only baseline.  NO nu-nu0 beat.
================================================================================

Geometry (choice (i): literal far-field two-source, consistent with the prior
doppler_interference_analysis):
    - lambda0 = 1, c = 1, baseline W = 5 (lambda0 units), source separation.
    - baseline d PERPENDICULAR to velocity beta  => path diff = W sin(theta_screen)
      (no Lorentz contraction).
    - velocity beta = 0.6 at angle theta_v = 30 deg from the screen normal z,
      same plane.   gamma = 1/sqrt(1-beta^2) = 1.25.
    - theta_screen = lab-frame observation angle (screen position).
    - x-axis phase  phi = 2 pi W sin(theta_screen) / lambda0   (rest-frame geometric
      phase of the fundamental).  Hence  sin(theta_screen) = phi / (2 pi W).
      For phi in [-360,360] deg this is sin(theta_screen) in [-0.2, 0.2]
      (theta_screen in [-11.54, +11.54] deg): a near-paraxial window, so the
      asymmetry is small but EXACT (not dropped).

Direction-dependent observed wavelength of mode n (rest wavelength lambda0/n):
    lambda_obs,n(theta_screen) = (lambda0/n) * gamma * (1 - beta cos(theta_screen - theta_v))
Two-source phase of mode n at theta_screen:
    dphi_n = (2 pi / lambda_obs,n) * W sin(theta_screen)
           = n * phi / L(theta_screen),    L = gamma (1 - beta cos(theta_screen - theta_v))
Combined two-source field (sources at common offset +/- alpha, alpha = 15 deg):
    psi(phi) = sum_n [cos(dphi_n - alpha) + cos(dphi_n + alpha)]
             = 2 cos(alpha) sum_n cos( n * phi / L(theta_screen(phi)) )
    I(phi)   = psi(phi)^2   (normalized to its own max = 1).

Baseline = same with L == 1 (no Doppler): symmetric, peaks at 0,+/-180,+/-360.
Because L(theta_screen - theta_v) is asymmetric in theta_screen for theta_v != 0,
the EXACT Doppler pattern is left-right ASYMMETRIC: on the approach (+) side
L < center => more compression => more/tighter peaks; on the (-) side fewer/wider.

Author: N. Kihara (peer-review dialogue note), 2026.  EXACT Doppler, beat OFF.
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
THETA_V_DEG = 30.0                 # velocity angle from screen normal
BETA      = 0.6
LAM0      = 1.0
W         = 5.0
PHI_MAX   = 360.0
NPHI      = 1200001                # fine: the +side gets compressed
OUTDIR    = os.path.dirname(os.path.abspath(__file__))

alpha   = np.deg2rad(ALPHA_DEG)
theta_v = np.deg2rad(THETA_V_DEG)
gamma   = 1.0 / np.sqrt(1.0 - BETA**2)

# ------------------------------- geometry ----------------------------------- #
phi_deg = np.linspace(-PHI_MAX, PHI_MAX, NPHI)
phi = np.deg2rad(phi_deg)                                  # rest geometric phase
sin_ts = phi / (2.0 * np.pi * W)                           # = sin(theta_screen)
sin_ts = np.clip(sin_ts, -1.0, 1.0)
theta_s = np.arcsin(sin_ts)                                # observation angle
L = gamma * (1.0 - BETA * np.cos(theta_s - theta_v))       # lambda_obs of fundamental
#   (L = lambda_obs,1 ; lambda_obs,n = L/n ; dphi_n = n*phi/L)

# ------------------------------- fields ------------------------------------- #
def two_source(arg_per_mode):
    """sum_n [cos(n*g - a) + cos(n*g + a)] = 2cos(a) sum_n cos(n*g);
    arg_per_mode is the array g(phi)."""
    out = np.zeros_like(arg_per_mode)
    for n in MODES:
        out += np.cos(n * arg_per_mode - alpha) + np.cos(n * arg_per_mode + alpha)
    return out

g0 = phi                       # baseline: dphi_n = n*phi  (L = 1)
gD = phi / L                   # exact Doppler: dphi_n = n*phi/L(theta_screen)

I0 = two_source(g0) ** 2
I0 /= I0.max()
ID = two_source(gD) ** 2
ID /= ID.max()

# ----------------------- locate main peaks (I > 0.5) ------------------------ #
def main_peaks(y, xdeg, thr=0.5):
    up = (y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:]) & (y[1:-1] > thr)
    idx = np.where(up)[0] + 1
    return [round(float(xdeg[i]), 1) for i in idx]

pk0 = main_peaks(I0, phi_deg)
pkD = main_peaks(ID, phi_deg)

# ------------------------------- figure ------------------------------------- #
fig, ax = plt.subplots(figsize=(13.6, 6.9))
ax.axvspan(-90, 90, color="#2ca02c", alpha=0.10,
           label=r"Central closed-system region ($|\varphi|\leq 90^\circ$)")
ax.plot(phi_deg, I0, color="#1f77b4", lw=1.0,
        label=r"baseline (no Doppler, symmetric)  $I_0(\varphi)$")
ax.plot(phi_deg, ID, color="#d62728", lw=1.0, alpha=0.92,
        label=r"EXACT Doppler ($\beta=0.6,\ \theta_v=30^\circ,\ W=5$)  $I_D(\varphi)$")
ax.plot(pk0, [1.0]*len(pk0), "v", color="#1f77b4", ms=6, zorder=6)
ax.plot(pkD, [1.0]*len(pkD), "v", color="#d62728", ms=6, zorder=6)

ax.set_xlim(-PHI_MAX, PHI_MAX)
ax.set_ylim(-0.02, 1.08)
ax.set_xticks(range(-360, 361, 90))
ax.set_xlabel(r"Phase  $\varphi = 2\pi W \sin\theta_{\rm screen}$  (degrees)  "
              r"[$\sin\theta_{\rm screen}=\varphi/(2\pi W)$]")
ax.set_ylabel(r"Normalized Intensity  $I(\varphi)/\max I$")
ax.grid(alpha=0.25)
ax.legend(loc="upper right", fontsize=8.6, framealpha=0.95)

# L (=lambda_obs of fundamental) at the three reference points
def Lval(phid):
    s = np.deg2rad(phid)/(2*np.pi*W); t = np.arcsin(np.clip(s,-1,1))
    return gamma*(1-BETA*np.cos(t-theta_v))
info = (
    "EXACT direction-dependent Doppler (choice i, far field).\n"
    "lambda0=1, W=5, beta=0.6, gamma=1.25, theta_v=30 deg, baseline _|_ velocity.\n"
    "sin(theta_screen)=phi/(2 pi W) -> theta_screen in [-11.54, +11.54] deg over +/-360.\n"
    "lambda_obs,1 = L(theta_screen) = gamma(1 - beta cos(theta_screen - theta_v)):\n"
    "   phi=-360: L=%.3f   phi=0: L=%.3f   phi=+360: L=%.3f\n"
    "dphi_n = n*phi / L  =>  asymmetric chirp (approach + side more compressed).\n"
    "baseline peaks : %s\n"
    "Doppler  peaks : %s\n"
    "(NO nu-nu0 beat: pattern is static; left-right ASYMMETRIC because theta_v!=0.)"
    % (Lval(-360), Lval(0), Lval(360), pk0, pkD)
)
ax.text(0.012, 0.975, info, transform=ax.transAxes, va="top", ha="left",
        fontsize=7.0, family="monospace",
        bbox=dict(boxstyle="round", fc="#fffbe6", ec="0.7", alpha=0.95))

fig.suptitle("Combined two-source interference (n=1,3,5,7,9): baseline vs "
             "EXACT relativistic Doppler (c=1, tilt 30 deg, no beat) -> "
             "left-right ASYMMETRIC", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
base = os.path.join(OUTDIR, "relphase_combined_doppler_exact")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)

# ------------------------------- report ------------------------------------- #
print("EXACT Doppler:  beta=%.2f gamma=%.4f theta_v=%.0f W=%.1f" %
      (BETA, gamma, THETA_V_DEG, W))
print("L=lambda_obs,1 :  phi=-360 -> %.4f | phi=0 -> %.4f | phi=+360 -> %.4f"
      % (Lval(-360), Lval(0), Lval(360)))
print("baseline peaks (deg):", pk0)
print("Doppler  peaks (deg):", pkD)
# asymmetry quantification: |+ side first peak| vs |- side first peak|
posD = [p for p in pkD if p > 1]; negD = [p for p in pkD if p < -1]
if posD and negD:
    print("first +peak = %.1f deg ; first -peak = %.1f deg  -> asymmetry %.1f deg"
          % (min(posD), max(negD), min(posD) + max(negD)))
print("n(+side peaks)=%d  n(-side peaks)=%d" % (len(posD), len(negD)))
print("Outputs ->", base + ".png / .svg")
