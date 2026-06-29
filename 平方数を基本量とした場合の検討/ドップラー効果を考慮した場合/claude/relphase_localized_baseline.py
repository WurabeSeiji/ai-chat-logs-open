#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Baseline localized peak wave of the relative-phase-contrast model
 (c = 1, optical-axis tilt theta = 30 deg, far field L -> inf,
  NO observer-system (nu, nu0) beat / Doppler fluctuation)
================================================================================

This is STEP 1 of connecting the odd-harmonic localized wave to a real two-source
(two-slit) interference geometry.  Here we compute only the BASELINE pattern:
the state with NO nu-nu0 beat (Doppler) fluctuation.  In this baseline the
inter-source relative phase is purely geometric and symmetric, so the pattern is
the localized peak wave centered at phase phi = 0.

Model (matching paper_relative_phase_contrast_ja_v0_1):
    Each source emits the constant-amplitude odd-harmonic sum on the half-
    wavelength phase interval phi in [-pi/2, pi/2]:

        S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi)
                 = sin((N+1) phi) / (2 sin phi)        (closed form; (N+1)/2 at 0)

    Observable of that paper = the (non-negative) SQUARED amplitude
        I_N(phi) = |S_N(phi)|^2 .
    Figure 1 of that paper shows the amplitude S_N; Figure 2 shows I_N.
    We therefore plot BOTH (top: normalized amplitude, bottom: normalized I_N).

Geometric frame (documented; it does NOT change the baseline shape, it sets the
scale for the NEXT step where the nu-nu0 beat is switched on):
    - c = 1 (units).
    - far field: screen distance L >> baseline W, so the two-source pattern is
      the clean Fraunhofer interference and the position <-> phase map is linear.
    - optical-axis tilt theta = 30 deg.  Its role appears only once the beat is
      added: a temporal beat phase phi_beat shifts the fringe CENTER by
          Delta X_center = (w / 2pi) * sin(theta) * phi_beat ,   w = lambda L / W,
      i.e. the transverse projection factor is sin(theta) = sin 30 deg = 0.5.
      At theta = 0 the projection vanishes (no center shift); see the c=1 endpoint
      degeneracy (forward whiteout / backward blackout).  Here, with the beat OFF,
      the center sits at phi = 0 and nothing jitters yet.

Author: N. Kihara (peer-review dialogue note), 2026.  Baseline, beat OFF.
================================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------- parameters --------------------------------- #
C        = 1.0                     # speed of light (units)
THETA_DEG = 30.0                   # optical-axis tilt (deg); sin(theta)=0.5
W        = 5.0                     # source separation (lambda = 1 units)
NLIST    = [1, 3, 5, 7, 9]         # highest odd-harmonic order N (odd)
NPHI     = 200001                  # phase samples on [-pi/2, pi/2]
OUTDIR   = os.path.dirname(os.path.abspath(__file__))

SIN_THETA = np.sin(np.deg2rad(THETA_DEG))   # = 0.5 for 30 deg

# ------------------------------- the wave ----------------------------------- #
def S_direct(phi, N):
    """Direct constant-amplitude odd-harmonic sum S_N(phi)."""
    out = np.zeros_like(phi)
    for m in range(0, (N + 1) // 2):
        out += np.cos((2 * m + 1) * phi)
    return out

def S_closed(phi, N):
    """Closed form sin((N+1) phi)/(2 sin phi), with the phi->0 limit (N+1)/2."""
    num = np.sin((N + 1) * phi)
    den = 2.0 * np.sin(phi)
    small = np.abs(phi) < 1e-9
    out = np.empty_like(phi)
    out[~small] = num[~small] / den[~small]
    out[small] = (N + 1) / 2.0
    return out

# ------------------------------- compute ------------------------------------ #
phi = np.linspace(-np.pi / 2, np.pi / 2, NPHI)

curves = {}
max_abs_err = 0.0
for N in NLIST:
    s_d = S_direct(phi, N)
    s_c = S_closed(phi, N)
    err = np.max(np.abs(s_d - s_c))
    max_abs_err = max(max_abs_err, err)
    peak = (N + 1) / 2.0                       # S_N(0)
    s_hat = s_d / peak                         # normalized amplitude  (peak = 1)
    i_hat = s_hat ** 2                          # normalized squared amplitude
    curves[N] = dict(s_hat=s_hat, i_hat=i_hat, peak=peak, err=err)

# ------------------------------- figure ------------------------------------- #
colors = ["#0b3d91", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
xdeg = np.rad2deg(phi)                          # x-axis in degrees of phase

fig, (axA, axI) = plt.subplots(2, 1, figsize=(9.2, 8.0), sharex=True)

for k, N in enumerate(NLIST):
    axA.plot(xdeg, curves[N]["s_hat"], color=colors[k], lw=1.5,
             label=r"$N=%d$" % N)
    axI.plot(xdeg, curves[N]["i_hat"], color=colors[k], lw=1.5,
             label=r"$N=%d$" % N)

axA.axhline(0.0, color="0.6", lw=0.7)
for ax in (axA, axI):
    ax.axvline(0.0, color="0.85", lw=0.8, zorder=0)
    ax.set_xlim(-90, 90)
    ax.grid(alpha=0.25)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])

axA.set_ylabel(r"normalized amplitude  $\hat S_N=S_N/S_N(0)$")
axA.set_title("Top: amplitude (paper Fig.1 convention)")
axA.legend(loc="upper right", fontsize=9, ncol=5, columnspacing=1.0)
axA.set_ylim(-0.35, 1.05)

axI.set_ylabel(r"normalized squared amplitude  $\hat I_N=|S_N|^2/|S_N(0)|^2$")
axI.set_title("Bottom: squared amplitude = observable (paper Fig.2 convention)")
axI.set_xlabel(r"phase  $\varphi$  (degrees;  half-wavelength interval $[-90^\circ,90^\circ]$)")
axI.set_ylim(-0.03, 1.05)

fig.suptitle(
    "Baseline localized peak wave  (c=1, tilt $\\theta=%g^\\circ$, far field, "
    "NO $\\nu\\!-\\!\\nu_0$ beat)\n"
    "two equal odd-harmonic sources, $N=1,3,5,7,9$  "
    "($\\sin\\theta=%.2f$ = projection scale for the future beat jitter)"
    % (THETA_DEG, SIN_THETA), y=0.99, fontsize=11)

fig.tight_layout(rect=[0, 0, 1, 0.96])
base = os.path.join(OUTDIR, "relphase_localized_baseline")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)

# ------------------------------- report ------------------------------------- #
print("Baseline (beat OFF):  c=%.0f,  theta=%.0f deg,  sin(theta)=%.4f,  W=%.1f" %
      (C, THETA_DEG, SIN_THETA, W))
print("direct-sum vs closed-form  max |abs err| = %.2e  (machine precision)" % max_abs_err)
print("%-4s  %-10s  %-22s  %-22s" % ("N", "peak S_N(0)", "amp half-width(deg @1/2)",
                                     "I half-width(deg @1/2)"))
for N in NLIST:
    s = curves[N]["s_hat"]; i = curves[N]["i_hat"]
    # half-maximum widths (first crossing from the centre)
    def half_width(y, level):
        idx = np.argmin(np.abs(phi))            # centre index
        j = idx
        while j < len(y) and y[j] >= level:
            j += 1
        return np.rad2deg(phi[j]) if j < len(y) else 90.0
    aw = half_width(s, 0.5)
    iw = half_width(i, 0.5)
    print("%-4d  %-10.1f  %-22.3f  %-22.3f" % (N, curves[N]["peak"], aw, iw))
print("Outputs ->", base + ".png / .svg")
