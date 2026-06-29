#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 nu-nu0 BEAT on the combined odd-harmonic localized wave.
 Base frequency nu = 1 (lambda=1, c=1).  Observer frequency nu0 = pi/2.
 The beat shifts the localized-peak CENTRE.  x-axis = phase, y-axis = I = |psi|^2.
 (Standalone; NOT yet combined with the exact-Doppler figure.)
================================================================================

Localized wave (combined odd modes n=1,3,5,7,9):
    S(phi) = sum_n cos(n phi) ,    I(phi) = S(phi)^2  (period 180 deg).

Beat:
    nu = 1, nu0 = pi/2  =>  ratio nu/nu0 = 2/pi  (IRRATIONAL) => Weyl equidistribution.
    beat phase   phi_beat(t) = 2 pi (nu - nu0) t        (mod 2 pi)
    centre shift phi_centre   = sin(theta) * phi_beat ,  sin(30 deg) = 0.5
    => over one beat cycle phi_beat in [0,2pi), the centre sweeps [0, pi) = one
       full intensity period (180 deg): one fringe of jitter.

Two things are shown:
  (top)    snapshots I(phi - phi_centre) at several beat phases: the peak JITTERS.
           Plus the long-time TIME-AVERAGE of the intensity -> nearly FLAT
           (the jitter washes the fringe out).  Time-average is NOT |psi|^2.
  (bottom) histogram of phi_centre over discrete measurement times t_n = n with
           nu0 = pi/2: it is FLAT -> the irrational ratio gives equidistribution,
           confirming pi/2 is a sound choice.
           (The Born-weighted SINGLE-SHOT histogram would instead reproduce
            |S(phi)|^2; that is the conjecture's separate claim, not shown here.)

Author: N. Kihara (peer-review dialogue note), 2026.  Beat ON, nu0 = pi/2.
================================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------- parameters --------------------------------- #
MODES     = [1, 3, 5, 7, 9]
NU        = 1.0
NU0       = np.pi / 2.0             # observer frequency (irrational ratio 2/pi)
THETA_DEG = 30.0
NPHI      = 400001
NMEAS     = 20000                   # discrete measurement times t_n = n
OUTDIR    = os.path.dirname(os.path.abspath(__file__))

sin_theta = np.sin(np.deg2rad(THETA_DEG))      # 0.5
dnu = NU - NU0                                  # beat frequency = 1 - pi/2

def S(phi):
    out = np.zeros_like(phi)
    for n in MODES:
        out += np.cos(n * phi)
    return out

# ------------------------------- snapshots ---------------------------------- #
phi_deg = np.linspace(-90, 270, NPHI)
phi = np.deg2rad(phi_deg)
peak0 = float(len(MODES))                       # S(0) = number of modes = 5
I_unit = (S(phi - 0.0) ** 2) / peak0**2          # for shape ref

# beat-phase snapshots (evenly spaced over one beat cycle)
beat_phases = np.linspace(0, 2*np.pi, 5, endpoint=False)   # 0, 72,144,216,288 deg
centre_shifts = sin_theta * beat_phases                    # 0,36,72,108,144 deg(rad)

# long-time average over the discrete measurement times with nu0=pi/2
n_idx = np.arange(1, NMEAS + 1)
phi_beat_n = np.mod(2*np.pi * dnu * n_idx, 2*np.pi)         # equidistributes (2/pi irrational)
phi_centre_n = sin_theta * phi_beat_n                       # in [0, pi)
# time-average intensity on a coarser grid (cost)
phig = np.deg2rad(np.linspace(-90, 270, 4001))
acc = np.zeros_like(phig)
for pc in phi_centre_n:
    acc += S(phig - pc) ** 2
acc /= NMEAS
acc /= acc.max()

# ------------------------------- figure ------------------------------------- #
fig, (axT, axB) = plt.subplots(2, 1, figsize=(12.6, 8.4),
                               gridspec_kw=dict(height_ratios=[2.2, 1.0]))

cmap = plt.get_cmap("viridis")
for k, (bp, cs) in enumerate(zip(beat_phases, centre_shifts)):
    I = (S(phi - cs) ** 2) / peak0**2
    axT.plot(phi_deg, I, color=cmap(k/len(beat_phases)), lw=1.4,
             label=r"$\varphi_{\rm beat}=%3.0f^\circ \Rightarrow \varphi_{\rm centre}=%3.0f^\circ$"
                   % (np.rad2deg(bp), np.rad2deg(cs)))
axT.plot(np.rad2deg(phig), acc, color="#d62728", lw=2.0, ls="--",
         label="long-time average (washes out -> nearly flat; NOT $|\\psi|^2$)")
axT.set_xlim(-90, 270)
axT.set_ylim(-0.03, 1.06)
axT.set_xticks(range(-90, 271, 45))
axT.set_xlabel(r"phase  $\varphi$  (degrees)")
axT.set_ylabel(r"squared amplitude  $I=|S|^2/|S(0)|^2$")
axT.set_title(r"Beat-induced centre jitter:  $\nu=1,\ \nu_0=\pi/2$  "
              r"($\nu/\nu_0=2/\pi$ irrational),  $\sin\theta=0.5$  "
              r"$\Rightarrow$ centre sweeps one fringe (0–180°)")
axT.legend(loc="upper right", fontsize=8.0)
axT.grid(alpha=0.25)

# bottom: equidistribution of the centre over discrete times t_n with nu0=pi/2
axB.hist(np.rad2deg(phi_centre_n), bins=90, range=(0, 180),
         color="#1f77b4", alpha=0.8, density=True)
axB.axhline(1.0/180.0, color="k", lw=1.0, ls="--", label="uniform $1/180$")
axB.set_xlim(0, 180)
axB.set_xlabel(r"centre position  $\varphi_{\rm centre}=\sin\theta\,\varphi_{\rm beat}$  (deg)")
axB.set_ylabel("density")
axB.set_title(r"Equidistribution of the centre over $t_n=n$ with $\nu_0=\pi/2$ "
              r"(%d shots): FLAT $\Rightarrow$ irrational ratio works" % NMEAS)
axB.legend(loc="upper right", fontsize=8.5)
axB.grid(alpha=0.25)

fig.tight_layout()
base = os.path.join(OUTDIR, "relphase_beat_nu0")
fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
plt.close(fig)

# ------------------------------- report ------------------------------------- #
# star-discrepancy-ish check: max deviation of the empirical CDF from uniform
xs = np.sort(phi_centre_n / np.pi)              # in [0,1)
emp = np.arange(1, len(xs)+1) / len(xs)
disc = np.max(np.abs(emp - xs))
print("nu=%.4f  nu0=pi/2=%.4f  beat dnu=%.4f  ratio nu/nu0=2/pi=%.6f" %
      (NU, NU0, dnu, NU/NU0))
print("sin(theta)=%.3f  centre-jitter range = [0, %.0f] deg (one fringe)" %
      (sin_theta, np.rad2deg(sin_theta*2*np.pi)))
print("equidistribution star-discrepancy (%d shots) = %.4e  (-> 0 = uniform)" %
      (NMEAS, disc))
print("Outputs ->", base + ".png / .svg")
