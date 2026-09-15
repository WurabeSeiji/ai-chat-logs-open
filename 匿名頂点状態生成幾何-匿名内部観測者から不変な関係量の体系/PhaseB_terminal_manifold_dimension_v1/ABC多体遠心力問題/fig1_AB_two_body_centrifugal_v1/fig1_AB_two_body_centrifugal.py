#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1: Two-body centrifugal motion (A-B)

Series : ABC多体遠心力問題 (ABC many-body centrifugal problem)
Step 1 : two-body A-B configuration.

Setup (as specified)
--------------------
- A sits at the origin O, mass M_A.
- B sits on a circle of radius R around A, mass M_B.
- B rotates around A with angular velocity omega.
- The gravitational acceleration g_ab (of B toward A) and the
  centrifugal acceleration c_ab (of B, outward) balance:

      G * M_A / R^2  =  omega^2 * R          (balance condition)

Assumptions recorded in the figure
----------------------------------
(1) Only the two bodies A, B exist in spacetime.
(2) Strictly, the pair should rotate about the A-B barycenter;
    this is ignored for now: A is held fixed at the origin O.
(3) B rotates with angular velocity omega, but since there is no
    counterpart against which to measure the phase, the phase is
    frozen at theta = 0 (B drawn at rest on the +x axis).

Units (deliberate choice for this series)
-----------------------------------------
c = 1, G = 1, and the Planck constant h = 1.
h = 1 is set deliberately even though h does not appear at this
Newtonian step.  Formulas are first written with the constants
(c, G, h) kept explicit; at this step only G actually enters
(c and h will enter in later steps of the series).  The reduced
natural-unit form is then obtained by setting c = G = h = 1.

    Full form   : g_ab = G * M_A / R^2 ,  c_ab = omega^2 * R
    Natural unit: g_ab = M_A / R^2     ,  c_ab = omega^2 * R
    Balance     : M_A / R^2 = omega^2 * R  =>  omega = sqrt(M_A / R^3)

This is an analytic figure: every plotted element (positions, orbit,
arrow lengths) is computed from the formulas above, not drawn by hand.

Outputs
-------
  fig1_AB_two_body_centrifugal.svg
  fig1_AB_two_body_centrifugal.png
  results.txt
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# ----------------------------------------------------------------------
# Constants: natural units  c = G = h = 1  (h deliberately set to 1)
# ----------------------------------------------------------------------
c_light = 1.0   # speed of light          c = 1
G_grav  = 1.0   # gravitational constant  G = 1
h_planck = 1.0  # Planck constant         h = 1 (deliberate; unused at
                #                          this Newtonian step)

# ----------------------------------------------------------------------
# Parameters of the A-B two-body configuration (natural units)
# ----------------------------------------------------------------------
M_A = 1.0       # mass of A (at origin O)
M_B = 0.1       # mass of B (on the circle; does not enter the balance)
R   = 1.0       # orbital radius of B around A
theta = 0.0     # phase frozen at 0 (no counterpart to measure phase)

# ----------------------------------------------------------------------
# Analytic content
#   Full form (constants explicit):  g_ab = G M_A / R^2,  c_ab = omega^2 R
#   Balance:  G M_A / R^2 = omega^2 R  =>  omega = sqrt(G M_A / R^3)
#   Natural units (c = G = h = 1):   g_ab = M_A / R^2,  omega = sqrt(M_A/R^3)
# ----------------------------------------------------------------------
omega = np.sqrt(G_grav * M_A / R**3)   # angular velocity from balance
g_ab  = G_grav * M_A / R**2            # gravitational acceleration on B
c_ab  = omega**2 * R                   # centrifugal acceleration on B
residual = abs(g_ab - c_ab)

# Position of B (phase frozen at theta = 0)
xB, yB = R * np.cos(theta), R * np.sin(theta)

lines = [
    "Figure 1: Two-body centrifugal motion (A-B)  --  computed values",
    f"  units      : c = {c_light}, G = {G_grav}, h = {h_planck} (natural units, h deliberate)",
    f"  M_A        = {M_A}",
    f"  M_B        = {M_B}   (does not enter the balance)",
    f"  R          = {R}",
    f"  theta      = {theta}  (phase frozen)",
    f"  omega      = sqrt(G M_A / R^3) = {omega:.15g}",
    f"  g_ab       = G M_A / R^2       = {g_ab:.15g}",
    f"  c_ab       = omega^2 R         = {c_ab:.15g}",
    f"  |g_ab-c_ab|= {residual:.15g}   (balance check)",
]
results_text = "\n".join(lines)
print(results_text)
with open("results.txt", "w") as f:
    f.write(results_text + "\n")

# ----------------------------------------------------------------------
# Figure (all labels in English)
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11.5, 7.0))
ax.set_aspect("equal")

# Orbit circle (radius R, centered at O)
tt = np.linspace(0.0, 2.0 * np.pi, 720)
ax.plot(R * np.cos(tt), R * np.sin(tt), color="0.55", lw=1.0, ls="--",
        zorder=1)

# Axes through O
ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)
ax.axvline(0.0, color="0.85", lw=0.8, zorder=0)

# Body A at the origin O
ax.plot(0.0, 0.0, "o", ms=14, color="black", zorder=5)
ax.annotate("A  (mass $M_A$)", xy=(0.0, 0.0), xytext=(-0.05, -0.17),
            ha="right", va="top", fontsize=12)
ax.annotate("O (origin)", xy=(0.0, 0.0), xytext=(-0.05, 0.10),
            ha="right", va="bottom", fontsize=10, color="0.35")

# Body B on the circle at theta = 0
ax.plot(xB, yB, "o", ms=10, color="tab:green", zorder=5)
ax.annotate("B  (mass $M_B$)", xy=(xB, yB), xytext=(xB + 0.04, yB - 0.17),
            ha="left", va="top", fontsize=12)
ax.annotate(r"$\theta = 0$ (phase frozen)", xy=(xB, yB),
            xytext=(xB - 0.02, yB - 0.32), ha="center", va="top",
            fontsize=9, color="0.35")

# Separation R (line O-B), label below the line
ax.plot([0.0, xB], [0.0, yB], color="0.3", lw=1.2, zorder=2)
ax.annotate("$R$", xy=(0.35 * xB, 0.0), xytext=(0.35 * xB, -0.09),
            ha="center", va="top", fontsize=13)

# Acceleration arrows at B, lengths proportional to computed magnitudes
s = 0.45  # common scale factor [figure length per acceleration unit]
ax.add_patch(FancyArrowPatch((xB, yB), (xB - s * g_ab, yB),
                             arrowstyle="-|>", mutation_scale=22,
                             color="tab:blue", lw=2.4, zorder=6))
ax.annotate(r"$g_{ab} = \dfrac{G M_A}{R^2}$",
            xy=(xB - s * g_ab, yB), xytext=(xB - s * g_ab - 0.10, yB + 0.12),
            ha="left", va="bottom", fontsize=13, color="tab:blue")

ax.add_patch(FancyArrowPatch((xB, yB), (xB + s * c_ab, yB),
                             arrowstyle="-|>", mutation_scale=22,
                             color="tab:red", lw=2.4, zorder=6))
ax.annotate(r"$c_{ab} = \omega^2 R$",
            xy=(xB + s * c_ab, yB), xytext=(xB + s * c_ab + 0.06, yB + 0.12),
            ha="right", va="bottom", fontsize=13, color="tab:red")

# Rotation sense omega (arc arrow on the orbit), drawn but frozen
a0, a1 = np.deg2rad(38.0), np.deg2rad(72.0)
r_arc = 1.13 * R
arc_t = np.linspace(a0, a1, 60)
ax.plot(r_arc * np.cos(arc_t), r_arc * np.sin(arc_t), color="0.4", lw=1.4,
        zorder=3)
ax.add_patch(FancyArrowPatch(
    (r_arc * np.cos(a1 - 0.04), r_arc * np.sin(a1 - 0.04)),
    (r_arc * np.cos(a1), r_arc * np.sin(a1)),
    arrowstyle="-|>", mutation_scale=16, color="0.4", lw=1.4, zorder=3))
ax.annotate(r"$\omega$ (angular velocity;" "\n"
            r"phase frozen at $\theta=0$)",
            xy=(r_arc * np.cos(0.5 * (a0 + a1)),
                r_arc * np.sin(0.5 * (a0 + a1))),
            xytext=(0.30, 1.32), ha="left", va="bottom", fontsize=11,
            color="0.25")

# --- Equation box (full form with constants, then natural units) -------
eq_text = (
    "Full form (constants explicit):\n"
    r"  $g_{ab} = \dfrac{G\,M_A}{R^2}$" "\n"
    r"  $c_{ab} = \omega^2 R$" "\n"
    r"  balance:  $\dfrac{G\,M_A}{R^2} = \omega^2 R$"
    r"  $\Rightarrow\ \omega = \sqrt{\dfrac{G\,M_A}{R^3}}$" "\n"
    "\n"
    "Natural units  (c = 1, G = 1, h = 1):\n"
    r"  $g_{ab} = \dfrac{M_A}{R^2}$,   $c_{ab} = \omega^2 R$,"
    r"   $\omega = \sqrt{\dfrac{M_A}{R^3}}$"
)
ax.text(1.82, 1.58, eq_text, fontsize=11, ha="left", va="top",
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", fc="#f4f7ff", ec="tab:blue",
                  lw=1.0))

# --- Units / assumptions box ------------------------------------------
note_text = (
    "Units:  c = 1,  G = 1,  and Planck constant h = 1\n"
    "  (h = 1 is deliberate; only G enters at this Newtonian step,\n"
    "   c and h will enter in later steps of the series)\n"
    "\n"
    "Assumptions:\n"
    "  (1) only two bodies A, B exist in spacetime\n"
    "  (2) barycenter rotation ignored: A held fixed at O\n"
    "  (3) no counterpart to measure phase:\n"
    r"        B frozen at $\theta = 0$" "\n"
    "\n"
    f"Computed ($M_A$={M_A}, $R$={R}):\n"
    f"  $\\omega$ = {omega:g},  $g_{{ab}}$ = {g_ab:g},  $c_{{ab}}$ = {c_ab:g},"
    f"  $|g_{{ab}}-c_{{ab}}|$ = {residual:g}"
)
ax.text(1.82, -0.38, note_text, fontsize=10, ha="left", va="top",
        linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.5", fc="#fffaf0", ec="0.5",
                  lw=0.8))

ax.set_title("Figure 1: Two-body centrifugal motion (A-B)", fontsize=15,
             pad=14)
ax.set_xlabel("$x$", fontsize=12)
ax.set_ylabel("$y$", fontsize=12)
ax.set_xlim(-1.55, 4.55)
ax.set_ylim(-1.85, 1.80)

fig.tight_layout()
fig.savefig("fig1_AB_two_body_centrifugal.svg")
fig.savefig("fig1_AB_two_body_centrifugal.png", dpi=200)
print("saved: fig1_AB_two_body_centrifugal.svg / .png / results.txt")
