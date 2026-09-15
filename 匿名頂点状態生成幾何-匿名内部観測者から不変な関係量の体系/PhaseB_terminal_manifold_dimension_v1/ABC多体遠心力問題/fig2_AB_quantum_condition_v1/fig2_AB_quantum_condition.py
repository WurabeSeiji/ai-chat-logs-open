#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2: Quantum condition on the two-body system (A-B)

Series : ABC多体遠心力問題 (ABC many-body centrifugal problem)
Step 2 : quantum condition added to the Fig.1 configuration.
The drawing is a faithful copy of Fig.1 (same geometry); only the
equation boxes are replaced.

Quantum condition (as specified)
--------------------------------
Frequencies nu_a, nu_b and wavelengths lambda_a, lambda_b are assigned
to A and B; winding numbers k_a, k_b are defined by

    2 pi R = k_a * (1/2) lambda_a ,   2 pi R = k_b * (1/2) lambda_b .

Derived relations
-----------------
    lambda_a = 4 pi R / k_a ,        lambda_b = 4 pi R / k_b
    nu_a     = k_a omega / (4 pi) ,  nu_b     = k_b omega / (4 pi)
    c = nu * lambda = omega R   =>   R = c / omega

Mass ratio
----------
    E = h nu ,  E = M c^2   =>   M = h nu / c^2
    M_a / M_b = nu_a / nu_b = k_a / k_b

Units: c = 1, G = 1, and Planck constant h = 1 (deliberate).
Formulas are written with the constants explicit first, then in the
reduced natural-unit form:
    nu = k omega/(4 pi),  lambda = 4 pi/(k omega),  R = 1/omega.

Outputs
-------
  fig2_AB_quantum_condition.svg
  fig2_AB_quantum_condition.png
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
c_light  = 1.0
G_grav   = 1.0
h_planck = 1.0

# ----------------------------------------------------------------------
# Parameters (same as Fig.1)
# ----------------------------------------------------------------------
M_A = 1.0
M_B = 0.1
R   = 1.0
theta = 0.0

omega = np.sqrt(G_grav * M_A / R**3)
g_ab  = G_grav * M_A / R**2
c_ab  = omega**2 * R
xB, yB = R * np.cos(theta), R * np.sin(theta)

lines = [
    "Figure 2: Quantum condition on the two-body system (A-B)",
    f"  units : c = {c_light}, G = {G_grav}, h = {h_planck} (natural units)",
    "  quantum condition : 2 pi R = k_a (1/2) lambda_a = k_b (1/2) lambda_b",
    "  derived           : lambda = 4 pi R / k",
    "                      nu     = k omega / (4 pi)",
    "                      c = nu lambda = omega R  =>  R = c / omega",
    "  natural units     : nu = k omega/(4 pi), lambda = 4 pi/(k omega), R = 1/omega",
    "  mass ratio        : M = h nu / c^2  =>  M_a / M_b = nu_a / nu_b = k_a / k_b",
]
results_text = "\n".join(lines)
print(results_text)
with open("results.txt", "w") as f:
    f.write(results_text + "\n")

# ----------------------------------------------------------------------
# Figure (all labels in English) -- same geometry as Fig.1
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

# Acceleration arrows at B (state of Fig.1, labels shortened to symbols)
s = 0.45
ax.add_patch(FancyArrowPatch((xB, yB), (xB - s * g_ab, yB),
                             arrowstyle="-|>", mutation_scale=22,
                             color="tab:blue", lw=2.4, zorder=6))
ax.annotate(r"$g_{ab}$", xy=(xB - s * g_ab, yB),
            xytext=(xB - s * g_ab - 0.02, yB + 0.10),
            ha="left", va="bottom", fontsize=13, color="tab:blue")

ax.add_patch(FancyArrowPatch((xB, yB), (xB + s * c_ab, yB),
                             arrowstyle="-|>", mutation_scale=22,
                             color="tab:red", lw=2.4, zorder=6))
ax.annotate(r"$c_{ab}$", xy=(xB + s * c_ab, yB),
            xytext=(xB + s * c_ab + 0.02, yB + 0.10),
            ha="right", va="bottom", fontsize=13, color="tab:red")

# Rotation sense omega (arc arrow on the orbit)
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

# --- Quantum condition and derived relations ---------------------------
eq_text = (
    "Quantum condition:\n"
    r"  $2\pi R = k_a\,\dfrac{\lambda_a}{2}$,"
    r"   $2\pi R = k_b\,\dfrac{\lambda_b}{2}$" "\n"
    "\n"
    "Derived (constants explicit):\n"
    r"  $\lambda_a = \dfrac{4\pi R}{k_a}$,"
    r"   $\lambda_b = \dfrac{4\pi R}{k_b}$" "\n"
    r"  $\nu_a = \dfrac{k_a\,\omega}{4\pi}$,"
    r"   $\nu_b = \dfrac{k_b\,\omega}{4\pi}$" "\n"
    r"  $c = \nu\lambda = \omega R\ \Rightarrow\ R = c/\omega$" "\n"
    "\n"
    "Natural units  (c = 1, G = 1, h = 1):\n"
    r"  $\nu = \dfrac{k\omega}{4\pi}$,   $\lambda = \dfrac{4\pi}{k\omega}$,"
    r"   $R = \dfrac{1}{\omega}$"
)
ax.text(1.82, 1.58, eq_text, fontsize=11, ha="left", va="top",
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", fc="#f4f7ff", ec="tab:blue",
                  lw=1.0))

# --- Mass ratio --------------------------------------------------------
mass_text = (
    "Mass ratio:\n"
    r"  $E = h\nu$,   $E = Mc^2$"
    r"  $\Rightarrow\ M = \dfrac{h\nu}{c^2}$" "\n"
    "\n"
    r"  $\dfrac{M_a}{M_b} = \dfrac{\nu_a}{\nu_b} = \dfrac{k_a}{k_b}$"
)
ax.text(1.82, -0.72, mass_text, fontsize=11, ha="left", va="top",
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", fc="#fffaf0", ec="0.5",
                  lw=0.8))

ax.set_title("Figure 2: Quantum condition on the two-body system (A-B)",
             fontsize=15, pad=14)
ax.set_xlabel("$x$", fontsize=12)
ax.set_ylabel("$y$", fontsize=12)
ax.set_xlim(-1.55, 4.55)
ax.set_ylim(-1.85, 1.80)

fig.tight_layout()
fig.savefig("fig2_AB_quantum_condition.svg")
fig.savefig("fig2_AB_quantum_condition.png", dpi=200)
print("saved: fig2_AB_quantum_condition.svg / .png / results.txt")
