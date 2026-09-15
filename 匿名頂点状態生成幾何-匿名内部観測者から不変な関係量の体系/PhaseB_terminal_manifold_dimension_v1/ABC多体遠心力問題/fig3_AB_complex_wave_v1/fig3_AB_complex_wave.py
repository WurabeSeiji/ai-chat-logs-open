#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3: Complex wave representation of A and B

Series : ABC多体遠心力問題 (ABC many-body centrifugal problem)
Step 3 : complex wave representation and amplitudes.
The drawing is a faithful copy of Fig.2 (same geometry, no waves are
drawn); only the equation boxes are replaced.

Complex wave representation (the central result)
------------------------------------------------
    z_a = A_a * exp(i k_a theta / 2)
    z_b = A_b * exp(i k_b theta / 2)
theta: angle on the circle.  No omega, nu, lambda, t appear: the wave
is the winding phase itself.  |z| is constant, so there is no wave
shape to draw on the circle.

Amplitudes (the postponed barycenter assumption is restored)
------------------------------------------------------------
    M_a A_a = M_b A_b ,   A_a + A_b = R
    =>  A_a = k_b R / (k_a + k_b) ,  A_b = k_a R / (k_a + k_b)
    A_a / A_b = k_b / k_a = M_b / M_a

Eliminated quantities (aliases of R and k)
------------------------------------------
    omega = c / R ,  nu = k c / (4 pi R) ,  lambda = 4 pi R / k

Invariants
----------
    A_a k_a = A_b k_b ,   A_a M_a = A_b M_b = mu R
    M lambda = h / c      (natural units: M lambda = 1)

Units: c = 1, G = 1, and Planck constant h = 1 (deliberate).

Outputs
-------
  fig3_AB_complex_wave.svg
  fig3_AB_complex_wave.png
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
# Parameters (same as Fig.1 / Fig.2)
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
    "Figure 3: Complex wave representation of A and B",
    f"  units : c = {c_light}, G = {G_grav}, h = {h_planck} (natural units)",
    "  complex waves : z_a = A_a exp(i k_a theta/2), z_b = A_b exp(i k_b theta/2)",
    "  amplitudes    : M_a A_a = M_b A_b, A_a + A_b = R",
    "                  => A_a = k_b R/(k_a+k_b), A_b = k_a R/(k_a+k_b)",
    "                  A_a/A_b = k_b/k_a = M_b/M_a",
    "  eliminated    : omega = c/R, nu = k c/(4 pi R), lambda = 4 pi R/k",
    "  invariants    : A_a k_a = A_b k_b, A_a M_a = A_b M_b = mu R",
    "                  M lambda = h/c (natural units: M lambda = 1)",
    "  expansion     : g_ab = G M_a/R^2 = G h k_a/(4 pi c R^3)",
    "                  c_ab = omega^2 A_b = c^2 k_a/(R(k_a+k_b))",
    "                  (for A: g_ba, c_ba by k_a -> k_b)",
    "                  balance (B) = balance (A) => k_a + k_b = 4 pi c^3 R^2/(G h)",
    "                  natural units: g_ab = k_a/(4 pi R^3), c_ab = k_a/(R(k_a+k_b))",
    "                  => k_a + k_b = 4 pi R^2,  M_a + M_b = R",
]
results_text = "\n".join(lines)
print(results_text)
with open("results.txt", "w") as f:
    f.write(results_text + "\n")

# ----------------------------------------------------------------------
# Figure (all labels in English) -- same geometry as Fig.2, no waves
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11.5, 8.9))
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
ax.annotate(r"$\theta = 0$ (drawing convention)", xy=(xB, yB),
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

# Rotation sense (arc arrow on the orbit)
a0, a1 = np.deg2rad(38.0), np.deg2rad(72.0)
r_arc = 1.13 * R
arc_t = np.linspace(a0, a1, 60)
ax.plot(r_arc * np.cos(arc_t), r_arc * np.sin(arc_t), color="0.4", lw=1.4,
        zorder=3)
ax.add_patch(FancyArrowPatch(
    (r_arc * np.cos(a1 - 0.04), r_arc * np.sin(a1 - 0.04)),
    (r_arc * np.cos(a1), r_arc * np.sin(a1)),
    arrowstyle="-|>", mutation_scale=16, color="0.4", lw=1.4, zorder=3))
ax.annotate(r"$\theta$ (angle on the circle)",
            xy=(r_arc * np.cos(0.5 * (a0 + a1)),
                r_arc * np.sin(0.5 * (a0 + a1))),
            xytext=(0.30, 1.32), ha="left", va="bottom", fontsize=11,
            color="0.25")

# --- Complex wave representation and amplitudes ------------------------
eq_text = (
    "Complex wave representation:\n"
    r"  $z_a = A_a\,e^{\,i k_a \theta/2}$,"
    r"   $z_b = A_b\,e^{\,i k_b \theta/2}$" "\n"
    "\n"
    "Amplitudes (barycenter restored):\n"
    r"  $M_a A_a = M_b A_b$,   $A_a + A_b = R$" "\n"
    r"  $\Rightarrow\ A_a = \dfrac{k_b}{k_a+k_b}\,R$,"
    r"   $A_b = \dfrac{k_a}{k_a+k_b}\,R$" "\n"
    r"  $\dfrac{A_a}{A_b} = \dfrac{k_b}{k_a} = \dfrac{M_b}{M_a}$"
)
ax.text(1.82, 1.58, eq_text, fontsize=11, ha="left", va="top",
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", fc="#f4f7ff", ec="tab:blue",
                  lw=1.0))

# --- Eliminated quantities and invariants ------------------------------
inv_text = (
    "Eliminated (aliases of $R$ and $k$):\n"
    r"  $\omega = \dfrac{c}{R}$,   $\nu = \dfrac{k\,c}{4\pi R}$,"
    r"   $\lambda = \dfrac{4\pi R}{k}$" "\n"
    "\n"
    "Invariants  (c = 1, G = 1, h = 1):\n"
    r"  $A_a k_a = A_b k_b$,   $A_a M_a = A_b M_b = \mu R$" "\n"
    r"  $M\lambda = h/c$   (natural units: $M\lambda = 1$)"
)
ax.text(1.82, -0.35, inv_text, fontsize=11, ha="left", va="top",
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", fc="#fffaf0", ec="0.5",
                  lw=0.8))

# --- Expansion of g_ab and c_ab in (R, k) ------------------------------
exp_text = (
    r"Expansion of $g_{ab}$, $c_{ab}$ (barycenter restored):" "\n"
    r"  $g_{ab} = G M_a/R^2 = G h k_a/(4\pi c R^3)$" "\n"
    r"  $c_{ab} = \omega^2 A_b = c^2 k_a\,/\,[R\,(k_a+k_b)]$" "\n"
    r"  (for A: $g_{ba}$, $c_{ba}$ by $k_a \to k_b$)" "\n"
    r"  balance for B $=$ balance for A $\Rightarrow$" "\n"
    r"  $k_a + k_b = (4\pi c^3/\,G h)\,R^2$" "\n"
    "\n"
    "Natural units (c = 1, G = 1, h = 1):\n"
    r"  $g_{ab} = k_a/(4\pi R^3)$,   $c_{ab} = k_a/[R\,(k_a+k_b)]$" "\n"
    r"  $\Rightarrow\ k_a + k_b = 4\pi R^2$,   $M_a + M_b = R$"
)
ax.text(-1.42, -1.22, exp_text, fontsize=11, ha="left", va="top",
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", fc="#f2fbf2", ec="tab:green",
                  lw=1.0))

ax.set_title("Figure 3: Complex wave representation of A and B",
             fontsize=15, pad=14)
ax.set_xlabel("$x$", fontsize=12)
ax.set_ylabel("$y$", fontsize=12)
ax.set_xlim(-1.55, 4.55)
ax.set_ylim(-2.95, 1.80)

fig.tight_layout()
fig.savefig("fig3_AB_complex_wave.svg")
fig.savefig("fig3_AB_complex_wave.png", dpi=200)
print("saved: fig3_AB_complex_wave.svg / .png / results.txt")
