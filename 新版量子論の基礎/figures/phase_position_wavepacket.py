"""
Position as a complex phase with uncertainty.

Left panel:  complex phase representation
  - unit circle
  - central phase arrow at angle θ_x
  - sector (fan) of half-width Δθ_x  (an angular spread, NOT a radial one)
  - boundary lines at θ_x ± Δθ_x

Right panel: phase-space image of the wave packet
  - horizontal axis: phase θ in degrees (0° … 360°)
  - vertical axis : wave amplitude (wave height)
  - Re ψ(θ) = A(θ) cos(k (θ − θ_x) + θ_x)
  - A(θ)    = exp(−(θ − θ_x)² / (2 Δθ_x²))
  - solid Re ψ, dashed envelope ±A(θ)
  - θ_x marker and Δθ_x spread arrow

Convention used in this figure:
  • the central phase corresponds to θ_x = 180°
  • the phase spread (half-width) is Δθ_x = 60°
  • the carrier wavenumber is k = 1.0  (so phase θ = k · position)

Output: phase_position_wavepacket.png (300 dpi, 14 × 6 inch)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, FancyArrowPatch, Arc

# ---------------------------------------------------------------------------
# parameters (modify here)
# ---------------------------------------------------------------------------
theta_x_deg        = 180.0   # central phase                         [deg]
delta_theta_x_deg  = 60.0    # phase spread (HALF-WIDTH)             [deg]
k                  = 1.0     # carrier wavenumber  (θ = k · x)       [-]
theta_min_deg      = 0.0     # plot range, left edge                  [deg]
theta_max_deg      = 360.0   # plot range, right edge                 [deg]
out_path           = "phase_position_wavepacket.png"

# derived
theta_x       = np.deg2rad(theta_x_deg)
delta_theta_x = np.deg2rad(delta_theta_x_deg)

# ---------------------------------------------------------------------------
# figure layout
# ---------------------------------------------------------------------------
fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 6), facecolor="white")
fig.suptitle("Position as a complex phase with uncertainty",
             fontsize=16, fontweight="bold", y=0.98)

# ===========================================================================
# LEFT PANEL: complex phase representation
# ===========================================================================
ax_l.set_aspect("equal")
ax_l.set_xlim(-1.5, 1.5)
ax_l.set_ylim(-1.5, 1.5)
ax_l.set_xlabel("Re", fontsize=12)
ax_l.set_ylabel("Im", fontsize=12)
ax_l.set_title("Complex phase representation", fontsize=13, pad=10)

# axes through origin
ax_l.axhline(0, color="0.6", lw=0.6)
ax_l.axvline(0, color="0.6", lw=0.6)

# unit circle
phi = np.linspace(0, 2 * np.pi, 400)
ax_l.plot(np.cos(phi), np.sin(phi),
          color="0.55", lw=1.0, alpha=0.7)

# angle reference ticks on unit circle (every 90°)
for a_deg in (0, 90, 180, 270):
    a = np.deg2rad(a_deg)
    ax_l.plot([0.96 * np.cos(a), 1.04 * np.cos(a)],
              [0.96 * np.sin(a), 1.04 * np.sin(a)],
              color="0.45", lw=0.9)
    # push 180° label slightly off the horizontal so it doesn't collide
    # with the central-phase arrow tip
    if a_deg == 180:
        ax_l.text(1.20 * np.cos(a), 0.10,
                  f"{a_deg}°", fontsize=9, color="0.35",
                  ha="center", va="center")
    else:
        ax_l.text(1.13 * np.cos(a), 1.13 * np.sin(a),
                  f"{a_deg}°", fontsize=9, color="0.35",
                  ha="center", va="center")

# semi-transparent sector for ±Δθ_x
sector = Wedge(center=(0, 0), r=1.0,
               theta1=theta_x_deg - delta_theta_x_deg,
               theta2=theta_x_deg + delta_theta_x_deg,
               facecolor="tab:orange", edgecolor="none",
               alpha=0.28, zorder=1)
ax_l.add_patch(sector)

# two boundary lines θ_x ± Δθ_x
for sign in (-1, +1):
    a = theta_x + sign * delta_theta_x
    ax_l.plot([0, np.cos(a)], [0, np.sin(a)],
              color="tab:orange", lw=1.4, linestyle="--",
              alpha=0.9, zorder=2)

# central phase arrow (bold)
arrow_tip = (np.cos(theta_x), np.sin(theta_x))
arrow = FancyArrowPatch((0, 0), arrow_tip,
                        arrowstyle="-|>", mutation_scale=20,
                        lw=2.6, color="tab:blue", zorder=4)
ax_l.add_patch(arrow)

# label for the central phase arrow — placed in the upper-right area
# with a leader line pointing to the middle of the arrow shaft
shaft_mid = (0.55 * np.cos(theta_x), 0.55 * np.sin(theta_x))
ax_l.annotate(r"central phase $\theta_x = 180°$",
              xy=shaft_mid, xytext=(0.95, 1.05),
              fontsize=12, color="tab:blue",
              ha="center", va="bottom",
              arrowprops=dict(arrowstyle="->", color="tab:blue",
                              lw=1.0, connectionstyle="arc3,rad=0.25"))

# arc for θ_x (from positive Re axis around to the central direction)
arc_r = 0.32
theta_arc = Arc((0, 0), 2 * arc_r, 2 * arc_r,
                angle=0, theta1=0, theta2=theta_x_deg,
                color="tab:blue", lw=1.4)
ax_l.add_patch(theta_arc)
# place θ_x label near the top of that arc (where it crosses the +Im axis)
ax_l.text(0.0, arc_r + 0.10,
          r"$\theta_x$", color="tab:blue", fontsize=13,
          ha="center", va="bottom")

# outer arc showing that Δθ_x is an angular width (not a radial spread)
arc2_r = 1.22
ang_arc = Arc((0, 0), 2 * arc2_r, 2 * arc2_r, angle=0,
              theta1=theta_x_deg - delta_theta_x_deg,
              theta2=theta_x_deg + delta_theta_x_deg,
              color="tab:orange", lw=1.8)
ax_l.add_patch(ang_arc)

# tick marks at the two ends of the angular-width arc
for sign in (-1, +1):
    a = theta_x + sign * delta_theta_x
    ax_l.plot([arc2_r * np.cos(a) * 0.96, arc2_r * np.cos(a) * 1.04],
              [arc2_r * np.sin(a) * 0.96, arc2_r * np.sin(a) * 1.04],
              color="tab:orange", lw=1.8)

# label for the angular spread — annotation pointing to the outer arc midpoint
mid_outer = theta_x
arc_mid_pt = (arc2_r * np.cos(mid_outer), arc2_r * np.sin(mid_outer))
ax_l.annotate(r"phase spread $\pm\Delta\theta_x = \pm 60°$",
              xy=arc_mid_pt, xytext=(0.95, -1.05),
              fontsize=11.5, color="tab:orange",
              ha="center", va="top",
              arrowprops=dict(arrowstyle="->", color="tab:orange",
                              lw=1.0, connectionstyle="arc3,rad=-0.3"))

# tick labels and tidy
ax_l.set_xticks([-1, 0, 1])
ax_l.set_yticks([-1, 0, 1])
ax_l.grid(alpha=0.15)

# clarifying note
ax_l.text(0.0, -1.42,
          r"$\Delta\theta_x$ is an angular width, not a radial spread",
          fontsize=10.5, ha="center", va="center",
          style="italic", color="0.25")

# ===========================================================================
# RIGHT PANEL: phase-space image of the wave packet
# ===========================================================================
theta_deg = np.linspace(theta_min_deg, theta_max_deg, 1800)
theta_rad = np.deg2rad(theta_deg)

envelope = np.exp(-(theta_rad - theta_x) ** 2 / (2 * delta_theta_x ** 2))
psi_re   = envelope * np.cos(k * (theta_rad - theta_x) + theta_x)

ax_r.set_title("Phase-space image (wave packet)", fontsize=13, pad=10)
ax_r.set_xlabel(r"phase $\theta$ (degrees)", fontsize=12)
ax_r.set_ylabel("wave amplitude (wave height)", fontsize=12)
ax_r.set_xlim(theta_min_deg, theta_max_deg)
ax_r.set_ylim(-1.40, 1.65)
ax_r.axhline(0, color="0.6", lw=0.6)

# x ticks every 45°
ax_r.set_xticks(np.arange(0, 361, 45))

# envelope (dashed)
ax_r.plot(theta_deg,  envelope,
          color="tab:orange", lw=1.6, linestyle="--",
          label=r"envelope $+A(\theta)$")
ax_r.plot(theta_deg, -envelope,
          color="tab:orange", lw=1.6, linestyle="--",
          label=r"envelope $-A(\theta)$")

# Re ψ (solid)
ax_r.plot(theta_deg, psi_re, color="tab:blue", lw=2.2,
          label=r"Re $\psi(\theta)=A(\theta)\cos(k(\theta-\theta_x)+\theta_x)$")

# centre line at θ_x = 180°
ax_r.axvline(theta_x_deg, color="0.35", lw=1.0, linestyle=":")
ax_r.text(theta_x_deg - 5, 1.18,
          r"central phase $\theta_x = 180°$",
          color="0.20", fontsize=11, ha="right", va="center")

# spread arrow Δθ_x at a clear y-level
y_arrow = -1.20
ax_r.annotate("", xy=(theta_x_deg + delta_theta_x_deg, y_arrow),
              xytext=(theta_x_deg - delta_theta_x_deg, y_arrow),
              arrowprops=dict(arrowstyle="<->",
                              color="tab:orange", lw=1.8))
ax_r.text(theta_x_deg, y_arrow + 0.11,
          r"phase spread $\Delta\theta_x = 60°$",
          color="tab:orange", fontsize=11.5,
          ha="center", va="bottom")

# carrier-phase callout near a visible point of the wave
cx_deg = theta_x_deg + 30          # 30° to the right of centre
cx_rad = np.deg2rad(cx_deg)
cy = (np.exp(-(cx_rad - theta_x) ** 2 / (2 * delta_theta_x ** 2))
      * np.cos(k * (cx_rad - theta_x) + theta_x))
ax_r.annotate(r"carrier phase $\theta_x$",
              xy=(cx_deg, cy),
              xytext=(cx_deg + 35, cy + 0.55),
              fontsize=11, color="tab:blue",
              arrowprops=dict(arrowstyle="->", color="tab:blue",
                              lw=1.2,
                              connectionstyle="arc3,rad=-0.2"))

# parameter / correspondence box (lower-right corner of right panel)
param_text = ("parameters:\n"
              r"  $k = 1.0$,  $\theta_x = 180°$,  $\Delta\theta_x = 60°$"
              "\n"
              "correspondence:\n"
              r"  $\theta_x = k\,x_0$,  $\Delta\theta_x = k\,\Delta x$")
ax_r.text(0.985, -1.32, param_text,
          transform=ax_r.transData,
          fontsize=10, ha="right", va="bottom",
          color="0.20",
          bbox=dict(boxstyle="round,pad=0.4",
                    fc="white", ec="0.7", alpha=0.92))

ax_r.grid(alpha=0.15)
ax_r.legend(loc="upper right", fontsize=9.5, framealpha=0.9)

# ---------------------------------------------------------------------------
# save
# ---------------------------------------------------------------------------
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
print(f"saved: {out_path}")
