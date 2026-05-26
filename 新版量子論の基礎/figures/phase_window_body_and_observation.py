"""
Rectangular phase-window body and its low-pass observation image,
shown as actual Fourier-series partial sums.

This figure visualises §7 of 思考実験(9):
  Body (本体)        : R_x(θ) = 1 for |θ−θ_x| ≤ Δθ_x/2, 0 otherwise (period 2π)
  Observation (観測像): P_x^obs(θ) = L_Λ[R_x(θ)]
                       = truncated Fourier series up to harmonic N

Convention (§7):
  Δθ_x is the FULL window width (so the half-width is Δθ_x/2).

For R_x(θ) centred at θ_x = π, full width Δθ_x = 2π/3, period 2π:

  R_x(θ)  =  a_0/2  +  Σ_{n≥1} a_n cos(n(θ − θ_x))
  a_0     = 2 · (Δθ_x / 2π)        =  2/3
  a_n     = (2 / (nπ)) · sin(n · Δθ_x / 2)
          = (2 / (nπ)) · sin(nπ/3)

  (sin(nπ/3) = 0 at n = 3, 6, 9, … so those harmonics vanish.)

Layout: 1 horizontal PNG, two panels side by side:
  Left  : body R_x(θ) — sharp rectangle
  Right : body (faint) + Fourier partial sums at N = 1, 3, 7

Output: phase_window_body_and_observation.png (300 dpi, 14 × 6 inch)
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# parameters (modify here)
# ---------------------------------------------------------------------------
theta_x_deg        = 180.0     # central phase                        [deg]
delta_theta_x_deg  = 120.0     # FULL window width  (per §7)          [deg]
theta_min_deg      = 0.0
theta_max_deg      = 360.0
N_list             = [1, 3, 7] # harmonics to show in partial sums
out_path           = "phase_window_body_and_observation.png"

# ---------------------------------------------------------------------------
# derived
# ---------------------------------------------------------------------------
theta_x       = np.deg2rad(theta_x_deg)
delta_theta_x = np.deg2rad(delta_theta_x_deg)
half_width    = delta_theta_x / 2.0
theta_left    = theta_x - half_width
theta_right   = theta_x + half_width

theta_rad  = np.linspace(0.0, 2 * np.pi, 4000)
theta_deg  = np.rad2deg(theta_rad)

# Body : rectangular window (period 2π)
R_x = np.where(np.abs(theta_rad - theta_x) <= half_width, 1.0, 0.0)


# ---------------------------------------------------------------------------
# Fourier series: R_x(θ) = a0/2 + Σ a_n cos(n(θ − θ_x))
#   (only cosine terms in the body-centred expansion, since R_x is
#    even about θ_x)
# ---------------------------------------------------------------------------
def fourier_coeff(n):
    """Cosine coefficient a_n for the rectangle of full width Δθ_x
    centred at θ_x on a period-2π domain."""
    if n == 0:
        return 2.0 * (delta_theta_x / (2 * np.pi))  # = 2 · duty cycle
    return (2.0 / (n * np.pi)) * np.sin(n * delta_theta_x / 2.0)


def partial_sum(theta_rad, N):
    """Sum a_0/2 + Σ_{n=1..N} a_n cos(n(θ − θ_x))."""
    out = np.full_like(theta_rad, fourier_coeff(0) / 2.0)
    for n in range(1, N + 1):
        out = out + fourier_coeff(n) * np.cos(n * (theta_rad - theta_x))
    return out


partial = {N: partial_sum(theta_rad, N) for N in N_list}

# ---------------------------------------------------------------------------
# figure layout
# ---------------------------------------------------------------------------
fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 6), facecolor="white")
fig.suptitle("Rectangular phase-window body and its Fourier partial sums",
             fontsize=16, fontweight="bold", y=0.98)

common_xlim   = (theta_min_deg, theta_max_deg)
common_ylim   = (-0.45, 1.45)
common_xticks = np.arange(0, 361, 45)

# =========================================================================
# LEFT PANEL : Body — rectangular window
# =========================================================================
ax_l.set_title("Body:  $R_x(\\theta)$  (rectangular phase window)",
               fontsize=13, pad=10)
ax_l.fill_between(theta_deg, R_x, color="tab:blue", alpha=0.35,
                  label=r"$R_x(\theta)$")
ax_l.plot(theta_deg, R_x, color="tab:blue", lw=2.4)

ax_l.set_xlim(*common_xlim)
ax_l.set_ylim(*common_ylim)
ax_l.set_xticks(common_xticks)
ax_l.set_xlabel(r"phase $\theta$ (degrees)", fontsize=12)
ax_l.set_ylabel(r"$R_x(\theta)$", fontsize=12)
ax_l.axhline(0, color="0.6", lw=0.6)
ax_l.grid(alpha=0.2)

# centre line and width arrow
ax_l.axvline(theta_x_deg, color="0.35", lw=1.0, linestyle=":")
ax_l.text(theta_x_deg + 4, 1.20, r"$\theta_x = 180°$",
          color="0.20", fontsize=11, ha="left", va="center")

ax_l.annotate("", xy=(theta_x_deg + delta_theta_x_deg / 2, -0.22),
              xytext=(theta_x_deg - delta_theta_x_deg / 2, -0.22),
              arrowprops=dict(arrowstyle="<->", color="tab:orange", lw=1.8))
ax_l.text(theta_x_deg, -0.32,
          r"window width $\Delta\theta_x = 120°$ (full width)",
          color="tab:orange", fontsize=11, ha="center", va="top")

# defining condition
ax_l.text(theta_min_deg + 6, 0.62,
          r"$R_x(\theta) = 1$  if  $|\theta - \theta_x| \leq \Delta\theta_x/2$",
          fontsize=10.5, ha="left", va="center", color="0.20",
          bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.95))
ax_l.text(theta_min_deg + 6, 0.42,
          r"$R_x(\theta) = 0$  otherwise",
          fontsize=10.5, ha="left", va="center", color="0.20",
          bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", alpha=0.95))

# =========================================================================
# RIGHT PANEL : Fourier partial sums
# =========================================================================
ax_r.set_title("Observation image:  Fourier partial sums  "
               r"$S_N(\theta) = \frac{a_0}{2} + \sum_{n=1}^{N} a_n \cos(n(\theta-\theta_x))$",
               fontsize=12.5, pad=10)

# ghost of the rectangular body
ax_r.fill_between(theta_deg, R_x, color="tab:blue", alpha=0.10,
                  label=r"$R_x(\theta)$  (body, ghost)")
ax_r.plot(theta_deg, R_x, color="tab:blue", lw=1.2, alpha=0.4)

# partial sums
colors    = ["tab:red", "tab:green", "tab:purple"]
linestyles = ["-", "-", "-"]
for color, ls, N in zip(colors, linestyles, N_list):
    ax_r.plot(theta_deg, partial[N], color=color, lw=2.2, linestyle=ls,
              label=fr"$S_{{{N}}}(\theta)$  (up to $n={N}$)")

ax_r.set_xlim(*common_xlim)
ax_r.set_ylim(*common_ylim)
ax_r.set_xticks(common_xticks)
ax_r.set_xlabel(r"phase $\theta$ (degrees)", fontsize=12)
ax_r.set_ylabel("amplitude", fontsize=12)
ax_r.axhline(0, color="0.6", lw=0.6)
ax_r.grid(alpha=0.2)

# centre line
ax_r.axvline(theta_x_deg, color="0.35", lw=1.0, linestyle=":")

# coefficients table — show the actual computed values
coeff_lines = [r"Fourier coefficients (centred at $\theta_x$):"]
coeff_lines.append(r"  $a_0/2 = 1/3$  (DC, = duty cycle)")
for n in [1, 2, 3, 4, 5, 6, 7]:
    a = fourier_coeff(n)
    coeff_lines.append(fr"  $a_{n} = {a:+.3f}$"
                        + (r"   ($\sin(n\pi/3)=0$)" if abs(a) < 1e-9 else ""))
ax_r.text(0.02, 0.04, "\n".join(coeff_lines),
          transform=ax_r.transAxes, fontsize=9, ha="left", va="bottom",
          color="0.20",
          bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7", alpha=0.92))

ax_r.legend(loc="upper right", fontsize=9.5, framealpha=0.92)

# ---------------------------------------------------------------------------
# transformation note between the two panels
# ---------------------------------------------------------------------------
fig.text(0.5, 0.03,
         r"$L_\Lambda$ : truncation of the Fourier series at harmonic $N$.   "
         r"Partial sums $S_1, S_3, S_7$ converge to $R_x$ with Gibbs-type ripples at the edges.",
         fontsize=11, ha="center", color="0.20", style="italic")

# ---------------------------------------------------------------------------
# save
# ---------------------------------------------------------------------------
fig.tight_layout(rect=[0, 0.06, 1, 0.94])
fig.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
print(f"saved: {out_path}")
print("First non-vanishing coefficients:")
for n in range(0, 9):
    a = fourier_coeff(n)
    print(f"  a_{n} = {a:+.4f}")
