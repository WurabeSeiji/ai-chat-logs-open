#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Doppler-aware two-source interference of finite wave trains  (RIGOROUS v2)
 角度依存ドップラーを織り込んだ二光源干渉  ―  λ_obs(θ) 版
================================================================================

正規化: lam0 = 1, c = 1, omega0 = 2 pi, source separation W0 = 5 lam0,
far-field screen.  Each source emits an N-cycle rectangular train.

GEOMETRY (made explicit after peer review):
    - baseline  d  (|d| = W0)  is PERPENDICULAR to the velocity  beta
      => no Lorentz contraction of the baseline, path difference = W0 sin(theta).
    - velocity makes angle theta_v with the screen normal z, in the same plane.
    - theta is the LAB-FRAME diffraction angle (screen position), s = sin(theta).

RELATIVISTIC DOPPLER is DIRECTION DEPENDENT (this is the v2 correction).
The light reaching screen angle theta has observed wavelength

    lam_obs(theta) = lam0 * gamma * ( 1 - beta cos(theta - theta_v) )

Writing it in s = sin(theta)  (cos(theta) = +sqrt(1-s^2), theta in (-90,90)):

    cos(theta - theta_v) = sqrt(1-s^2) cos(theta_v) + s sin(theta_v)
    lam_obs(s) = gamma ( 1 - beta ( sqrt(1-s^2) cos(theta_v) + s sin(theta_v) ) )

At s = 0 (screen centre) this reduces to the OLD single value
    lam_obs(0) = gamma (1 - beta cos(theta_v))          <-- paraxial approximation.
The v1 script applied that central value over the whole screen; that is only
valid for |theta| << (1-beta cos th_v)/(beta sin th_v).  v2 keeps the full s
dependence, so the fringes are chirped (unequal spacing) and left/right
ASYMMETRIC, and a single number lam_obs no longer describes the screen.

Path difference in observed wavelengths:
    delta(s) = W0 sin(theta) / lam_obs(s) = W0 s / lam_obs(s)
Intensity (finite N-cycle trains, triangular coherence envelope):
    I(s)/I0 = 1 + V(delta) cos(2 pi delta),   V = max(0, 1 - |delta|/N)
Bright fringe  <=>  delta(s) = integer m  with  |m| < N  (coherence)
                    and a real s in (-1,1) exists (geometry).
No closed form: delta(s)=m is TRANSCENDENTAL and is solved numerically.
delta(s) can be NON-monotonic (fold-back near grazing), giving extra fringes
beyond the naive integer-in-range count; all such genuine maxima are counted.

Author: Iris (peer-review dialogue) for N. Kihara, 2026.  Rigorous v2.
================================================================================
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ------------------------------ parameters ---------------------------------- #
LAM0   = 1.0
C      = 1.0
W0     = 5.0
BETA   = 0.6                       # emitter speed (representative relativistic)
ANGLES = [30, 60, 90, 120, 150]    # theta_v in degrees (0=approach, 180=recede)
NLIST  = list(range(1, 11))        # N = 1..10
NGRID  = 30000                     # fine grid: chirp + fold-back near grazing
S_MAX  = 0.99995                   # avoid exact grazing s = +/-1
OUTDIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUTDIR, exist_ok=True)

GAMMA = 1.0 / np.sqrt(1.0 - BETA**2)

# ------------------------------- physics ------------------------------------ #
def lam_obs_s(s, theta_v_deg, beta=BETA):
    """Direction-dependent observed wavelength as a function of s = sin(theta)."""
    tv = np.deg2rad(theta_v_deg)
    cosk = np.sqrt(np.maximum(0.0, 1.0 - s * s)) * np.cos(tv) + s * np.sin(tv)
    return GAMMA * (1.0 - beta * cosk)

def lam_obs_par(theta_v_deg, beta=BETA):
    """Paraxial (screen-centre) value used by the v1 script: s = 0."""
    return GAMMA * (1.0 - beta * np.cos(np.deg2rad(theta_v_deg)))

def regime(theta_v_deg):
    c = np.cos(np.deg2rad(theta_v_deg))
    if c > 1e-9:
        return "approach (blueshift at centre)"
    if c < -1e-9:
        return "recede (redshift at centre)"
    return "transverse (relativistic redshift at centre)"

def visibility(delta, N):
    return np.maximum(0.0, 1.0 - np.abs(delta) / float(N))

def delta_rig(s, theta_v):
    return W0 * s / lam_obs_s(s, theta_v)

def delta_par(s, theta_v):
    return (W0 / lam_obs_par(theta_v)) * s

def intensity(delta, N):
    return 1.0 + visibility(delta, N) * np.cos(2.0 * np.pi * delta)

def bright_orders(d, s, N):
    """Every screen point where delta(s)=m (integer, |m|<N) is a bright maximum,
    since cos(2 pi m)=+1.  Found as sign changes of delta-m (handles fold-back:
    a non-monotonic delta can cross the same m twice => two genuine fringes)."""
    out = []
    for m in range(-(N - 1), N):          # |m| < N  (coherence gate)
        g = d - m
        prod = g[:-1] * g[1:]
        cross = np.where(prod < 0.0)[0]
        for i in cross:
            s0 = s[i] - g[i] * (s[i + 1] - s[i]) / (g[i + 1] - g[i])
            out.append((m, float(s0)))
    return out

def count_rig(N, theta_v, s=None, d=None):
    if s is None:
        s = np.linspace(-S_MAX, S_MAX, NGRID)
        d = delta_rig(s, theta_v)
    return len(bright_orders(d, s, N))

def count_par(N, theta_v):
    weff = W0 / lam_obs_par(theta_v)
    return 2 * int(min(N - 1, np.floor(weff + 1e-9))) + 1

# --------------------------- one detailed figure ---------------------------- #
def make_figure(N, theta_v):
    s = np.linspace(-S_MAX, S_MAX, NGRID)
    lam_r = lam_obs_s(s, theta_v)
    d_r = W0 * s / lam_r
    V_r = visibility(d_r, N)
    I_r = 1.0 + V_r * np.cos(2.0 * np.pi * d_r)

    lam_p = lam_obs_par(theta_v)
    d_p = (W0 / lam_p) * s
    I_p = 1.0 + visibility(d_p, N) * np.cos(2.0 * np.pi * d_p)

    orders = bright_orders(d_r, s, N)
    n_r = len(orders)
    n_p = count_par(N, theta_v)

    lam_c = lam_obs_par(theta_v)            # centre
    lam_pl = lam_obs_s(S_MAX, theta_v)      # +grazing (velocity side)
    lam_mn = lam_obs_s(-S_MAX, theta_v)     # -grazing

    fig = plt.figure(figsize=(9.8, 6.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.6], hspace=0.42)

    ax0 = fig.add_subplot(gs[0])
    ax0.imshow(np.tile(I_r, (60, 1)), extent=[-1, 1, 0, 1], aspect="auto",
               cmap="inferno", vmin=0.0, vmax=2.0)
    ax0.set_yticks([])
    ax0.set_xlabel("sin θ   (lab-frame screen position)")
    ax0.set_title("Rigorous fringes   θv = %d°   "
                  "λ_obs:  −graz %.3f  |  centre %.3f  |  +graz %.3f"
                  % (theta_v, lam_mn, lam_c, lam_pl))

    ax1 = fig.add_subplot(gs[1])
    ax1.plot(s, I_r, color="#1f77b4", lw=1.3,
             label="rigorous  I/I₀ = 1 + V·cos(2πδ(s)),  δ=W·s/λ_obs(s)")
    ax1.plot(s, I_p, color="#7f7f7f", lw=1.0, ls="-", alpha=0.55,
             label="paraxial (v1: λ_obs=const=centre)")
    ax1.plot(s,  1 + V_r, color="#d62728", lw=1.0, ls="--", label="rigorous envelope 1 ± V(s)")
    ax1.plot(s,  1 - V_r, color="#d62728", lw=1.0, ls="--")
    for (m, s0) in orders:
        ax1.plot([s0], [2.0 - abs(m) / float(N)], "o", color="black", ms=4.0, zorder=5)
    # rigorous coherent region: |delta(s)| < N
    coh = np.abs(d_r) < N
    if np.any(coh):
        ax1.fill_between(s, 0, 2.2, where=coh, color="gold", alpha=0.10,
                         label="coherent region |δ| < N")
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(0, 2.2)
    ax1.set_xlabel("sin θ")
    ax1.set_ylabel("I / I₀")
    ax1.legend(loc="upper right", fontsize=7.3, framealpha=0.92)
    info = ("rigorous fringes = %d\nparaxial (v1) = %d\n"
            "λ_obs  −graz/centre/+graz\n  %.3f / %.3f / %.3f\n%s"
            % (n_r, n_p, lam_mn, lam_c, lam_pl, regime(theta_v)))
    ax1.text(0.015, 0.96, info, transform=ax1.transAxes, va="top", ha="left",
             fontsize=8.3, bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.92))
    ax1.grid(alpha=0.25)

    fig.suptitle("Doppler interference β=%.1f (RIGOROUS, baseline ⊥ v):  "
                 "N = %d,  θv = %d°  →  fringes = %d  (paraxial %d)"
                 % (BETA, N, theta_v, n_r, n_p), y=0.995, fontsize=10.5)

    base = os.path.join(OUTDIR, "dop_a%03d_N%02d" % (theta_v, N))
    fig.savefig(base + ".png", dpi=140, bbox_inches="tight")
    fig.savefig(base + ".svg", bbox_inches="tight")
    plt.close(fig)
    return dict(N=N, theta_v=theta_v, lam_c=lam_c, lam_mn=lam_mn, lam_pl=lam_pl,
                n_r=n_r, n_p=n_p, orders=[m for (m, _) in orders],
                regime=regime(theta_v))

# ------------------------------ overview grid ------------------------------- #
def make_grid(rows):
    nN, nA = len(NLIST), len(ANGLES)
    fig = plt.figure(figsize=(2.7 * nA, 1.05 * nN + 1.2))
    gs = GridSpec(nN, nA, figure=fig, hspace=0.18, wspace=0.08,
                  left=0.07, right=0.995, top=0.92, bottom=0.04)
    s = np.linspace(-S_MAX, S_MAX, 4000)
    by = {(r["N"], r["theta_v"]): r for r in rows}
    for i, N in enumerate(NLIST):
        for j, th in enumerate(ANGLES):
            d = delta_rig(s, th)
            I = intensity(d, N)
            ax = fig.add_subplot(gs[i, j])
            ax.imshow(np.tile(I, (10, 1)), extent=[-1, 1, 0, 1], aspect="auto",
                      cmap="inferno", vmin=0.0, vmax=2.0)
            ax.set_xticks([]); ax.set_yticks([])
            cnt = by[(N, th)]["n_r"]
            ax.text(0.5, 0.5, str(cnt), color="white", fontsize=9, fontweight="bold",
                    ha="center", va="center", transform=ax.transAxes)
            if i == 0:
                ax.set_title("θv=%d°\nλc=%.2f" % (th, lam_obs_par(th)), fontsize=9)
            if j == 0:
                ax.set_ylabel("N=%d" % N, fontsize=9, rotation=0, ha="right",
                              va="center", labelpad=14)
    fig.suptitle("Rigorous Doppler fringe map (β=%.1f, baseline ⊥ v): rows N=1..10, "
                 "cols θv  (number = rigorous fringe count)" % BETA, fontsize=12, y=0.985)
    fig.savefig(os.path.join(OUTDIR, "dop_grid_strips.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, "dop_grid_strips.svg"), bbox_inches="tight")
    plt.close(fig)

# ----------------------------- summary counts ------------------------------- #
def make_summary(rows):
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    colors = ["#0b3d91", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    for k, th in enumerate(ANGLES):
        yr = [r["n_r"] for r in rows if r["theta_v"] == th]
        yp = [r["n_p"] for r in rows if r["theta_v"] == th]
        ax.plot(NLIST, yr, "o-", color=colors[k], lw=1.8, ms=6,
                label="θv=%d° rigorous" % th)
        ax.plot(NLIST, yp, "x--", color=colors[k], lw=1.0, ms=5, alpha=0.6,
                label="θv=%d° paraxial" % th)
    ax.set_xlabel("vibration interval  N  (cycles)")
    ax.set_ylabel("number of visible fringes")
    ax.set_title("Doppler effect on fringe count (β=%.1f): rigorous (solid) vs "
                 "paraxial v1 (dashed)" % BETA)
    ax.set_xticks(NLIST)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7.5, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "dop_summary_counts.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, "dop_summary_counts.svg"), bbox_inches="tight")
    plt.close(fig)

# --------------- where the uniform (paraxial) approximation breaks ----------- #
def make_lambda_profiles():
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    colors = ["#0b3d91", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    s = np.linspace(-S_MAX, S_MAX, 2000)
    for k, th in enumerate(ANGLES):
        ax.plot(s, lam_obs_s(s, th), color=colors[k], lw=1.8,
                label="θv=%d° rigorous λ_obs(s)" % th)
        ax.axhline(lam_obs_par(th), color=colors[k], ls=":", lw=1.0, alpha=0.7)
    ax.axhline(1.0, color="0.3", lw=0.8)
    ax.text(-0.98, 1.02, "λ0 = 1 (no shift)", fontsize=8, color="0.3")
    ax.set_xlabel("sin θ   (screen position)")
    ax.set_ylabel("observed wavelength  λ_obs")
    ax.set_title("Where the uniform approximation breaks (β=%.1f):\n"
                 "solid = direction-dependent λ_obs(s),  dotted = paraxial centre value"
                 % BETA)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.0)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "dop_lambda_profiles.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, "dop_lambda_profiles.svg"), bbox_inches="tight")
    plt.close(fig)

# --------------------------------- main ------------------------------------- #
def main():
    rows = []
    for th in ANGLES:
        for N in NLIST:
            rows.append(make_figure(N, th))
    make_grid(rows)
    make_summary(rows)
    make_lambda_profiles()

    with open(os.path.join(OUTDIR, "dop_fringe_counts.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["beta", "theta_v_deg", "regime",
                    "lambda_obs_minus_graze", "lambda_obs_centre", "lambda_obs_plus_graze",
                    "N", "fringe_count_rigorous", "fringe_count_paraxial",
                    "visible_orders_m"])
        for r in rows:
            w.writerow([BETA, r["theta_v"], r["regime"],
                        round(r["lam_mn"], 4), round(r["lam_c"], 4), round(r["lam_pl"], 4),
                        r["N"], r["n_r"], r["n_p"],
                        " ".join(str(m) for m in sorted(set(r["orders"])))])

    by = {(r["N"], r["theta_v"]): r for r in rows}
    print("beta = %.2f,  gamma = %.4f,  baseline PERP v" % (BETA, GAMMA))
    print("theta_v        :" + "".join("%9d" % a for a in ANGLES))
    print("lam_obs centre :" + "".join("%9.3f" % lam_obs_par(a) for a in ANGLES))
    print("lam_obs -graze :" + "".join("%9.3f" % lam_obs_s(-S_MAX, a) for a in ANGLES))
    print("lam_obs +graze :" + "".join("%9.3f" % lam_obs_s(S_MAX, a) for a in ANGLES))
    print("-" * 62)
    print("RIGOROUS fringe count   N \\ theta_v")
    print("        " + "".join("%9d" % a for a in ANGLES))
    for N in NLIST:
        print("  N=%2d  " % N + "".join("%9d" % by[(N, a)]["n_r"] for a in ANGLES))
    print("-" * 62)
    print("PARAXIAL (v1) fringe count for comparison")
    print("        " + "".join("%9d" % a for a in ANGLES))
    for N in NLIST:
        print("  N=%2d  " % N + "".join("%9d" % by[(N, a)]["n_p"] for a in ANGLES))
    print("Outputs ->", OUTDIR)

if __name__ == "__main__":
    main()
