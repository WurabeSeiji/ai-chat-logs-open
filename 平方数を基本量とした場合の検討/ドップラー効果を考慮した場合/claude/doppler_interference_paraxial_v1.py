#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Doppler-aware two-source interference of finite wave trains
 速度方向と視線のなす角 theta_v によるドップラー効果込みの干渉縞
================================================================================

正規化: lam0 = 1, c = 1, omega0 = 2 pi, source separation W0 = 5 lam0,
far-field screen (L >> W0). Each source emits an N-cycle rectangular train.

Source/emitter moves at speed beta along a direction making angle theta_v with
the line of sight (theta_v = 0 : head-on approach ; 180 : head-on recede).
RELATIVISTIC DOPPLER shifts the OBSERVED wavelength:

    lam_obs = lam0 * (1 - beta cos theta_v) / sqrt(1 - beta^2)
    W_eff   = W0 / lam_obs        (effective separation in OBSERVED wavelengths)

NOTE on beta=1 (a true photon): lam_obs -> infinity for any theta_v > 0
(relativistic beaming -> radiation only exactly forward), so nothing is seen
off-axis.  To visualise the angular dependence at finite observed frequency we
use a representative relativistic emitter speed beta (set below; adjustable).

Fringe pattern (path difference delta_lam = W_eff sin theta, in observed lam):
    I(sin theta)/I0 = 1 + V cos(2 pi delta_lam),  V = max(0, 1 - |delta_lam|/N)
Bright orders m :  sin theta = m / W_eff ,  |m| <= floor(W_eff), |m| < N.
Fringe count = 2 min(N-1, floor(W_eff)) + 1.

Author: Iris (peer-review dialogue) for N. Kihara, 2026.
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
LAM0  = 1.0
C     = 1.0
W0    = 5.0
BETA  = 0.6                       # emitter speed (representative relativistic; ADJUSTABLE)
ANGLES = [30, 60, 90, 120, 150]   # theta_v in degrees (0=approach, 180=recede)
NLIST  = list(range(1, 11))       # N = 1..10
NGRID  = 6000
OUTDIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUTDIR, exist_ok=True)

GAMMA = 1.0 / np.sqrt(1.0 - BETA**2)

# ------------------------------- physics ------------------------------------ #
def lambda_obs(theta_v_deg, beta=BETA):
    """Relativistic Doppler-shifted observed wavelength (lam0 = 1)."""
    c = np.cos(np.deg2rad(theta_v_deg))
    return (1.0 - beta * c) / np.sqrt(1.0 - beta**2)

def w_eff(theta_v_deg, beta=BETA):
    """Effective source separation in OBSERVED wavelengths."""
    return W0 / lambda_obs(theta_v_deg, beta)

def regime(theta_v_deg):
    c = np.cos(np.deg2rad(theta_v_deg))
    if c > 1e-9:
        return "approach (blueshift)"
    if c < -1e-9:
        return "recede (redshift)"
    return "transverse (relativistic redshift)"

def visibility(delta_lam, N):
    return np.maximum(0.0, 1.0 - np.abs(delta_lam) / float(N))

def intensity(sin_theta, N, weff):
    dl = weff * sin_theta
    return 1.0 + visibility(dl, N) * np.cos(2.0 * np.pi * dl)

def visible_orders(N, weff):
    mmax = int(np.floor(weff + 1e-9))
    return [m for m in range(-mmax, mmax + 1) if (1.0 - abs(m) / float(N)) > 1e-9]

def fringe_count_analytic(N, weff):
    return 2 * int(min(N - 1, np.floor(weff + 1e-9))) + 1

def fringe_count_numeric(N, weff, sin_theta, I):
    """Lattice-confirmation cross-check: each predicted order m must coincide
    with a genuine local maximum of I within +/- 0.4 fringe spacing."""
    half = 0.4 / weff
    confirmed = 0
    for m in visible_orders(N, weff):
        s0 = m / weff
        mask = (sin_theta >= s0 - half) & (sin_theta <= s0 + half)
        if not np.any(mask):
            continue
        idx = np.where(mask)[0]
        j = idx[int(np.argmax(I[mask]))]
        is_max = True if not (0 < j < len(I) - 1) else (I[j] >= I[j-1] and I[j] >= I[j+1])
        if is_max and I[j] > 1.0 + 1e-3:
            confirmed += 1
    return confirmed

# --------------------------- one detailed figure ---------------------------- #
def make_figure(N, theta_v):
    weff = w_eff(theta_v)
    lam = lambda_obs(theta_v)
    sin_theta = np.linspace(-1.0, 1.0, NGRID)
    I = intensity(sin_theta, N, weff)
    Venv = visibility(weff * sin_theta, N)
    orders = visible_orders(N, weff)
    n_an = fringe_count_analytic(N, weff)
    n_num = fringe_count_numeric(N, weff, sin_theta, I)
    assert n_an == n_num, "mismatch N=%d th=%d (%d vs %d)" % (N, theta_v, n_an, n_num)

    fig = plt.figure(figsize=(9.6, 6.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.5], hspace=0.42)

    ax0 = fig.add_subplot(gs[0])
    ax0.imshow(np.tile(I, (60, 1)), extent=[-1, 1, 0, 1], aspect="auto",
               cmap="inferno", vmin=0.0, vmax=2.0)
    ax0.set_yticks([])
    ax0.set_xlabel("sin θ")
    ax0.set_title("Fringes on screen   θv = %d°,  λ_obs = %.3f,  W_eff = %.2f"
                  % (theta_v, lam, weff))

    ax1 = fig.add_subplot(gs[1])
    ax1.plot(sin_theta, I, color="#1f77b4", lw=1.2, label="I/I₀ = 1 + V·cos(2πδλ)")
    ax1.plot(sin_theta, 1 + Venv, color="#d62728", lw=1.0, ls="--", label="envelope 1 ± V")
    ax1.plot(sin_theta, 1 - Venv, color="#d62728", lw=1.0, ls="--")
    for m in orders:
        s = m / weff
        ax1.plot([s], [intensity(np.array([s]), N, weff)[0]], "o",
                 color="black", ms=4.0, zorder=5)
    sc = min(1.0, N / weff)
    ax1.axvspan(-sc, sc, color="gold", alpha=0.10, label="coherent region")
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(0, 2.2)
    ax1.set_xlabel("sin θ")
    ax1.set_ylabel("I / I₀")
    ax1.legend(loc="upper right", fontsize=7.5, framealpha=0.92)
    info = ("fringes = %d\nλ_obs = %.3f  (×%.2f vs λ0)\nW_eff = W/λ_obs = %.2f\n%s"
            % (n_an, lam, lam, weff, regime(theta_v)))
    ax1.text(0.015, 0.96, info, transform=ax1.transAxes, va="top", ha="left",
             fontsize=8.5, bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.92))
    ax1.grid(alpha=0.25)

    fig.suptitle("Doppler interference (β=%.1f):  N = %d,  θv = %d°  →  fringe count = %d"
                 % (BETA, N, theta_v, n_an), y=0.995, fontsize=11)

    base = os.path.join(OUTDIR, "dop_a%03d_N%02d" % (theta_v, N))
    fig.savefig(base + ".png", dpi=140, bbox_inches="tight")
    fig.savefig(base + ".svg", bbox_inches="tight")
    plt.close(fig)
    return dict(N=N, theta_v=theta_v, lam=lam, weff=weff, count=n_an, num=n_num,
                orders=orders, regime=regime(theta_v))

# ------------------------------ overview grid ------------------------------- #
def make_grid(rows_by):
    nN, nA = len(NLIST), len(ANGLES)
    fig = plt.figure(figsize=(2.6 * nA, 1.05 * nN + 1.2))
    gs = GridSpec(nN, nA, figure=fig, hspace=0.18, wspace=0.08,
                  left=0.06, right=0.995, top=0.93, bottom=0.04)
    sin_theta = np.linspace(-1.0, 1.0, 2000)
    for i, N in enumerate(NLIST):
        for j, th in enumerate(ANGLES):
            weff = w_eff(th)
            I = intensity(sin_theta, N, weff)
            ax = fig.add_subplot(gs[i, j])
            ax.imshow(np.tile(I, (10, 1)), extent=[-1, 1, 0, 1], aspect="auto",
                      cmap="inferno", vmin=0.0, vmax=2.0)
            ax.set_xticks([]); ax.set_yticks([])
            cnt = fringe_count_analytic(N, weff)
            ax.text(0.5, 0.5, str(cnt), color="white", fontsize=9, fontweight="bold",
                    ha="center", va="center", transform=ax.transAxes,
                    path_effects=None)
            if i == 0:
                ax.set_title("θv=%d°\nλ=%.2f, W=%.1f" % (th, lambda_obs(th), weff),
                             fontsize=9)
            if j == 0:
                ax.set_ylabel("N=%d" % N, fontsize=9, rotation=0, ha="right", va="center",
                              labelpad=14)
    fig.suptitle("Doppler interference fringe map (β=%.1f): rows N=1..10, cols θv "
                 "(number = fringe count)" % BETA, fontsize=12, y=0.985)
    fig.savefig(os.path.join(OUTDIR, "dop_grid_strips.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, "dop_grid_strips.svg"), bbox_inches="tight")
    plt.close(fig)

# ----------------------------- summary counts ------------------------------- #
def make_summary(rows):
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    colors = ["#0b3d91", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    for k, th in enumerate(ANGLES):
        ys = [r["count"] for r in rows if r["theta_v"] == th]
        weff = w_eff(th)
        ax.plot(NLIST, ys, "o-", color=colors[k], lw=1.7, ms=6,
                label="θv=%d°  (λ_obs=%.2f, ⌊W_eff⌋=%d)"
                      % (th, lambda_obs(th), int(np.floor(weff + 1e-9))))
        ax.axhline(2 * int(np.floor(weff + 1e-9)) + 1, color=colors[k], ls=":", lw=0.9, alpha=0.6)
    ax.set_xlabel("vibration interval  N  (cycles)")
    ax.set_ylabel("number of visible fringes")
    ax.set_title("Doppler effect on fringe count (β=%.1f):  approach → more, recede → fewer"
                 % BETA)
    ax.set_xticks(NLIST)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5, title="dotted = geometric saturation")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "dop_summary_counts.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, "dop_summary_counts.svg"), bbox_inches="tight")
    plt.close(fig)

# --------------------------------- main ------------------------------------- #
def main():
    rows = []
    for th in ANGLES:
        for N in NLIST:
            rows.append(make_figure(N, th))
    make_grid(rows)
    make_summary(rows)

    with open(os.path.join(OUTDIR, "dop_fringe_counts.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["beta", "theta_v_deg", "regime", "lambda_obs", "W_eff",
                    "floor_Weff", "N", "fringe_count", "fringe_count_numeric",
                    "visible_orders_m"])
        for r in rows:
            w.writerow([BETA, r["theta_v"], r["regime"], round(r["lam"], 4),
                        round(r["weff"], 4), int(np.floor(r["weff"] + 1e-9)),
                        r["N"], r["count"], r["num"],
                        " ".join(str(m) for m in r["orders"])])

    # console summary table : rows N, cols angle
    print("beta = %.2f,  gamma = %.4f" % (BETA, GAMMA))
    print("theta_v :   " + "".join("%8d°" % a for a in ANGLES))
    print("lam_obs :   " + "".join("%9.3f" % lambda_obs(a) for a in ANGLES))
    print("W_eff   :   " + "".join("%9.2f" % w_eff(a) for a in ANGLES))
    print("floorW  :   " + "".join("%9d" % int(np.floor(w_eff(a) + 1e-9)) for a in ANGLES))
    print("-" * 60)
    print(" N \\ th  " + "".join("%9d" % a for a in ANGLES))
    for N in NLIST:
        line = "  %2d     " % N
        for a in ANGLES:
            line += "%9d" % fringe_count_analytic(N, w_eff(a))
        print(line)
    print("Outputs ->", OUTDIR)

if __name__ == "__main__":
    main()
