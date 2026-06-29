#!/usr/bin/env python3
"""
EXACT two-slit interference of a localized odd-harmonic wave -- parameterized.

Pass L, W, lambda0 (cage = fundamental wavelength) and N on the command line:
    python3 fig_oddharm_interference.py --L 240 --W 140 --lam0 250 --N 99
Optional: --halfdeg <x>  (x-axis half range in lambda0-based phase; default 720).

Physics (nothing dropped):
  source = equal-amplitude odd harmonics, lambda_n = lambda0 / n, n=1,3,...,N
    (n=1 is the fundamental: lambda_1 = lambda0)
  source centred at y=0, distance L from the slit plane, slits at y = +-W/2
  source->slit distance r_k = sqrt(L^2 + (W/2)^2)  (same for both slits)
  Phi_k^(n) = (2 pi n / lambda0)(r_k - y_slit,k s)
  psi(s) = sum_n [exp(i Phi1^(n)) + exp(i Phi2^(n))],   I = |psi|^2

Horizontal axis = lambda0-based phase theta = pi W s / lambda0 (deg). When aligned
psi = 2 S_N(theta): localized fringes sit at theta = 0, +-180, +-360 deg; one fringe
fills the basic domain theta in [-90, 90] deg.

Single central fringe <=> r_k / lambda0 = integer or half-integer
(the slit sits on a spike of the localized wave). Otherwise the pattern scrambles.

Output (filename carries the parameters):
    fig_oddharm_interference_L<L>_W<W>_lam<lam0>_N<N>.png / .svg
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})


def fmt(x):
    """Compact value for filenames: 240 -> '240', 27.5 -> '27p5'."""
    return str(int(x)) if float(x).is_integer() else str(x).replace(".", "p")


def decomposition_mode(L, W, lam0, N, n_list, ysl1, ysl2, dlam, ngrid, halfdeg):
    """Source-position fluctuation (Paper-1 fig_decomposition_static, generalized
    to the localized odd-harmonic source). Normalisation: case-1 = by the y=0 peak.

    Given dlam (full width), internal half-range = dlam/2. Source positions
    y = x*(dlam/2), x in [-0.8, 0.8] (Paper-1 grid), weight P(y) = cos^2(pi x/2).
    Each exact localized two-slit intensity I(theta;y) (all N odd harmonics,
    off-axis r1(y) != r2(y), full source phase) is scaled by P(y) and divided by
    the y=0 peak. Fringe peaks are compared with the cos^2 envelope (push-forward
    of the source distribution). For N=1 this reduces EXACTLY to Paper 1.

    x-axis: two-slit reference phase Phi0 = 2 pi W s / lam0 (deg), as in Paper 1.
    """
    half_src = dlam / 2.0
    x_norm = np.linspace(-0.8, 0.8, ngrid)
    y_src = x_norm * half_src
    w = np.cos(np.pi * x_norm / 2.0)**2

    Phi0_deg = np.linspace(-halfdeg, halfdeg, 16001)
    s = np.radians(Phi0_deg) / (2.0 * np.pi * W / lam0)        # Phi0 = 2 pi W s / lam0

    def loc_I(yi):
        r1 = np.sqrt(L**2 + (yi - ysl1)**2)
        r2 = np.sqrt(L**2 + (yi - ysl2)**2)
        psi = np.zeros_like(s, dtype=complex)
        for nn in n_list:
            k = 2.0 * np.pi * nn / lam0
            psi += np.exp(1j * k * (r1 - ysl1 * s)) + np.exp(1j * k * (r2 - ysl2 * s))
        return (psi * np.conj(psi)).real

    peak0 = loc_I(0.0).max()                                   # case-1 reference

    fig, ax = plt.subplots(figsize=(13.0, 5.8))
    env_x = np.linspace(-90.0, 90.0, 600)
    env_y = np.cos(np.pi * env_x / 180.0)**2
    ax.plot(env_x, env_y, color="#f0a500", lw=3.2, zorder=5,
            label=r"$\cos^2$ envelope = push-forward of source distribution")

    hdr = "  x_i     peak Phi0 (exact)   peak height   cos^2 at that Phi0   dev"
    print(hdr)
    rows = [hdr]
    peaks_x, peaks_y = [], []
    for x_i, y_i, w_i in zip(x_norm, y_src, w):
        curve = w_i * loc_I(y_i) / peak0                      # case-1: /peak0
        if abs(x_i) < 1e-12:
            ax.plot(Phi0_deg, curve, color="#1f5fbf", lw=2.3, zorder=4,
                    label=r"$x=0$ (centred source, peak $1$)")
        else:
            ax.plot(Phi0_deg, curve, color="#2ca02c", lw=1.0, alpha=0.85, zorder=2)
        m = np.abs(Phi0_deg) < 180.0
        idx = np.argmax(curve[m])
        px, py = Phi0_deg[m][idx], curve[m][idx]
        peaks_x.append(px); peaks_y.append(py)
        env_at = np.cos(np.pi * px / 180.0)**2
        line = (f"{x_i:+.1f}   {px:+8.2f}            {py:.4f}        {env_at:.4f}"
                f"            {py - env_at:+.4f}")
        print(line); rows.append(line)

    ax.plot(peaks_x, peaks_y, 'o', color="#d11f2d", ms=5.5, zorder=6,
            label="fringe peaks (exact)")
    ax.plot([], [], color="#2ca02c", lw=1.2,
            label=r"each source $x_i$ (green), height $P(y_i)\cdot I/I_{y=0}^{\rm peak}$")

    for k in range(-2, 3):
        ax.axvline(360 * k, color="0.92", lw=0.8, zorder=0)
    ax.axhline(0, color="0.92", lw=0.8, zorder=0)
    ax.set_xlim(-halfdeg, halfdeg)
    ax.set_ylim(-0.03, 1.08)
    step = next((c for c in [180, 360, 720, 1170] if halfdeg / c <= 4), 2340)
    ax.set_xticks(np.arange(-halfdeg, halfdeg + 1e-6, step))
    ax.set_xlabel(r"Two-slit reference phase $\Phi_0 = 2\pi W s/\lambda_0$  (degrees)")
    ax.set_ylabel(r"Intensity $|\psi|^2$  (normalised to $y{=}0$ peak)")
    ax.legend(loc="upper right", fontsize=9.0, framealpha=0.96)
    cap = (rf"Source-position fluctuation. $L={fmt(L)}$, $W={fmt(W)}$, "
           rf"$\lambda_0={fmt(lam0)}$, $N={N}$ ({len(n_list)} harmonics); "
           rf"$\Delta\lambda={fmt(dlam)}$ (half-range ${fmt(half_src)}$). Each source "
           r"$y=x\,\Delta\lambda/2$ ($x\in[-0.8,0.8]$) gives its EXACT localized two-slit "
           r"intensity (all harmonics, off-axis $r_1\neq r_2$), scaled by $P(y)=\cos^2(\pi x/2)$ "
           r"and divided by the $y{=}0$ peak (case-1). Peaks (red) vs $\cos^2$ push-forward "
           r"(yellow). For $N=1$ this is exactly Paper 1.")
    fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.6, wrap=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.17)

    outdir = os.path.dirname(os.path.abspath(__file__))
    base = (f"fig_oddharm_fluct_L{fmt(L)}_W{fmt(W)}_lam{fmt(lam0)}_N{N}_dlam{fmt(dlam)}")
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{base}.{ext}"), dpi=200)
    print(f"\npeak0(y=0) = {peak0:.4f}   saved {base}.png / .svg")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=240.0, help="source-to-slit-plane distance")
    ap.add_argument("--W", type=float, default=140.0, help="slit separation")
    ap.add_argument("--lam0", type=float, default=250.0, help="cage = fundamental wavelength (=lambda at n=1)")
    ap.add_argument("--N", type=int, default=99, help="highest odd harmonic order")
    ap.add_argument("--halfdeg", type=float, default=720.0, help="x half-range in lambda0-based phase (deg)")
    ap.add_argument("--dlam", type=float, default=0.0,
                    help="source-position fluctuation FULL width; internal half-range = dlam/2 "
                         "(dlam=0: fixed source; dlam=lam0 reproduces Paper-1 +-lambda/2)")
    ap.add_argument("--ngrid", type=int, default=9, help="number of source-position grid points (dlam>0)")
    args = ap.parse_args()

    L, W, lam0, N, half = args.L, args.W, args.lam0, args.N, args.halfdeg
    ysl1, ysl2 = +W / 2.0, -W / 2.0
    y = 0.0
    n_list = np.arange(1, N + 1, 2)

    if args.dlam > 0:
        decomposition_mode(L, W, lam0, N, n_list, ysl1, ysl2, args.dlam, args.ngrid, half)
        return

    r1 = np.sqrt(L**2 + (y - ysl1)**2)
    r2 = np.sqrt(L**2 + (y - ysl2)**2)
    ratio = r1 / lam0
    align = np.abs(np.sum(np.exp(1j * 2 * np.pi * n_list * r1 / lam0)))**2 / len(n_list)**2
    # aligned within the practical tolerance band ~1/(2N) around the nearest
    # integer/half-integer of r_k/lambda0 (band narrows as N grows).
    dev = abs(ratio - round(2 * ratio) / 2.0)
    tol = 1.0 / (2.0 * N)
    aligned = dev < tol

    # horizontal axis: lambda0-based phase theta = pi W s / lam0  (degrees)
    theta_deg = np.linspace(-half, half, 80001)
    theta = np.radians(theta_deg)
    s = theta * lam0 / (np.pi * W)                          # invert theta = pi W s / lam0

    psi = np.zeros_like(s, dtype=complex)
    for n in n_list:
        k = 2.0 * np.pi * n / lam0
        psi += np.exp(1j * k * (r1 - ysl1 * s)) + np.exp(1j * k * (r2 - ysl2 * s))
    I = (psi * np.conj(psi)).real
    I_norm = I / I.max()

    # ideal aligned reference: 4 S_N(theta)^2
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.where(np.abs(np.sin(theta)) < 1e-12, (N + 1) / 2.0,
                     np.sin((N + 1) * theta) / (2.0 * np.sin(theta)))
    I_ideal = (S / S.max())**2

    n_waves = 2 * len(n_list)
    n_pairs = n_waves * (n_waves - 1) // 2
    loc = (I_norm[1:-1] > I_norm[:-2]) & (I_norm[1:-1] > I_norm[2:]) & (I_norm[1:-1] >= 0.5)
    n_peaks = int(loc.sum())
    ic = I_norm[np.argmin(np.abs(theta_deg))]

    fig, ax = plt.subplots(figsize=(13.5, 5.6))
    ax.plot(theta_deg, I_norm, color="#1f5fbf", lw=1.5, zorder=3,
            label=rf"{len(n_list)} odd harmonics $\times$ 2 slits $=$ {n_waves} waves "
                  rf"(exact, {n_pairs} pairs)")
    if aligned:
        ax.plot(theta_deg, I_ideal, color="#f0a500", lw=1.0, ls="--", zorder=4,
                label=r"ideal aligned $4\,S_N(\theta)^2$")
    # basic-domain boundaries (odd multiples of 90 deg) and fringe centres (180 deg)
    for x0 in np.arange(-half, half + 1, 180.0):
        ax.axvline(x0, color="0.90", lw=0.8, zorder=0)
    for x0 in np.arange(-half + 90, half, 180.0):
        ax.axvline(x0, color="0.95", lw=0.7, ls=":", zorder=0)
    ax.axhline(0, color="0.92", lw=0.8, zorder=0)
    ax.set_xlim(-half, half)
    ax.set_ylim(-0.03, 1.08)
    step = next((c for c in [5, 10, 15, 30, 45, 60, 90, 180, 360, 720, 1170]
                 if half / c <= 4), 2340)
    ax.set_xticks(np.arange(-half, half + 1e-6, step))
    ax.set_xlabel(r"$\lambda_0$-based phase  $\theta=\pi W s/\lambda_0$  (degrees;  "
                  r"basic domain $[-90^\circ,90^\circ]$, fringes every $180^\circ$)")
    ax.set_ylabel(r"Intensity $|\psi|^2$  (peak-normalised)")
    ax.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
    state = (rf"aligned (single fringe per domain): $r_k/\lambda_0={ratio:.4f}$, "
             rf"dev$={dev:.4f}<1/(2N)={tol:.4f}$"
             if aligned else
             rf"NOT aligned (scrambled): $r_k/\lambda_0={ratio:.4f}$, "
             rf"dev$={dev:.4f}\geq1/(2N)={tol:.4f}$")
    cap = (rf"EXACT two-slit interference. $L={fmt(L)}$, $W={fmt(W)}$, "
           rf"$\lambda_0={fmt(lam0)}$ ($=\lambda$ at $n{{=}}1$), $N={N}$ "
           rf"({len(n_list)} odd harmonics, $\lambda_n=\lambda_0/n$). "
           rf"$r_k=\sqrt{{L^2+(W/2)^2}}={r1:.3f}$; source phase, $L$ and slit offset all kept. "
           rf"{state}; alignment$={align:.3f}$, $I(0)={ic:.3f}$, {n_peaks} peak(s)$\geq0.5$.")
    fig.text(0.5, 0.005, cap, ha="center", va="bottom", fontsize=8.4, wrap=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.20)

    outdir = os.path.dirname(os.path.abspath(__file__))
    base = f"fig_oddharm_interference_L{fmt(L)}_W{fmt(W)}_lam{fmt(lam0)}_N{N}_pm{fmt(half)}"
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{base}.{ext}"), dpi=200)

    print(f"L={fmt(L)} W={fmt(W)} lam0={fmt(lam0)} N={N}")
    print(f"r_k={r1:.4f}  r_k/lam0={ratio:.4f}  aligned={aligned}  alignment={align:.4f}")
    print(f"I(0)={ic:.4f}  peaks>=0.5: {n_peaks}  theta-range=+-{half:.0f} deg")
    print(f"saved {base}.png / .svg")


if __name__ == "__main__":
    main()
