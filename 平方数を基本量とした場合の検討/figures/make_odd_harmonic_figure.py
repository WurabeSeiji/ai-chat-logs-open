#!/usr/bin/env python3
"""Generate the figure for the odd-harmonic localization observation paper.

Plots the normalized squared amplitude

    I_N(phi) / I_N(0) = [ sin((N+1) phi) / ((N+1) sin phi) ]^2

of the equal-amplitude odd-harmonic sum

    S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi)
             = sin((N+1) phi) / (2 sin phi)

for N = 99, 999, 9999, on a phase axis expressed as a percentage of the
full period [-pi, pi]. Labels are in English so the same figure can be
reused in an English manuscript. Writes PNG and SVG.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent


def normalized_intensity(N, phi):
    """I_N(phi)/I_N(0) via the closed-form odd-harmonic Dirichlet-type kernel.

    The kernel sin((N+1)phi)/((N+1) sin phi) is a removable 0/0 wherever
    sin(phi) = 0, i.e. at phi = 0 and phi = +/- pi. All of these are equal-
    height peaks of the normalized intensity (value 1), so they are set
    explicitly to avoid floating-point blow-up at the +/- pi endpoints.
    """
    phi = np.asarray(phi, dtype=float)
    out = np.empty_like(phi)
    small = np.abs(np.sin(phi)) < 1e-9
    out[small] = 1.0
    p = phi[~small]
    out[~small] = (np.sin((N + 1) * p) / ((N + 1) * np.sin(p))) ** 2
    return out


def main():
    # Phase axis in percent of the full period [-pi, pi]: 100% <-> 2*pi.
    pct = np.linspace(-50.0, 50.0, 600001)
    phi = pct / 100.0 * (2.0 * np.pi)

    cases = [
        (99, "tab:red", "up to the 99th harmonic (N = 99)"),
        (999, "tab:blue", "up to the 999th harmonic (N = 999)"),
        (9999, "tab:green", "up to the 9999th harmonic (N = 9999)"),
    ]

    fig, ax = plt.subplots(figsize=(11, 5.4))
    for N, color, label in cases:
        ax.plot(pct, normalized_intensity(N, phi), color=color, lw=1.1, label=label)

    ax.set_xlim(-50, 50)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("phase  (% of full period $[-\\pi,\\pi]$)")
    ax.set_ylabel("squared amplitude  (normalized to peak = 1)")
    ax.set_title(
        "Figure 1.  Localization of the equal-amplitude odd-harmonic sum\n"
        "$S_N(\\varphi)=\\sum_{m=0}^{(N-1)/2}\\cos((2m+1)\\varphi)"
        "=\\dfrac{\\sin((N+1)\\varphi)}{2\\sin\\varphi}$",
        fontsize=12,
    )
    ax.grid(True, which="both", color="0.85", lw=0.6)
    ax.legend(loc="upper right", framealpha=0.95)

    for ext in ("png", "svg"):
        fig.savefig(OUT / f"fig01_odd_harmonic_localization.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
