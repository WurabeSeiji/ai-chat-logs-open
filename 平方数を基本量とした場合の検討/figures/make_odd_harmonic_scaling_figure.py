#!/usr/bin/env python3
"""Three-panel figure demonstrating the 1/N localization -- HONEST overlay.

ALL THREE curves (N = 99 red, N = 999 blue, N = 9999 green) are drawn together
in EVERY panel; only the horizontal zoom changes between panels:

    panel 1: phase +/-10 %
    panel 2: phase +/- 1 %
    panel 3: phase +/-0.1 %

The plotted quantity is the normalized squared amplitude

    I_N(phi) / max I_N = | S_N(phi) |^2 / ((N+1)/2)^2 ,
    S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi) = sin((N+1) phi)/(2 sin phi)

with phase in percent, 100% = 180 deg (half-wavelength). Because the width
scales as 1/(N+1) with N+1 = 100, 1000, 10000, in each panel exactly ONE curve
has its main lobe fill the window while the others are 10x wider (a broad,
nearly flat top) or 10x narrower (a thin central spike). Overlaying them is the
honest demonstration of W ~ 1/N: NOT three separately-tuned plots that
trivially look identical, but the same three curves seen at three zooms.

The vertical scale is the same in every panel (0..1.02, peak = 1.0); curve
parts outside the panel's x-range are clipped. The amplitude is evaluated and
squared honestly (no shortcut); the closed form is verified against the direct
term-by-term sum to machine precision first.

Labels are English so the figure drops into an English manuscript. PNG + SVG.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).resolve().parent


def amplitude_closed_form(N, phi):
    phi = np.asarray(phi, dtype=float)
    s = np.sin(phi)
    out = np.empty_like(phi)
    small = np.abs(s) < 1e-12
    out[small] = (N + 1) / 2.0
    out[~small] = np.sin((N + 1) * phi[~small]) / (2.0 * s[~small])
    return out


def amplitude_direct_sum(N, phi):
    phi = np.asarray(phi, dtype=float)
    out = np.zeros_like(phi)
    for m in range((N - 1) // 2 + 1):
        out += np.cos((2 * m + 1) * phi)
    return out


def verify_closed_form(N):
    phi = np.linspace(-np.pi / 2, np.pi / 2, 2001)
    err = np.max(np.abs(amplitude_closed_form(N, phi) - amplitude_direct_sum(N, phi)))
    assert err < 1e-8, f"closed form vs direct sum mismatch for N={N}: {err}"
    return err


def normalized_squared(N, phi):
    A = amplitude_closed_form(N, phi)
    I = A ** 2
    return I / I.max()


def widths_percent(N):
    """Return (first-zero half-width, 1%-level half-width) in percent.

    The 1%-level half-width is the FIRST descent through 1% inside the main
    lobe (phi_1% ~ 2.85/(N+1), eq. (4.2) of the paper), consistent with the
    paper's definition of W(N). Phase percent uses 100% = pi (= 180 deg)."""
    first_zero = 100.0 / (N + 1)  # phi = pi/(N+1) rad -> /pi*100 %
    # first crossing of the 1% level, on a fine grid out to the first zero
    pct = np.linspace(0.0, first_zero, 400001)
    phi = pct / 100.0 * np.pi
    y = normalized_squared(N, phi)
    below = np.where(y <= 0.01)[0]
    first_001 = pct[below[0]] if below.size else first_zero
    return first_zero, first_001


def main():
    # The three curves -- drawn together in EVERY panel.
    curves = [
        (99, "tab:red", "$N=99$"),
        (999, "tab:blue", "$N=999$"),
        (9999, "tab:green", "$N=9999$"),
    ]
    xranges = [10.0, 1.0, 0.1]

    # verify + report widths once
    for N, _, _ in curves:
        err = verify_closed_form(N)
        fz, w01 = widths_percent(N)
        print(f"N={N:5d}: err={err:.1e}, first-zero=+/-{fz:.4g}%, 1%-level(first)=+/-{w01:.4g}%")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0), sharey=True)

    for ax, xr in zip(axes, xranges):
        pct = np.linspace(-xr, xr, 200001)
        phi = pct / 100.0 * np.pi
        for N, color, label in curves:
            y = normalized_squared(N, phi)
            ax.plot(pct, y, color=color, lw=1.3, label=label)
        ax.set_xlim(-xr, xr)
        ax.set_ylim(0, 1.02)
        ax.grid(True, which="both", color="0.85", lw=0.6)
        ax.set_title(f"phase window $\\pm{xr:g}\\%$", fontsize=12)
        ax.set_xlabel("phase $\\varphi$  (%)")
        ax.legend(loc="upper right", fontsize=9, framealpha=0.95)

    axes[0].set_ylabel("squared amplitude $|S_N(\\varphi)|^2$\n(normalized to peak $=1.0$)")

    fig.suptitle(
        "Figure 2.  The same three curves ($N=99$ red, $999$ blue, $9999$ green) overlaid at three "
        "zooms ($\\pm10\\%,\\ \\pm1\\%,\\ \\pm0.1\\%$).\n"
        "In each window exactly one main lobe fits; the others are $\\times10$ wider or narrower, so the "
        "localization width scales as $W(N)\\propto 1/N$.",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))

    for ext in ("png", "svg"):
        fig.savefig(OUT / f"fig02_odd_harmonic_scaling.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
