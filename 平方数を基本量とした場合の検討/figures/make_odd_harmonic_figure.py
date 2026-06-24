#!/usr/bin/env python3
"""Generate the figure for the odd-harmonic localization observation paper.

The wave is the equal-amplitude odd-harmonic sum

    S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi)

and the observable is the squared amplitude (squared norm of the real wave)

    I_N(phi) = | S_N(phi) |^2 = S_N(phi)^2 ,

normalized so that its maximum equals 1.0:

    I_N(phi) / max_phi I_N(phi) .

The vertical axis is computed honestly: the amplitude S_N is evaluated, then
squared, then divided by the array maximum (the peak at phi = 0). It is NOT a
shortcut. The closed form S_N = sin((N+1)phi)/(2 sin phi) is used for speed but
is verified against the DIRECT term-by-term sum to machine precision before
plotting (see verify_closed_form), so squaring the closed form and squaring the
direct sum give the same curve.

Domain / axis convention:
  * Fundamental domain = the HALF-WAVELENGTH interval phi in [-pi/2, pi/2]
    (width pi = 180 deg). Every odd harmonic vanishes there at the ends
    (cos((2m+1) pi/2) = 0), so S_N is maximal at the center and zero at both
    ends: a single isolated pulse, with NO endpoint peak even after squaring.
  * Horizontal axis = phase in percent with 100% = 180 deg, i.e.
    +/-50% <-> +/- 90 deg <-> +/- pi/2.
  * Increasing N narrows the central peak as 1/N.

Labels are Japanese (Hiragino Sans). Writes PNG and SVG.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

# All labels are in English so the same figure can be dropped into an English
# manuscript as-is.
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).resolve().parent


def amplitude_closed_form(N, phi):
    """S_N(phi) = sin((N+1) phi) / (2 sin phi), with the removable value
    S_N -> (N+1)/2 at sin(phi) = 0 (here only phi = 0 inside [-pi/2, pi/2])."""
    phi = np.asarray(phi, dtype=float)
    s = np.sin(phi)
    out = np.empty_like(phi)
    small = np.abs(s) < 1e-12
    out[small] = (N + 1) / 2.0
    out[~small] = np.sin((N + 1) * phi[~small]) / (2.0 * s[~small])
    return out


def amplitude_direct_sum(N, phi):
    """S_N(phi) by the literal definition: sum_{m=0}^{(N-1)/2} cos((2m+1) phi).

    No closed form, no shortcut. Used to verify amplitude_closed_form."""
    phi = np.asarray(phi, dtype=float)
    out = np.zeros_like(phi)
    for m in range((N - 1) // 2 + 1):
        out += np.cos((2 * m + 1) * phi)
    return out


def verify_closed_form(N):
    """Assert the closed form equals the direct sum to machine precision."""
    phi = np.linspace(-np.pi / 2, np.pi / 2, 2001)
    err = np.max(np.abs(amplitude_closed_form(N, phi) - amplitude_direct_sum(N, phi)))
    assert err < 1e-8, f"closed form vs direct sum mismatch for N={N}: {err}"
    return err


def normalized_squared_amplitude(N, phi):
    """I_N(phi)/max I_N = S_N(phi)^2 / max(S_N^2), peak normalized to exactly 1."""
    A = amplitude_closed_form(N, phi)
    I = A ** 2
    return I / I.max()


def main():
    # Half-wavelength domain phi in [-pi/2, pi/2]. Phase axis in percent with
    # 100% = 180 deg, so +/-50% <-> +/- 90 deg <-> +/- pi/2.
    pct = np.linspace(-50.0, 50.0, 600001)
    phi = pct / 100.0 * np.pi  # 100% <-> pi (= 180 deg)

    cases = [
        (99, "tab:red", "up to the 99th harmonic ($N=99$)"),
        (999, "tab:blue", "up to the 999th harmonic ($N=999$)"),
        (9999, "tab:green", "up to the 9999th harmonic ($N=9999$)"),
    ]

    fig, ax = plt.subplots(figsize=(11, 5.6))
    for N, color, label in cases:
        err = verify_closed_form(N)
        y = normalized_squared_amplitude(N, phi)
        print(
            f"N={N:5d}: closed-vs-direct max err={err:.2e}, "
            f"peak={y.max():.6f}, value at +/-50% (=+/-90deg)={y[0]:.2e}"
        )
        ax.plot(pct, y, color=color, lw=1.2, label=label)

    # Zoom in on the central pulse: outside +/-10% the squared amplitude is
    # essentially zero, so plotting the full +/-50% just shows flat baseline.
    ax.set_xlim(-10, 10)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("phase $\\varphi$  (%,  $100\\%=180^\\circ$,  $\\pm50\\%\\leftrightarrow\\pm90^\\circ=\\pm\\pi/2$)")
    ax.set_ylabel("squared amplitude $|S_N(\\varphi)|^2$  (normalized to peak $=1.0$)")
    ax.set_title(
        "Figure 1.  Isolated peak wave: the equal-amplitude odd-harmonic sum on the half-wavelength interval\n"
        "$S_N(\\varphi)=\\sum_{m=0}^{(N-1)/2}\\cos((2m+1)\\varphi)"
        "=\\dfrac{\\sin((N+1)\\varphi)}{2\\sin\\varphi}$",
        fontsize=12,
    )
    ax.grid(True, which="both", color="0.85", lw=0.6)
    ax.legend(loc="upper right", framealpha=0.95)

    # Specification box (matches the intended setup).
    spec = (
        "Phase interval: $\\varphi\\in[-\\pi/2,\\,\\pi/2]$  (half-wavelength, $180^\\circ$)\n"
        "Vertical: squared amplitude $|S_N(\\varphi)|^2$, normalized to peak $1.0$\n"
        "Horizontal: phase $\\varphi$  (%,  $100\\%=180^\\circ$)"
    )
    ax.text(
        0.015, 0.97, spec, transform=ax.transAxes, va="top", ha="left",
        fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.95),
    )

    for ext in ("png", "svg"):
        fig.savefig(OUT / f"fig01_odd_harmonic_localization.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
