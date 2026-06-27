"""
Figure 1 -- single source, odd-harmonic localized peak on the HALF-WAVE closed
system (the correct Paper-2 construction).

Odd harmonics k = 1,3,5,...,17 (9 terms, central harmonic 9). The squared-amplitude
uses the REAL cosine sum (Paper 2's I_N), NOT the analytic |sum e^{ikPhi}|^2:

    S(Phi) = sum_{k=1,3,...,17} cos(k Phi) = sin(18 Phi) / (2 sin Phi)
    P1(Phi) = |S(Phi)|^2 = sin^2(18 Phi) / (4 sin^2 Phi)      (Phi in degrees)

Closed system (fundamental domain) = HALF-WAVE 180 deg = [-90, +90]:
the odd-harmonic sum is anti-periodic, S(Phi+180) = -S(Phi), so |S|^2 has period
180 deg; on [-90,90] there is exactly ONE central peak and the boundaries +/-90 deg
are zeros (cos(odd*90)=0).

Checkpoints:
  peak at Phi=0  (S=9, |S|^2=81 before normalisation)
  first zero at Phi = +/-10 deg   (sin(18 Phi)=0 -> 18 Phi = 180 -> Phi = 10)
  zeros at +/-10, +/-20, ..., +/-90 deg
The horizontal ruler runs -180..+180 deg, but the curve is drawn only on [-90,90].

Normalised to PEAK height 1.0 (divide by the value at Phi=0).

Outputs (this folder):
    odd17_halfwave_single.png / .svg
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
K = np.arange(1, 18, 2)              # 1,3,...,17  (central harmonic 9)
_trapz = getattr(np, "trapezoid", None) or np.trapz


def main():
    phi_deg = np.linspace(-90.0, 90.0, 180001)   # plotted domain (closed system)
    phi = np.deg2rad(phi_deg)

    S = np.cos(np.outer(K, phi)).sum(axis=0)      # real cosine sum
    P1 = S ** 2

    # closed-form check (away from Phi=0)
    cf = np.sin(18 * phi) ** 2 / (4 * np.sin(phi) ** 2)
    good = np.abs(phi_deg) > 1e-4
    print("max|direct-closed| =", np.max(np.abs(P1[good] - cf[good])))

    P1n = P1 / P1.max()                 # normalize to peak height 1.0
    print("peak P1n =", P1n.max(), " (should be 1.0)")

    fig, ax = plt.subplots(figsize=(11.0, 5.0))
    ax.plot(phi_deg, P1n, color="#b00050", lw=1.3)
    ax.axhline(0, color="0.6", lw=0.8)

    # closed-system boundaries at +/-90 deg
    for b in (-90, 90):
        ax.axvline(b, color="0.4", lw=1.0, ls="--")
    ax.text(0, P1n.max() * 1.02, "closed system  [-90$^\\circ$, +90$^\\circ$]  (half-wave)",
            ha="center", va="bottom", fontsize=9, color="0.3")

    # mark first zeros at +/-10 deg
    for z in (-10, 10):
        ax.plot([z], [0], "o", color="#444", ms=5, zorder=5)
    ax.annotate(r"first zero $\pm10^\circ$", (10, 0), textcoords="offset points",
                xytext=(8, 14), fontsize=10, color="#222",
                arrowprops=dict(arrowstyle="->", color="#222", lw=0.8))

    ax.set_xlabel(r"phase  $\Phi$  (deg)", fontsize=12)
    ax.set_ylabel(r"normalized squared amplitude  $P_1$  (peak $=1.0$)", fontsize=11)
    ax.set_title(
        r"(1) Single source: odd-harmonic localized peak on the half-wave system "
        r"($k=1,3,\dots,17$)", fontsize=12.5)
    ax.set_xlim(-180, 180)                  # ruler -180..180
    ax.set_ylim(0, P1n.max() * 1.12)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.grid(True, ls=":", alpha=0.45)

    formula = (
        r"$S=\sum_{k=1,3,\dots,17}\cos(k\Phi)=\dfrac{\sin 18\Phi}{2\sin\Phi}$" "\n"
        r"$P_1=|S|^2=\dfrac{\sin^2 18\Phi}{4\sin^2\Phi}$" "\n"
        r"peak at $0$; zeros at $\pm10^\circ,\pm20^\circ,\dots,\pm90^\circ$"
    )
    ax.text(0.985, 0.97, formula, transform=ax.transAxes, va="top", ha="right",
            fontsize=10.5, bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "odd17_halfwave_single.png"), dpi=300)
    fig.savefig(os.path.join(HERE, "odd17_halfwave_single.svg"))
    print("wrote odd17_halfwave_single.png/.svg")


if __name__ == "__main__":
    main()
