"""
Interference vs. relative phase (phase on the horizontal axis, L-free).

Reproduces the final figure of the thought-experiment conversation:
a single unit-amplitude rotation z(theta) = cos(theta) - i sin(theta) = e^{-i theta}
has |z|^2 = 1 everywhere (flat), yet superposing two such waves with a relative
phase phi gives a genuine interference pattern in the squared amplitude:

    z_total = e^{-i theta} + e^{-i(theta + phi)} = e^{-i theta} (1 + e^{-i phi})
    |z_total|^2 = |1 + e^{-i phi}|^2 = 2 + 2 cos(phi) = 4 cos^2(phi/2)
    I(phi) = cos^2(phi/2)   (normalized to a maximum of 1)

The horizontal axis is the relative phase phi alone (no wavelength, slit spacing,
or screen distance needed): phi = (2 pi / lambda) * (path difference) ties it back
to a spatial position when desired.

Outputs (this folder):
    interference_phase_comparison.png  (300 dpi)
    interference_phase_comparison.svg  (vector)
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    phi = np.linspace(-np.pi, np.pi, 4001)

    # single unit-amplitude wave: |z|^2 is always 1
    I_single = np.ones_like(phi)
    # two equal-amplitude waves, relative phase phi
    I_interf = np.cos(phi / 2.0) ** 2

    # sanity check of the closed form against the direct two-wave sum
    theta = 0.37  # arbitrary; |z_total|^2 is independent of theta
    z_total = np.exp(-1j * theta) + np.exp(-1j * (theta + phi))
    I_direct = np.abs(z_total) ** 2 / 4.0  # divide by max (=4) to normalize
    max_err = np.max(np.abs(I_direct - I_interf))
    print(f"max |direct - closed form| = {max_err:.2e}")

    fig, ax = plt.subplots(figsize=(8.2, 5.0))

    ax.plot(
        phi, I_single,
        color="#1f77b4", lw=2.2, ls="--",
        label=r"single wave  $|z|^2=\cos^2\theta+\sin^2\theta=1$",
    )
    ax.plot(
        phi, I_interf,
        color="#d62728", lw=2.6,
        label=r"two-wave interference  $I(\phi)=\cos^2(\phi/2)$",
    )

    # mark the canonical points
    for x, y, txt in [
        (0.0, 1.0, r"$\phi=0,\ I=1$"),
        (np.pi, 0.0, r"$\phi=\pi,\ I=0$"),
        (-np.pi, 0.0, None),
        (np.pi / 2, 0.5, r"$\phi=\pi/2,\ I=0.5$"),
        (-np.pi / 2, 0.5, None),
    ]:
        ax.plot([x], [y], "o", color="#d62728", ms=6, zorder=5)
        if txt:
            ax.annotate(
                txt, (x, y), textcoords="offset points", xytext=(6, 8),
                fontsize=10, color="#7a0d0d",
            )

    ax.set_xlabel(r"relative phase  $\phi$  (rad)", fontsize=12)
    ax.set_ylabel(r"normalized squared amplitude  $I$", fontsize=12)
    ax.set_title(
        "Each unit-amplitude wave has $|z|^2=1$, "
        "yet their superposition interferes",
        fontsize=12.5,
    )

    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-0.04, 1.08)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(loc="center right", fontsize=10, framealpha=0.95)

    fig.tight_layout()

    png = os.path.join(HERE, "interference_phase_comparison.png")
    svg = os.path.join(HERE, "interference_phase_comparison.svg")
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    print("wrote", png)
    print("wrote", svg)


if __name__ == "__main__":
    main()
