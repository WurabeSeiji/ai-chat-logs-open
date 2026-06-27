"""
Figure 3 -- the 9 odd harmonics k = 1,3,5,...,17 of Figure 2, SYNTHESISED
(coherently summed) into one interference wave through the double slit.

Same conditions as Figure 2 (nothing changed): two sources separated by
W = 120 deg (= 20 x 6), horizontal phase axis Phi from -180 to +180 deg (no clip),
peak normalised to 1.0.

Coherent synthesis: sum the two-source AMPLITUDES of all harmonics first, then
square once (Born), exactly like Figure 1 ("amplitudes add, square at the end"):

    psi_k(Phi)    = 2 cos(delta_k/2),   delta_k = k W Phi / 180 deg
    psi_total(Phi)= sum_{k=1,3,...,17} 2 cos( k W Phi / 360 deg )
    I_total(Phi)  = | psi_total |^2

Closed form for W = 120 deg (delta_k/2 = k Phi/3):
    psi_total = sin(6 Phi) / sin(Phi/3)
    I_total   = sin^2(6 Phi) / sin^2(Phi/3)     (then normalised to peak 1.0)

Central peak at Phi=0 with interference fringes (sin^2 6Phi, 30 deg pitch).

Outputs (this folder):
    odd17_halfwave_synth_W120.png / .svg
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
K = np.arange(1, 18, 2)        # 1,3,...,17
W = 120.0                      # source separation (phase, degrees) = 20 x 6


def main():
    phi_deg = np.linspace(-180.0, 180.0, 360001)

    # coherent amplitude sum, then square (Born)
    arg = np.deg2rad(np.outer(K, W * phi_deg / 360.0))     # k W Phi / 360 (rad)
    psi_total = (2.0 * np.cos(arg)).sum(axis=0)            # sum_k 2 cos(k W Phi/360)
    I_total = psi_total ** 2

    # closed-form check away from Phi=0: sin^2(6Phi)/sin^2(Phi/3)
    cf = np.sin(np.deg2rad(6 * phi_deg)) ** 2 / np.sin(np.deg2rad(phi_deg / 3.0)) ** 2
    good = np.abs(phi_deg) > 1e-3
    print("max|direct-closed| =", np.max(np.abs(I_total[good] - cf[good])))

    In = I_total / I_total.max()                           # peak -> 1.0
    print("peak =", In.max(), "(should be 1.0)")

    fig, ax = plt.subplots(figsize=(13.0, 5.2))
    ax.plot(phi_deg, In, color="#5a2d82", lw=1.1)
    ax.axhline(0, color="0.6", lw=0.8)

    ax.set_xlabel(r"observation phase  $\Phi$  (deg)", fontsize=12)
    ax.set_ylabel(r"normalized intensity  $I_{\mathrm{total}}$  (peak $=1.0$)", fontsize=11)
    ax.set_title(
        r"(3) Double slit ($W=120^\circ$): synthesised (coherent) interference wave "
        r"of $k=1,3,\dots,17$", fontsize=12.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 1.12)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.grid(True, ls=":", alpha=0.45)

    cond = (
        "Conditions (same as Fig. 2)\n"
        r"$k=1,3,\dots,17$ coherently summed (synthesised)" "\n"
        r"two sources $W=120^\circ=20^\circ\times6$;  axis $\pm180^\circ$, no clip"
    )
    ax.text(0.012, 0.97, cond, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.5, bbox=dict(boxstyle="round", fc="#eef6ff", ec="#4a7fb5", alpha=0.95))

    formula = (
        r"$\psi_{\mathrm{total}}=\sum_k 2\cos(kW\Phi/360^\circ)$" "\n"
        r"$I_{\mathrm{total}}=|\psi_{\mathrm{total}}|^2=\dfrac{\sin^2 6\Phi}{\sin^2(\Phi/3)}$" "\n"
        r"(coherent: add amplitudes, square once; peak $=1.0$)"
    )
    ax.text(0.988, 0.97, formula, transform=ax.transAxes, va="top", ha="right",
            fontsize=10, bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "odd17_halfwave_synth_W120.png"), dpi=300)
    fig.savefig(os.path.join(HERE, "odd17_halfwave_synth_W120.svg"))
    print("wrote odd17_halfwave_synth_W120.png/.svg")


if __name__ == "__main__":
    main()
