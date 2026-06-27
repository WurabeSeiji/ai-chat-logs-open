"""
Figure 2 -- two-source (double-slit) interference of the odd harmonics
k = 1,3,5,...,17, each interfered INDEPENDENTLY (not synthesised / not summed),
overlaid on one horizontal phase axis Phi from -180 to +180 deg (NO clipping).

The two sources are separated by a phase  W = 120 deg ( = 20 deg x 6 ).
For source-separation W, the relative phase of harmonic k at observation phase Phi
grows linearly (so the pattern spreads across the full +/-180 deg, no clip):

    delta_k(Phi) = k * W * Phi / 180 deg
    I_k(Phi)     = | e^{+i delta_k/2} + e^{-i delta_k/2} |^2
                 = 4 cos^2( delta_k / 2 )
                 = 4 cos^2( k W Phi / 360 deg )                  (Born, two equal sources)

With W = 120 deg this is  I_k(Phi) = 4 cos^2( k Phi / 3 ).

Fringe spacing of harmonic k = 540 deg / k ; number of fringes over +/-180 deg
= 2*floor(k/3)+1 :  k=1 ->1, 3->3, 5->3, 7->5, 9->7, 11->7, 13->9, 15->11, 17->11.
(The carrier k=9 gives 7 fringes: peaks at 0, +/-60, +/-120, +/-180 deg.)

Outputs (this folder):
    odd17_halfwave_double_W120.png / .svg
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

    fig, ax = plt.subplots(figsize=(13.0, 5.6))
    cmap = plt.cm.viridis

    for i, k in enumerate(K):
        delta_k = k * W * phi_deg / 180.0          # relative phase (deg)
        Ik = np.cos(np.deg2rad(delta_k / 2.0)) ** 2   # 4cos^2(.)/4 -> peak normalized to 1.0
        n_fr = 2 * int(k // 3) + 1
        ax.plot(phi_deg, Ik, color=cmap(i / (len(K) - 1)), lw=0.9, alpha=0.85,
                label=f"$k={k}$  ({n_fr} fr.)")

    ax.set_xlabel(r"observation phase  $\Phi$  (deg)", fontsize=12)
    ax.set_ylabel(r"normalized intensity  $I_k$  (peak $=1.0$)", fontsize=12)
    ax.set_title(
        r"(2) Double slit ($W=120^\circ=20^\circ\!\times\!6$): odd harmonics "
        r"$k=1,3,\dots,17$ interfered independently", fontsize=12.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 1.12)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.grid(True, ls=":", alpha=0.45)
    ax.legend(loc="upper center", ncol=9, fontsize=7.5, framealpha=0.95,
              columnspacing=0.8, handlelength=1.0)

    cond = (
        "Conditions\n"
        r"odd harmonics $k=1,3,\dots,17$ (independent, not summed)" "\n"
        r"two sources $W=120^\circ=20^\circ\times6$;  axis $\pm180^\circ$, no clip"
    )
    ax.text(0.012, 0.97, cond, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.5, bbox=dict(boxstyle="round", fc="#eef6ff", ec="#4a7fb5", alpha=0.95))

    formula = (
        r"$\delta_k(\Phi)=k\,W\,\Phi/180^\circ$" "\n"
        r"$I_k=4\cos^2(\delta_k/2)=4\cos^2(kW\Phi/360^\circ)$" "\n"
        r"$W=120^\circ:\ I_k=4\cos^2(k\Phi/3)$"
    )
    ax.text(0.988, 0.04, formula, transform=ax.transAxes, va="bottom", ha="right",
            fontsize=9.5, bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "odd17_halfwave_double_W120.png"), dpi=300)
    fig.savefig(os.path.join(HERE, "odd17_halfwave_double_W120.svg"))
    print("wrote odd17_halfwave_double_W120.png/.svg")
    for k in K:
        print(f"  k={k}: fringes = {2*int(k//3)+1}")


if __name__ == "__main__":
    main()
