"""
Figure 2 with W = 10 deg -- two-source (double-slit) interference of the odd
harmonics k = 1,3,5,...,17, each interfered INDEPENDENTLY (not summed), overlaid
on the phase axis Phi from -180 to +180 deg (no clip).

Exactly the same equation as the W=120 case, only W changed to 10 deg:
    delta_k(Phi) = k * W * Phi / 180 deg
    I_k(Phi)     = 4 cos^2(delta_k/2) = 4 cos^2( k W Phi / 360 deg )
peak-normalised to 1.0:  I_k = cos^2( k W Phi / 360 deg ).

Threshold for a harmonic to show a complete dark fringe (full destructive
interference): W >= lambda_k/2 = (360/k)/2 = 180/k deg, i.e. k >= 180/W.
For W = 10 deg this is k >= 18, but the highest harmonic is 17, so NONE of the
nine harmonics reach the threshold: every curve is a single broad central lobe
that decreases toward the edges WITHOUT reaching zero (no dark fringe).
Edge value at Phi=+/-180 deg: I_k = cos^2(k*W/2 deg) = cos^2(5k deg).

Outputs (this folder):
    odd17_halfwave_double_W10.png / .svg
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
K = np.arange(1, 18, 2)        # 1,3,...,17
W = 10.0                       # source separation (phase, degrees)


def main():
    phi_deg = np.linspace(-180.0, 180.0, 360001)

    fig, ax = plt.subplots(figsize=(13.0, 5.6))
    cmap = plt.cm.viridis

    print(f"W = {W} deg;  fringe threshold k >= 180/W = {180.0/W:.1f}")
    for i, k in enumerate(K):
        delta_k = k * W * phi_deg / 180.0          # relative phase (deg)
        Ik = np.cos(np.deg2rad(delta_k / 2.0)) ** 2   # peak-normalized to 1.0
        n_max = 2 * int(k * W / 360.0) + 1          # number of maxima in +/-180
        edge = np.cos(np.deg2rad(k * W / 2.0)) ** 2  # value at Phi=+/-180
        has_fringe = "dark fringe" if k >= 180.0 / W else "NO dark fringe"
        print(f"  k={k:2d}: maxima={n_max}, edge I={edge:.4f}, {has_fringe}")
        ax.plot(phi_deg, Ik, color=cmap(i / (len(K) - 1)), lw=1.0, alpha=0.85,
                label=f"$k={k}$")

    ax.set_xlabel(r"observation phase  $\Phi$  (deg)", fontsize=12)
    ax.set_ylabel(r"normalized intensity  $I_k$  (peak $=1.0$)", fontsize=12)
    ax.set_title(
        r"Figure 2 with $W=10^\circ$: odd harmonics $k=1,3,\dots,17$ "
        r"interfered independently", fontsize=12.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 1.12)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.grid(True, ls=":", alpha=0.45)
    ax.legend(loc="lower center", ncol=9, fontsize=8, framealpha=0.95,
              columnspacing=0.8, handlelength=1.0)

    cond = (
        "Conditions\n"
        r"$k=1,3,\dots,17$ independent (not summed)" "\n"
        r"two sources $W=10^\circ$;  axis $\pm180^\circ$, no clip" "\n"
        r"fringe threshold $k\geq180/W=18$  $\Rightarrow$  no $k\leq17$ reaches it"
    )
    ax.text(0.012, 0.97, cond, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.5, bbox=dict(boxstyle="round", fc="#eef6ff", ec="#4a7fb5", alpha=0.95))

    formula = (
        r"$\delta_k(\Phi)=k\,W\,\Phi/180^\circ$" "\n"
        r"$I_k=4\cos^2(\delta_k/2)=4\cos^2(kW\Phi/360^\circ)$" "\n"
        r"$W=10^\circ:$ edge $I_k(\pm180^\circ)=\cos^2(5k^\circ)>0$"
    )
    ax.text(0.988, 0.97, formula, transform=ax.transAxes, va="top", ha="right",
            fontsize=9.5, bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "odd17_halfwave_double_W10.png"), dpi=300)
    fig.savefig(os.path.join(HERE, "odd17_halfwave_double_W10.svg"))
    print("wrote odd17_halfwave_double_W10.png/.svg")


if __name__ == "__main__":
    main()
