"""
Double-slit interference for several wavelengths, shown SEPARATELY (not summed),
in the phase representation, each Born-normalised (int P dPhi = 1).

Two unit-amplitude point sources at a FIXED physical separation d = 6 (the distance
that equals 6 wavelengths = 6*360 deg for lambda = 1). The reference phase Phi is the
phase the lambda = 1 wave accumulates over the path difference; Phi spans the full
angular range sin(alpha) in [-1, 1], i.e. Phi in [-6*360, +6*360] deg = [-12pi, 12pi].

For a wave of wavelength lambda = 1/m (m = 1, 3, 5, 9) the actual relative phase over
the same physical path difference is m times larger:
    phi_lambda = m * Phi
Born rule (squared amplitude = amplitude times complex conjugate):
    P_raw = |psi|^2 = psi* psi = |1 + e^{i phi_lambda}|^2 = 4 cos^2(phi_lambda/2)
          = 4 cos^2(m Phi / 2)
Rigorous normalisation (m integer => integral is the same for every wavelength):
    Integral_{-12pi}^{+12pi} 4 cos^2(m Phi/2) dPhi = 48 pi
    P_norm(Phi) = cos^2(m Phi/2) / (12 pi),   peak 1/(12pi) ~ 0.02653
Number of fringes = 12 m + 1  ->  lambda=1:13, 1/3:37, 1/5:61, 1/9:109.

The four curves share the same peak height and the same unit area; only the fringe
DENSITY differs (3 : 5 : 9 relative to the black lambda=1 reference).

Outputs (this folder):
    double_slit_phase_multilambda.png  (300 dpi)
    double_slit_phase_multilambda.svg  (vector)
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

W = 6  # fixed physical spacing in lambda=1 wavelengths -> Phi half-range = 6*360 deg

# (m = 1/lambda, label, colour, zorder)
CURVES = [
    (9, r"$\lambda=1/9$", "#2ca02c", 1),  # green, densest -> draw behind
    (5, r"$\lambda=1/5$", "#1f77b4", 2),  # blue
    (3, r"$\lambda=1/3$", "#d62728", 3),  # red
    (1, r"$\lambda=1$ (reference)", "black", 4),  # black, sparsest -> on top
]


def main():
    phi_max_deg = W * 360.0
    Phi_deg = np.linspace(-phi_max_deg, phi_max_deg, 400001)
    Phi = np.deg2rad(Phi_deg)

    _trapz = getattr(np, "trapezoid", None) or np.trapz
    peak = 1.0 / (12 * np.pi)

    fig, ax = plt.subplots(figsize=(13.0, 5.4))

    for m, label, colour, z in CURVES:
        P = np.cos(m * Phi / 2.0) ** 2 / (12 * np.pi)   # normalized Born density
        total = _trapz(P, Phi)
        n_fringe = 12 * m + 1
        lw = {1: 1.0, 3: 0.8, 5: 0.65, 9: 0.5}[m]
        alpha = {1: 0.95, 3: 0.85, 5: 0.8, 9: 0.7}[m]
        ax.plot(Phi_deg, P, color=colour, lw=lw, alpha=alpha, zorder=z,
                label=f"{label}  ($m={m}$, {n_fringe} fringes)")
        print(f"m={m}: int P dPhi = {total:.6f}, fringes = {n_fringe}")

    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel(r"reference phase  $\Phi$  (deg, the $\lambda=1$ phase; "
                  r"$\sin\alpha\in[-1,1]$)", fontsize=12)
    ax.set_ylabel(r"normalized probability density  $P$" "\n"
                  r"($\int P\,d\Phi=1$, $\Phi$ in rad)", fontsize=11)
    ax.set_title(
        r"Double-slit interference, separate wavelengths (not summed), "
        r"Born $P=|\psi|^2=\psi^{*}\psi$, each normalized $\int P\,d\Phi=1$",
        fontsize=12.5)

    ax.set_xlim(-phi_max_deg, phi_max_deg)
    ax.set_ylim(0, peak * 1.5)
    ax.set_xticks(np.arange(-phi_max_deg, phi_max_deg + 1, 360))
    ax.set_xticklabels([f"{int(t)}" for t in np.arange(-phi_max_deg, phi_max_deg + 1, 360)],
                       fontsize=8, rotation=45)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.95, ncol=1)

    # ---- conditions box ----
    cond = (
        "Conditions\n"
        r"fixed physical spacing $d=6$ ($=6\lambda$ at $\lambda=1$)" "\n"
        r"unit-amplitude sources, Born $\int P\,d\Phi=1$" "\n"
        r"$\Phi\in[-6\!\cdot\!360^\circ,+6\!\cdot\!360^\circ]=[-12\pi,12\pi]$"
    )
    ax.text(0.30, 0.97, cond, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, bbox=dict(boxstyle="round", fc="#eef6ff", ec="#4a7fb5", alpha=0.95))

    # ---- formula box ----
    formula = (
        r"$\lambda=1/m,\quad \varphi_\lambda=m\,\Phi$" "\n"
        r"$P=|\psi|^2=\psi^{*}\psi=|1+e^{i\varphi_\lambda}|^2=4\cos^2(m\Phi/2)$" "\n"
        r"$\int_{-12\pi}^{12\pi}4\cos^2(m\Phi/2)\,d\Phi=48\pi$" "\n"
        r"$P_{\mathrm{norm}}=\dfrac{\cos^2(m\Phi/2)}{12\pi},\ \ "
        r"P_{\max}=\dfrac{1}{12\pi}$"
    )
    ax.text(0.988, 0.97, formula, transform=ax.transAxes, va="top", ha="right",
            fontsize=9.5, bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))

    fig.tight_layout()
    png = os.path.join(HERE, "double_slit_phase_multilambda.png")
    svg = os.path.join(HERE, "double_slit_phase_multilambda.svg")
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    print("wrote", png)
    print("wrote", svg)


if __name__ == "__main__":
    main()
