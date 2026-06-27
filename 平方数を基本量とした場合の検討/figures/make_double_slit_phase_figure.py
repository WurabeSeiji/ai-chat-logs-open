"""
Double-slit interference in the phase representation, Born-normalised.

Two unit-amplitude point sources (slits) emit the conjugate-complex wave
    Z(theta) = cos(theta) + i sin(theta) = e^{i theta}.
Wavelength lambda = 1, slit spacing W = 6 (in wavelengths), so the maximum
path-difference phase is W * 360 deg = 6 * 360 deg = 2160 deg = 12 pi rad.

Probability amplitude (superposition principle), relative phase phi between paths:
    psi = psi_1 + psi_2 = e^{i theta} (1 + e^{i phi})
Born rule (squared amplitude = amplitude times its complex conjugate):
    P_raw = |psi|^2 = psi* psi = (1 + e^{-i phi})(1 + e^{i phi})
          = 2 + 2 cos(phi) = 4 cos^2(phi/2)
The carrier e^{i theta} cancels exactly, so P depends only on the relative phase phi.

Rigorous probability interpretation: normalise so that the total probability over
the whole domain is 1.
    Integral_{-12pi}^{+12pi} 4 cos^2(phi/2) dphi = 48 pi
    P_norm(phi) = 4 cos^2(phi/2) / (48 pi) = cos^2(phi/2) / (12 pi)     [phi in rad]
Peak value 1/(12 pi) ~ 0.02653 at phi = m*360 deg (m = -6..+6, 13 fringes).

Outputs (this folder):
    double_slit_phase_normalized.png  (300 dpi)
    double_slit_phase_normalized.svg  (vector)
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

W = 6          # slit spacing in wavelengths -> phase half-range = W*360 deg
LAMBDA = 1.0   # wavelength
AMP = 1.0      # source amplitude


def main():
    phi_max_deg = W * 360.0                 # 2160
    phi_deg = np.linspace(-phi_max_deg, phi_max_deg, 200001)
    phi = np.deg2rad(phi_deg)               # radians

    # Born squared amplitude (psi* psi) and its rigorous normalisation
    P_raw = 4.0 * np.cos(phi / 2.0) ** 2     # = 2 + 2 cos phi, max 4
    norm = 48.0 * np.pi                      # int_{-12pi}^{12pi} P_raw dphi
    P = P_raw / norm                         # = cos^2(phi/2)/(12 pi)

    # numerical check that the total probability is 1
    _trapz = getattr(np, "trapezoid", None) or np.trapz
    total = _trapz(P, phi)
    print(f"integral P dphi (should be 1) = {total:.8f}")
    print(f"peak value 1/(12 pi)          = {1.0/(12*np.pi):.6f}")

    fig, ax = plt.subplots(figsize=(12.5, 5.2))

    ax.plot(phi_deg, P, color="#d62728", lw=1.6,
            label=r"$P_{\mathrm{norm}}(\varphi)=\cos^2(\varphi/2)/(12\pi)$")
    ax.axhline(0, color="0.6", lw=0.8)

    # mark the 13 principal maxima
    m = np.arange(-W, W + 1)
    peak_x = m * 360.0
    peak_y = np.full_like(peak_x, 1.0 / (12 * np.pi), dtype=float)
    ax.plot(peak_x, peak_y, "o", color="#7a0d0d", ms=4, zorder=5)

    ax.set_xlabel(r"relative phase  $\varphi$  (deg)", fontsize=12)
    ax.set_ylabel(r"normalized probability density  $P$" "\n"
                  r"($\int P\,d\varphi=1$, $\varphi$ in rad)", fontsize=11)
    ax.set_title(
        r"Double-slit interference (phase representation), Born rule "
        r"$P=|\psi|^2=\psi^{*}\psi$, normalized $\int P\,d\varphi=1$",
        fontsize=12.5)

    ax.set_xlim(-phi_max_deg, phi_max_deg)
    ax.set_ylim(0, 1.0 / (12 * np.pi) * 1.55)
    ax.set_xticks(np.arange(-phi_max_deg, phi_max_deg + 1, 360))
    ax.set_xticklabels([f"{int(t)}" for t in np.arange(-phi_max_deg, phi_max_deg + 1, 360)],
                       fontsize=8, rotation=45)
    ax.grid(True, ls=":", alpha=0.45)

    # ---- conditions box (upper left) ----
    cond = (
        "Conditions\n"
        rf"$\lambda=1$,  slit spacing $W={W}$ ($={W}\lambda$)" "\n"
        r"unit-amplitude sources $A=1$" "\n"
        rf"phase range $\pm W\!\cdot\!360^\circ=\pm{int(phi_max_deg)}^\circ\,(=\pm12\pi)$" "\n"
        rf"orders $m=-{W}\dots+{W}$  (13 fringes)"
    )
    ax.text(0.012, 0.97, cond, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, family="DejaVu Sans",
            bbox=dict(boxstyle="round", fc="#eef6ff", ec="#4a7fb5", alpha=0.95))

    # ---- formula box (upper right) ----
    formula = (
        r"$\psi=\psi_1+\psi_2=e^{i\theta}\,(1+e^{i\varphi})$" "\n"
        r"$P=|\psi|^2=\psi^{*}\psi=(1+e^{-i\varphi})(1+e^{i\varphi})$" "\n"
        r"$\quad=2+2\cos\varphi=4\cos^2(\varphi/2)$" "\n"
        r"$\int_{-12\pi}^{12\pi}\!4\cos^2(\varphi/2)\,d\varphi=48\pi$" "\n"
        r"$P_{\mathrm{norm}}=\dfrac{\cos^2(\varphi/2)}{12\pi},\ \ "
        r"P_{\max}=\dfrac{1}{12\pi}$"
    )
    ax.text(0.988, 0.97, formula, transform=ax.transAxes, va="top", ha="right",
            fontsize=9.5,
            bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))

    fig.tight_layout()

    png = os.path.join(HERE, "double_slit_phase_normalized.png")
    svg = os.path.join(HERE, "double_slit_phase_normalized.svg")
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    print("wrote", png)
    print("wrote", svg)


if __name__ == "__main__":
    main()
