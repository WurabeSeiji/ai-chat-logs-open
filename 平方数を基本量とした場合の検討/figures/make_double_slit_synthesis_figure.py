"""
Synthesis of four wavelengths through the double slit: the two interpretations.

Wavelengths lambda = 1, 1/3, 1/5, 1/9  (m = 1/lambda = 1, 3, 5, 9), equal amplitude,
fixed physical slit spacing d = 6, reference phase Phi in [-12pi, +12pi] (= +/-6*360 deg).
With the symmetric slit convention each harmonic's two-slit amplitude is real:
    psi_m(Phi) = e^{+i m Phi/2} + e^{-i m Phi/2} = 2 cos(m Phi / 2).

(1) COHERENT  (square ONCE at the end; the four harmonics are phase-locked at t = 0,
    i.e. the localized packet/envelope is formed on the source side, then it interferes):
        A(Phi) = sum_m 2 cos(m Phi / 2)
        P_coh(Phi) = |A(Phi)|^2          <-- all cross-wavelength terms present
(2) INCOHERENT  (square each harmonic's slit interference, then add; this is the
    time-averaged limit, since cross terms beat at (m-m') omega and average to 0):
        P_inc(Phi) = sum_m |psi_m|^2 = sum_m 4 cos^2(m Phi / 2)   <-- no cross terms

Both have the same total integral over the domain:
    Integral_{-12pi}^{+12pi} P_coh dPhi = Integral P_inc dPhi = 192 pi
(the cross terms integrate to zero), so after normalising int P dPhi = 1 the two
share the same unit area but differ in SHAPE: coherent concentrates probability into
sharp localized peaks, incoherent spreads it out. Their difference is exactly the
cross-wavelength interference that finite-time observation removes.

Outputs (this folder):
    synthesis_coherent.png / .svg
    synthesis_incoherent.png / .svg
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
W = 6
M = np.array([1, 3, 5, 9])   # m = 1/lambda for lambda = 1, 1/3, 1/5, 1/9


def main():
    phi_max_deg = W * 360.0
    Phi_deg = np.linspace(-phi_max_deg, phi_max_deg, 600001)
    Phi = np.deg2rad(Phi_deg)
    _trapz = getattr(np, "trapezoid", None) or np.trapz

    # per-harmonic real two-slit amplitude 2 cos(m Phi/2)
    amps = np.array([2.0 * np.cos(m * Phi / 2.0) for m in M])  # shape (4, N)

    P_coh_raw = np.abs(amps.sum(axis=0)) ** 2          # |sum psi|^2
    P_inc_raw = (amps ** 2).sum(axis=0)                # sum |psi|^2

    Z_coh = _trapz(P_coh_raw, Phi)
    Z_inc = _trapz(P_inc_raw, Phi)
    P_coh = P_coh_raw / Z_coh
    P_inc = P_inc_raw / Z_inc
    print(f"Z_coh = {Z_coh:.4f}  (analytic 192 pi = {192*np.pi:.4f})")
    print(f"Z_inc = {Z_inc:.4f}  (analytic 192 pi = {192*np.pi:.4f})")
    print(f"int P_coh = {_trapz(P_coh, Phi):.6f}   int P_inc = {_trapz(P_inc, Phi):.6f}")
    print(f"peak P_coh = {P_coh.max():.5f}   peak P_inc = {P_inc.max():.5f}")

    ymax = max(P_coh.max(), P_inc.max()) * 1.10  # shared y-scale for fair comparison

    cond = (
        "Conditions\n"
        r"$\lambda=1,\,1/3,\,1/5,\,1/9$  ($m=1,3,5,9$), equal amplitude" "\n"
        r"fixed spacing $d=6$,  $\Phi\in[-12\pi,12\pi]$,  Born $\int P\,d\Phi=1$"
    )

    def make(P, title, formula, colour, fname):
        fig, ax = plt.subplots(figsize=(13.0, 5.2))
        ax.plot(Phi_deg, P, color=colour, lw=0.9)
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_xlabel(r"reference phase  $\Phi$  (deg, the $\lambda=1$ phase)", fontsize=12)
        ax.set_ylabel(r"normalized probability density  $P$  ($\int P\,d\Phi=1$)", fontsize=11)
        ax.set_title(title, fontsize=12.5)
        ax.set_xlim(-phi_max_deg, phi_max_deg)
        ax.set_ylim(0, ymax)
        ax.set_xticks(np.arange(-phi_max_deg, phi_max_deg + 1, 360))
        ax.set_xticklabels([f"{int(t)}" for t in np.arange(-phi_max_deg, phi_max_deg + 1, 360)],
                           fontsize=8, rotation=45)
        ax.grid(True, ls=":", alpha=0.4)
        ax.text(0.012, 0.97, cond, transform=ax.transAxes, va="top", ha="left",
                fontsize=9, bbox=dict(boxstyle="round", fc="#eef6ff", ec="#4a7fb5", alpha=0.95))
        ax.text(0.988, 0.97, formula, transform=ax.transAxes, va="top", ha="right",
                fontsize=10, bbox=dict(boxstyle="round", fc="#fff7ee", ec="#c9892f", alpha=0.95))
        fig.tight_layout()
        png = os.path.join(HERE, fname + ".png")
        svg = os.path.join(HERE, fname + ".svg")
        fig.savefig(png, dpi=300)
        fig.savefig(svg)
        plt.close(fig)
        print("wrote", png)
        print("wrote", svg)

    make(
        P_coh,
        r"(1) COHERENT: square ONCE at the end  $P=\left|\sum_m\psi_m\right|^2$"
        r"  (cross-terms present; source-side packet)",
        r"$\psi_m=2\cos(m\Phi/2)$" "\n"
        r"$A=\sum_m 2\cos(m\Phi/2)$" "\n"
        r"$P_{\rm coh}=|A|^2/(192\pi)$" "\n"
        r"(all cross terms $m\neq m'$ kept; $t=0$)",
        "#7a0d8c", "synthesis_coherent",
    )
    make(
        P_inc,
        r"(2) INCOHERENT: square each, then add  $P=\sum_m|\psi_m|^2$"
        r"  (no cross-terms; time-averaged limit)",
        r"$\psi_m=2\cos(m\Phi/2)$" "\n"
        r"$P_{\rm inc}=\sum_m 4\cos^2(m\Phi/2)/(192\pi)$" "\n"
        r"(cross terms $m\neq m'$ averaged to 0)",
        "#0d7a6b", "synthesis_incoherent",
    )


if __name__ == "__main__":
    main()
