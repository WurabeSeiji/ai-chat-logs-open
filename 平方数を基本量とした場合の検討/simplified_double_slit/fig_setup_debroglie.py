#!/usr/bin/env python3
"""
Experimental-setup schematic for the de Broglie two-slit thought experiment,
in the SAME style as fig_setup_double_slit.py but at the PHYSICAL scale.

Because the three characteristic lengths differ by huge factors
    lambda0 ~ 1e-10 m (0.1 nm) ,  W = 5e-6 m (5 um) ,  L = 5e-2 m (5 cm)
    -> W/lambda0 ~ 5e4 ,  L/W ~ 1e4 ,  L/lambda0 ~ 5e8
a single true-scale drawing is impossible. Here EACH of lambda0, W, L is drawn
at its OWN independent visual scale, and the TRUE values are annotated on the
figure (nothing is to scale; only the topology source->slits->screen is real).

lambda0 is the (rough) de Broglie wavelength at the chosen voltage V and also
sets the source-position fluctuation width Delta = lambda0 (Paper-1 convention).

Usage:
    python3 fig_setup_debroglie.py                  # V=150 V (lambda0 ~ 0.1 nm)
    python3 fig_setup_debroglie.py --V 50 --L 0.05 --W 5e-6 --D 1.0
Outputs: fig_setup_debroglie_V<V>.png / .svg
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "cm"})

h  = 6.62607015e-34
me = 9.1093837015e-31
e  = 1.602176634e-19


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=0.05, help="source-slit distance [m]")
    ap.add_argument("--W", type=float, default=5e-6, help="slit separation [m]")
    ap.add_argument("--V", type=float, default=150.0, help="accelerating voltage [V]")
    ap.add_argument("--D", type=float, default=1.0, help="slit-screen distance [m]")
    args = ap.parse_args()
    L, W, V, D = args.L, args.W, args.V, args.D

    lam0 = h / np.sqrt(2 * me * e * V)          # rough de Broglie wavelength
    theta = np.degrees(np.arctan2(W / 2.0, L))  # true half-angle (tiny)
    dX = lam0 * D / W                           # fringe spacing on the screen

    # ---- visual layout (each dimension at its OWN independent scale) --------
    xS, xB, xScr = 0.0, 6.0, 11.6               # source, barrier, screen (visual x)
    a = 1.05                                    # visual half-height of source width lambda0
    w = 1.75                                    # visual half slit separation W

    fig, ax = plt.subplots(figsize=(11.4, 6.2))

    # optical axis
    ax.plot([xS - 1.6, xScr + 0.3], [0, 0], ls="--", color="0.55", lw=1.0, zorder=1)

    # ---- vibrating source (a short transverse standing wave) + fluctuation --
    yy = np.linspace(-a, a, 260)
    ax.plot(xS + 0.34 * np.sin(2 * np.pi * 2.5 * (yy + a) / (2 * a)), yy,
            color="#c0392b", lw=1.7, zorder=5)
    ax.plot(xS, 0, "o", color="black", ms=9, zorder=6)
    ax.annotate("Vibrating source $S$\n(localized odd-harmonic wave;\n"
                r"emits $\lambda_0$, speed $c$)",
                (xS, -a), textcoords="offset points", xytext=(-4, -30),
                ha="center", va="top", fontsize=10.5)

    # lambda0 dimension (its OWN scale)
    xL0 = xS - 1.15
    ax.annotate("", (xL0, a), (xL0, -a), arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
    ax.text(xL0 - 0.15, 0, r"$\lambda_0$", ha="right", va="center", color="#c0392b", fontsize=13)

    # ---- rays: source -> slits ---------------------------------------------
    ax.plot([xS, xB], [0, w], color="#1f5fbf", lw=1.6, zorder=3)
    ax.plot([xS, xB], [0, -w], color="#1f5fbf", lw=1.6, zorder=3)

    # exaggerated half-angle arc (true value annotated)
    arc = Arc((xS, 0), 3.2, 3.2, angle=0, theta1=0,
              theta2=np.degrees(np.arctan2(w, xB)), color="#7f4fbf", lw=1.5, zorder=4)
    ax.add_patch(arc)
    ax.annotate(r"$\theta=\arctan\dfrac{W/2}{L}\approx %.1e^\circ$" % theta,
                (xS + 1.85, 0.55), color="#7f4fbf", fontsize=10.5)

    # ---- barrier with two slits --------------------------------------------
    gap = 0.34
    for (y0, y1) in [(-3.3, -w - gap), (-w + gap, w - gap), (w + gap, 3.3)]:
        ax.plot([xB, xB], [y0, y1], color="black", lw=5, solid_capstyle="butt", zorder=4)
    for s, lab, dy in [(w, r"Slit $1$", 11), (-w, r"Slit $2$", -21)]:
        ax.plot(xB, s, "o", mfc="white", mec="black", ms=7, zorder=5)
        ax.annotate(lab, (xB, s), textcoords="offset points", xytext=(9, dy), fontsize=10.5)

    # W dimension (its OWN scale)
    xWd = xB + 0.85
    ax.annotate("", (xWd, w), (xWd, -w), arrowprops=dict(arrowstyle="<->", color="black", lw=1.1))
    ax.text(xWd + 0.15, 0, r"$W$", ha="left", va="center", fontsize=13)

    # ---- rays: slits -> screen (paraxial, near-parallel) with a scale break -
    for s in [w, -w]:
        ax.plot([xB, 8.4], [s, s], color="#1f5fbf", lw=1.1, alpha=0.6, zorder=2)
        ax.plot([9.4, xScr], [s, s], color="#1f5fbf", lw=1.1, alpha=0.6, zorder=2)
    for xb in (8.75, 9.05):
        ax.plot([xb - 0.12, xb + 0.12], [-0.3, 0.3], color="0.4", lw=1.2, zorder=3)

    # ---- screen + fringe spacing -------------------------------------------
    ax.plot([xScr, xScr], [-3.4, 3.4], color="black", lw=3.5, zorder=4)
    # little fringe pattern on the screen
    yf = np.linspace(-3.0, 3.0, 400)
    ax.plot(xScr + 0.30 * (0.5 + 0.5 * np.cos(2 * np.pi * yf / 0.8)), yf,
            color="#1f5fbf", lw=1.0, zorder=5)
    ax.annotate("Observation screen\n" r"($D$, fringe spacing $\Delta X=\lambda_0 D/W$)",
                (xScr, -3.4), textcoords="offset points", xytext=(0, -10),
                ha="center", va="top", fontsize=10.5)

    # ---- L dimension (its OWN scale) ---------------------------------------
    yLd = -3.9
    ax.annotate("", (xS, yLd), (xB, yLd), arrowprops=dict(arrowstyle="<->", color="black", lw=1.1))
    ax.text((xS + xB) / 2, yLd - 0.28, r"$L$", ha="center", va="top", fontsize=13)

    # ---- true-values box + not-to-scale note -------------------------------
    vals = (r"$\bf{True\ values}$ (each drawn at an independent scale):" "\n"
            r"$\lambda_0 = %.4f$ nm  $= h/\sqrt{2m_e eV}$  ($V=%g$ V)" "\n"
            r"$W = %g\ \mu$m,   $L = %g$ cm,   $D = %g$ m" "\n"
            r"$\Delta X = \lambda_0 D/W = %.1f\ \mu$m" "\n"
            r"ratios: $W/\lambda_0\approx%.0e$, $L/W\approx%.0e$, $L/\lambda_0\approx%.0e$"
            % (lam0 * 1e9, V, W * 1e6, L * 100, D, dX * 1e6,
               W / lam0, L / W, L / lam0))
    ax.text(xS - 2.25, 3.95, vals, ha="left", va="top", fontsize=8.7,
            bbox=dict(boxstyle="round", fc="#f7f7f7", ec="0.6", alpha=0.95))

    ax.set_xlim(xS - 2.4, xScr + 2.3)
    ax.set_ylim(-5.4, 4.15)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    cap = (r"de Broglie two-slit configuration (physical scale). A vibrating source $S$ emits "
           r"the localized odd-harmonic wave of (rough) wavelength $\lambda_0=h/\sqrt{2m_e eV}$; "
           r"two slits separated by $W$ sit at distance $L$; the screen is at $D$. "
           r"$\bf{Not\ to\ scale}$: $\lambda_0$, $W$ and $L$ are each drawn at an INDEPENDENT visual "
           r"scale because $W/\lambda_0\sim5\times10^4$ and $L/W\sim10^4$ make a common scale "
           r"unreadable; the true magnitudes are annotated. Only the topology "
           r"source$\to$slits$\to$screen and the paraxiality ($\theta\sim10^{-3\,\circ}$) are physical.")
    fig.text(0.5, 0.015, cap, ha="center", va="bottom", fontsize=8.8, wrap=True)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.99, bottom=0.17)

    outdir = os.path.dirname(os.path.abspath(__file__))
    base = "fig_setup_debroglie_V%03d" % int(V)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "%s.%s" % (base, ext)), dpi=200)
    print("saved %s.png / .svg   (lam0=%.4f nm, theta=%.2e deg, dX=%.1f um)"
          % (base, lam0 * 1e9, theta, dX * 1e6))


if __name__ == "__main__":
    main()
