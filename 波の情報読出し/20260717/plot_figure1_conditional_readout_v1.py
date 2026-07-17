#!/usr/bin/env python3
"""Generate Figure 1: the conditional N(R) readout under H2 and D4.

The empirical principal candidates R1 and R2 are inputs from the independent
full-range recurrence scan.  The inverse images of the two external diagnostic
values are computed from the same conditional formula and are therefore shown
as diagnostic positions, not as independent predictions.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, LinearLocator


R1 = 0.697177902556148
R2 = 0.688363902556148
ALPHA_INV_ZERO = 137.035999177
ALPHA_INV_MZ = 128.946

CURVE_COLOR = "#1455cc"
R1_COLOR = "#d62728"
R2_COLOR = "#e69f00"
TARGET_COLOR = "#6a3d9a"


def n_readout(r_value: float | np.ndarray) -> float | np.ndarray:
    """Conditional readout N(R)=4*pi/(1-R)^2 obtained under H2 and D4."""
    return 4.0 * math.pi / (1.0 - r_value) ** 2


def inverse_readout(n_value: float) -> float:
    """Inverse image within the same conditional readout formula."""
    return 1.0 - math.sqrt(4.0 * math.pi / n_value)


def create_figure(output_path: Path) -> None:
    r1_n = float(n_readout(R1))
    r2_n = float(n_readout(R2))
    r_zero = inverse_readout(ALPHA_INV_ZERO)
    r_mz = inverse_readout(ALPHA_INV_MZ)

    figure = plt.figure(figsize=(18.0, 10.0))
    grid = figure.add_gridspec(2, 2, width_ratios=(1.18, 1.0))
    full_axis = figure.add_subplot(grid[:, 0])
    low_axis = figure.add_subplot(grid[0, 1])
    high_axis = figure.add_subplot(grid[1, 1])

    # Panel A: full conditional curve and all four logically distinct points.
    r_full = np.linspace(0.686, 0.700, 2000)
    full_axis.plot(
        r_full,
        n_readout(r_full),
        color=CURVE_COLOR,
        linewidth=3.0,
        label=r"conditional readout $N(R)=4\pi/(1-R)^2$",
    )
    full_axis.axhline(
        ALPHA_INV_ZERO,
        color=R1_COLOR,
        linestyle="--",
        linewidth=1.8,
        alpha=0.65,
        label=r"diagnostic $\alpha^{-1}(0)=137.035999177$",
    )
    full_axis.axhline(
        ALPHA_INV_MZ,
        color=TARGET_COLOR,
        linestyle="--",
        linewidth=1.8,
        alpha=0.70,
        label=r"diagnostic $\alpha^{-1}(M_Z^2)=128.946$",
    )
    full_axis.scatter(
        [R1],
        [r1_n],
        s=150,
        facecolors="white",
        edgecolors=R1_COLOR,
        linewidths=3.0,
        zorder=5,
        label=rf"empirical $R_1$; $N(R_1)={r1_n:.6f}$",
    )
    full_axis.scatter(
        [R2],
        [r2_n],
        s=150,
        facecolors="white",
        edgecolors=R2_COLOR,
        linewidths=3.0,
        zorder=5,
        label=rf"empirical $R_2$; $N(R_2)={r2_n:.6f}$",
    )
    full_axis.scatter(
        [r_zero, r_mz],
        [ALPHA_INV_ZERO, ALPHA_INV_MZ],
        marker="D",
        s=85,
        color=TARGET_COLOR,
        edgecolors="white",
        linewidths=1.0,
        zorder=6,
        label="conditional inverse images of diagnostics",
    )
    full_axis.set_xlim(0.6857, 0.7003)
    full_axis.set_ylim(126.8, 140.1)
    full_axis.set_xlabel(r"reflection coefficient $R$", fontsize=14)
    full_axis.set_ylabel(r"conditional readout $N(R)\equiv\alpha^{-1}$", fontsize=14)
    full_axis.set_title("A. Full conditional readout", fontsize=16, weight="bold")
    full_axis.grid(True, alpha=0.25)
    full_axis.legend(loc="upper left", fontsize=10.5, framealpha=0.95)
    full_axis.text(
        0.03,
        0.03,
        r"Used: hypothesis H2 $G_R=e$; definition D4 $N:=\alpha^{-1}$",
        transform=full_axis.transAxes,
        fontsize=11,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.90},
    )

    # Panel B: resolve the very small difference between R1 and the low-energy
    # diagnostic inverse image.
    low_center = 0.5 * (R1 + r_zero)
    low_half_width = 1.0e-7
    r_low = np.linspace(low_center - low_half_width, low_center + low_half_width, 1000)
    low_axis.plot(r_low, n_readout(r_low), color=CURVE_COLOR, linewidth=3.0)
    low_axis.axhline(ALPHA_INV_ZERO, color=TARGET_COLOR, linestyle="--", linewidth=1.8)
    low_axis.axvline(r_zero, color=TARGET_COLOR, linestyle=":", linewidth=1.8)
    low_axis.scatter(
        [R1], [r1_n], s=165, facecolors="white", edgecolors=R1_COLOR,
        linewidths=3.0, zorder=5, label=r"empirical $R_1$",
    )
    low_axis.scatter(
        [r_zero], [ALPHA_INV_ZERO], marker="D", s=95, color=TARGET_COLOR,
        edgecolors="white", linewidths=1.0, zorder=6,
        label=r"inverse image of $\alpha^{-1}(0)$",
    )
    low_axis.set_xlim(low_center - low_half_width, low_center + low_half_width)
    low_values = n_readout(np.array([low_center - low_half_width, low_center + low_half_width]))
    low_axis.set_ylim(float(low_values.min()) - 1.5e-5, float(low_values.max()) + 1.5e-5)
    low_axis.xaxis.set_major_locator(LinearLocator(5))
    low_axis.xaxis.set_major_formatter(FormatStrFormatter("%.9f"))
    low_axis.yaxis.set_major_locator(LinearLocator(6))
    low_axis.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
    low_axis.set_xlabel(r"$R$", fontsize=13)
    low_axis.set_ylabel(r"$N(R)$", fontsize=13)
    low_axis.set_title("B. Low-energy diagnostic", fontsize=16, weight="bold")
    low_axis.grid(True, alpha=0.25)
    low_axis.legend(loc="upper left", fontsize=10)
    low_axis.text(
        0.98,
        0.05,
        "$\\Delta R=+2.33\\times10^{-8}$\n"
        "$\\Delta N=+2.11\\times10^{-5}$\n"
        "relative difference $=1.54\\times10^{-5}\\%$",
        transform=low_axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=10.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92},
    )

    # Panel C: show explicitly that R2 does not map to the M_Z diagnostic value.
    r_high = np.linspace(r_mz - 0.00020, R2 + 0.00020, 1200)
    high_axis.plot(r_high, n_readout(r_high), color=CURVE_COLOR, linewidth=3.0)
    high_axis.axhline(ALPHA_INV_MZ, color=TARGET_COLOR, linestyle="--", linewidth=1.8)
    high_axis.axvline(r_mz, color=TARGET_COLOR, linestyle=":", linewidth=1.8)
    high_axis.scatter(
        [R2], [r2_n], s=165, facecolors="white", edgecolors=R2_COLOR,
        linewidths=3.0, zorder=5, label=rf"empirical $R_2$: {r2_n:.6f}",
    )
    high_axis.scatter(
        [r_mz], [ALPHA_INV_MZ], marker="D", s=95, color=TARGET_COLOR,
        edgecolors="white", linewidths=1.0, zorder=6,
        label=rf"diagnostic inverse image: {ALPHA_INV_MZ:.3f}",
    )
    high_axis.set_xlim(r_mz - 0.00020, R2 + 0.00020)
    high_values = n_readout(np.array([r_mz - 0.00020, R2 + 0.00020]))
    high_axis.set_ylim(float(high_values.min()) - 0.08, float(high_values.max()) + 0.08)
    high_axis.ticklabel_format(axis="x", style="plain", useOffset=False)
    high_axis.set_xlabel(r"$R$", fontsize=13)
    high_axis.set_ylabel(r"$N(R)$", fontsize=13)
    high_axis.set_title("C. High-energy diagnostic", fontsize=16, weight="bold")
    high_axis.grid(True, alpha=0.25)
    high_axis.legend(loc="upper left", fontsize=10)
    high_axis.text(
        0.98,
        0.05,
        "$\\Delta R=+5.41\\times10^{-4}$\n"
        "$\\Delta N=+0.448063$\n"
        "relative difference $=0.3475\\%$",
        transform=high_axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=10.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92},
    )

    for axis in (full_axis, low_axis, high_axis):
        axis.tick_params(labelsize=11)

    figure.subplots_adjust(
        left=0.070,
        right=0.985,
        bottom=0.095,
        top=0.875,
        wspace=0.230,
        hspace=0.340,
    )
    figure.suptitle(
        r"Conditional coupling readout $N(R)$ under H2 and D4",
        fontsize=21,
        weight="bold",
        y=0.965,
    )
    figure.text(
        0.5,
        0.025,
        "Diagnostic inverse images use the same conditional formula; they are not independent derivations.",
        ha="center",
        fontsize=11.5,
        style="italic",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def main() -> None:
    default_output = Path(__file__).resolve().parent / "figures" / "fig1_N_R_derivation.png"
    parser = argparse.ArgumentParser(description="Generate the corrected Figure 1")
    parser.add_argument("--output", type=Path, default=default_output)
    args = parser.parse_args()
    create_figure(args.output)
    print(f"Saved corrected Figure 1 to {args.output}")


if __name__ == "__main__":
    main()
