#!/usr/bin/env python3
"""Phase 5 v2: exact eigenphase and finite-period resonance analysis.

The exchange map is known analytically, so its spectrum must be computed from
the exact 2x2 complex matrix rather than reconstructed by least squares from a
single trajectory.  This script verifies the neutral unit-modulus spectrum and
finds the nearest low-order roots of unity to the two empirical R candidates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
R1_CANDIDATE = 0.697177902556148
R2_CANDIDATE = 0.688363902556148
DEFAULT_MAX_ORDER = 256
EMPIRICAL_R_GRID_STEP = 1.0e-7


@dataclass(frozen=True)
class ResonanceMatch:
    label: str
    candidate_R: float
    period_n: int
    mode_m: int
    root_index_k: int
    resonant_R: float
    delta_R: float
    relative_delta_R: float
    candidate_residual: float
    exact_residual: float
    matrix_return_error: float
    eigenvalue_real: float
    eigenvalue_imag: float
    eigenvalue_modulus: float
    eigenphase_cycles: float


def theta_from_R(reflection_rate: float) -> float:
    if not 0.0 <= reflection_rate <= 1.0:
        raise ValueError(f"R must be in [0, 1]: {reflection_rate}")
    return math.asin(math.sqrt(reflection_rate))


def scattering_coefficients(reflection_rate: float) -> tuple[complex, complex]:
    theta = theta_from_R(reflection_rate)
    phase = complex(math.cos(theta), math.sin(theta))
    transmission = phase * math.cos(theta)
    reflection = -1j * phase * math.sin(theta)
    return complex(transmission), complex(reflection)


def scattering_matrix(reflection_rate: float) -> np.ndarray:
    transmission, reflection = scattering_coefficients(reflection_rate)
    return np.array(
        [[reflection, transmission], [transmission, reflection]],
        dtype=np.complex128,
    )


def exact_eigenvalues(reflection_rate: float) -> tuple[complex, complex]:
    """Return the symmetric and antisymmetric eigenvalues exactly."""
    theta = theta_from_R(reflection_rate)
    lambda_symmetric = 1.0 + 0.0j
    lambda_antisymmetric = -complex(math.cos(2.0 * theta), math.sin(2.0 * theta))
    return lambda_symmetric, lambda_antisymmetric


def resonant_R(period_n: int, mode_m: int) -> float:
    """R for which lambda_a is an exact n-th root of unity.

    The fundamental-order representation uses 1 <= m < n/2 and gcd(m,n)=1.
    The eigenphase root index is k=n-m.
    """
    if period_n < 2:
        raise ValueError("period_n must be at least 2")
    if not 1 <= mode_m < period_n / 2:
        raise ValueError("mode_m must satisfy 1 <= m < n/2")
    return math.cos(math.pi * mode_m / period_n) ** 2


def resonance_residual(reflection_rate: float, period_n: int) -> float:
    _, lambda_a = exact_eigenvalues(reflection_rate)
    return float(abs(lambda_a**period_n - 1.0))


def matrix_return_error(reflection_rate: float, period_n: int) -> float:
    matrix = scattering_matrix(reflection_rate)
    returned = np.linalg.matrix_power(matrix, period_n)
    return float(np.linalg.norm(returned - np.eye(2), ord="fro"))


def find_nearest_low_order_resonance(
    label: str,
    candidate_R: float,
    max_order: int = DEFAULT_MAX_ORDER,
) -> ResonanceMatch:
    best: tuple[float, int, int, float] | None = None
    for period_n in range(3, max_order + 1):
        for mode_m in range(1, (period_n - 1) // 2 + 1):
            if math.gcd(mode_m, period_n) != 1:
                continue
            exact_R = resonant_R(period_n, mode_m)
            error = abs(exact_R - candidate_R)
            candidate = (error, period_n, mode_m, exact_R)
            if best is None or candidate < best:
                best = candidate
    if best is None:
        raise RuntimeError("no resonance candidates were generated")

    _, period_n, mode_m, exact_R = best
    _, lambda_a = exact_eigenvalues(candidate_R)
    eigenphase_cycles = (math.atan2(lambda_a.imag, lambda_a.real) / (2.0 * math.pi)) % 1.0
    return ResonanceMatch(
        label=label,
        candidate_R=candidate_R,
        period_n=period_n,
        mode_m=mode_m,
        root_index_k=period_n - mode_m,
        resonant_R=exact_R,
        delta_R=candidate_R - exact_R,
        relative_delta_R=(candidate_R - exact_R) / exact_R,
        candidate_residual=resonance_residual(candidate_R, period_n),
        exact_residual=resonance_residual(exact_R, period_n),
        matrix_return_error=matrix_return_error(exact_R, period_n),
        eigenvalue_real=float(lambda_a.real),
        eigenvalue_imag=float(lambda_a.imag),
        eigenvalue_modulus=float(abs(lambda_a)),
        eigenphase_cycles=eigenphase_cycles,
    )


def build_resonance_curves(
    matches: list[ResonanceMatch],
    half_width: float = 1.0e-6,
    points: int = 4001,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for match in matches:
        r_values = np.linspace(match.resonant_R - half_width, match.resonant_R + half_width, points)
        theta = np.arcsin(np.sqrt(r_values))
        lambda_a = -np.exp(2j * theta)
        residuals = np.abs(lambda_a**match.period_n - 1.0)
        for reflection_rate, residual in zip(r_values, residuals):
            rows.append(
                {
                    "label": match.label,
                    "period_n": match.period_n,
                    "mode_m": match.mode_m,
                    "R": float(reflection_rate),
                    "delta_R_micro": float((reflection_rate - match.resonant_R) * 1.0e6),
                    "residual": float(residual),
                }
            )
    return rows


def plot_resonance_residuals(
    matches: list[ResonanceMatch],
    curves: list[dict[str, float | int | str]],
    output: Path,
) -> None:
    colors = {"R1": "#d94b4b", "R2": "#3478bf"}
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.1))
    for axis, match in zip(axes, matches):
        selected = [row for row in curves if row["label"] == match.label]
        x = np.array([float(row["delta_R_micro"]) for row in selected])
        y = np.maximum(np.array([float(row["residual"]) for row in selected]), 1.0e-15)
        color = colors[match.label]
        candidate_x = match.delta_R * 1.0e6
        half_grid_bin_x = 0.5 * EMPIRICAL_R_GRID_STEP * 1.0e6
        axis.semilogy(x, y, color=color, linewidth=2.0)
        axis.axvspan(
            candidate_x - half_grid_bin_x,
            candidate_x + half_grid_bin_x,
            color=color,
            alpha=0.16,
            linewidth=0.0,
            label=r"empirical grid bin ($\Delta R=10^{-7}$)",
        )
        axis.axvline(0.0, color="#222222", linestyle="-", linewidth=1.1, label="analytic root")
        axis.axvline(candidate_x, color=color, linestyle="--", linewidth=1.4, label="empirical candidate")
        axis.scatter([candidate_x], [match.candidate_residual], color=color, edgecolor="black", zorder=5)
        axis.set_xlabel(r"$(R-R_{n,m})\times10^6$")
        axis.set_ylabel(r"$d_n(R)=|\lambda_a(R)^n-1|$")
        axis.set_title(
            rf"{match.label}: $n={match.period_n}$, $m={match.mode_m}$" "\n"
            rf"$R_{{n,m}}={match.resonant_R:.15f}$"
        )
        axis.grid(alpha=0.25, which="both")
        axis.legend(loc="upper center")
        axis.text(
            0.03,
            0.04,
            rf"$R_j-R_{{n,m}}={match.delta_R:+.2e}$" "\n"
            rf"$d_n(R_j)={match.candidate_residual:.2e}$",
            transform=axis.transAxes,
            fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.90},
        )
    fig.suptitle("Eigenphase resonance residuals and empirical scan resolution", fontsize=14)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_eigenphase_orbits(matches: list[ResonanceMatch], output: Path) -> None:
    colors = {"R1": "#d94b4b", "R2": "#3478bf"}
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.7))
    angle = np.linspace(0.0, 2.0 * np.pi, 1000)
    for axis, match in zip(axes, matches):
        _, lambda_root = exact_eigenvalues(match.resonant_R)
        powers = lambda_root ** np.arange(match.period_n)
        color_values = np.arange(match.period_n)
        axis.plot(np.cos(angle), np.sin(angle), color="#777777", linestyle="--", linewidth=1.0)
        scatter = axis.scatter(
            powers.real,
            powers.imag,
            c=color_values,
            cmap="viridis",
            s=22,
            edgecolors="none",
        )
        axis.scatter([1.0], [0.0], marker="*", s=180, color="#f2b134", edgecolor="black", zorder=5)
        axis.scatter(
            [lambda_root.real],
            [lambda_root.imag],
            s=80,
            color=colors[match.label],
            edgecolor="black",
            zorder=5,
        )
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(-1.12, 1.12)
        axis.set_ylim(-1.12, 1.12)
        axis.axhline(0.0, color="#bbbbbb", linewidth=0.7)
        axis.axvline(0.0, color="#bbbbbb", linewidth=0.7)
        axis.set_xlabel("Re")
        axis.set_ylabel("Im")
        axis.set_title(
            rf"{match.label}-associated root: {match.period_n} steps" "\n"
            rf"$\lambda_a=e^{{2\pi i\,{match.root_index_k}/{match.period_n}}}$"
        )
        colorbar = fig.colorbar(scatter, ax=axis, fraction=0.046, pad=0.04)
        colorbar.set_label("power $j$")
    fig.suptitle("Closed eigenphase orbits of the exact exchange map", fontsize=14)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_quarter_period_structure(match: ResonanceMatch, output: Path) -> None:
    """Plot the exact 31-step stroboscopic phases of the 124-step root.

    The old Figure 5 plotted a period-31 sine wave, which returned to the same
    scalar value at every marked checkpoint.  That does not represent the
    derived sequence 1 -> i -> -1 -> -i -> 1.  This figure evaluates the exact
    root of unity and displays both complex components and the stroboscopic
    orbit explicitly.
    """
    if (match.period_n, match.mode_m) != (124, 23):
        raise ValueError("Figure 5 requires the R1-associated (n,m)=(124,23) root")

    steps = np.arange(match.period_n + 1)
    lambda_root = np.exp(2j * np.pi * match.root_index_k / match.period_n)
    powers = lambda_root**steps
    checkpoints = np.array([0, 31, 62, 93, 124], dtype=int)
    checkpoint_values = lambda_root**checkpoints

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.6))

    time_axis = axes[0]
    time_axis.plot(steps, powers.real, color="#d84b4b", linewidth=1.6, label=r"$\Re(\lambda_a^n)$")
    time_axis.plot(steps, powers.imag, color="#3478bf", linewidth=1.6, label=r"$\Im(\lambda_a^n)$")
    for step in checkpoints:
        time_axis.axvline(step, color="#777777", linestyle="--", linewidth=0.8, alpha=0.65)
    time_axis.scatter(
        checkpoints,
        checkpoint_values.real,
        color="#d84b4b",
        edgecolor="black",
        s=48,
        zorder=5,
    )
    time_axis.scatter(
        checkpoints,
        checkpoint_values.imag,
        color="#3478bf",
        edgecolor="black",
        s=48,
        zorder=5,
    )
    time_axis.set_xlim(0, match.period_n)
    time_axis.set_ylim(-1.15, 1.15)
    time_axis.set_xticks(checkpoints)
    time_axis.set_xlabel(r"iteration step $n$")
    time_axis.set_ylabel("eigenphase component")
    time_axis.set_title("Exact 124-step eigenphase components")
    time_axis.grid(alpha=0.22)
    time_axis.legend(loc="upper right")

    phase_axis = axes[1]
    phase_axis.plot(np.cos(np.linspace(0.0, 2.0 * np.pi, 600)),
                    np.sin(np.linspace(0.0, 2.0 * np.pi, 600)),
                    color="#aaaaaa", linestyle="--", linewidth=1.0)
    unique_values = checkpoint_values[:4]
    closed_values = np.append(unique_values, unique_values[0])
    phase_axis.plot(closed_values.real, closed_values.imag, color="#444444", linewidth=1.3)
    phase_axis.scatter(
        unique_values.real,
        unique_values.imag,
        c=["#f2b134", "#3478bf", "#555555", "#d84b4b"],
        edgecolor="black",
        s=125,
        zorder=5,
    )
    phase_labels = [
        (1.0, 0.0, r"$n=0,124:\ 1$", (8, 8)),
        (0.0, 1.0, r"$n=31:\ i$", (8, -20)),
        (-1.0, 0.0, r"$n=62:\ -1$", (8, 8)),
        (0.0, -1.0, r"$n=93:\ -i$", (8, 8)),
    ]
    for x, y, label, offset in phase_labels:
        phase_axis.annotate(label, (x, y), xytext=offset, textcoords="offset points", fontsize=11)
    phase_axis.axhline(0.0, color="#bbbbbb", linewidth=0.7)
    phase_axis.axvline(0.0, color="#bbbbbb", linewidth=0.7)
    phase_axis.set_aspect("equal", adjustable="box")
    phase_axis.set_xlim(-1.28, 1.28)
    phase_axis.set_ylim(-1.28, 1.28)
    phase_axis.set_xlabel(r"$\Re(\lambda_a^n)$")
    phase_axis.set_ylabel(r"$\Im(\lambda_a^n)$")
    phase_axis.set_title("31-step stroboscopic four-phase cycle")
    phase_axis.grid(alpha=0.18)

    fig.suptitle(
        r"Exact $(n,m)=(124,23)$ resonance: "
        r"$1\rightarrow i\rightarrow-1\rightarrow-i\rightarrow1$",
        fontsize=15,
    )
    fig.text(
        0.5,
        0.012,
        rf"exact root $R_{{124,23}}={match.resonant_R:.15f}$; "
        rf"empirical $R_1$ residual $|\lambda_a(R_1)^{{124}}-1|={match.candidate_residual:.3e}$",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 0.94))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_curve_csv(rows: list[dict[str, float | int | str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=DEFAULT_MAX_ORDER)
    parser.add_argument("--output-dir", type=Path, default=HERE / "phase5_resonance_v2_results")
    parser.add_argument("--figure-dir", type=Path, default=HERE / "figures")
    args = parser.parse_args()

    matches = [
        find_nearest_low_order_resonance("R1", R1_CANDIDATE, args.max_order),
        find_nearest_low_order_resonance("R2", R2_CANDIDATE, args.max_order),
    ]
    curves = build_resonance_curves(matches)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "phase5_resonance_summary_v2.json").open("w") as handle:
        json.dump(
            {
                "model": "exact exchange-matrix eigenphase resonance",
                "max_order": args.max_order,
                "empirical_R_grid_step": EMPIRICAL_R_GRID_STEP,
                "analytic_spectrum": {
                    "lambda_s": "1",
                    "lambda_a": "-exp(2 i theta)",
                    "theta": "asin(sqrt(R))",
                },
                "search_scope": "coprime 1 <= m < n/2, 3 <= n <= max_order",
                "derived_phase_relations": {
                    "R1_associated_124_root": {
                        "lambda_a^31": "i",
                        "lambda_a^62": "-1",
                        "lambda_a^93": "-i",
                        "lambda_a^124": "1",
                    },
                    "R2_associated_122_root": {
                        "lambda_a^61": "-1",
                        "lambda_a^122": "1",
                    },
                },
                "matches": [asdict(match) for match in matches],
            },
            handle,
            indent=2,
        )
    write_curve_csv(curves, args.output_dir / "phase5_resonance_curves_v2.csv")
    plot_resonance_residuals(matches, curves, args.figure_dir / "fig2_floquet_stability.png")
    plot_eigenphase_orbits(matches, args.figure_dir / "fig3_eigenvalue_spectrum.png")
    plot_quarter_period_structure(matches[0], args.figure_dir / "fig5_31_period_structure.png")

    print("Phase 5 v2 complete")
    for match in matches:
        print(
            f"{match.label}: R={match.candidate_R:.15f} -> "
            f"cos^2({match.mode_m} pi/{match.period_n})={match.resonant_R:.15f}, "
            f"delta_R={match.delta_R:+.3e}, "
            f"|lambda_a^{match.period_n}-1|={match.candidate_residual:.3e}"
        )


if __name__ == "__main__":
    main()
