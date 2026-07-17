#!/usr/bin/env python3
"""Independent 4D geometric-running experiment for the 1+8+128 model.

This program does not import or reuse the Phase 4--6 scattering programs.
It tests the following separate model:

1. Enumerate all integer-centred unit 4-cubes fully contained in the
   four-ball of radius 3.  The expected decomposition is 1 + 8 + 128.
2. Treat the centre-plus-inner layer (9 cells) and the outer layer
   (128 cells) as two finite 4D scattering sets.
3. Calculate their isotropically averaged mutual coherence exactly from
   pair-distance histograms.  In four dimensions,

       <exp(i q n_hat . r)>_{n_hat in S^3} = 2 J_1(q r) / (q r).

4. Use the squared normalized mutual coherence as the visibility of the
   inner 9 cells as external channels.
5. Weight the continuous correction by the squared form factor of the
   unit four-ball,

       F_ball(q) = 8 J_2(q) / q^2.

6. Solve, at every dimensionless wave number q, the self-consistency
   equation

       1/alpha(q) = 128 + 9 W_core(q)
                    + (pi^2/2) alpha(q) W_ball(q).

The cell count is an exact consequence of the Paper 7 containment rule.
Steps 4--6 are explicit model hypotheses to be tested, not derivations.
Reference alpha values are used only after the curve is constructed, as
diagnostic crossing levels.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j1, jv


GEOMETRIC_RADIUS = 3.0
CELL_WIDTH = 1.0
HALF_WIDTH = CELL_WIDTH / 2.0
V4_UNIT = math.pi**2 / 2.0

DEFAULT_TARGETS = (
    137.035999177,
    129.394062925467,
    128.946,
)

Cell = Tuple[int, int, int, int]
Histogram = Dict[int, int]


def containment_value(cell: Cell) -> float:
    """Maximum squared radius attained by a unit 4-cube centred at cell."""
    return sum((abs(coordinate) + HALF_WIDTH) ** 2 for coordinate in cell)


def absolute_pattern(cell: Cell) -> Tuple[int, int, int, int]:
    return tuple(sorted(abs(coordinate) for coordinate in cell))


def enumerate_paper7_cells() -> List[Cell]:
    """Enumerate the Paper 7 cells from the containment inequality alone."""
    coordinate_limit = math.ceil(GEOMETRIC_RADIUS)
    cells = []
    for cell in product(range(-coordinate_limit, coordinate_limit + 1), repeat=4):
        if containment_value(cell) <= GEOMETRIC_RADIUS**2 + 1.0e-12:
            cells.append(tuple(int(value) for value in cell))
    cells.sort()
    return cells


def layer_name(cell: Cell) -> str:
    pattern = absolute_pattern(cell)
    if pattern == (0, 0, 0, 0):
        return "center_1"
    if pattern == (0, 0, 0, 1):
        return "inner_8"
    return "outer_128"


def validate_paper7_decomposition(cells: Sequence[Cell]) -> Dict[str, object]:
    layer_counts = Counter(layer_name(cell) for cell in cells)
    pattern_counts = Counter(absolute_pattern(cell) for cell in cells)
    outer_cells = [cell for cell in cells if layer_name(cell) == "outer_128"]
    outer_set = set(outer_cells)

    expected_layers = {
        "center_1": 1,
        "inner_8": 8,
        "outer_128": 128,
    }
    expected_patterns = {
        (0, 0, 0, 0): 1,
        (0, 0, 0, 1): 8,
        (0, 0, 1, 1): 24,
        (0, 1, 1, 1): 32,
        (1, 1, 1, 1): 16,
        (0, 0, 0, 2): 8,
        (0, 0, 1, 2): 48,
    }

    if len(cells) != 137:
        raise AssertionError(f"Paper 7 enumeration produced {len(cells)} cells, not 137")
    if dict(layer_counts) != expected_layers:
        raise AssertionError(f"Unexpected 1+8+128 decomposition: {dict(layer_counts)}")
    if dict(pattern_counts) != expected_patterns:
        raise AssertionError(f"Unexpected absolute-value patterns: {dict(pattern_counts)}")
    if any(tuple(-value for value in cell) not in outer_set for cell in outer_cells):
        raise AssertionError("Outer layer is not closed under central inversion")

    antipodal_representatives = [
        cell
        for cell in outer_cells
        if cell > tuple(-value for value in cell)
    ]
    if len(antipodal_representatives) != 64:
        raise AssertionError(
            "Outer layer did not decompose into 64 antipodal pairs"
        )

    return {
        "total": len(cells),
        "layers": expected_layers,
        "outer_antipodal_decomposition": {
            "pair_count": len(antipodal_representatives),
            "oriented_halves": [len(antipodal_representatives)] * 2,
        },
        "patterns": {
            ",".join(str(value) for value in pattern): count
            for pattern, count in sorted(pattern_counts.items())
        },
    }


def squared_distance(a: Cell, b: Cell) -> int:
    return sum((left - right) ** 2 for left, right in zip(a, b))


def pair_distance_histogram(left: Sequence[Cell], right: Sequence[Cell]) -> Histogram:
    histogram: Counter[int] = Counter()
    for left_cell in left:
        for right_cell in right:
            histogram[squared_distance(left_cell, right_cell)] += 1
    return dict(sorted(histogram.items()))


def four_dimensional_powder_kernel(q_values: np.ndarray, distance: float) -> np.ndarray:
    """Return 2 J1(q r)/(q r), with the continuous value 1 at q r = 0."""
    argument = q_values * distance
    values = np.ones_like(argument, dtype=float)
    nonzero = np.abs(argument) > 1.0e-12
    values[nonzero] = 2.0 * j1(argument[nonzero]) / argument[nonzero]
    return values


def powder_pair_average(
    q_values: np.ndarray,
    histogram: Histogram,
    pair_count: int,
) -> np.ndarray:
    result = np.zeros_like(q_values, dtype=float)
    for distance_squared, count in histogram.items():
        distance = math.sqrt(distance_squared)
        result += count * four_dimensional_powder_kernel(q_values, distance)
    return result / pair_count


def unit_four_ball_amplitude(q_values: np.ndarray) -> np.ndarray:
    """Normalized Fourier amplitude 8 J2(q)/q^2 of the unit four-ball."""
    result = np.ones_like(q_values, dtype=float)
    nonzero = np.abs(q_values) > 1.0e-12
    q_nonzero = q_values[nonzero]
    result[nonzero] = 8.0 * jv(2, q_nonzero) / q_nonzero**2
    return result


def solve_self_consistent_alpha(
    outer_count: int,
    core_count: int,
    core_visibility: np.ndarray,
    continuous_visibility: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    discrete_count = outer_count + core_count * core_visibility
    continuous_coefficient = V4_UNIT * continuous_visibility
    discriminant = discrete_count**2 + 4.0 * continuous_coefficient

    if np.any(discriminant <= 0.0):
        raise ValueError("Self-consistency equation has no positive real solution")

    # Algebraically equivalent to the positive quadratic root, but stable
    # when the continuous coefficient approaches zero.
    inverse_alpha = 0.5 * (discrete_count + np.sqrt(discriminant))
    alpha = 1.0 / inverse_alpha
    return alpha, inverse_alpha, discrete_count


def find_crossings(q_values: np.ndarray, values: np.ndarray, target: float) -> List[float]:
    crossings: List[float] = []
    shifted = values - target
    for index in range(len(q_values) - 1):
        left = shifted[index]
        right = shifted[index + 1]
        if left == 0.0:
            crossings.append(float(q_values[index]))
            continue
        if left * right < 0.0:
            fraction = -left / (right - left)
            crossings.append(float(q_values[index] + fraction * (q_values[index + 1] - q_values[index])))
    if shifted[-1] == 0.0:
        crossings.append(float(q_values[-1]))
    return crossings


def nearest_point(q_values: np.ndarray, values: np.ndarray, target: float) -> Dict[str, float]:
    index = int(np.argmin(np.abs(values - target)))
    return {
        "q": float(q_values[index]),
        "inverse_alpha": float(values[index]),
        "difference": float(values[index] - target),
    }


def save_cell_catalog(path: Path, cells: Sequence[Cell]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x1", "x2", "x3", "x4", "layer", "absolute_pattern", "containment_value"])
        for cell in cells:
            writer.writerow(
                [
                    *cell,
                    layer_name(cell),
                    " ".join(str(value) for value in absolute_pattern(cell)),
                    f"{containment_value(cell):.12g}",
                ]
            )


def save_distance_histograms(
    path: Path,
    histograms: Iterable[Tuple[str, Histogram, int]],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pair_type", "distance_squared", "distance", "pair_count", "pair_fraction"])
        for pair_type, histogram, total_pairs in histograms:
            for distance_squared, count in histogram.items():
                writer.writerow(
                    [
                        pair_type,
                        distance_squared,
                        f"{math.sqrt(distance_squared):.12g}",
                        count,
                        f"{count / total_pairs:.16g}",
                    ]
                )


def save_curve(
    path: Path,
    q_values: np.ndarray,
    core_intensity: np.ndarray,
    outer_intensity: np.ndarray,
    cross_amplitude: np.ndarray,
    normalized_coherence: np.ndarray,
    core_visibility: np.ndarray,
    ball_amplitude: np.ndarray,
    ball_visibility: np.ndarray,
    discrete_count: np.ndarray,
    alpha: np.ndarray,
    inverse_alpha: np.ndarray,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "q",
                "core_intensity",
                "outer_intensity",
                "cross_amplitude",
                "normalized_coherence",
                "core_visibility",
                "ball_amplitude",
                "ball_visibility",
                "effective_discrete_count",
                "alpha",
                "inverse_alpha",
            ]
        )
        for row in zip(
            q_values,
            core_intensity,
            outer_intensity,
            cross_amplitude,
            normalized_coherence,
            core_visibility,
            ball_amplitude,
            ball_visibility,
            discrete_count,
            alpha,
            inverse_alpha,
        ):
            writer.writerow([f"{float(value):.16g}" for value in row])


def plot_results(
    path: Path,
    q_values: np.ndarray,
    inverse_alpha: np.ndarray,
    core_visibility: np.ndarray,
    ball_visibility: np.ndarray,
    targets: Sequence[float],
) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

    axes[0].plot(q_values, inverse_alpha, color="black", linewidth=1.5, label=r"$\alpha^{-1}(q)$")
    target_colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for index, target in enumerate(targets):
        axes[0].axhline(
            target,
            color=target_colors[index % len(target_colors)],
            linewidth=1.0,
            linestyle="--",
            alpha=0.8,
            label=f"diagnostic {target:.9g}",
        )
    axes[0].set_ylabel(r"$\alpha^{-1}$")
    axes[0].set_title("Independent 4D geometric running: 1 + 8 + 128")
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="best", fontsize=8)

    axes[1].plot(q_values, core_visibility, label=r"$W_{\rm core}$", color="tab:purple")
    axes[1].plot(q_values, ball_visibility, label=r"$W_{\rm ball}$", color="tab:cyan")
    axes[1].set_xlabel("dimensionless wave number q")
    axes[1].set_ylabel("visibility")
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].grid(alpha=0.25)
    axes[1].legend(loc="best")

    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def parse_targets(raw_targets: str) -> Tuple[float, ...]:
    return tuple(float(item.strip()) for item in raw_targets.split(",") if item.strip())


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Independent Paper 7 geometric-running experiment for 137 -> 128"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "results_geometric_running_137_to_128_v1",
    )
    parser.add_argument("--q-max", type=float, default=40.0)
    parser.add_argument("--q-points", type=int, default=4001)
    parser.add_argument(
        "--targets",
        default=",".join(str(value) for value in DEFAULT_TARGETS),
        help="Comma-separated post-construction diagnostic levels",
    )
    args = parser.parse_args()

    if args.q_max <= 0.0:
        raise ValueError("--q-max must be positive")
    if args.q_points < 2:
        raise ValueError("--q-points must be at least 2")

    targets = parse_targets(args.targets)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cells = enumerate_paper7_cells()
    decomposition = validate_paper7_decomposition(cells)
    core_cells = [cell for cell in cells if layer_name(cell) != "outer_128"]
    outer_cells = [cell for cell in cells if layer_name(cell) == "outer_128"]

    core_histogram = pair_distance_histogram(core_cells, core_cells)
    outer_histogram = pair_distance_histogram(outer_cells, outer_cells)
    cross_histogram = pair_distance_histogram(core_cells, outer_cells)

    q_values = np.linspace(0.0, args.q_max, args.q_points)
    core_intensity = powder_pair_average(
        q_values, core_histogram, len(core_cells) * len(core_cells)
    )
    outer_intensity = powder_pair_average(
        q_values, outer_histogram, len(outer_cells) * len(outer_cells)
    )
    cross_amplitude = powder_pair_average(
        q_values, cross_histogram, len(core_cells) * len(outer_cells)
    )

    denominator = np.sqrt(np.maximum(core_intensity * outer_intensity, 1.0e-30))
    normalized_coherence = cross_amplitude / denominator
    core_visibility = np.clip(normalized_coherence**2, 0.0, 1.0)

    ball_amplitude = unit_four_ball_amplitude(q_values)
    ball_visibility = ball_amplitude**2

    alpha, inverse_alpha, discrete_count = solve_self_consistent_alpha(
        outer_count=len(outer_cells),
        core_count=len(core_cells),
        core_visibility=core_visibility,
        continuous_visibility=ball_visibility,
    )

    diagnostics = {}
    for target in targets:
        crossings = find_crossings(q_values, inverse_alpha, target)
        diagnostics[f"{target:.15g}"] = {
            "crossings": crossings,
            "first_crossing": crossings[0] if crossings else None,
            "nearest_grid_point": nearest_point(q_values, inverse_alpha, target),
        }

    low_energy_expected = 0.5 * (137.0 + math.sqrt(137.0**2 + 4.0 * V4_UNIT))
    minimum_index = int(np.argmin(inverse_alpha))
    maximum_index = int(np.argmax(inverse_alpha))
    cross_amplitude_zeroes = find_crossings(q_values, cross_amplitude, 0.0)
    summary = {
        "model": "independent_geometric_running_137_to_128_v1",
        "independent_of": [
            "phase4_basin_of_attraction_v1.py",
            "phase5_eigenvalue_fixed_point_v1.py",
            "phase6_N_R_derivation_v1.py",
        ],
        "geometry": {
            "dimension": 4,
            "radius": GEOMETRIC_RADIUS,
            "cell_width": CELL_WIDTH,
            "decomposition": decomposition,
        },
        "explicit_model_hypotheses": {
            "core_visibility": "squared normalized isotropic mutual coherence of core_9 and outer_128",
            "continuous_visibility": "squared normalized Fourier amplitude of the unit four-ball",
            "self_consistency": "1/alpha = 128 + 9 W_core + (pi^2/2) alpha W_ball",
        },
        "scan": {
            "q_min": float(q_values[0]),
            "q_max": float(q_values[-1]),
            "q_points": len(q_values),
        },
        "limits_and_checks": {
            "inverse_alpha_q0": float(inverse_alpha[0]),
            "paper7_inverse_alpha_expected": low_energy_expected,
            "q0_difference": float(inverse_alpha[0] - low_energy_expected),
            "core_visibility_q0": float(core_visibility[0]),
            "ball_visibility_q0": float(ball_visibility[0]),
            "inverse_alpha_qmax": float(inverse_alpha[-1]),
            "core_visibility_qmax": float(core_visibility[-1]),
            "ball_visibility_qmax": float(ball_visibility[-1]),
            "minimum_inverse_alpha": float(np.min(inverse_alpha)),
            "q_at_minimum_inverse_alpha": float(q_values[minimum_index]),
            "maximum_inverse_alpha": float(np.max(inverse_alpha)),
            "q_at_maximum_inverse_alpha": float(q_values[maximum_index]),
            "monotonic_nonincreasing_on_full_scan": bool(
                np.all(np.diff(inverse_alpha) <= 1.0e-12)
            ),
            "cross_amplitude_zeroes": cross_amplitude_zeroes,
        },
        "diagnostic_targets_not_used_to_construct_curve": diagnostics,
        "unfixed_physical_mapping": "q = Q * ell / (hbar c); ell is not set by this experiment",
    }

    save_cell_catalog(args.output_dir / "cell_catalog_v1.csv", cells)
    save_distance_histograms(
        args.output_dir / "distance_histograms_v1.csv",
        [
            ("core_core", core_histogram, len(core_cells) * len(core_cells)),
            ("outer_outer", outer_histogram, len(outer_cells) * len(outer_cells)),
            ("core_outer", cross_histogram, len(core_cells) * len(outer_cells)),
        ],
    )
    save_curve(
        args.output_dir / "geometric_running_curve_v1.csv",
        q_values,
        core_intensity,
        outer_intensity,
        cross_amplitude,
        normalized_coherence,
        core_visibility,
        ball_amplitude,
        ball_visibility,
        discrete_count,
        alpha,
        inverse_alpha,
    )
    with (args.output_dir / "geometric_running_summary_v1.json").open("w") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    plot_results(
        args.output_dir / "geometric_running_curve_v1.png",
        q_values,
        inverse_alpha,
        core_visibility,
        ball_visibility,
        targets,
    )

    print("Independent 4D geometric-running experiment")
    print(f"  cells: {len(cells)} = {len(core_cells)} + {len(outer_cells)}")
    print(f"  inverse alpha at q=0: {inverse_alpha[0]:.12f}")
    print(f"  inverse alpha at q={q_values[-1]:g}: {inverse_alpha[-1]:.12f}")
    for target, result in diagnostics.items():
        first_crossing = result["first_crossing"]
        if first_crossing is None:
            print(f"  target {target}: no crossing in scan")
        else:
            print(f"  target {target}: first crossing q={first_crossing:.12f}")
    print(f"  outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
