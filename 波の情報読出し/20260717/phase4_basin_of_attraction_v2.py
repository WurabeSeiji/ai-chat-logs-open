#!/usr/bin/env python3
"""Phase 4 v2: principal-candidate selection and conditional basin analysis.

The fixed-R exchange map used in Phase 4 v1 is unitary and does not evolve R.
It therefore cannot, by itself, possess attracting basins in R-space.  This
revision separates two operations:

1. Under the principal-pair hypothesis H7, extract the seven recurrence bands
   from the 2026-07-15 sweep and select ranks 1 and 2 by the largest adjacent
   gap in ranked gray depth.
2. Under the explicitly additional effective-relaxation hypothesis H8, evolve
   R (in a reduced coordinate x) in a damped empirical double-well and measure
   the basins generated inside that conditional model.

The lower five bands are retained in the output as finite-window alias
candidates; they are not silently deleted and are not treated as independent
attractors by the reduced model.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap


HERE = Path(__file__).resolve().parent
DEFAULT_CANDIDATE_CSV = (
    HERE.parent
    / "20260715"
    / "minimal_system_B_gray_bugcheck_result_v1"
    / "direct_depth_probe_v5_sweep_control"
    / "high_to_ext_full_delta1e-7_candidates_v5.csv"
)


@dataclass(frozen=True)
class CandidateBand:
    rank: int
    r_start: float
    r_end: float
    r_peak: float
    depth_peak: float
    error_peak: float
    best_step: int
    sample_count: int
    classification: str


def _read_candidate_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "R": float(raw["R"]),
                    "depth": float(raw["best_prefix_gray_depth_no_phase"]),
                    "error": float(raw["best_prefix_gray_error_no_phase"]),
                    "best_step": float(raw["best_step"]),
                }
            )
    if not rows:
        raise ValueError(f"candidate CSV is empty: {path}")
    return sorted(rows, key=lambda row: row["R"])


def _split_contiguous_bands(rows: list[dict[str, float]]) -> tuple[list[list[dict[str, float]]], float]:
    r_values = np.array([row["R"] for row in rows], dtype=float)
    positive_gaps = np.diff(r_values)
    positive_gaps = positive_gaps[positive_gaps > 0]
    if positive_gaps.size == 0:
        raise ValueError("candidate CSV must contain at least two distinct R values")

    grid_step = float(np.median(positive_gaps))
    gap_threshold = 5.0 * grid_step
    bands: list[list[dict[str, float]]] = []
    current = [rows[0]]
    for previous, row in zip(rows, rows[1:]):
        if row["R"] - previous["R"] > gap_threshold:
            bands.append(current)
            current = []
        current.append(row)
    bands.append(current)
    return bands, gap_threshold


def extract_candidate_bands(path: Path) -> tuple[list[CandidateBand], dict[str, float | int]]:
    rows = _read_candidate_rows(path)
    raw_bands, gap_threshold = _split_contiguous_bands(rows)

    summaries: list[dict[str, float | int]] = []
    for band in raw_bands:
        peak = max(band, key=lambda row: row["depth"])
        summaries.append(
            {
                "r_start": band[0]["R"],
                "r_end": band[-1]["R"],
                "r_peak": peak["R"],
                "depth_peak": peak["depth"],
                "error_peak": peak["error"],
                "best_step": int(peak["best_step"]),
                "sample_count": len(band),
            }
        )

    summaries.sort(key=lambda item: float(item["depth_peak"]), reverse=True)
    ranked_depths = np.array([float(item["depth_peak"]) for item in summaries])
    depth_gaps = ranked_depths[:-1] - ranked_depths[1:]
    if depth_gaps.size == 0:
        raise ValueError("at least two recurrence bands are required")
    principal_count = int(np.argmax(depth_gaps)) + 1
    largest_depth_gap = float(depth_gaps[principal_count - 1])

    bands = [
        CandidateBand(
            rank=index,
            classification="principal" if index <= principal_count else "finite-window alias candidate",
            **summary,
        )
        for index, summary in enumerate(summaries, start=1)
    ]

    diagnostics: dict[str, float | int] = {
        "source_row_count": len(rows),
        "band_count": len(bands),
        "R_grid_gap_threshold": gap_threshold,
        "principal_count_from_largest_depth_gap": principal_count,
        "largest_adjacent_depth_gap": largest_depth_gap,
        "equivalent_error_ratio": 10.0**largest_depth_gap,
    }
    return bands, diagnostics


def identify_principal_pair(bands: Iterable[CandidateBand]) -> tuple[CandidateBand, CandidateBand]:
    principal = [band for band in bands if band.classification == "principal"]
    if len(principal) != 2:
        raise ValueError(
            "H8 is a two-well closure, but the empirical H7 depth-gap rule selected "
            f"{len(principal)} principal bands"
        )
    high_r, low_r = sorted(principal, key=lambda band: band.r_peak, reverse=True)
    return high_r, low_r


def h8_parameters(r1: CandidateBand, r2: CandidateBand) -> dict[str, float]:
    r_mid = 0.5 * (r1.r_peak + r2.r_peak)
    r_half_separation = 0.5 * (r1.r_peak - r2.r_peak)
    if r_half_separation <= 0:
        raise ValueError("R1 must be the higher-R member of the principal pair")

    # The cubic Hermite tilt h(x)=(3x-x^3)/2 changes the two well depths while
    # preserving stationary points at x=+1 and x=-1 exactly.
    epsilon = (r1.depth_peak - r2.depth_peak) / (r1.depth_peak + r2.depth_peak)
    return {
        "R1": r1.r_peak,
        "R2": r2.r_peak,
        "R_mid": r_mid,
        "R_half_separation": r_half_separation,
        "depth_R1": r1.depth_peak,
        "depth_R2": r2.depth_peak,
        "epsilon": epsilon,
    }


def effective_depth(x: np.ndarray, epsilon: float) -> np.ndarray:
    """H8 effective depth Phi(x), with stable maxima at x=+1 and x=-1."""
    return -(x * x - 1.0) ** 2 + 0.5 * epsilon * (3.0 * x - x**3)


def effective_depth_gradient(x: np.ndarray, epsilon: float) -> np.ndarray:
    """Analytic derivative dPhi/dx."""
    return (1.0 - x * x) * (4.0 * x + 1.5 * epsilon)


def run_h8_basin_sweep(
    parameters: dict[str, float],
    r_points: int = 150,
    velocity_points: int = 75,
    steps: int = 512,
    x_limit: float = 1.35,
    velocity_limit: float = 0.08,
    eta: float = 0.02,
    damping: float = 0.82,
    tolerance_x: float = 1.0e-3,
    tolerance_velocity: float = 1.0e-4,
) -> tuple[list[dict[str, float | int | str]], dict[str, float | int]]:
    """Measure basins of the explicitly assumed H8 damped R dynamics."""
    if not (0.0 <= damping < 1.0):
        raise ValueError("damping must be in [0, 1)")
    if eta <= 0.0:
        raise ValueError("eta must be positive")

    x_axis = np.linspace(-x_limit, x_limit, r_points)
    velocity_axis = np.linspace(-velocity_limit, velocity_limit, velocity_points)
    x0_grid, u0_grid = np.meshgrid(x_axis, velocity_axis)
    x = x0_grid.ravel().copy()
    velocity = u0_grid.ravel().copy()
    first_entry_step = np.full(x.shape, -1, dtype=int)
    epsilon = parameters["epsilon"]

    for step in range(1, steps + 1):
        velocity = damping * velocity + eta * effective_depth_gradient(x, epsilon)
        x = np.clip(x + velocity, -2.5, 2.5)
        distance = np.minimum(np.abs(x - 1.0), np.abs(x + 1.0))
        newly_entered = (
            (first_entry_step < 0)
            & (distance <= tolerance_x)
            & (np.abs(velocity) <= tolerance_velocity)
        )
        first_entry_step[newly_entered] = step

    near_r1 = (np.abs(x - 1.0) <= tolerance_x) & (np.abs(velocity) <= tolerance_velocity)
    near_r2 = (np.abs(x + 1.0) <= tolerance_x) & (np.abs(velocity) <= tolerance_velocity)
    labels = np.full(x.shape, "ambiguous", dtype=object)
    labels[near_r1] = "R1"
    labels[near_r2] = "R2"

    r0 = parameters["R_mid"] + parameters["R_half_separation"] * x0_grid.ravel()
    final_r = parameters["R_mid"] + parameters["R_half_separation"] * x
    records: list[dict[str, float | int | str]] = []
    for index in range(x.size):
        records.append(
            {
                "R0": float(r0[index]),
                "x0": float(x0_grid.ravel()[index]),
                "u0": float(u0_grid.ravel()[index]),
                "attractor": str(labels[index]),
                "final_R": float(final_r[index]),
                "final_x": float(x[index]),
                "final_u": float(velocity[index]),
                "first_entry_step": int(first_entry_step[index]),
            }
        )

    label_grid = labels.reshape(velocity_points, r_points)
    mirrored = np.flip(np.flip(label_grid, axis=0), axis=1)
    opposite = np.full(mirrored.shape, "ambiguous", dtype=object)
    opposite[mirrored == "R1"] = "R2"
    opposite[mirrored == "R2"] = "R1"
    paired_duality_fraction = float(np.mean(label_grid == opposite))

    counts = {name: int(np.sum(labels == name)) for name in ("R1", "R2", "ambiguous")}
    total = int(labels.size)
    converged_steps = first_entry_step[first_entry_step >= 0]
    statistics: dict[str, float | int] = {
        "total_points": total,
        "count_R1": counts["R1"],
        "count_R2": counts["R2"],
        "count_ambiguous": counts["ambiguous"],
        "fraction_R1": counts["R1"] / total,
        "fraction_R2": counts["R2"] / total,
        "fraction_ambiguous": counts["ambiguous"] / total,
        "paired_duality_fraction": paired_duality_fraction,
        "median_first_entry_step": float(np.median(converged_steps)),
        "max_first_entry_step": int(np.max(converged_steps)),
        "r_points": r_points,
        "velocity_points": velocity_points,
        "steps": steps,
        "eta": eta,
        "damping": damping,
        "x_limit": x_limit,
        "velocity_limit": velocity_limit,
        "tolerance_x": tolerance_x,
        "tolerance_velocity": tolerance_velocity,
    }
    return records, statistics


def plot_basin(records: list[dict[str, float | int | str]], parameters: dict[str, float], output: Path) -> None:
    r_values = sorted({float(row["R0"]) for row in records})
    u_values = sorted({float(row["u0"]) for row in records})
    r_index = {value: index for index, value in enumerate(r_values)}
    u_index = {value: index for index, value in enumerate(u_values)}
    label_grid = np.full((len(u_values), len(r_values)), 2, dtype=int)
    step_grid = np.full((len(u_values), len(r_values)), np.nan)

    for row in records:
        i = u_index[float(row["u0"])]
        j = r_index[float(row["R0"])]
        label_grid[i, j] = {"R2": 0, "R1": 1, "ambiguous": 2}[str(row["attractor"])]
        if int(row["first_entry_step"]) >= 0:
            step_grid[i, j] = int(row["first_entry_step"])

    extent = [r_values[0], r_values[-1], u_values[0], u_values[-1]]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
    cmap = ListedColormap(["#3478bf", "#d84b4b", "#b8b8b8"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    axes[0].imshow(label_grid, origin="lower", aspect="auto", extent=extent, cmap=cmap, norm=norm)
    axes[0].axvline(parameters["R1"], color="#9d1c1c", linewidth=1.2, linestyle="--")
    axes[0].axvline(parameters["R2"], color="#165a9b", linewidth=1.2, linestyle="--")
    axes[0].set_xlabel(r"initial exchange coefficient $R_0$")
    axes[0].set_ylabel(r"initial relaxation velocity $u_0$")
    axes[0].set_title("Basins generated by the H8 model\nred: $R_1$, blue: $R_2$, gray: unresolved")

    image = axes[1].imshow(step_grid, origin="lower", aspect="auto", extent=extent, cmap="viridis")
    axes[1].set_xlabel(r"initial exchange coefficient $R_0$")
    axes[1].set_ylabel(r"initial relaxation velocity $u_0$")
    axes[1].set_title("First entry into the H8 convergence tolerance")
    colorbar = fig.colorbar(image, ax=axes[1])
    colorbar.set_label("relaxation step")

    fig.suptitle("Phase 4 v2: basins generated by the explicitly assumed H8 relaxation", fontsize=13)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_selection_and_statistics(
    bands: list[CandidateBand],
    statistics: dict[str, float | int],
    symmetric_control: dict[str, float | int],
    output: Path,
) -> None:
    ranks = np.array([band.rank for band in bands])
    depths = np.array([band.depth_peak for band in bands])
    colors = ["#d84b4b" if band.classification == "principal" else "#8d939a" for band in bands]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    axes[0].bar(ranks, depths / depths[0], color=colors)
    axes[0].axhline(depths[2] / depths[0], color="#666666", linewidth=0.8, linestyle=":")
    axes[0].set_xticks(ranks)
    axes[0].set_xlabel("candidate rank")
    axes[0].set_ylabel("normalized gray depth")
    axes[0].set_title("H7 depth-gap rule selects ranks 1 and 2")
    axes[0].set_ylim(0.0, 1.08)

    names = ["R1", "R2", "unresolved"]
    counts = [
        int(statistics["count_R1"]),
        int(statistics["count_R2"]),
        int(statistics["count_ambiguous"]),
    ]
    bars = axes[1].bar(names, counts, color=["#d84b4b", "#3478bf", "#b8b8b8"])
    axes[1].set_ylabel("initial-condition count")
    axes[1].set_title("Initial-condition counts in the H8 model")
    axes[1].bar_label(bars, padding=3)
    axes[1].set_ylim(0, max(counts) * 1.15)
    axes[1].text(
        0.97,
        0.95,
        "symmetric control $\\varepsilon=0$:\n"
        f"$R_1$={int(symmetric_control['count_R1'])}, "
        f"$R_2$={int(symmetric_control['count_R2'])}",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.92},
    )

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_csv(records: list[dict[str, float | int | str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def portable_source_path(path: Path) -> str:
    """Avoid embedding a machine-specific absolute path in public results."""
    try:
        return str(path.resolve().relative_to(HERE.parents[1]))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", type=Path, default=DEFAULT_CANDIDATE_CSV)
    parser.add_argument("--output-dir", type=Path, default=HERE / "phase4_basin_v2_results")
    parser.add_argument("--figure-dir", type=Path, default=HERE / "figures")
    parser.add_argument("--r-points", type=int, default=150)
    parser.add_argument("--velocity-points", type=int, default=75)
    parser.add_argument("--steps", type=int, default=512)
    parser.add_argument("--eta", type=float, default=0.02)
    parser.add_argument("--damping", type=float, default=0.82)
    args = parser.parse_args()

    bands, selection_diagnostics = extract_candidate_bands(args.candidate_csv)
    r1, r2 = identify_principal_pair(bands)
    parameters = h8_parameters(r1, r2)
    records, statistics = run_h8_basin_sweep(
        parameters,
        r_points=args.r_points,
        velocity_points=args.velocity_points,
        steps=args.steps,
        eta=args.eta,
        damping=args.damping,
    )
    symmetric_parameters = dict(parameters)
    symmetric_parameters["epsilon"] = 0.0
    _, symmetric_control = run_h8_basin_sweep(
        symmetric_parameters,
        r_points=args.r_points,
        velocity_points=args.velocity_points,
        steps=args.steps,
        eta=args.eta,
        damping=args.damping,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(records, args.output_dir / "phase4_basin_sweep_v2.csv")
    with (args.output_dir / "phase4_candidate_bands_v2.json").open("w") as handle:
        json.dump(
            {
                "source": portable_source_path(args.candidate_csv),
                "selection_diagnostics": selection_diagnostics,
                "bands": [asdict(band) for band in bands],
            },
            handle,
            indent=2,
        )
    with (args.output_dir / "phase4_basin_statistics_v2.json").open("w") as handle:
        json.dump(
            {
                "logical_status": (
                    "Conditional output of the explicitly assumed H8 effective model; "
                    "not an independent derivation of physical R-space attractors."
                ),
                "H8_parameters": parameters,
                "statistics": statistics,
                "symmetric_control_epsilon_zero": symmetric_control,
            },
            handle,
            indent=2,
        )

    plot_basin(records, parameters, args.figure_dir / "fig4_basin_of_attraction.png")
    plot_selection_and_statistics(
        bands,
        statistics,
        symmetric_control,
        args.figure_dir / "fig4_basin_statistics.png",
    )

    print("Phase 4 v2 complete")
    print(
        "Principal pair: "
        f"R1={parameters['R1']:.15f} (depth={parameters['depth_R1']:.6f}), "
        f"R2={parameters['R2']:.15f} (depth={parameters['depth_R2']:.6f})"
    )
    print(
        "Depth gap after rank 2: "
        f"{selection_diagnostics['largest_adjacent_depth_gap']:.6f} "
        f"({selection_diagnostics['equivalent_error_ratio']:.1f}x in error)"
    )
    print(
        f"Basins: R1={statistics['count_R1']}, R2={statistics['count_R2']}, "
        f"unresolved={statistics['count_ambiguous']}"
    )
    print(f"Paired duality fraction: {statistics['paired_duality_fraction']:.6f}")
    print(
        "Symmetric H8 control (epsilon=0): "
        f"R1={symmetric_control['count_R1']}, R2={symmetric_control['count_R2']}, "
        f"unresolved={symmetric_control['count_ambiguous']}"
    )


if __name__ == "__main__":
    main()
