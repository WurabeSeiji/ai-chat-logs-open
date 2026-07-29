#!/usr/bin/env python3
"""AB 回転不変量 theta トイモデル v1 の結果を図化する。"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
RUNNER_PATH = HERE / "run_ab_invariant_theta_toy_v1.py"
DEFAULT_RESULT_DIR = HERE / "result_v1"
DEFAULT_FIGURE_DIR = DEFAULT_RESULT_DIR / "figures_v1"


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location("ab_invariant_theta_toy_for_figures_v1", RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load toy runner: {RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy_module()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_COLORS = {
    "fundamental_control": "#4c78a8",
    "even_boson_control_B62": "#54a24b",
    "odd_fermion_candidate_B63": "#e45756",
}

CASE_LABELS = {
    "fundamental_control": "fundamental control",
    "even_boson_control_B62": "even-harmonic control B62",
    "odd_fermion_candidate_B63": "odd-harmonic candidate B63",
}


def read_rows(path: Path) -> dict[str, list[dict[str, float]]]:
    grouped: dict[str, list[dict[str, float]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for raw in csv.DictReader(handle):
            case = str(raw["case"])
            row: dict[str, float] = {}
            for key, value in raw.items():
                if key != "case":
                    row[key] = float(value)
            grouped.setdefault(case, []).append(row)
    for rows in grouped.values():
        rows.sort(key=lambda item: item["collision"])
    return grouped


def values(rows: list[dict[str, float]], key: str) -> np.ndarray:
    return np.asarray([row[key] for row in rows], dtype=float)


def save_figure(fig: Any, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def make_overview(grouped: dict[str, list[dict[str, float]]], output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2), constrained_layout=True)
    ax_theta, ax_rate, ax_transfer, ax_error = axes.ravel()

    for case, rows in grouped.items():
        collision = values(rows, "collision")
        color = CASE_COLORS.get(case)
        label = CASE_LABELS.get(case, case)
        ax_theta.plot(
            collision,
            values(rows, "theta_generated"),
            marker="o",
            markevery=max(1, len(rows) // 8),
            linewidth=2.0,
            color=color,
            label=label,
        )
        ax_rate.plot(
            collision,
            values(rows, "R_generated"),
            marker="o",
            markevery=max(1, len(rows) // 8),
            linewidth=2.0,
            color=color,
            label=label,
        )

    odd_rows = grouped["odd_fermion_candidate_B63"]
    collision = values(odd_rows, "collision")
    ax_transfer.plot(
        collision,
        values(odd_rows, "origin_B_in_A"),
        color="#e45756",
        linewidth=2.2,
        label=r"$B_0$ content in $A$",
    )
    ax_transfer.plot(
        collision,
        values(odd_rows, "origin_A_in_B"),
        color="#4c78a8",
        linewidth=2.2,
        linestyle="--",
        label=r"$A_0$ content in $B$",
    )
    ax_transfer.scatter(
        [2],
        [values(odd_rows, "origin_B_in_A")[2]],
        s=70,
        color="#f2cf5b",
        edgecolor="black",
        zorder=5,
        label="near-complete exchange at collision 2",
    )

    numerical_floor = 1.0e-18
    error_series = (
        ("theta drift", "theta_drift", "#e45756"),
        ("pair-norm drift", "pair_norm_drift", "#4c78a8"),
        ("zero-closure drift", "closure_drift_abs", "#54a24b"),
        ("combined-spectrum drift", "combined_spectrum_max_drift", "#b279a2"),
    )
    for label, key, color in error_series:
        error = np.maximum(values(odd_rows, key), numerical_floor)
        ax_error.semilogy(collision, error, linewidth=1.9, color=color, label=label)

    ax_theta.set_title(r"Internally generated rotation angle $\theta$")
    ax_theta.set_ylabel(r"$\theta$ [rad]")
    ax_theta.set_ylim(bottom=-0.03)
    ax_theta.text(
        16,
        0.025,
        "both controls overlap at 0",
        ha="center",
        color="#3b6f3b",
        fontsize=9,
    )
    ax_theta.legend(frameon=False, fontsize=9)

    ax_rate.set_title(r"Internally generated reflection rate $R=\sin^2\theta$")
    ax_rate.set_ylabel(r"$R$")
    ax_rate.set_ylim(-0.03, 0.53)
    ax_rate.text(
        16,
        0.018,
        "both controls overlap at 0",
        ha="center",
        color="#3b6f3b",
        fontsize=9,
    )
    ax_rate.legend(frameon=False, fontsize=9)

    ax_transfer.set_title("Repeated AB exchange: odd-harmonic candidate B63")
    ax_transfer.set_ylabel("origin projection weight")
    ax_transfer.set_ylim(-0.04, 1.04)
    ax_transfer.legend(frameon=False, fontsize=9)

    ax_error.set_title("Invariant drift: odd-harmonic candidate B63")
    ax_error.set_ylabel("absolute numerical drift")
    ax_error.set_ylim(5.0e-19, 2.0e-14)
    ax_error.legend(frameon=False, fontsize=9)

    for axis in axes.ravel():
        axis.set_xlabel("collision count")
        axis.grid(alpha=0.25)

    fig.suptitle(
        "AB-invariant scattering angle toy model v1\n"
        r"no external $R$ or $\theta$; $A=(1)$, fermion candidate $B=(1,3,\ldots,63)$",
        fontsize=15,
    )
    save_figure(fig, output_dir, "fig1_theta_R_transfer_invariants_v1")


def chi_density(vector: np.ndarray, source_params: Any) -> np.ndarray:
    array = vector.reshape(source_params.chi_grid_n, source_params.eta_grid_n)
    density = np.sum(np.abs(array) ** 2, axis=1)
    total = float(np.sum(density))
    return density / total if total > 0.0 else density


def capture_odd_snapshots() -> tuple[Any, Any, dict[int, tuple[np.ndarray, np.ndarray]]]:
    params = toy.base.Params(high_n=63, recursive_collision_count=2)
    source_params = toy.base.build_source_params(params)
    metric_context = toy.base.MetricContext(source_params)
    case = next(case for case in toy.build_cases(63) if case.name == "odd_fermion_candidate_B63")
    base_case = toy.to_base_case(case)
    a = toy.base.make_case_state(source_params, base_case, "A", hair_enabled=True)
    b = toy.base.make_case_state(source_params, base_case, "B", hair_enabled=True)
    snapshots: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for collision in range(3):
        snapshots[collision] = (a.copy(), b.copy())
        if collision < 2:
            readout = toy.theta_from_ab(a, b, source_params)
            a, b = toy.rotate_ab(a, b, readout.theta)
    return source_params, metric_context, snapshots


def make_waveform_exchange(output_dir: Path) -> None:
    source_params, metric_context, snapshots = capture_odd_snapshots()
    chi = np.linspace(-math.pi, math.pi, source_params.chi_grid_n, endpoint=False)
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.2, 7.8),
        sharey="row",
        constrained_layout=True,
    )

    for column, collision in enumerate((0, 1, 2)):
        a, b = snapshots[collision]
        density_a = chi_density(a, source_params)
        density_b = chi_density(b, source_params)
        axis_density = axes[0, column]
        axis_density.plot(chi, density_a, color="#4c78a8", linewidth=1.8, label="A channel")
        axis_density.plot(chi, density_b, color="#e45756", linewidth=1.8, label="B channel")
        axis_density.set_title(f"collision {collision}")
        axis_density.set_xlabel(r"$\chi$")
        axis_density.set_ylabel(r"normalized $\chi$ density")
        axis_density.set_xlim(-math.pi, math.pi)
        axis_density.grid(alpha=0.22)
        if column == 0:
            axis_density.legend(frameon=False)

        distribution_a = metric_context.harmonic_distribution(a)
        distribution_b = metric_context.harmonic_distribution(b)
        harmonic_max = max(max(distribution_a), max(distribution_b))
        harmonic = np.arange(harmonic_max + 1)
        power_a = np.asarray([distribution_a.get(int(n), 0.0) for n in harmonic])
        power_b = np.asarray([distribution_b.get(int(n), 0.0) for n in harmonic])
        axis_harmonic = axes[1, column]
        axis_harmonic.plot(harmonic, power_a, color="#4c78a8", linewidth=1.6, label="A channel")
        axis_harmonic.plot(harmonic, power_b, color="#e45756", linewidth=1.6, label="B channel")
        axis_harmonic.fill_between(harmonic, power_a, color="#4c78a8", alpha=0.15)
        axis_harmonic.fill_between(harmonic, power_b, color="#e45756", alpha=0.12)
        axis_harmonic.set_xlabel(r"raw harmonic bin $|k|$")
        axis_harmonic.set_ylabel("normalized harmonic power")
        axis_harmonic.set_xlim(0, 66)
        axis_harmonic.grid(alpha=0.22)

    fig.suptitle(
        "Waveform and harmonic exchange generated by the AB-invariant angle\n"
        r"$\theta=0.761952071831$, $R=61/128$; near-complete channel exchange after two collisions",
        fontsize=15,
    )
    save_figure(fig, output_dir, "fig2_odd_B63_waveform_harmonic_exchange_v1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot AB-invariant theta toy-model results")
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows_path = result_dir / "ab_invariant_theta_toy_rows_v1.csv"
    grouped = read_rows(rows_path)
    required = set(CASE_LABELS)
    missing = sorted(required - set(grouped))
    if missing:
        raise ValueError(f"missing cases in result CSV: {', '.join(missing)}")
    make_overview(grouped, output_dir)
    make_waveform_exchange(output_dir)
    print(f"saved figures: {output_dir}")


if __name__ == "__main__":
    main()
