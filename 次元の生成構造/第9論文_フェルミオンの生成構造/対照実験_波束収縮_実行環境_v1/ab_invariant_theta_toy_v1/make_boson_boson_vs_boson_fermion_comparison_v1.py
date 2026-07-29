#!/usr/bin/env python3
"""広域ボゾンAと局在ボゾン／フェルミオンBの長時間反応を並列図化する。"""

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
DEFAULT_RESULT_DIR = HERE / "result_longrun_v1"
DEFAULT_FIGURE_DIR = DEFAULT_RESULT_DIR / "comparison_figures_v1"
MAX_COLLISION = 256
TAIL_START = 193

CASE_BOSON = "even_boson_control_B62"
CASE_FERMION = "odd_fermion_candidate_B63"

COLOR_A = "#4c78a8"
COLOR_B = "#e45756"
COLOR_MEAN = "#f2cf5b"


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location("ab_theta_toy_for_pair_comparison_v1", RUNNER_PATH)
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


def read_rows(path: Path) -> dict[str, list[dict[str, float]]]:
    grouped: dict[str, list[dict[str, float]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for raw in csv.DictReader(handle):
            case = str(raw["case"])
            row = {key: float(value) for key, value in raw.items() if key != "case"}
            grouped.setdefault(case, []).append(row)
    for rows in grouped.values():
        rows.sort(key=lambda item: item["collision"])
    return grouped


def values(rows: list[dict[str, float]], key: str) -> np.ndarray:
    return np.asarray([row[key] for row in rows], dtype=float)


def running_mean(data: np.ndarray) -> np.ndarray:
    result = np.full_like(data, np.nan, dtype=float)
    if data.size > 1:
        result[1:] = np.cumsum(data[1:]) / np.arange(1, data.size)
    return result


def save_figure(fig: Any, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def make_dynamics_figure(
    grouped: dict[str, list[dict[str, float]]],
    output_dir: Path,
) -> dict[str, float]:
    cases = (
        (
            CASE_BOSON,
            r"broad boson $A$ + localized boson $B$",
            r"$A=(1),\ B=(1,2,4,\ldots,62),\ \theta=0$",
        ),
        (
            CASE_FERMION,
            r"broad boson $A$ + localized fermion candidate $B$",
            r"$A=(1),\ B=(1,3,5,\ldots,63),\ \theta=0.761952071831$",
        ),
    )
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(16.0, 9.0),
        sharex=True,
        sharey="col",
        constrained_layout=True,
    )
    tail_stats: dict[str, float] = {}

    for row_index, (case, row_title, formula) in enumerate(cases):
        rows = grouped[case]
        collision = values(rows, "collision")
        transfer = values(rows, "origin_B_in_A")
        l_a = values(rows, "L_A")
        l_b = values(rows, "L_B")
        n_a = values(rows, "N_eff_A")
        n_b = values(rows, "N_eff_B")
        transfer_mean = running_mean(transfer)
        tail_mask = collision >= TAIL_START
        tail_min = float(np.min(transfer[tail_mask]))
        tail_max = float(np.max(transfer[tail_mask]))
        tail_mean = float(np.mean(transfer[tail_mask]))
        tail_stats[f"{case}_tail_min"] = tail_min
        tail_stats[f"{case}_tail_max"] = tail_max
        tail_stats[f"{case}_tail_mean"] = tail_mean
        tail_stats[f"{case}_full_running_mean"] = float(transfer_mean[-1])

        ax_transfer, ax_localization, ax_harmonic = axes[row_index]
        ax_transfer.plot(
            collision,
            transfer,
            color=COLOR_B,
            linewidth=1.45,
            alpha=0.86,
            label=r"instantaneous $B_0$ content in $A$",
        )
        ax_transfer.plot(
            collision,
            transfer_mean,
            color="#8c6d1f",
            linewidth=2.3,
            label="running mean",
        )
        ax_transfer.axvspan(TAIL_START, MAX_COLLISION, color=COLOR_MEAN, alpha=0.13)
        ax_transfer.set_ylim(-0.04, 1.04)
        if case == CASE_BOSON:
            ax_transfer.text(
                128,
                0.13,
                "stationary from collision 0\n(no exchange)",
                ha="center",
                color="#435f43",
                fontsize=10,
            )
        else:
            ax_transfer.text(
                224,
                0.10,
                "tail instantaneous range\n"
                f"{tail_min:.4f} to {tail_max:.4f}\n"
                "(no pointwise settling)",
                ha="center",
                color="#7d3838",
                fontsize=9,
            )
            ax_transfer.text(
                150,
                0.60,
                f"running mean at 256\n{transfer_mean[-1]:.4f}",
                ha="center",
                color="#6d5819",
                fontsize=9,
            )

        ax_localization.plot(collision, l_a, color=COLOR_A, linewidth=1.55, label=r"$L_A$")
        ax_localization.plot(collision, l_b, color=COLOR_B, linewidth=1.55, label=r"$L_B$")
        ax_localization.set_ylim(-0.0002, 0.0055)

        ax_harmonic.plot(collision, n_a, color=COLOR_A, linewidth=1.55, label=r"$N_{\rm eff,A}$")
        ax_harmonic.plot(collision, n_b, color=COLOR_B, linewidth=1.55, label=r"$N_{\rm eff,B}$")
        ax_harmonic.set_ylim(-0.8, 33.0)

        ax_transfer.set_ylabel(f"{row_title}\norigin projection weight")
        ax_transfer.text(
            0.02,
            0.95,
            formula,
            transform=ax_transfer.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none"},
        )
        ax_localization.set_ylabel("localization $L$")
        ax_harmonic.set_ylabel(r"effective harmonic $N_{\rm eff}$")

        for axis in axes[row_index]:
            axis.set_xlabel("collision count")
            axis.grid(alpha=0.22)
            axis.legend(frameon=False, fontsize=9, loc="upper right")

    axes[0, 0].set_title("Wave-origin transfer and time average")
    axes[0, 1].set_title("Localization response")
    axes[0, 2].set_title("Effective-harmonic response")
    fig.suptitle(
        "Broad boson A colliding with localized boson B versus localized fermion-candidate B\n"
        "same initial A and collision axis; yellow region is the final 64-collision window",
        fontsize=15,
    )
    save_figure(fig, output_dir, "fig3_boson_boson_vs_boson_fermion_dynamics_v1")
    return tail_stats


def chi_density(vector: np.ndarray, source_params: Any) -> np.ndarray:
    array = vector.reshape(source_params.chi_grid_n, source_params.eta_grid_n)
    density = np.sum(np.abs(array) ** 2, axis=1)
    total = float(np.sum(density))
    return density / total if total > 0.0 else density


def evolve_case_densities(
    case_name: str,
) -> tuple[np.ndarray, dict[int, tuple[np.ndarray, np.ndarray]], tuple[np.ndarray, np.ndarray]]:
    params = toy.base.Params(high_n=63, recursive_collision_count=MAX_COLLISION)
    source_params = toy.base.build_source_params(params)
    case = next(case for case in toy.build_cases(63) if case.name == case_name)
    base_case = toy.to_base_case(case)
    a = toy.base.make_case_state(source_params, base_case, "A", hair_enabled=True)
    b = toy.base.make_case_state(source_params, base_case, "B", hair_enabled=True)
    chi = np.linspace(-math.pi, math.pi, source_params.chi_grid_n, endpoint=False)
    snapshots: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    tail_a: list[np.ndarray] = []
    tail_b: list[np.ndarray] = []
    for collision in range(MAX_COLLISION + 1):
        density_a = chi_density(a, source_params)
        density_b = chi_density(b, source_params)
        if collision in (0, 1, 2):
            snapshots[collision] = (density_a.copy(), density_b.copy())
        if collision >= TAIL_START:
            tail_a.append(density_a)
            tail_b.append(density_b)
        if collision < MAX_COLLISION:
            readout = toy.theta_from_ab(a, b, source_params)
            a, b = toy.rotate_ab(a, b, readout.theta)
    tail_mean = (np.mean(tail_a, axis=0), np.mean(tail_b, axis=0))
    return chi, snapshots, tail_mean


def make_waveform_figure(output_dir: Path) -> None:
    cases = (
        (
            CASE_BOSON,
            r"broad boson $A$ + localized boson $B$",
            r"$\theta=0$: unchanged from collision 0",
        ),
        (
            CASE_FERMION,
            r"broad boson $A$ + localized fermion candidate $B$",
            r"$\theta=0.761952071831$: persistent exchange; no instantaneous steady state",
        ),
    )
    fig, axes = plt.subplots(
        2,
        4,
        figsize=(18.0, 8.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    column_titles = ("before collision", "after 1 collision", "after 2 collisions", "tail time average")

    for row_index, (case, row_title, row_note) in enumerate(cases):
        chi, snapshots, tail_mean = evolve_case_densities(case)
        series = (
            snapshots[0],
            snapshots[1],
            snapshots[2],
            tail_mean,
        )
        for column_index, ((density_a, density_b), column_title) in enumerate(
            zip(series, column_titles)
        ):
            axis = axes[row_index, column_index]
            axis.plot(chi, density_a, color=COLOR_A, linewidth=1.8, label="A channel")
            axis.plot(chi, density_b, color=COLOR_B, linewidth=1.8, label="B channel")
            axis.set_xlim(-math.pi, math.pi)
            axis.set_ylim(-0.004, 0.132)
            axis.set_xlabel(r"$\chi$")
            axis.grid(alpha=0.22)
            if row_index == 0:
                axis.set_title(column_title)
            if column_index == 0:
                axis.set_ylabel(f"{row_title}\nnormalized $\\chi$ density")
                axis.legend(frameon=False, fontsize=9)
            if column_index == 3:
                axis.axvspan(-math.pi, math.pi, color=COLOR_MEAN, alpha=0.07)
                axis.text(
                    0.5,
                    0.93,
                    "mean over collisions 193-256",
                    transform=axis.transAxes,
                    ha="center",
                    va="top",
                    fontsize=9,
                    color="#6d5819",
                )
                if row_index == 1:
                    axis.text(
                        0.5,
                        0.82,
                        "A/B tail-mean curves overlap",
                        transform=axis.transAxes,
                        ha="center",
                        va="top",
                        fontsize=9,
                        color="#6d5819",
                    )
        axes[row_index, 1].text(
            0.5,
            0.94,
            row_note,
            transform=axes[row_index, 1].transAxes,
            ha="center",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.80, "edgecolor": "none"},
        )

    fig.suptitle(
        "Waveform comparison from pre-collision state to long-run behavior\n"
        "rightmost column is a time average; it is not an instantaneous converged state",
        fontsize=15,
    )
    save_figure(fig, output_dir, "fig4_boson_boson_vs_boson_fermion_waveforms_v1")


def write_summary(path: Path, tail_stats: dict[str, float]) -> None:
    lines = [
        "# 広域ボゾンAと局在ボゾン／フェルミオンBの長時間比較 v1",
        "",
        "- 衝突回数: 0–256",
        f"- 尾部平均区間: {TAIL_START}–{MAX_COLLISION}",
        "",
        "## 瞬時状態",
        "",
        "- 広域ボゾンA＋局在ボゾンB: theta=0 のため衝突0から不変。",
        "- 広域ボゾンA＋局在フェルミオン候補B: 256衝突後も交換振動が持続し、点ごとの定常状態には収束しない。",
        "",
        "## 時間平均",
        "",
        (
            "- ボゾン・フェルミオン条件の尾部64衝突における "
            f"B0→A瞬時移乗範囲: {tail_stats[f'{CASE_FERMION}_tail_min']:.12g}–"
            f"{tail_stats[f'{CASE_FERMION}_tail_max']:.12g}"
        ),
        (
            "- 同じ尾部平均: "
            f"{tail_stats[f'{CASE_FERMION}_tail_mean']:.12g}"
        ),
        (
            "- 衝突1–256の走行平均の終値: "
            f"{tail_stats[f'{CASE_FERMION}_full_running_mean']:.12g}"
        ),
        "",
        "したがって、現行の損失なし実直交回転では、ボゾン・フェルミオン条件に",
        "瞬時の散逸的定常化は生じない。定常化するのは約1/2へ向かう時間平均である。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare boson-boson and boson-fermion long-run responses")
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped = read_rows(result_dir / "ab_invariant_theta_toy_rows_v1.csv")
    missing = sorted({CASE_BOSON, CASE_FERMION} - set(grouped))
    if missing:
        raise ValueError(f"missing comparison cases: {', '.join(missing)}")
    tail_stats = make_dynamics_figure(grouped, output_dir)
    make_waveform_figure(output_dir)
    write_summary(output_dir / "boson_boson_vs_boson_fermion_summary_v1.md", tail_stats)
    print(f"saved comparison figures: {output_dir}")


if __name__ == "__main__":
    main()
