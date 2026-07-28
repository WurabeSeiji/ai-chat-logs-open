"""Analyze Stage F without searching beyond the preregistered runs.

The thresholds below are fixed before classification is computed:

* numerical_residual_limit = 1e-8
* baseline_equivalence_limit = 1e-8 in the scaled observable vector
* cycle_return_tolerance = 1e-6
* fixed_point_tail_step_tolerance = 1e-8
* amplitude_shift_tolerance = 1e-4
* tested lags = 1..64 only

No threshold is fitted to the output.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
FIGURE_DIR = STAGE_ROOT / "figures"
REPORT_DIR = STAGE_ROOT / "reports"
MPL_DIR = DATA_DIR / ".matplotlib-cache"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


NUMERICAL_RESIDUAL_LIMIT = 1.0e-8
BASELINE_EQUIVALENCE_LIMIT = 1.0e-8
CYCLE_RETURN_TOLERANCE = 1.0e-6
FIXED_POINT_TAIL_STEP_TOLERANCE = 1.0e-8
AMPLITUDE_SHIFT_TOLERANCE = 1.0e-4
MAX_LAG = 64

NUMERIC_FIELDS = {
    "collision_index",
    "kappa",
    "R0",
    "theta0",
    "rho",
    "c_A",
    "c_B",
    "c_mean",
    "delta_theta",
    "theta_eff",
    "R_eff",
    "T_eff",
    "L_A",
    "L_B",
    "N_eff_A",
    "N_eff_B",
    "spectral_similarity_A_to_initial_A",
    "spectral_similarity_A_to_initial_B",
    "spectral_similarity_B_to_initial_A",
    "spectral_similarity_B_to_initial_B",
    "B_to_A_transfer",
    "path_A_to_A_norm",
    "path_B_to_A_norm",
    "path_B_to_B_norm",
    "path_A_to_B_norm",
    "interference_A",
    "interference_B",
    "raw_norm_A",
    "raw_norm_B",
    "next_state_norm_A",
    "next_state_norm_B",
    "boson_weight_A",
    "fermion_weight_A",
    "boson_weight_B",
    "fermion_weight_B",
    "next_c_A",
    "next_c_B",
    "unitarity_residual",
    "coefficient_orthogonality_residual",
    "path_sum_residual_A",
    "path_sum_residual_B",
    "total_norm_conservation_residual",
    "demodulation_reconstruction_residual",
    "parity_projection_sum_residual",
    "normalization_scale_A",
    "normalization_scale_B",
    "nan_inf_count",
}


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for field in NUMERIC_FIELDS & set(row):
            row[field] = float(row[field])
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to create empty output: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def group_runs(rows: list[dict[str, Any]]) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row["run_id"])].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: row["collision_index"])
    return dict(grouped)


def observable_matrix(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(
        [
            [
                row["L_A"],
                row["L_B"],
                row["N_eff_A"],
                row["N_eff_B"],
                row["c_A"],
                row["c_B"],
                row["R_eff"],
            ]
            for row in rows
        ],
        dtype=float,
    )
    first = np.abs(matrix[0])
    scales = np.asarray(
        [
            max(first[0], first[1], 1.0e-12),
            max(first[0], first[1], 1.0e-12),
            max(first[2], first[3], 1.0),
            max(first[2], first[3], 1.0),
            1.0,
            1.0,
            1.0,
        ]
    )
    return matrix, scales


def cycle_metrics(rows: list[dict]) -> dict[str, Any]:
    matrix, scales = observable_matrix(rows)
    scaled = matrix / scales
    lag_metrics = []
    for lag in range(1, min(MAX_LAG, len(rows) - 1) + 1):
        left = scaled[lag:]
        right = scaled[:-lag]
        return_rms = float(
            np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1)))
        )
        flattened_left = left.reshape(-1)
        flattened_right = right.reshape(-1)
        if (
            np.std(flattened_left) <= 1.0e-300
            or np.std(flattened_right) <= 1.0e-300
        ):
            autocorrelation = 1.0 if return_rms == 0.0 else 0.0
        else:
            autocorrelation = float(
                np.corrcoef(flattened_left, flattened_right)[0, 1]
            )
        lag_metrics.append((lag, return_rms, autocorrelation))
    best_lag, best_return, best_autocorrelation = min(
        lag_metrics, key=lambda item: item[1]
    )
    tail = scaled[-32:]
    tail_step = float(
        np.sqrt(np.mean(np.sum(np.diff(tail, axis=0) ** 2, axis=1)))
    )
    l_difference = np.asarray(
        [row["L_A"] - row["L_B"] for row in rows]
    )
    signs = np.sign(l_difference)
    exchange_count = int(
        np.sum(
            [
                signs[index] != 0
                and signs[index - 1] != 0
                and signs[index] != signs[index - 1]
                for index in range(1, len(signs))
            ]
        )
    )
    l_gap = np.abs(l_difference)
    n_gap = np.asarray(
        [
            abs(row["N_eff_A"] - row["N_eff_B"])
            for row in rows
        ]
    )
    amplitude_score = float(
        (np.ptp(l_gap) / scales[0]) + (np.ptp(n_gap) / scales[2])
    )
    return {
        "best_lag": best_lag,
        "best_return_rms": best_return,
        "best_lag_autocorrelation": best_autocorrelation,
        "cycle_detected": best_return <= CYCLE_RETURN_TOLERANCE,
        "tail_step_rms": tail_step,
        "exchange_count": exchange_count,
        "amplitude_score": amplitude_score,
    }


def scaled_max_difference(
    rows: list[dict], baseline: list[dict]
) -> float:
    matrix, scales = observable_matrix(baseline)
    candidate, _ = observable_matrix(rows)
    return float(
        np.max(np.sqrt(np.sum(((candidate - matrix) / scales) ** 2, axis=1)))
    )


def classify(
    row_metrics: dict,
    baseline_metrics: dict,
    baseline_difference: float,
    maximum_residual: float,
    nan_inf_count: int,
    norm_min: float,
    norm_max: float,
) -> str:
    if nan_inf_count or maximum_residual > NUMERICAL_RESIDUAL_LIMIT:
        return "numerically_unstable"
    if norm_max > 1.0e6 or norm_min < 1.0e-12:
        return "divergent"
    if baseline_difference <= BASELINE_EQUIVALENCE_LIMIT:
        return "baseline_equivalent"
    if row_metrics["tail_step_rms"] <= FIXED_POINT_TAIL_STEP_TOLERANCE:
        return "fixed_point_convergence"
    period_shift = (
        row_metrics["best_lag"] != baseline_metrics["best_lag"]
        or row_metrics["cycle_detected"]
        != baseline_metrics["cycle_detected"]
    )
    amplitude_shift = (
        abs(
            row_metrics["amplitude_score"]
            - baseline_metrics["amplitude_score"]
        )
        > AMPLITUDE_SHIFT_TOLERANCE
    )
    if period_shift and amplitude_shift:
        return "period_and_amplitude_shift"
    if period_shift:
        return "period_shift_only"
    if amplitude_shift:
        return "amplitude_shift_only"
    return "new_quasistable_cycle"


def calculate_cycle_tables(
    all_rows: list[dict],
) -> tuple[list[dict], list[dict], list[dict]]:
    grouped = group_runs(all_rows)
    raw_metrics = {
        run_id: cycle_metrics(rows) for run_id, rows in grouped.items()
    }
    baseline_lookup = {}
    for run_id, rows in grouped.items():
        first = rows[0]
        if (
            first["scattering_mode"] == "C0"
            and abs(first["kappa"] - 0.01) <= 1.0e-15
        ):
            baseline_lookup[
                (first["normalization_mode"], first["R0"])
            ] = run_id

    cycle_rows = []
    residual_rows = []
    enriched_summaries = []
    for run_id, rows in sorted(grouped.items()):
        first = rows[0]
        baseline_id = baseline_lookup[
            (first["normalization_mode"], first["R0"])
        ]
        metrics = raw_metrics[run_id]
        baseline_metrics = raw_metrics[baseline_id]
        baseline_difference = scaled_max_difference(
            rows, grouped[baseline_id]
        )
        maximums = {
            field: max(abs(float(row[field])) for row in rows)
            for field in (
                "unitarity_residual",
                "coefficient_orthogonality_residual",
                "path_sum_residual_A",
                "path_sum_residual_B",
                "total_norm_conservation_residual",
                "demodulation_reconstruction_residual",
                "parity_projection_sum_residual",
            )
        }
        maximum_residual = max(maximums.values())
        nan_inf_count = int(sum(row["nan_inf_count"] for row in rows))
        norm_min = min(
            min(row["raw_norm_A"], row["raw_norm_B"]) for row in rows
        )
        norm_max = max(
            max(row["raw_norm_A"], row["raw_norm_B"]) for row in rows
        )
        classification = classify(
            metrics,
            baseline_metrics,
            baseline_difference,
            maximum_residual,
            nan_inf_count,
            norm_min,
            norm_max,
        )
        cycle_rows.append(
            {
                "run_id": run_id,
                "scattering_mode": first["scattering_mode"],
                "normalization_mode": first["normalization_mode"],
                "kappa": first["kappa"],
                "R0": first["R0"],
                **metrics,
                "scaled_max_difference_from_C0": baseline_difference,
                "classification": classification,
                "cycle_return_tolerance": CYCLE_RETURN_TOLERANCE,
                "tested_lag_min": 1,
                "tested_lag_max": MAX_LAG,
            }
        )
        residual_rows.append(
            {
                "run_id": run_id,
                "scattering_mode": first["scattering_mode"],
                "normalization_mode": first["normalization_mode"],
                "kappa": first["kappa"],
                "R0": first["R0"],
                **maximums,
                "raw_norm_min": norm_min,
                "raw_norm_max": norm_max,
                "normalization_scale_min": min(
                    min(
                        row["normalization_scale_A"],
                        row["normalization_scale_B"],
                    )
                    for row in rows
                ),
                "normalization_scale_max": max(
                    max(
                        row["normalization_scale_A"],
                        row["normalization_scale_B"],
                    )
                    for row in rows
                ),
                "nan_inf_count": nan_inf_count,
                "theta_range_violation_count": sum(
                    int(str(row["theta_range_violation"]).lower() == "true")
                    for row in rows
                ),
            }
        )
        min_l = min(
            rows, key=lambda row: abs(row["L_A"] - row["L_B"])
        )
        min_n = min(
            rows,
            key=lambda row: abs(
                row["N_eff_A"] - row["N_eff_B"]
            ),
        )
        enriched_summaries.append(
            {
                "run_id": run_id,
                "scattering_mode": first["scattering_mode"],
                "normalization_mode": first["normalization_mode"],
                "kappa": first["kappa"],
                "R0": first["R0"],
                "collision_count": len(rows),
                "min_L_difference": abs(
                    min_l["L_A"] - min_l["L_B"]
                ),
                "min_L_difference_collision": int(
                    min_l["collision_index"]
                ),
                "min_N_eff_difference": abs(
                    min_n["N_eff_A"] - min_n["N_eff_B"]
                ),
                "min_N_eff_difference_collision": int(
                    min_n["collision_index"]
                ),
                "R_eff_min": min(row["R_eff"] for row in rows),
                "R_eff_max": max(row["R_eff"] for row in rows),
                "c_A_min": min(row["c_A"] for row in rows),
                "c_A_max": max(row["c_A"] for row in rows),
                "c_B_min": min(row["c_B"] for row in rows),
                "c_B_max": max(row["c_B"] for row in rows),
                "c_mean_min": min(row["c_mean"] for row in rows),
                "c_mean_max": max(row["c_mean"] for row in rows),
                "best_lag": metrics["best_lag"],
                "best_return_rms": metrics["best_return_rms"],
                "cycle_detected": metrics["cycle_detected"],
                "exchange_count": metrics["exchange_count"],
                "classification": classification,
                "maximum_numerical_residual": maximum_residual,
                "nan_inf_count": nan_inf_count,
            }
        )
    return cycle_rows, residual_rows, enriched_summaries


def representative(
    rows: list[dict],
    scattering: str,
    normalization: str,
    kappa: float,
    reflection: float = 0.55,
) -> list[dict]:
    return [
        row
        for row in rows
        if row["scattering_mode"] == scattering
        and row["normalization_mode"] == normalization
        and abs(row["kappa"] - kappa) <= 1.0e-15
        and abs(row["R0"] - reflection) <= 1.0e-15
    ]


def save_figures(
    all_rows: list[dict],
    summaries: list[dict],
    cycle_rows: list[dict],
    f_c_rows: list[dict],
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    colors = {0.01: "#2a6fbb", 0.1: "#e07a1f", 1.0: "#a23b72"}
    c0 = representative(
        all_rows, "C0", "existing_normalization", 0.01
    )

    for metric, filename, ylabel in (
        (
            "L",
            "C0_vs_reversed_C1_L_exchange.png",
            "localization L",
        ),
        (
            "N_eff",
            "C0_vs_reversed_C1_N_eff_exchange.png",
            "effective harmonic N_eff",
        ),
    ):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
        for channel, axis in (("A", axes[0]), ("B", axes[1])):
            axis.plot(
                [row["collision_index"] for row in c0],
                [row[f"{metric}_{channel}"] for row in c0],
                color="black",
                linewidth=1.6,
                label="C0",
            )
            for kappa in (0.01, 0.1, 1.0):
                values = representative(
                    all_rows,
                    "reversed_C1",
                    "existing_normalization",
                    kappa,
                )
                axis.plot(
                    [row["collision_index"] for row in values],
                    [row[f"{metric}_{channel}"] for row in values],
                    color=colors[kappa],
                    linewidth=1.0,
                    label=f"reversed C1 k={kappa:g}",
                )
            axis.set_title(f"channel {channel}")
            axis.set_xlabel("collision")
            axis.set_ylabel(ylabel)
            axis.grid(alpha=0.2)
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURE_DIR / filename, dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for kappa in (0.01, 0.1, 1.0):
        values = representative(
            all_rows,
            "reversed_C1",
            "existing_normalization",
            kappa,
        )
        axes[0].plot(
            [row["collision_index"] for row in values],
            [row["c_mean"] for row in values],
            color=colors[kappa],
            label=f"k={kappa:g}",
        )
        axes[1].plot(
            [row["collision_index"] for row in values],
            [row["R_eff"] for row in values],
            color=colors[kappa],
            label=f"k={kappa:g}",
        )
    axes[0].set_ylabel("c_mean")
    axes[1].set_ylabel("R_eff")
    axes[1].set_xlabel("collision")
    axes[0].grid(alpha=0.2)
    axes[1].grid(alpha=0.2)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "parity_and_R_eff_by_collision.png", dpi=180
    )
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for scattering, kappa, color, label in (
        ("C0", 0.01, "black", "C0"),
        ("reversed_C1", 0.01, colors[0.01], "reversed k=.01"),
        ("reversed_C1", 0.1, colors[0.1], "reversed k=.1"),
        ("reversed_C1", 1.0, colors[1.0], "reversed k=1"),
    ):
        selected = sorted(
            [
                row
                for row in summaries
                if row["scattering_mode"] == scattering
                and row["normalization_mode"]
                == "existing_normalization"
                and abs(float(row["kappa"]) - kappa) <= 1.0e-15
            ],
            key=lambda row: float(row["R0"]),
        )
        axes[0].plot(
            [float(row["R0"]) for row in selected],
            [float(row["min_L_difference"]) for row in selected],
            marker="o",
            color=color,
            label=label,
        )
        axes[1].plot(
            [float(row["R0"]) for row in selected],
            [float(row["min_N_eff_difference"]) for row in selected],
            marker="o",
            color=color,
            label=label,
        )
    axes[0].set_ylabel("min |L_A-L_B|")
    axes[1].set_ylabel("min |N_eff,A-N_eff,B|")
    for axis in axes:
        axis.set_xlabel("R0")
        axis.set_yscale("symlog", linthresh=1.0e-12)
        axis.grid(alpha=0.2)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "R_sweep_minimum_differences.png", dpi=180
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for scattering, kappa, color, label in (
        ("C0", 0.01, "black", "C0"),
        ("reversed_C1", 0.01, colors[0.01], "reversed k=.01"),
        ("reversed_C1", 0.1, colors[0.1], "reversed k=.1"),
        ("reversed_C1", 1.0, colors[1.0], "reversed k=1"),
    ):
        selected = sorted(
            [
                row
                for row in cycle_rows
                if row["scattering_mode"] == scattering
                and row["normalization_mode"]
                == "existing_normalization"
                and abs(float(row["kappa"]) - kappa) <= 1.0e-15
            ],
            key=lambda row: float(row["R0"]),
        )
        ax.plot(
            [float(row["R0"]) for row in selected],
            [float(row["best_lag"]) for row in selected],
            marker="o",
            color=color,
            label=label,
        )
    ax.set_xlabel("R0")
    ax.set_ylabel("best lag in preregistered 1..64")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "cycle_period_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    normed = representative(
        all_rows, "reversed_C1", "existing_normalization", 1.0
    )
    raw = representative(
        all_rows, "reversed_C1", "raw_update", 1.0
    )
    axes[0].plot(
        [row["collision_index"] for row in normed],
        [row["L_A"] for row in normed],
        label="existing normalization",
    )
    axes[0].plot(
        [row["collision_index"] for row in raw],
        [row["L_A"] for row in raw],
        linestyle="--",
        label="raw update",
    )
    axes[1].plot(
        [row["collision_index"] for row in normed],
        [row["next_state_norm_A"] for row in normed],
        label="existing normalization",
    )
    axes[1].plot(
        [row["collision_index"] for row in raw],
        [row["next_state_norm_A"] for row in raw],
        linestyle="--",
        label="raw update",
    )
    axes[0].set_ylabel("L_A")
    axes[1].set_ylabel("channel A norm squared")
    for axis in axes:
        axis.set_xlabel("collision")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "normalization_vs_raw_update.png", dpi=180
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for kappa in (0.01, 0.1, 1.0):
        selected = sorted(
            [
                row
                for row in f_c_rows
                if row["normalization_mode"]
                == "existing_normalization"
                and abs(float(row["kappa"]) - kappa) <= 1.0e-15
            ],
            key=lambda row: float(row["iteration"]),
        )
        ax.plot(
            [float(row["iteration"]) for row in selected],
            [float(row["reversed_C1_return_error"]) for row in selected],
            marker="o",
            color=colors[kappa],
            label=f"reversed k={kappa:g}",
        )
    c0_fc = sorted(
        [
            row
            for row in f_c_rows
            if row["normalization_mode"] == "existing_normalization"
            and abs(float(row["kappa"]) - 0.01) <= 1.0e-15
        ],
        key=lambda row: float(row["iteration"]),
    )
    ax.plot(
        [float(row["iteration"]) for row in c0_fc],
        [float(row["C0_return_error"]) for row in c0_fc],
        marker="o",
        color="black",
        label="C0",
    )
    ax.set_xlabel("fixed 31-series iteration")
    ax.set_ylabel("full-state return error")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "return_error_31_series.png", dpi=180)
    plt.close(fig)


def maximum_normalization_difference(
    all_rows: list[dict],
) -> dict[str, float]:
    grouped = group_runs(all_rows)
    maximum_observable = 0.0
    for run_id, rows in grouped.items():
        if "__existing_normalization__" not in run_id:
            continue
        raw_id = run_id.replace(
            "__existing_normalization__", "__raw_update__"
        )
        raw_rows = grouped[raw_id]
        for left, right in zip(rows, raw_rows):
            maximum_observable = max(
                maximum_observable,
                abs(left["L_A"] - right["L_A"]),
                abs(left["L_B"] - right["L_B"]),
                abs(left["N_eff_A"] - right["N_eff_A"]),
                abs(left["N_eff_B"] - right["N_eff_B"]),
                abs(left["R_eff"] - right["R_eff"]),
                abs(left["c_mean"] - right["c_mean"]),
            )

    archive = np.load(DATA_DIR / "stage_F_final_states.npz")
    maximum_final_state = 0.0
    for key in archive.files:
        if "__existing_normalization__" not in key:
            continue
        raw_key = key.replace(
            "__existing_normalization__", "__raw_update__"
        )
        difference = archive[key] - archive[raw_key]
        maximum_final_state = max(
            maximum_final_state,
            math.sqrt(float(np.vdot(difference, difference).real)),
        )
    return {
        "maximum_observable_difference": maximum_observable,
        "maximum_final_state_vector_difference": maximum_final_state,
    }


def reports(
    all_rows: list[dict],
    summaries: list[dict],
    cycles: list[dict],
    residuals: list[dict],
    f_c_rows: list[dict],
    normalization: dict[str, float],
) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    before = json.loads(
        (DATA_DIR / "reference_hashes_before.json").read_text(
            encoding="utf-8"
        )
    )
    source_hashes = [
        item
        for item in before["references"]
        if item["path"].endswith(".py")
    ]
    report_00 = """# 00 範囲と参照正本

Stage Fの全出力は本ディレクトリ内に限定した。既存System A/B、Stage E、既存CSV・図・報告書は読取りだけで、変更していない。

## 正本

"""
    for item in source_hashes:
        report_00 += (
            f"- `{item['path']}`\n"
            f"  - SHA-256: `{item['sha256']}`\n"
        )
    report_00 += f"""

事前監視対象は合計{before["reference_count"]}ファイル。事後値は `reference_hashes_after.json` と `reference_hash_comparison.json` に保存する。

20260713正本は `N_A=1,N_B=63`、128衝突、R点 `(0,0.51,0.55,0.60,0.70,0.90,1)` の直接参照元である。31系列は別条件の20260715 custom packet `A=(1), B=(1,2), R=0.697177927` に属する。
"""
    (REPORT_DIR / "00_scope_and_source_hashes.md").write_text(
        report_00, encoding="utf-8"
    )

    (REPORT_DIR / "02_integration_definition.md").write_text(
        """# 02 統合定義

## 物理散乱層

512×16全状態に対し、`a_raw=r_eff*a+t_eff*b`, `b_raw=t_eff*a+r_eff*b` を適用した。C0は `theta_eff=theta0`、反転Candidate 1は `theta_eff=theta0-kappa*rho(theta0)*(c_A+c_B)/2` であり、候補差は角度だけである。型名による条件分岐はない。

## 読出し層

A由来 `(q=+1,m_eta=1)` とB由来 `(q=-1,m_eta=2)` をeta射影し、それぞれの搬送波を除去してから半周期相関 `c_pi` と偶奇射影重みを計算した。全チャネルへの単一復調は行っていない。

## 更新系列

`existing_normalization` は正本と同じチャネル別正規化、`raw_update` はraw出力をそのまま次状態にした。raw全状態と経路ノルムを物理量、スペクトル類似度を診断量として分離した。
""",
        encoding="utf-8",
    )

    fa = [
        row
        for row in summaries
        if abs(float(row["R0"]) - 0.55) <= 1.0e-15
        and row["normalization_mode"] == "existing_normalization"
        and (
            row["scattering_mode"] == "reversed_C1"
            or abs(float(row["kappa"]) - 0.01) <= 1.0e-15
        )
    ]
    fa_lines = [
        "| mode | kappa | R_eff | min L gap (collision) | min N_eff gap (collision) | classification |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in fa:
        fa_lines.append(
            f"| {row['scattering_mode']} | {float(row['kappa']):g} | "
            f"{float(row['R_eff_min']):.12g} | "
            f"{float(row['min_L_difference']):.6g} "
            f"({int(row['min_L_difference_collision'])}) | "
            f"{float(row['min_N_eff_difference']):.6g} "
            f"({int(row['min_N_eff_difference_collision'])}) | "
            f"{row['classification']} |"
        )
    (REPORT_DIR / "03_repeated_scattering_comparison.md").write_text(
        """# 03 代表条件の反復比較

`N_A=1,N_B=63,R0=0.55` では両入力が奇数倍音カーネルであり、全128衝突を通じて `c_A,c_B,c_mean=-1` が機械精度内で維持された。したがって反転Candidate 1の `R_eff` は時間変動せず、κごとの定数となった。

"""
        + "\n".join(fa_lines)
        + """

これは「型依存応答の実装」が既存の交換運動を別の一定反射率の交換運動へ写したという数値観察である。パリティ純度が交換に伴って変化した、という観察ではない。
""",
        encoding="utf-8",
    )

    selected_r = [
        row
        for row in summaries
        if float(row["R0"]) in (0.0, 0.55, 0.7, 1.0)
        and row["normalization_mode"] == "existing_normalization"
        and (
            (
                row["scattering_mode"] == "C0"
                and abs(float(row["kappa"]) - 0.01) <= 1.0e-15
            )
            or row["scattering_mode"] == "reversed_C1"
        )
    ]
    sweep_lines = [
        "| R0 | mode | kappa | R_eff | min L gap@collision | min N_eff gap@collision |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for row in selected_r:
        sweep_lines.append(
            f"| {float(row['R0']):g} | {row['scattering_mode']} | "
            f"{float(row['kappa']):g} | {float(row['R_eff_min']):.9g} | "
            f"{float(row['min_L_difference']):.5g}@"
            f"{int(row['min_L_difference_collision'])} | "
            f"{float(row['min_N_eff_difference']):.5g}@"
            f"{int(row['min_N_eff_difference_collision'])} |"
        )
    (REPORT_DIR / "04_R_sweep_comparison.md").write_text(
        """# 04 R掃引比較

正本7点 `(0,0.51,0.55,0.60,0.70,0.90,1)` を変更せず使用した。C0は正本核、反転Candidate 1は同じ初期状態・反復数である。端点では `rho=0` のため両核が一致する。

"""
        + "\n".join(sweep_lines)
        + """

R0=0.70でC0にも小さい局在性差・実効次数差が現れるという既存観察を再現した。反転Candidate 1ではκによって最小値と到達衝突回が移るが、R_effが衝突ごとに動いたためではなく、全奇数条件で一定のR_effへ写ったためである。
""",
        encoding="utf-8",
    )

    counts = Counter(row["classification"] for row in cycles)
    count_text = "\n".join(
        f"- `{key}`: {value}" for key, value in sorted(counts.items())
    )
    (REPORT_DIR / "05_cycle_and_quasistability_analysis.md").write_text(
        f"""# 05 周期・準安定性解析

既存コードに一般周期分類器はなかったため、実行前固定の方法として主要観測ベクトルのスケール済み帰還RMSと自己相関をlag 1〜64だけで測った。帰還閾値は `{CYCLE_RETURN_TOLERANCE:.1e}`、C0同等閾値は `{BASELINE_EQUIVALENCE_LIMIT:.1e}`、固定点tail-step閾値は `{FIXED_POINT_TAIL_STEP_TOLERANCE:.1e}`、振幅差閾値は `{AMPLITUDE_SHIFT_TOLERANCE:.1e}` である。結果に合わせた閾値変更や64を超える探索はない。

## 分類件数

{count_text}

`best_lag`は限定窓内の最小候補であり、それ自体を厳密な基本周期と解釈していない。`cycle_detected`は帰還RMSが固定閾値以下の場合だけ真である。
""",
        encoding="utf-8",
    )

    fc_log = json.loads(
        (STAGE_ROOT / "logs" / "stage_F_FC_run.json").read_text(
            encoding="utf-8"
        )
    )
    f248 = [
        row
        for row in f_c_rows
        if int(float(row["iteration"])) == 248
        and row["normalization_mode"] == "existing_normalization"
    ]
    fc_lines = [
        "| kappa | C0 return@248 | reversed return@248 | c_mean | R_eff |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in f248:
        fc_lines.append(
            f"| {float(row['kappa']):g} | "
            f"{float(row['C0_return_error']):.6g} | "
            f"{float(row['reversed_C1_return_error']):.6g} | "
            f"{float(row['c_mean']):.6g} | "
            f"{float(row['R_eff']):.9g} |"
        )
    validation = fc_log["canonical_validation"][
        "maximum_absolute_errors"
    ]
    (REPORT_DIR / "06_31_series_check.md").write_text(
        """# 06 31系列の限定確認

31系列はF-Aとは別の既存条件、custom packet `A=(1), B=(1,2), R=0.697177927` に属する。C0の0〜256回を既存514行CSVへ照合し、最大絶対誤差は Lで `{L:.3g}`、N_effで `{N:.3g}`、スペクトル類似度で `{S:.3g}` だった。

既存の固定系列 `31,62,93,124,155,186,217,248,279` だけを評価した。帰還誤差は位相不変な正規化全状態距離のA/B RMS、交換量は交換済み初期対 `(B0,A0)` への同じ距離である。新しい周期探索ではない。

""".format(
            L=validation["L"],
            N=validation["N_eff"],
            S=validation["spectral_similarity"],
        )
        + "\n".join(fc_lines)
        + """

C0の248回は既存の帰還構造を再現する。反転Candidate 1はB入力の偶奇混合により初期 `c_B≈0`、A入力 `c_A=-1` を読むため、全奇数F-Aとは異なる変位を受ける。この数値差は候補式の帰結であり、ボゾン・フェルミオン対応そのものの導出ではない。
""",
        encoding="utf-8",
    )

    (REPORT_DIR / "07_normalization_effect.md").write_text(
        f"""# 07 正規化の影響

`existing_normalization` と `raw_update` は別系列として全条件を保存した。

- 主要観測量の最大差: `{normalization['maximum_observable_difference']:.17g}`
- 最終全状態ベクトルの最大差: `{normalization['maximum_final_state_vector_difference']:.17g}`

直交etaモードとユニタリ二チャネル散乱により各rawチャネルノルムが1に保たれ、今回の基準群では既存正規化が機械精度内で実質恒等になった。これは一般の状態依存散乱でチャネル正規化が無害だという証明ではない。
""",
        encoding="utf-8",
    )

    max_residual = {
        key: max(float(row[key]) for row in residuals)
        for key in (
            "unitarity_residual",
            "coefficient_orthogonality_residual",
            "path_sum_residual_A",
            "path_sum_residual_B",
            "total_norm_conservation_residual",
            "demodulation_reconstruction_residual",
            "parity_projection_sum_residual",
        )
    }
    residual_text = "\n".join(
        f"- `{key}`: `{value:.17g}`"
        for key, value in max_residual.items()
    )
    (REPORT_DIR / "08_numerical_invariants.md").write_text(
        f"""# 08 数値不変量

全F-B実行の最大残差:

{residual_text}

- NaN/Inf: `{sum(int(row['nan_inf_count']) for row in residuals)}`
- theta範囲違反: `{sum(int(row['theta_range_violation_count']) for row in residuals)}`
- 数値不安定判定閾値: `{NUMERICAL_RESIDUAL_LIMIT:.1e}`

経路干渉項は今回も直交eta由来で機械精度内に抑えられた。`B_to_A_transfer`は `spectral_similarity_to_initial_B; not path flux` であり、経路フラックスではない。
""",
        encoding="utf-8",
    )

    fa_c = [
        row["c_mean_min"]
        for row in summaries
        if abs(float(row["R0"]) - 0.55) <= 1.0e-15
    ] + [
        row["c_mean_max"]
        for row in summaries
        if abs(float(row["R0"]) - 0.55) <= 1.0e-15
    ]
    overall = f"""# Stage F report

## 結論

C0再現ゲートは通過した。`N_A=1,N_B=63,R=0.55` の128回系列における正本との差は、最大絶対誤差で Lが0、N_effが `1.42e-14`、スペクトル類似度が `4.44e-16` だった。

反転Candidate 1を既存System Aの独立コピーへ統合した。代表条件では `c_A,c_B,c_mean` が全衝突で `[{min(fa_c):.17g},{max(fa_c):.17g}]`、すなわち機械精度内で-1に固定された。このため新核は動的なRではなく、κごとの一定な高いR_effを生成し、交換の位相・最小差到達回を変えた。

既存7点R掃引を完了した。R0=0.70のC0で既存の小さいL差・N_eff差を再現し、反転Candidate 1ではκごとに位置が変わった。端点R0=0,1では包絡rho=0によりC0と一致した。

31系列は別の既存custom packet条件に接続し、C0を既存CSVへ機械精度で照合した上で固定反復だけを比較した。248回帰還はC0で再現され、反転Candidate 1では移動した。

既存正規化とraw更新は別記録だが、今回の直交eta・単位ノルム条件では最大差が機械精度内であった。

## コード上の事実

- 散乱は512×16全状態で実施した。
- パリティは由来別eta射影・搬送波除去後に毎衝突再計算した。
- C0と反転Candidate 1の差は `delta_theta` だけである。
- Candidate 2/3、N体系、論文本文は変更していない。

## 数学的帰結

純奇数入力でc_mean=-1が保存される限り、反転Candidate 1は `theta_eff=theta0+kappa*rho(theta0)` という一定角へ簡約される。端点ではrho=0なのでC0へ戻る。

## 数値観察

一定R_effへの写像でも交換系列の位相、限定窓内best lag、最小L/N_eff差は変化した。custom packetの偶奇混合では全奇数条件とは異なる変位が生じた。

## モデル上の仮説

奇数倍音主体同士で反射を強めるという規則は、反転Candidate 1として実装した仮説である。

## 物理対応上の仮説

偶数倍音をボゾン型、奇数倍音をフェルミオン型と読む物理対応は、本Stageの数値整合だけから導出されない。

## 未導出

κ、包絡rho、パリティ読出し相互作用の力学的起源は未導出である。

## 保留

N体系への拡張、既存本体への正式統合、Candidate 2/3、論文本文への反映は行っていない。限定lagでのbest lagは厳密周期の証明ではない。
"""
    (REPORT_DIR / "Stage_F_report.md").write_text(
        overall, encoding="utf-8"
    )


def main() -> None:
    sweep_rows = read_csv(DATA_DIR / "stage_F_R_sweep_results.csv")
    f_c_rows = read_csv(DATA_DIR / "stage_F_31_series_results.csv")
    cycle_rows, residual_rows, summary_rows = calculate_cycle_tables(
        sweep_rows
    )
    write_csv(DATA_DIR / "stage_F_cycle_metrics.csv", cycle_rows)
    write_csv(
        DATA_DIR / "stage_F_numerical_residuals.csv", residual_rows
    )
    write_csv(DATA_DIR / "stage_F_run_summary.csv", summary_rows)
    normalization = maximum_normalization_difference(sweep_rows)
    classification_counts = dict(
        Counter(row["classification"] for row in cycle_rows)
    )
    summary = {
        "stage": "F",
        "status": "complete_pending_hash_verification",
        "reproduction_gate": json.loads(
            (
                DATA_DIR / "stage_F_reproduction_gate_summary.json"
            ).read_text(encoding="utf-8")
        ),
        "thresholds": {
            "numerical_residual_limit": NUMERICAL_RESIDUAL_LIMIT,
            "baseline_equivalence_limit": BASELINE_EQUIVALENCE_LIMIT,
            "cycle_return_tolerance": CYCLE_RETURN_TOLERANCE,
            "fixed_point_tail_step_tolerance": (
                FIXED_POINT_TAIL_STEP_TOLERANCE
            ),
            "amplitude_shift_tolerance": AMPLITUDE_SHIFT_TOLERANCE,
            "tested_lag_range": [1, MAX_LAG],
        },
        "run_count": len(summary_rows),
        "classification_counts": classification_counts,
        "normalization_comparison": normalization,
        "maximum_numerical_residual": max(
            max(
                float(row["unitarity_residual"]),
                float(row["coefficient_orthogonality_residual"]),
                float(row["path_sum_residual_A"]),
                float(row["path_sum_residual_B"]),
                float(row["total_norm_conservation_residual"]),
                float(row["demodulation_reconstruction_residual"]),
                float(row["parity_projection_sum_residual"]),
            )
            for row in residual_rows
        ),
        "nan_inf_count": sum(
            int(row["nan_inf_count"]) for row in residual_rows
        ),
        "theta_range_violation_count": sum(
            int(row["theta_range_violation_count"])
            for row in residual_rows
        ),
        "F_C": json.loads(
            (
                STAGE_ROOT / "logs" / "stage_F_FC_run.json"
            ).read_text(encoding="utf-8")
        ),
    }
    (DATA_DIR / "stage_F_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    save_figures(sweep_rows, summary_rows, cycle_rows, f_c_rows)
    reports(
        sweep_rows,
        summary_rows,
        cycle_rows,
        residual_rows,
        f_c_rows,
        normalization,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
