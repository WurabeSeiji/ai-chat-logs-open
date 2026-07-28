"""Stage G-C limited rerun at the existing custom-packet 31-series points."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from parity_demodulation import norm2
from system_A_stage_G_copy import (
    KAPPA_VALUES,
    Params,
    custom_31_initial_state_pair,
    run_series,
)


STAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[5]
F_RESULTS = (
    REPO_ROOT
    / "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "codex_systemA_scattering_audit_v1/"
    "stage_F_original_systemA_integration/data/"
    "stage_F_31_series_results.csv"
)
OUTPUT = STAGE_ROOT / "data" / "stage_G_31_series_results.csv"
FINAL_STATES = STAGE_ROOT / "data" / "stage_G_final_states.npz"
LOG_OUTPUT = STAGE_ROOT / "logs" / "stage_G_31_series_run.json"
R_31 = 0.697177927
FIXED_ITERATIONS = (31, 62, 93, 124, 155, 186, 217, 247, 248, 279)
MAX_AUTOCORRELATION_LAG = 64
REGRESSION_TOLERANCE = 5.0e-10


def normalized_distance(left: np.ndarray, right: np.ndarray) -> float:
    left_norm = math.sqrt(max(norm2(left), 0.0))
    right_norm = math.sqrt(max(norm2(right), 0.0))
    if left_norm <= 0.0 or right_norm <= 0.0:
        return float("inf")
    overlap = abs(
        np.vdot(left / left_norm, right / right_norm)
    )
    return math.sqrt(max(2.0 - 2.0 * float(overlap), 0.0))


def pair_distance(
    a: np.ndarray,
    b: np.ndarray,
    target_a: np.ndarray,
    target_b: np.ndarray,
) -> float:
    return math.sqrt(
        (
            normalized_distance(a, target_a) ** 2
            + normalized_distance(b, target_b) ** 2
        )
        / 2.0
    )


def autocorrelation_peak(rows: list[dict]) -> tuple[int, float]:
    matrix = np.asarray(
        [
            [
                row["L_A"],
                row["L_B"],
                row["N_eff_A"],
                row["N_eff_B"],
                row["Gamma_AB"],
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
        ]
    )
    values = matrix / scales
    candidates = []
    for lag in range(1, MAX_AUTOCORRELATION_LAG + 1):
        left = values[lag:].reshape(-1)
        right = values[:-lag].reshape(-1)
        if np.std(left) <= 1.0e-300 or np.std(right) <= 1.0e-300:
            correlation = 1.0 if np.array_equal(left, right) else 0.0
        else:
            correlation = float(np.corrcoef(left, right)[0, 1])
        candidates.append((lag, correlation))
    return max(candidates, key=lambda item: item[1])


def stage_f_reference() -> list[dict]:
    with F_RESULTS.open(encoding="utf-8", newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle)
            if row["normalization_mode"] == "raw_update"
        ]


def validate_fixed_rows(rows: list[dict], references: list[dict]) -> dict:
    maximum = 0.0
    comparisons = 0
    for row in rows:
        mode = row["scattering_mode"]
        if mode not in ("C0", "reversed_C1"):
            continue
        matches = [
            ref
            for ref in references
            if abs(float(ref["kappa"]) - float(row["kappa"]))
            <= 1.0e-15
            and int(ref["iteration"]) == int(row["iteration"])
        ]
        if len(matches) != 1:
            raise ValueError("Stage F 31-series reference mismatch")
        ref = matches[0]
        prefix = "C0" if mode == "C0" else "reversed_C1"
        fields = (
            ("return_error", f"{prefix}_return_error"),
            ("exchange_measure", f"{prefix}_exchange_measure"),
            ("L_difference", f"{prefix}_L_difference"),
            ("N_eff_difference", f"{prefix}_N_eff_difference"),
        )
        for current_field, reference_field in fields:
            maximum = max(
                maximum,
                abs(
                    float(row[current_field])
                    - float(ref[reference_field])
                ),
            )
            comparisons += 1
        c_reference_field = (
            "C0_c_mean" if mode == "C0" else "c_mean"
        )
        r_reference_field = (
            "C0_R_eff" if mode == "C0" else "R_eff"
        )
        maximum = max(
            maximum,
            abs(
                float(row["c_mean"])
                - float(ref[c_reference_field])
            ),
            abs(
                float(row["R_eff"])
                - float(ref[r_reference_field])
            ),
        )
        comparisons += 2
    return {
        "passed": maximum <= REGRESSION_TOLERANCE,
        "comparison_count": comparisons,
        "maximum_absolute_error": maximum,
        "tolerance": REGRESSION_TOLERANCE,
    }


def main() -> None:
    params = Params()
    initial_a, initial_b = custom_31_initial_state_pair(params)
    rows: list[dict] = []
    summaries = []
    final_state_updates: dict[str, np.ndarray] = {}
    all_snapshot_indices = range(280)
    for scattering_mode in ("C0", "reversed_C1", "relational_C1"):
        for kappa in KAPPA_VALUES:
            result = run_series(
                params,
                reflection_baseline=R_31,
                scattering_mode=scattering_mode,
                normalization_mode="raw_update",
                kappa=kappa,
                collision_count=279,
                initial_state_override=(initial_a, initial_b),
                run_label="custom31_A1_B1plus2",
                snapshot_indices=all_snapshot_indices,
            )
            by_collision = {
                int(row["collision_index"]): row for row in result.rows
            }
            return_errors = {}
            exchange_measures = {}
            for iteration, (a, b) in result.snapshots.items():
                return_errors[iteration] = pair_distance(
                    a, b, initial_a, initial_b
                )
                exchange_measures[iteration] = pair_distance(
                    a, b, initial_b, initial_a
                )
            best_return_iteration = min(
                range(1, 280), key=lambda index: return_errors[index]
            )
            best_exchange_iteration = min(
                range(1, 280),
                key=lambda index: exchange_measures[index],
            )
            autocorrelation_lag, autocorrelation = (
                autocorrelation_peak(result.rows)
            )
            summaries.append(
                {
                    "scattering_mode": scattering_mode,
                    "kappa": kappa,
                    "minimum_return_error": return_errors[
                        best_return_iteration
                    ],
                    "minimum_return_error_iteration": (
                        best_return_iteration
                    ),
                    "minimum_exchange_measure": exchange_measures[
                        best_exchange_iteration
                    ],
                    "minimum_exchange_measure_iteration": (
                        best_exchange_iteration
                    ),
                    "deviation_from_32_exchange": abs(
                        best_exchange_iteration - 32
                    ),
                    "autocorrelation_peak_lag": autocorrelation_lag,
                    "autocorrelation_peak": autocorrelation,
                    "return_error_248": return_errors[248],
                    "Gamma_AB_min": min(
                        row["Gamma_AB"] for row in result.rows
                    ),
                    "Gamma_AB_max": max(
                        row["Gamma_AB"] for row in result.rows
                    ),
                    "Gamma_AB_range": max(
                        row["Gamma_AB"] for row in result.rows
                    )
                    - min(row["Gamma_AB"] for row in result.rows),
                    "R_eff_min": min(
                        row["R_eff"] for row in result.rows
                    ),
                    "R_eff_max": max(
                        row["R_eff"] for row in result.rows
                    ),
                    "R_eff_range": max(
                        row["R_eff"] for row in result.rows
                    )
                    - min(row["R_eff"] for row in result.rows),
                    "unitarity_residual_max": max(
                        row["unitarity_residual"]
                        for row in result.rows
                    ),
                    "orthogonality_residual_max": max(
                        row["orthogonality_residual"]
                        for row in result.rows
                    ),
                    "path_sum_residual_max": max(
                        max(
                            row["path_sum_residual_A"],
                            row["path_sum_residual_B"],
                        )
                        for row in result.rows
                    ),
                    "total_norm_residual_max": max(
                        row["total_norm_residual"]
                        for row in result.rows
                    ),
                    "demodulation_residual_max": max(
                        max(
                            row[
                                "demodulation_reconstruction_residual_A"
                            ],
                            row[
                                "demodulation_reconstruction_residual_B"
                            ],
                        )
                        for row in result.rows
                    ),
                    "gamma_range_violation_max": max(
                        max(
                            0.0,
                            -row["Gamma_AB"],
                            row["Gamma_AB"] - 1.0,
                        )
                        for row in result.rows
                    ),
                    "relation_wave_norm2_min": min(
                        min(
                            row["relation_wave_norm2_A"],
                            row["relation_wave_norm2_B"],
                        )
                        for row in result.rows
                    ),
                    "nan_inf_count": sum(
                        int(row["nan_inf_count"])
                        for row in result.rows
                    ),
                    "theta_range_violation_count": sum(
                        int(bool(row["theta_range_violation"]))
                        for row in result.rows
                    ),
                }
            )
            for iteration in FIXED_ITERATIONS:
                row = by_collision[iteration]
                rows.append(
                    {
                        "scattering_mode": scattering_mode,
                        "kappa": kappa,
                        "R0": R_31,
                        "iteration": iteration,
                        "return_error": return_errors[iteration],
                        "exchange_measure": exchange_measures[iteration],
                        "L_difference": abs(
                            row["L_A"] - row["L_B"]
                        ),
                        "N_eff_difference": abs(
                            row["N_eff_A"] - row["N_eff_B"]
                        ),
                        "c_A": row["c_A"],
                        "c_B": row["c_B"],
                        "c_mean": row["c_mean"],
                        "overlap_complex_real": row[
                            "overlap_complex_real"
                        ],
                        "overlap_complex_imag": row[
                            "overlap_complex_imag"
                        ],
                        "Gamma_AB": row["Gamma_AB"],
                        "candidate_response": row[
                            "candidate_response"
                        ],
                        "theta_eff": row["theta_eff"],
                        "R_eff": row["R_eff"],
                    }
                )
            key = (
                f"custom31__{scattering_mode}__raw_update__"
                f"kappa_{kappa:g}__R0_{R_31:g}"
            ).replace(".", "p")
            final_state_updates[f"{key}__A"] = result.final_a
            final_state_updates[f"{key}__B"] = result.final_b

    validation = validate_fixed_rows(rows, stage_f_reference())
    if not validation["passed"]:
        raise SystemExit(
            "Stage G-C C0/reversed regression failed; stop"
        )
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    existing_archive = np.load(FINAL_STATES)
    combined = {
        key: existing_archive[key].copy() for key in existing_archive.files
    }
    combined.update(final_state_updates)
    np.savez_compressed(FINAL_STATES, **combined)
    result = {
        "status": "complete",
        "condition": {
            "state_family": "custom_packet",
            "A_harmonics": [1],
            "B_harmonics": [1, 2],
            "R0": R_31,
            "collision_count": 279,
            "fixed_iterations": list(FIXED_ITERATIONS),
        },
        "autocorrelation_lag_range": [
            1,
            MAX_AUTOCORRELATION_LAG,
        ],
        "row_count": len(rows),
        "run_count": len(summaries),
        "stage_F_regression": validation,
        "run_summaries": summaries,
    }
    LOG_OUTPUT.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "run_count": result["run_count"],
                "row_count": result["row_count"],
                "stage_F_regression": validation,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
