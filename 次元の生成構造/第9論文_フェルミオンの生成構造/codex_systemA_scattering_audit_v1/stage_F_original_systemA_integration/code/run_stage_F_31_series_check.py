"""Limited Stage F-C check of the already observed custom-packet 31 series."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from parity_demodulation import channel_parity, norm2
from system_A_stage_F_copy import (
    KAPPA_VALUES,
    Params,
    custom_31_initial_state_pair,
    harmonic_distribution,
    localization,
    effective_n,
    distribution_similarity,
    make_grids,
    run_series,
)


STAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[5]
DATA_DIR = STAGE_ROOT / "data"
OUTPUT = DATA_DIR / "stage_F_31_series_results.csv"
LOG_OUTPUT = STAGE_ROOT / "logs" / "stage_F_FC_run.json"
CANONICAL_ROWS = (
    REPO_ROOT
    / "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/20260715/instrument_check_v1/"
    "base_run/system_A_custom_packet_A1_B1-2_R0-0p697177927_C256_rows_v1.csv"
)
SERIES_ITERATIONS = (31, 62, 93, 124, 155, 186, 217, 248, 279)
EVALUATION_ITERATIONS = (
    31,
    62,
    93,
    124,
    155,
    186,
    217,
    247,
    248,
    279,
)
R_31 = 0.697177927
VALIDATION_TOLERANCE = 5.0e-10


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


def initial_metric_row(
    params: Params, a: np.ndarray, b: np.ndarray
) -> dict[str, float]:
    h_a = harmonic_distribution(params, a)
    h_b = harmonic_distribution(params, b)
    return {
        "collision_index": 0,
        "L_A": localization(a),
        "L_B": localization(b),
        "N_eff_A": effective_n(h_a),
        "N_eff_B": effective_n(h_b),
        "spectral_similarity_A_to_initial_A": distribution_similarity(
            h_a, h_a
        ),
        "spectral_similarity_A_to_initial_B": distribution_similarity(
            h_a, h_b
        ),
        "spectral_similarity_B_to_initial_A": distribution_similarity(
            h_b, h_a
        ),
        "spectral_similarity_B_to_initial_B": distribution_similarity(
            h_b, h_b
        ),
    }


def validate_c0(
    params: Params,
    result_rows: list[dict],
    initial_a: np.ndarray,
    initial_b: np.ndarray,
) -> dict:
    integrated = {
        int(row["collision_index"]): row for row in result_rows
    }
    integrated[0] = initial_metric_row(params, initial_a, initial_b)
    with CANONICAL_ROWS.open(encoding="utf-8", newline="") as handle:
        source_rows = [
            row
            for row in csv.DictReader(handle)
            if abs(float(row["R"]) - R_31) <= 1.0e-12
        ]
    grouped: dict[int, dict[str, dict]] = {}
    for row in source_rows:
        grouped.setdefault(int(row["collision"]), {})[
            str(row["channel"])
        ] = row
    maximums = {
        "L": 0.0,
        "N_eff": 0.0,
        "spectral_similarity": 0.0,
    }
    for collision, channels in grouped.items():
        candidate = integrated[collision]
        source_a = channels["A_channel"]
        source_b = channels["B_channel"]
        maximums["L"] = max(
            maximums["L"],
            abs(float(source_a["L"]) - float(candidate["L_A"])),
            abs(float(source_b["L"]) - float(candidate["L_B"])),
        )
        maximums["N_eff"] = max(
            maximums["N_eff"],
            abs(
                float(source_a["N_eff"])
                - float(candidate["N_eff_A"])
            ),
            abs(
                float(source_b["N_eff"])
                - float(candidate["N_eff_B"])
            ),
        )
        maximums["spectral_similarity"] = max(
            maximums["spectral_similarity"],
            abs(
                float(source_a["sim_to_A0"])
                - float(
                    candidate[
                        "spectral_similarity_A_to_initial_A"
                    ]
                )
            ),
            abs(
                float(source_a["sim_to_B0"])
                - float(
                    candidate[
                        "spectral_similarity_A_to_initial_B"
                    ]
                )
            ),
            abs(
                float(source_b["sim_to_A0"])
                - float(
                    candidate[
                        "spectral_similarity_B_to_initial_A"
                    ]
                )
            ),
            abs(
                float(source_b["sim_to_B0"])
                - float(
                    candidate[
                        "spectral_similarity_B_to_initial_B"
                    ]
                )
            ),
        )
    return {
        "canonical_row_count": len(source_rows),
        "canonical_collision_count": len(grouped) - 1,
        "maximum_absolute_errors": maximums,
        "passed": max(maximums.values()) <= VALIDATION_TOLERANCE,
    }


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    params = Params()
    initial_a, initial_b = custom_31_initial_state_pair(params)
    snapshots = (0, 31, 62, 93, 124, 155, 186, 217, 247, 248, 279)
    c0_by_normalization = {}
    for normalization_mode in ("existing_normalization", "raw_update"):
        c0_by_normalization[normalization_mode] = run_series(
            params,
            reflection_baseline=R_31,
            scattering_mode="C0",
            normalization_mode=normalization_mode,
            kappa=0.01,
            collision_count=279,
            initial_state_override=(initial_a, initial_b),
            run_label="custom31_A1_B1plus2",
            snapshot_indices=snapshots,
        )

    validation = validate_c0(
        params,
        c0_by_normalization["existing_normalization"].rows,
        initial_a,
        initial_b,
    )
    if not validation["passed"]:
        result = {
            "status": "not_reproducible",
            "reason": "custom-packet C0 did not reproduce canonical rows",
            "validation": validation,
        }
        LOG_OUTPUT.write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        raise SystemExit(json.dumps(result, ensure_ascii=False))

    u, eta = make_grids(params)
    rows: list[dict] = []
    final_states: dict[str, np.ndarray] = {}
    for normalization_mode, c0_result in c0_by_normalization.items():
        c0_rows = {
            int(row["collision_index"]): row for row in c0_result.rows
        }
        for kappa in KAPPA_VALUES:
            reversed_result = run_series(
                params,
                reflection_baseline=R_31,
                scattering_mode="reversed_C1",
                normalization_mode=normalization_mode,
                kappa=kappa,
                collision_count=279,
                initial_state_override=(initial_a, initial_b),
                run_label="custom31_A1_B1plus2",
                snapshot_indices=snapshots,
            )
            reversed_rows = {
                int(row["collision_index"]): row
                for row in reversed_result.rows
            }
            for iteration in EVALUATION_ITERATIONS:
                c0_a, c0_b = c0_result.snapshots[iteration]
                new_a, new_b = reversed_result.snapshots[iteration]
                c0_row = c0_rows[iteration]
                new_row = reversed_rows[iteration]
                rows.append(
                    {
                        "normalization_mode": normalization_mode,
                        "kappa": kappa,
                        "iteration": iteration,
                        "probe_kind": (
                            "31_series"
                            if iteration in SERIES_ITERATIONS
                            else "247_near_return_probe"
                        ),
                        "C0_return_error": pair_distance(
                            c0_a, c0_b, initial_a, initial_b
                        ),
                        "reversed_C1_return_error": pair_distance(
                            new_a, new_b, initial_a, initial_b
                        ),
                        "C0_exchange_measure": pair_distance(
                            c0_a, c0_b, initial_b, initial_a
                        ),
                        "reversed_C1_exchange_measure": pair_distance(
                            new_a, new_b, initial_b, initial_a
                        ),
                        "C0_L_difference": abs(
                            c0_row["L_A"] - c0_row["L_B"]
                        ),
                        "reversed_C1_L_difference": abs(
                            new_row["L_A"] - new_row["L_B"]
                        ),
                        "C0_N_eff_difference": abs(
                            c0_row["N_eff_A"] - c0_row["N_eff_B"]
                        ),
                        "reversed_C1_N_eff_difference": abs(
                            new_row["N_eff_A"] - new_row["N_eff_B"]
                        ),
                        "C0_c_mean": c0_row["c_mean"],
                        "c_mean": new_row["c_mean"],
                        "C0_R_eff": c0_row["R_eff"],
                        "R_eff": new_row["R_eff"],
                    }
                )
            key = (
                f"{normalization_mode}__kappa_{kappa:g}"
            ).replace(".", "p")
            final_states[f"{key}__A"] = reversed_result.final_a
            final_states[f"{key}__B"] = reversed_result.final_b

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(
        STAGE_ROOT / "logs" / "stage_F_FC_final_states.npz",
        **final_states,
    )
    initial_parity_a = channel_parity(initial_a, u, eta, params.p0)
    initial_parity_b = channel_parity(initial_b, u, eta, params.p0)
    result = {
        "status": "complete",
        "experiment": "F-C",
        "condition": {
            "state_family": "custom_packet",
            "A_harmonics": [1],
            "B_harmonics": [1, 2],
            "R": R_31,
            "31_series_iterations": list(SERIES_ITERATIONS),
            "fixed_evaluation_iterations": list(
                EVALUATION_ITERATIONS
            ),
            "collision_count": 279,
        },
        "method_scope": (
            "This is the distinct existing custom-packet 31-series "
            "condition, not the N_A=1,N_B=63,R=0.55 F-A condition."
        ),
        "return_error_definition": (
            "RMS of phase-insensitive normalized full-state distances "
            "from (A0,B0)"
        ),
        "exchange_measure_definition": (
            "RMS of phase-insensitive normalized full-state distances "
            "from the swapped pair (B0,A0)"
        ),
        "initial_c_A": initial_parity_a.c_pi,
        "initial_c_B": initial_parity_b.c_pi,
        "canonical_validation": validation,
        "row_count": len(rows),
    }
    LOG_OUTPUT.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
