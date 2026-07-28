"""Deterministic tests for the Stage E repeated-scattering experiment copy."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from system_A_experimental_copy import (
    ExperimentConfig,
    SPEC_ORIGIN_A,
    SPEC_ORIGIN_B,
    boson_kernel,
    channel_parity,
    effective_angle,
    fermion_kernel,
    initial_state_pair,
    make_grids,
    modulate_kernel,
    norm2,
    run_repeated,
    scatter_once,
)


STAGE_E_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = STAGE_E_ROOT / "data" / "stage_E_test_results.json"
TOLERANCE = 5.0e-12


def _assert_below(name: str, value: float, limit: float = TOLERANCE) -> None:
    if value > limit:
        raise AssertionError(f"{name}: {value} > {limit}")


def run_tests() -> dict:
    config = ExperimentConfig()
    u, eta = make_grids(config)
    fermion = fermion_kernel(u, config.k_components)
    boson = boson_kernel(u, config.k_components)
    state_f = modulate_kernel(
        fermion, SPEC_ORIGIN_A, u, eta, config.p0
    )
    state_b = modulate_kernel(
        boson, SPEC_ORIGIN_B, u, eta, config.p0
    )
    metrics_f = channel_parity(state_f, u, eta, config.p0)
    metrics_b = channel_parity(state_b, u, eta, config.p0)
    _assert_below("fermion parity", abs(metrics_f.indicator + 1.0))
    _assert_below("boson parity", abs(metrics_b.indicator - 1.0))
    _assert_below(
        "initial reconstruction",
        max(
            metrics_f.reconstruction_residual,
            metrics_b.reconstruction_residual,
        ),
    )

    angle_ff = effective_angle(
        "C1_reversed",
        config.reflection_baseline,
        1.0,
        -1.0,
        -1.0,
    )
    angle_bb = effective_angle(
        "C1_reversed",
        config.reflection_baseline,
        1.0,
        1.0,
        1.0,
    )
    angle_fb = effective_angle(
        "C1_reversed",
        config.reflection_baseline,
        1.0,
        -1.0,
        1.0,
    )
    if not (
        angle_ff.reflection_probability > config.reflection_baseline
        and angle_bb.reflection_probability < config.reflection_baseline
    ):
        raise AssertionError("reversed Candidate 1 direction is incorrect")
    _assert_below(
        "unlike-pair baseline",
        abs(angle_fb.reflection_probability - config.reflection_baseline),
    )

    maximums = {
        "unitarity_residual": 0.0,
        "coefficient_orthogonality_residual": 0.0,
        "path_sum_residual": 0.0,
        "pair_norm_conservation_residual": 0.0,
        "lineage_reconstruction_residual": 0.0,
        "raw_channel_norm_error": 0.0,
        "pure_response_drift": 0.0,
        "unlike_response_from_baseline": 0.0,
        "unlike_parity_mean": 0.0,
        "C0_kappa_duplicate_difference": 0.0,
        "unlike_C0_vs_reversed_final_state_difference": 0.0,
    }
    all_runs = {}
    for kernel_name in ("C0", "C1_reversed"):
        for kappa in (0.01, 0.1, 1.0):
            for case_id in ("F_x_F", "B_x_B", "F_x_B"):
                rows, final_a, final_b = run_repeated(
                    case_id, kernel_name, kappa, config
                )
                if len(rows) != config.collision_count:
                    raise AssertionError("collision count mismatch")
                all_runs[(kernel_name, kappa, case_id)] = (
                    rows,
                    final_a,
                    final_b,
                )
                first_r = float(rows[0]["R_eff"])
                for row in rows:
                    maximums["unitarity_residual"] = max(
                        maximums["unitarity_residual"],
                        abs(float(row["unitarity_residual"])),
                    )
                    maximums["coefficient_orthogonality_residual"] = max(
                        maximums["coefficient_orthogonality_residual"],
                        abs(
                            float(
                                row[
                                    "coefficient_orthogonality_residual"
                                ]
                            )
                        ),
                    )
                    maximums["path_sum_residual"] = max(
                        maximums["path_sum_residual"],
                        abs(float(row["path_sum_residual_A"])),
                        abs(float(row["path_sum_residual_B"])),
                    )
                    maximums["pair_norm_conservation_residual"] = max(
                        maximums["pair_norm_conservation_residual"],
                        abs(
                            float(
                                row["pair_norm_conservation_residual"]
                            )
                        ),
                    )
                    maximums["lineage_reconstruction_residual"] = max(
                        maximums["lineage_reconstruction_residual"],
                        abs(
                            float(
                                row["input_reconstruction_residual_A"]
                            )
                        ),
                        abs(
                            float(
                                row["input_reconstruction_residual_B"]
                            )
                        ),
                        abs(
                            float(
                                row["raw_reconstruction_residual_A"]
                            )
                        ),
                        abs(
                            float(
                                row["raw_reconstruction_residual_B"]
                            )
                        ),
                    )
                    maximums["raw_channel_norm_error"] = max(
                        maximums["raw_channel_norm_error"],
                        abs(float(row["a_raw_norm2"]) - 1.0),
                        abs(float(row["b_raw_norm2"]) - 1.0),
                    )
                    if kernel_name == "C1_reversed":
                        if case_id in ("F_x_F", "B_x_B"):
                            maximums["pure_response_drift"] = max(
                                maximums["pure_response_drift"],
                                abs(float(row["R_eff"]) - first_r),
                            )
                        else:
                            maximums[
                                "unlike_response_from_baseline"
                            ] = max(
                                maximums[
                                    "unlike_response_from_baseline"
                                ],
                                abs(
                                    float(row["R_eff"])
                                    - config.reflection_baseline
                                ),
                            )
                            maximums["unlike_parity_mean"] = max(
                                maximums["unlike_parity_mean"],
                                abs(float(row["c_mean"])),
                            )
                    if row["channel_normalization_applied"]:
                        raise AssertionError(
                            "channel normalization was unexpectedly applied"
                        )

    for case_id in ("F_x_F", "B_x_B", "F_x_B"):
        reference_rows, reference_a, reference_b = all_runs[
            ("C0", 0.01, case_id)
        ]
        for kappa in (0.1, 1.0):
            rows, final_a, final_b = all_runs[("C0", kappa, case_id)]
            row_difference = max(
                abs(float(left["R_eff"]) - float(right["R_eff"]))
                for left, right in zip(reference_rows, rows)
            )
            state_difference = math.sqrt(
                norm2(reference_a - final_a)
                + norm2(reference_b - final_b)
            )
            maximums["C0_kappa_duplicate_difference"] = max(
                maximums["C0_kappa_duplicate_difference"],
                row_difference,
                state_difference,
            )

    for kappa in (0.01, 0.1, 1.0):
        _, c0_a, c0_b = all_runs[("C0", kappa, "F_x_B")]
        _, new_a, new_b = all_runs[
            ("C1_reversed", kappa, "F_x_B")
        ]
        difference = math.sqrt(
            norm2(c0_a - new_a) + norm2(c0_b - new_b)
        )
        maximums["unlike_C0_vs_reversed_final_state_difference"] = max(
            maximums["unlike_C0_vs_reversed_final_state_difference"],
            difference,
        )

    # Explicitly verify that raw output arrays, without normalization, are
    # passed to the next collision.
    initial_a, initial_b = initial_state_pair(
        "F_x_B", config, u, eta
    )
    first = scatter_once(
        initial_a,
        initial_b,
        kernel_name="C1_reversed",
        reflection_baseline=config.reflection_baseline,
        kappa=1.0,
        u=u,
        eta=eta,
        p0=config.p0,
    )
    second = scatter_once(
        first.raw_output_a,
        first.raw_output_b,
        kernel_name="C1_reversed",
        reflection_baseline=config.reflection_baseline,
        kappa=1.0,
        u=u,
        eta=eta,
        p0=config.p0,
    )
    _assert_below(
        "raw propagation A",
        math.sqrt(norm2(second.input_a - first.raw_output_a)),
    )
    _assert_below(
        "raw propagation B",
        math.sqrt(norm2(second.input_b - first.raw_output_b)),
    )

    for name, value in maximums.items():
        _assert_below(name, value)

    payload = {
        "schema": "stage_E_repeated_scattering_tests_v1",
        "status": "PASS",
        "tolerance": TOLERANCE,
        "configuration": {
            "grid": [config.u_grid_n, config.eta_grid_n],
            "K": config.k_components,
            "R_0": config.reflection_baseline,
            "collisions": config.collision_count,
            "kappa_values": [0.01, 0.1, 1.0],
        },
        "maximums": maximums,
        "verified_direction_at_kappa_1": {
            "F_x_F_R_eff": angle_ff.reflection_probability,
            "B_x_B_R_eff": angle_bb.reflection_probability,
            "F_x_B_R_eff": angle_fb.reflection_probability,
        },
        "raw_outputs_propagated_without_channel_normalization": True,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(run_tests(), indent=2, ensure_ascii=False))
