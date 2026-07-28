"""Deterministic Stage G unit checks, independent of experiment outputs."""

from __future__ import annotations

import math

import numpy as np

from parity_demodulation import (
    ORIGIN_A,
    ORIGIN_B,
    modulate_kernel,
    norm2,
)
from relational_scattering import (
    RelationUndefinedError,
    effective_angle,
    relation_readout,
    scatter_once,
)
from system_A_stage_G_copy import Params, make_grids, normalize


TOLERANCE = 5.0e-12


def normalized_kernel(values: np.ndarray) -> np.ndarray:
    return normalize(np.asarray(values, dtype=complex))


def test_row(
    case_id: str,
    quantity: str,
    actual: float,
    expected: float,
    tolerance: float = TOLERANCE,
) -> dict:
    error = abs(actual - expected)
    return {
        "case_id": case_id,
        "quantity": quantity,
        "expected": expected,
        "actual": actual,
        "absolute_error": error,
        "tolerance": tolerance,
        "status": "pass" if error <= tolerance else "fail",
    }


def relation_case(
    case_id: str,
    kernel_a: np.ndarray,
    kernel_b: np.ndarray,
    expected_c_a: float,
    expected_c_b: float,
    expected_gamma: float,
    expected_response: float,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> list[dict]:
    state_a = modulate_kernel(kernel_a, ORIGIN_A, u, eta, p0)
    state_b = modulate_kernel(kernel_b, ORIGIN_B, u, eta, p0)
    relation = relation_readout(state_a, state_b, u, eta, p0)
    angle = effective_angle(
        "relational_C1",
        0.55,
        1.0,
        relation.parity_a.c_pi,
        relation.parity_b.c_pi,
        relation.gamma_ab,
    )
    return [
        test_row(
            case_id,
            "c_A",
            relation.parity_a.c_pi,
            expected_c_a,
        ),
        test_row(
            case_id,
            "c_B",
            relation.parity_b.c_pi,
            expected_c_b,
        ),
        test_row(
            case_id,
            "Gamma_AB",
            relation.gamma_ab,
            expected_gamma,
        ),
        test_row(
            case_id,
            "candidate_response",
            angle.candidate_response,
            expected_response,
        ),
    ]


def run_all_tests() -> list[dict]:
    params = Params()
    u, eta = make_grids(params)
    f1 = normalized_kernel(np.cos(u))
    f3 = normalized_kernel(np.cos(3.0 * u))
    b2 = normalized_kernel(np.cos(2.0 * u))
    rows = []
    rows.extend(
        relation_case(
            "A1_identical_pure_odd",
            f1,
            f1,
            -1.0,
            -1.0,
            1.0,
            1.0,
            u,
            eta,
            params.p0,
        )
    )
    rows.extend(
        relation_case(
            "A2_orthogonal_pure_odd",
            f1,
            f3,
            -1.0,
            -1.0,
            0.0,
            0.0,
            u,
            eta,
            params.p0,
        )
    )
    rows.extend(
        relation_case(
            "A3_identical_pure_even",
            b2,
            b2,
            1.0,
            1.0,
            1.0,
            -1.0,
            u,
            eta,
            params.p0,
        )
    )
    rows.extend(
        relation_case(
            "A4_pure_odd_times_even",
            f1,
            b2,
            -1.0,
            1.0,
            0.0,
            0.0,
            u,
            eta,
            params.p0,
        )
    )
    for phase in (0.0, math.pi / 2.0, math.pi):
        rows.extend(
            relation_case(
                f"A5_phase_{phase:.12g}",
                f1,
                np.exp(1j * phase) * f1,
                -1.0,
                -1.0,
                1.0,
                1.0,
                u,
                eta,
                params.p0,
            )
        )

    state_a = modulate_kernel(f1, ORIGIN_A, u, eta, params.p0)
    state_b = modulate_kernel(f3, ORIGIN_B, u, eta, params.p0)
    collision = scatter_once(
        state_a,
        state_b,
        scattering_mode="relational_C1",
        reflection_baseline=0.55,
        kappa=1.0,
        u=u,
        eta=eta,
        p0=params.p0,
    )
    rows.extend(
        [
            test_row(
                "unitarity",
                "unitarity_residual",
                collision.angle.unitarity_residual,
                0.0,
            ),
            test_row(
                "unitarity",
                "orthogonality_residual",
                collision.angle.orthogonality_residual,
                0.0,
            ),
            test_row(
                "unitarity",
                "path_sum_residual_A",
                collision.path_sum_residual_a,
                0.0,
            ),
            test_row(
                "unitarity",
                "path_sum_residual_B",
                collision.path_sum_residual_b,
                0.0,
            ),
            test_row(
                "unitarity",
                "total_norm_residual",
                collision.total_norm_residual,
                0.0,
            ),
        ]
    )

    cancelling = modulate_kernel(
        f1, ORIGIN_A, u, eta, params.p0
    ) - modulate_kernel(f1, ORIGIN_B, u, eta, params.p0)
    zero_norm_rejected = False
    try:
        relation_readout(cancelling, state_b, u, eta, params.p0)
    except RelationUndefinedError:
        zero_norm_rejected = True
    rows.append(
        test_row(
            "zero_relation_norm",
            "undefined_relation_rejected",
            float(zero_norm_rejected),
            1.0,
            0.0,
        )
    )

    angle_c0 = effective_angle("C0", 0.55, 1.0, -1.0, -1.0, 0.37)
    angle_reversed = effective_angle(
        "reversed_C1", 0.55, 1.0, -1.0, -1.0, 0.37
    )
    angle_relational = effective_angle(
        "relational_C1", 0.55, 1.0, -1.0, -1.0, 0.37
    )
    rows.extend(
        [
            test_row(
                "mode_preservation",
                "C0_response",
                angle_c0.candidate_response,
                0.0,
            ),
            test_row(
                "mode_preservation",
                "reversed_C1_response",
                angle_reversed.candidate_response,
                1.0,
            ),
            test_row(
                "mode_preservation",
                "relational_C1_response",
                angle_relational.candidate_response,
                0.37,
            ),
        ]
    )
    if any(row["status"] != "pass" for row in rows):
        failures = [row for row in rows if row["status"] != "pass"]
        raise AssertionError(f"Stage G unit test failure: {failures}")
    return rows


if __name__ == "__main__":
    print(f"{len(run_all_tests())} Stage G checks passed")
