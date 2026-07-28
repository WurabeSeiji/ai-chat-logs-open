"""Deterministic verification suite for the independent Stage D implementation."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from candidate_responses import (
    ProductNormTooSmallError,
    candidate_response,
)
from parity_metrics import (
    ModulationSpec,
    demodulate_state,
    modulate_kernel,
    norm2,
    parity_measurement,
)
from scattering_api import scatter_wave_pair
from wave_generators import make_eta_grid, make_u_grid, wave_library


STAGE_D_ROOT = Path(__file__).resolve().parents[1]
TEST_OUTPUT = STAGE_D_ROOT / "data" / "stage_D_test_results.json"
CRITICAL_R = 0.6971778791282474
TOLERANCE = 5.0e-13
SPEC_A = ModulationSpec("System_A_channel_A", q=1.0, eta_mode=1)
SPEC_B = ModulationSpec("System_A_channel_B", q=-1.0, eta_mode=2)


def _assert_below(name: str, value: float, tolerance: float = TOLERANCE) -> None:
    if not value <= tolerance:
        raise AssertionError(f"{name}: {value} > {tolerance}")


def _assert_close(
    name: str,
    actual: float,
    expected: float,
    tolerance: float = TOLERANCE,
) -> None:
    _assert_below(name, abs(actual - expected), tolerance)


def run_tests() -> dict:
    u = make_u_grid(512)
    eta = make_eta_grid(16)
    maximums = {
        "demodulation_residual": 0.0,
        "roundtrip_reconstruction_residual": 0.0,
        "pure_parity_error": 0.0,
        "candidate_pure_response_error": 0.0,
        "unitarity_residual": 0.0,
        "orthogonality_residual": 0.0,
        "path_sum_residual": 0.0,
        "pair_norm_conservation_residual": 0.0,
        "half_shift_equivariance_residual": 0.0,
        "exchange_output_residual": 0.0,
        "exchange_path_residual": 0.0,
        "endpoint_baseline_error": 0.0,
        "kappa_zero_baseline_error": 0.0,
    }

    libraries = {K: wave_library(u, K) for K in (1, 4, 8, 16)}
    for library in libraries.values():
        for wave in library.values():
            expected_parity = {"F": -1.0, "B": 1.0}.get(wave.label, 0.0)
            parity = parity_measurement(wave.kernel)
            maximums["pure_parity_error"] = max(
                maximums["pure_parity_error"],
                abs(parity.indicator - expected_parity),
            )
            for spec in (SPEC_A, SPEC_B):
                state = modulate_kernel(wave.kernel, u, eta, spec)
                demodulated = demodulate_state(
                    state, wave.kernel, u, eta, spec
                )
                maximums["demodulation_residual"] = max(
                    maximums["demodulation_residual"],
                    demodulated.reference_residual,
                )
                maximums["roundtrip_reconstruction_residual"] = max(
                    maximums["roundtrip_reconstruction_residual"],
                    demodulated.roundtrip_residual,
                )

    library = libraries[4]
    pure_expected = {
        ("C1", "B", "B"): 1.0,
        ("C1", "F", "F"): -1.0,
        ("C1", "B", "F"): 0.0,
        ("C1", "F", "B"): 0.0,
        ("C3", "B", "B"): 1.0,
        ("C3", "F", "F"): 1.0,
        ("C3", "B", "F"): -1.0,
        ("C3", "F", "B"): -1.0,
    }
    for (candidate, label_a, label_b), expected in pure_expected.items():
        value = candidate_response(
            candidate,
            library[label_a].kernel,
            library[label_b].kernel,
        ).value
        maximums["candidate_pure_response_error"] = max(
            maximums["candidate_pure_response_error"],
            abs(value - expected),
        )

    kernel_a = library["F"].kernel
    kernel_b = library["B"].kernel
    state_a = modulate_kernel(kernel_a, u, eta, SPEC_A)
    state_b = modulate_kernel(kernel_b, u, eta, SPEC_B)

    baseline = scatter_wave_pair(
        state_a,
        state_b,
        expected_kernel_a=kernel_a,
        expected_kernel_b=kernel_b,
        spec_a=SPEC_A,
        spec_b=SPEC_B,
        u=u,
        eta=eta,
        reflection_parameter=CRITICAL_R,
        kappa=0.0,
        candidate="C0",
    )
    _assert_close(
        "Candidate 0 reflection probability",
        baseline.reflection_probability,
        CRITICAL_R,
    )
    _assert_close(
        "Candidate 0 transmission probability",
        baseline.transmission_probability,
        1.0 - CRITICAL_R,
    )

    for candidate in ("C0", "C1", "C3"):
        result = scatter_wave_pair(
            state_a,
            state_b,
            expected_kernel_a=kernel_a,
            expected_kernel_b=kernel_b,
            spec_a=SPEC_A,
            spec_b=SPEC_B,
            u=u,
            eta=eta,
            reflection_parameter=CRITICAL_R,
            kappa=1.0,
            candidate=candidate,
        )
        swapped = scatter_wave_pair(
            state_b,
            state_a,
            expected_kernel_a=kernel_b,
            expected_kernel_b=kernel_a,
            spec_a=SPEC_B,
            spec_b=SPEC_A,
            u=u,
            eta=eta,
            reflection_parameter=CRITICAL_R,
            kappa=1.0,
            candidate=candidate,
        )
        exchange_denominator = math.sqrt(
            norm2(result.raw_output_a) + norm2(result.raw_output_b)
        )
        exchange_output = math.sqrt(
            norm2(result.raw_output_a - swapped.raw_output_b)
            + norm2(result.raw_output_b - swapped.raw_output_a)
        ) / exchange_denominator
        exchange_path = max(
            math.sqrt(
                norm2(
                    result.path_a_to_a_amplitude
                    - swapped.path_b_to_b_amplitude
                )
            ),
            math.sqrt(
                norm2(
                    result.path_b_to_a_amplitude
                    - swapped.path_a_to_b_amplitude
                )
            ),
            math.sqrt(
                norm2(
                    result.path_b_to_b_amplitude
                    - swapped.path_a_to_a_amplitude
                )
            ),
            math.sqrt(
                norm2(
                    result.path_a_to_b_amplitude
                    - swapped.path_b_to_a_amplitude
                )
            ),
        )
        maximums["exchange_output_residual"] = max(
            maximums["exchange_output_residual"], exchange_output
        )
        maximums["exchange_path_residual"] = max(
            maximums["exchange_path_residual"], exchange_path
        )
        maximums["unitarity_residual"] = max(
            maximums["unitarity_residual"], result.unitarity_residual
        )
        maximums["orthogonality_residual"] = max(
            maximums["orthogonality_residual"],
            result.orthogonality_residual,
        )
        maximums["path_sum_residual"] = max(
            maximums["path_sum_residual"],
            result.path_sum_residual_a,
            result.path_sum_residual_b,
        )
        maximums["pair_norm_conservation_residual"] = max(
            maximums["pair_norm_conservation_residual"],
            result.pair_norm_conservation_residual,
        )
        maximums["half_shift_equivariance_residual"] = max(
            maximums["half_shift_equivariance_residual"],
            result.half_shift_equivariance_residual,
        )

        kappa_zero = scatter_wave_pair(
            state_a,
            state_b,
            expected_kernel_a=kernel_a,
            expected_kernel_b=kernel_b,
            spec_a=SPEC_A,
            spec_b=SPEC_B,
            u=u,
            eta=eta,
            reflection_parameter=CRITICAL_R,
            kappa=0.0,
            candidate=candidate,
        )
        maximums["kappa_zero_baseline_error"] = max(
            maximums["kappa_zero_baseline_error"],
            abs(
                kappa_zero.reflection_probability
                - baseline.reflection_probability
            ),
            abs(
                kappa_zero.transmission_probability
                - baseline.transmission_probability
            ),
        )

        for endpoint in (0.0, 1.0):
            endpoint_result = scatter_wave_pair(
                state_a,
                state_b,
                expected_kernel_a=kernel_a,
                expected_kernel_b=kernel_b,
                spec_a=SPEC_A,
                spec_b=SPEC_B,
                u=u,
                eta=eta,
                reflection_parameter=endpoint,
                kappa=1.0,
                candidate=candidate,
            )
            maximums["endpoint_baseline_error"] = max(
                maximums["endpoint_baseline_error"],
                abs(endpoint_result.reflection_probability - endpoint),
                abs(
                    endpoint_result.transmission_probability
                    - (1.0 - endpoint)
                ),
                abs(endpoint_result.rho),
            )

    threshold_was_enforced = False
    try:
        candidate_response(
            "C3",
            np.zeros_like(kernel_a),
            np.zeros_like(kernel_b),
        )
    except ProductNormTooSmallError:
        threshold_was_enforced = True
    if not threshold_was_enforced:
        raise AssertionError("Candidate 3 zero-product threshold was not enforced")

    for name, value in maximums.items():
        _assert_below(name, value)

    result = {
        "schema": "stage_D_candidate_tests_v1",
        "status": "PASS",
        "grid": {"u": 512, "eta": 16},
        "tolerance": TOLERANCE,
        "maximums": maximums,
        "candidate3_zero_product_threshold_enforced": threshold_was_enforced,
        "tests": [
            "modulation/demodulation roundtrip",
            "pure parity eigenvalues",
            "Candidate 1 and 3 pure-state response values",
            "Candidate 0 reflection/transmission baseline",
            "unitarity and coefficient orthogonality",
            "path-sum and pair-norm conservation",
            "literal input-swap symmetry",
            "lifted half-period equivariance",
            "endpoint envelope return",
            "kappa=0 baseline return",
            "Candidate 3 fixed zero-product threshold",
        ],
    }
    TEST_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    TEST_OUTPUT.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    outcome = run_tests()
    print(json.dumps(outcome, indent=2, ensure_ascii=False))
