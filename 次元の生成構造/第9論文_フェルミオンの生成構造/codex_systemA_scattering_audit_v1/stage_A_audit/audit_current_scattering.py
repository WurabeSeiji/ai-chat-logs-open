#!/usr/bin/env python3
"""Independent numerical audit of the current System A scattering formula.

No original program is imported or executed.  The formulas identified by the
static audit are transcribed here and evaluated on controlled vectors.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


GRID_N = 512
TOL = 5.0e-13


def normalize(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError("zero norm state")
    return vector / norm


def coefficients_from_reflection_rate(reflection_rate: float) -> tuple[complex, complex]:
    if reflection_rate < 0.0 or reflection_rate > 1.0:
        raise ValueError("reflection rate must be in [0, 1]")
    delta_f = 2.0 * math.asin(math.sqrt(reflection_rate))
    phase = np.exp(0.5j * delta_f)
    t = phase * math.cos(0.5 * delta_f)
    r = -1j * phase * math.sin(0.5 * delta_f)
    return complex(t), complex(r)


def norm2(vector: np.ndarray) -> float:
    return float(np.vdot(vector, vector).real)


def half_shift(vector: np.ndarray) -> np.ndarray:
    if vector.size % 2:
        raise ValueError("half-period roll requires an even grid")
    return np.roll(vector, vector.size // 2)


def parity_indicator(vector: np.ndarray) -> float:
    return float(np.vdot(vector, half_shift(vector)).real / norm2(vector))


def parity_weights(vector: np.ndarray) -> tuple[float, float]:
    shifted = half_shift(vector)
    even = 0.5 * (vector + shifted)
    odd = 0.5 * (vector - shifted)
    total = norm2(vector)
    return norm2(even) / total, norm2(odd) / total


def scatter_raw(a: np.ndarray, b: np.ndarray, reflection_rate: float) -> tuple[np.ndarray, np.ndarray]:
    t, r = coefficients_from_reflection_rate(reflection_rate)
    return r * a + t * b, t * a + r * b


def scatter_current(a: np.ndarray, b: np.ndarray, reflection_rate: float) -> tuple[np.ndarray, np.ndarray]:
    a_raw, b_raw = scatter_raw(a, b, reflection_rate)
    return normalize(a_raw), normalize(b_raw)


def pair_norm(a: np.ndarray, b: np.ndarray) -> float:
    return norm2(a) + norm2(b)


def commutator_response(
    a: np.ndarray,
    b: np.ndarray,
    reflection_rate: float,
    normalize_channels: bool,
) -> float:
    scatter = scatter_current if normalize_channels else scatter_raw
    lhs_a, lhs_b = scatter(half_shift(a), half_shift(b), reflection_rate)
    rhs_a, rhs_b = scatter(a, b, reflection_rate)
    delta = pair_norm(lhs_a - half_shift(rhs_a), lhs_b - half_shift(rhs_b))
    denominator = pair_norm(a, b)
    return math.sqrt(max(delta, 0.0) / denominator)


def complex_pair(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def case_record(
    name: str,
    a: np.ndarray,
    b: np.ndarray,
    reflection_rate: float,
) -> dict[str, Any]:
    t, r = coefficients_from_reflection_rate(reflection_rate)
    a_raw, b_raw = scatter_raw(a, b, reflection_rate)
    a_out, b_out = normalize(a_raw), normalize(b_raw)
    a_in_even, a_in_odd = parity_weights(a)
    b_in_even, b_in_odd = parity_weights(b)
    a_raw_even, a_raw_odd = parity_weights(a_raw)
    b_raw_even, b_raw_odd = parity_weights(b_raw)
    return {
        "case": name,
        "R_input": reflection_rate,
        "T": abs(t) ** 2,
        "R": abs(r) ** 2,
        "input_inner_product": complex_pair(complex(np.vdot(a, b))),
        "norm_before_a": norm2(a),
        "norm_before_b": norm2(b),
        "norm_before_combined": pair_norm(a, b),
        "norm_after_raw_a": norm2(a_raw),
        "norm_after_raw_b": norm2(b_raw),
        "norm_after_raw_combined": pair_norm(a_raw, b_raw),
        "norm_after_final_a": norm2(a_out),
        "norm_after_final_b": norm2(b_out),
        "norm_after_final_combined": pair_norm(a_out, b_out),
        "normalization_scale_a": 1.0 / math.sqrt(norm2(a_raw)),
        "normalization_scale_b": 1.0 / math.sqrt(norm2(b_raw)),
        "raw_combined_norm_error": pair_norm(a_raw, b_raw) - pair_norm(a, b),
        "a_parity_in_even": a_in_even,
        "a_parity_in_odd": a_in_odd,
        "b_parity_in_even": b_in_even,
        "b_parity_in_odd": b_in_odd,
        "a_parity_raw_even": a_raw_even,
        "a_parity_raw_odd": a_raw_odd,
        "b_parity_raw_even": b_raw_even,
        "b_parity_raw_odd": b_raw_odd,
        "commutator_raw": commutator_response(a, b, reflection_rate, False),
        "commutator_after_channel_normalization": commutator_response(a, b, reflection_rate, True),
    }


def main() -> int:
    u = np.linspace(-math.pi, math.pi, GRID_N, endpoint=False)
    even_2 = normalize(np.cos(2.0 * u).astype(complex))
    even_4 = normalize(np.cos(4.0 * u).astype(complex))
    odd_1 = normalize(np.cos(u).astype(complex))
    odd_3 = normalize(np.cos(3.0 * u).astype(complex))
    mixed_a = normalize(odd_1 + 0.6 * even_2)
    mixed_b = normalize(0.4j * odd_3 + even_4)
    relative_phase_a = odd_1
    relative_phase_b = 1j * odd_1

    test_r = 0.6971778791282474
    cases = [
        case_record("pure_even_pair", even_2, even_4, test_r),
        case_record("pure_odd_pair", odd_1, odd_3, test_r),
        case_record("even_odd_pair", even_2, odd_1, test_r),
        case_record("mixed_pair", mixed_a, mixed_b, test_r),
        case_record("identical_shape_relative_phase_i", relative_phase_a, relative_phase_b, test_r),
    ]

    coefficient_rows = []
    for reflection_rate in (0.0, 0.25, 0.5, test_r, 1.0):
        t, r = coefficients_from_reflection_rate(reflection_rate)
        matrix = np.asarray([[r, t], [t, r]], dtype=complex)
        unitarity_error = float(np.linalg.norm(matrix.conj().T @ matrix - np.eye(2)))
        coefficient_rows.append(
            {
                "R_input": reflection_rate,
                "t": complex_pair(t),
                "r": complex_pair(r),
                "T": abs(t) ** 2,
                "R": abs(r) ** 2,
                "T_plus_R_error": abs(t) ** 2 + abs(r) ** 2 - 1.0,
                "unitarity_error": unitarity_error,
            }
        )

    kernel_odd = normalize((np.cos(u) + np.cos(3.0 * u)).astype(complex))
    carrier = np.exp(1j * u)
    full_state = normalize(kernel_odd * carrier)
    inverse_shifted = normalize(full_state * np.exp(-1j * u))
    carrier_audit = {
        "kernel_odd_c_pi": parity_indicator(kernel_odd),
        "full_state_after_q_plus_1_carrier_c_pi": parity_indicator(full_state),
        "after_inverse_carrier_c_pi": parity_indicator(inverse_shifted),
        "interpretation": "the q=+1 carrier flips the half-period eigenvalue; kernel parity requires inverse carrier removal",
    }

    assertions = {
        "all_scattering_matrices_unitary": all(row["unitarity_error"] <= TOL for row in coefficient_rows),
        "all_raw_combined_norms_conserved": all(abs(row["raw_combined_norm_error"]) <= TOL for row in cases),
        "raw_scattering_commutes_with_half_shift": all(row["commutator_raw"] <= TOL for row in cases),
        "normalized_scattering_commutes_with_half_shift": all(
            row["commutator_after_channel_normalization"] <= TOL for row in cases
        ),
        "normalization_redundant_for_orthogonal_cases": all(
            abs(row["normalization_scale_a"] - 1.0) <= TOL
            and abs(row["normalization_scale_b"] - 1.0) <= TOL
            for row in cases[:4]
        ),
        "normalization_nontrivial_for_relative_phase_overlap": (
            abs(cases[-1]["normalization_scale_a"] - 1.0) > 1.0e-6
            or abs(cases[-1]["normalization_scale_b"] - 1.0) > 1.0e-6
        ),
        "carrier_flips_odd_kernel_to_even_full_state": (
            abs(carrier_audit["kernel_odd_c_pi"] + 1.0) <= TOL
            and abs(carrier_audit["full_state_after_q_plus_1_carrier_c_pi"] - 1.0) <= TOL
            and abs(carrier_audit["after_inverse_carrier_c_pi"] + 1.0) <= TOL
        ),
    }

    audit_root = Path(__file__).resolve().parents[1]
    log_dir = audit_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    json_path = log_dir / "current_scattering_diagnostic.json"
    csv_path = log_dir / "current_scattering_diagnostic.csv"
    payload = {
        "schema": "codex_systemA_current_scattering_diagnostic_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "independent transcription; no original module import or execution",
        "grid_n": GRID_N,
        "tolerance": TOL,
        "coefficient_rows": coefficient_rows,
        "cases": cases,
        "carrier_audit": carrier_audit,
        "assertions": assertions,
    }
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    csv_fields = [
        key
        for key, value in cases[0].items()
        if not isinstance(value, (list, dict))
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_fields)
        writer.writeheader()
        for row in cases:
            writer.writerow({key: row[key] for key in csv_fields})

    if not all(assertions.values()):
        raise SystemExit(f"diagnostic assertion failed: {assertions}")
    print(json_path)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
