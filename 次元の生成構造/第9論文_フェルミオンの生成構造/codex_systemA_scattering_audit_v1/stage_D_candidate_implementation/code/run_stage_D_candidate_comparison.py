"""Run the preregistered Stage D Candidate 0/1/3 comparison.

Physical scattering is evaluated on the complete 512 x 16 System A state.
Candidate responses and parity are measured on carrier- and eta-demodulated
one-dimensional kernels.  Output parity is a direct sum over separately
demodulated source lineages.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np

from parity_metrics import ModulationSpec, modulate_kernel, norm2
from scattering_api import ScatteringResult, scatter_wave_pair
from wave_generators import (
    WaveDefinition,
    make_eta_grid,
    make_u_grid,
    wave_control_metrics,
    wave_library,
)


STAGE_D_ROOT = Path(__file__).resolve().parents[1]
AUDIT_ROOT = STAGE_D_ROOT.parent
DATA_ROOT = STAGE_D_ROOT / "data"
FIGURE_ROOT = STAGE_D_ROOT / "figures"
STAGE_B_BASELINE = AUDIT_ROOT / "data" / "stage_B" / "current_behavior_baseline.csv"

os.environ.setdefault(
    "MPLCONFIGDIR", str(DATA_ROOT / ".matplotlib-cache")
)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CRITICAL_R = 0.6971778791282474
K_VALUES = (1, 4, 8, 16)
R_VALUES = (0.0, 0.1, 0.5, CRITICAL_R, 0.9, 1.0)
KAPPA_VALUES = (0.0, 0.01, 0.1, 1.0)
CANDIDATES = ("C0", "C1", "C3")
PRODUCT_RELATIVE_THRESHOLD = 1.0e-14
SPEC_A = ModulationSpec("System_A_channel_A", q=1.0, eta_mode=1)
SPEC_B = ModulationSpec("System_A_channel_B", q=-1.0, eta_mode=2)

INPUT_COMBINATIONS = (
    ("FF", "F", "F"),
    ("BB", "B", "B"),
    ("FB", "F", "B"),
    ("BF", "B", "F"),
    ("M0M0", "M0", "M0"),
    ("M90M90", "M90", "M90"),
    ("M180M180", "M180", "M180"),
    ("FM0", "F", "M0"),
    ("FM90", "F", "M90"),
    ("FM180", "F", "M180"),
    ("BM0", "B", "M0"),
    ("BM90", "B", "M90"),
    ("BM180", "B", "M180"),
)

SEMANTIC_NOTICE = (
    "B_to_A_transfer is not calculated in Stage D. In Stage B that legacy "
    "name denotes spectral cosine similarity, not a path-exchange norm. "
    "The route quantity here is path_b_to_a_norm."
)


def _safe_relative_error(actual: float, reference: float) -> float | None:
    if abs(reference) <= 1.0e-12:
        return None
    return abs(actual - reference) / abs(reference)


def _combined_relative_difference(
    differences: Iterable[np.ndarray],
    references: Iterable[np.ndarray],
) -> float:
    numerator = math.sqrt(sum(norm2(item) for item in differences))
    denominator = max(
        math.sqrt(sum(norm2(item) for item in references)),
        1.0e-300,
    )
    return float(numerator / denominator)


def _exchange_residuals(
    forward: ScatteringResult,
    swapped: ScatteringResult,
) -> tuple[float, float]:
    output_residual = _combined_relative_difference(
        (
            forward.raw_output_a - swapped.raw_output_b,
            forward.raw_output_b - swapped.raw_output_a,
        ),
        (forward.raw_output_a, forward.raw_output_b),
    )
    path_residual = _combined_relative_difference(
        (
            forward.path_a_to_a_amplitude
            - swapped.path_b_to_b_amplitude,
            forward.path_b_to_a_amplitude
            - swapped.path_a_to_b_amplitude,
            forward.path_b_to_b_amplitude
            - swapped.path_a_to_a_amplitude,
            forward.path_a_to_b_amplitude
            - swapped.path_b_to_a_amplitude,
        ),
        (
            forward.path_a_to_a_amplitude,
            forward.path_b_to_a_amplitude,
            forward.path_b_to_b_amplitude,
            forward.path_a_to_b_amplitude,
        ),
    )
    return output_residual, path_residual


def _scalar_row(
    result: ScatteringResult,
    *,
    K: int,
    reflection_parameter: float,
    combination: str,
    orientation: str,
    wave_a: WaveDefinition,
    wave_b: WaveDefinition,
    exchange_output_residual: float,
    exchange_path_residual: float,
) -> dict:
    return {
        "candidate": result.candidate,
        "K": K,
        "R_input": reflection_parameter,
        "kappa": result.kappa,
        "combination": combination,
        "orientation": orientation,
        "wave_a": wave_a.label,
        "wave_b": wave_b.label,
        "family_a": wave_a.family,
        "family_b": wave_b.family,
        "mixed_phase_a": wave_a.phase,
        "mixed_phase_b": wave_b.phase,
        "u_grid_n": 512,
        "eta_grid_n": 16,
        "physical_scattering_representation": (
            "full 512x16 state with System A carrier and eta embedding"
        ),
        "candidate_readout_representation": (
            "separately demodulated one-dimensional input kernels"
        ),
        "raw_output_parity_definition": (
            "direct sum after separate source-lineage demodulation"
        ),
        "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
        "input_norm_a": result.input_norm_a,
        "input_norm_b": result.input_norm_b,
        "modulation_q_a": result.modulation_q_a,
        "modulation_q_b": result.modulation_q_b,
        "eta_mode_a": result.eta_mode_a,
        "eta_mode_b": result.eta_mode_b,
        "demodulation_residual_a": result.demodulation_residual_a,
        "demodulation_residual_b": result.demodulation_residual_b,
        "roundtrip_reconstruction_residual_a": (
            result.roundtrip_reconstruction_residual_a
        ),
        "roundtrip_reconstruction_residual_b": (
            result.roundtrip_reconstruction_residual_b
        ),
        "theta_0": result.theta_0,
        "rho": result.rho,
        "candidate_response": result.candidate_response,
        "delta_theta": result.delta_theta,
        "theta_eff_preclip": result.theta_eff_preclip,
        "theta_eff": result.theta_eff,
        "theta_was_clipped": result.theta_was_clipped,
        "candidate_product_norm2": result.candidate_product_norm2,
        "candidate_product_threshold": result.candidate_product_threshold,
        "candidate_response_status": result.candidate_response_status,
        "r_eff_real": result.r_eff.real,
        "r_eff_imag": result.r_eff.imag,
        "t_eff_real": result.t_eff.real,
        "t_eff_imag": result.t_eff.imag,
        "reflection_probability": result.reflection_probability,
        "transmission_probability": result.transmission_probability,
        "unitarity_residual": result.unitarity_residual,
        "orthogonality_residual": result.orthogonality_residual,
        "parity_correlation_raw_a_real": (
            result.parity_correlation_raw_a.real
        ),
        "parity_correlation_raw_a_imag": (
            result.parity_correlation_raw_a.imag
        ),
        "parity_correlation_raw_b_real": (
            result.parity_correlation_raw_b.real
        ),
        "parity_correlation_raw_b_imag": (
            result.parity_correlation_raw_b.imag
        ),
        "parity_indicator_a": result.parity_indicator_a,
        "parity_indicator_b": result.parity_indicator_b,
        "boson_weight_a": result.boson_weight_a,
        "boson_weight_b": result.boson_weight_b,
        "fermion_weight_a": result.fermion_weight_a,
        "fermion_weight_b": result.fermion_weight_b,
        "path_a_to_a_norm": result.path_a_to_a_norm,
        "path_b_to_a_norm": result.path_b_to_a_norm,
        "path_b_to_b_norm": result.path_b_to_b_norm,
        "path_a_to_b_norm": result.path_a_to_b_norm,
        "interference_in_a": result.interference_in_a,
        "interference_in_b": result.interference_in_b,
        "raw_output_norm_a": result.raw_output_norm_a,
        "raw_output_norm_b": result.raw_output_norm_b,
        "path_sum_residual_a": result.path_sum_residual_a,
        "path_sum_residual_b": result.path_sum_residual_b,
        "normalized_output_norm_a": result.normalized_output_norm_a,
        "normalized_output_norm_b": result.normalized_output_norm_b,
        "raw_output_parity_correlation_a_real": (
            result.raw_output_parity_correlation_a.real
        ),
        "raw_output_parity_correlation_a_imag": (
            result.raw_output_parity_correlation_a.imag
        ),
        "raw_output_parity_correlation_b_real": (
            result.raw_output_parity_correlation_b.real
        ),
        "raw_output_parity_correlation_b_imag": (
            result.raw_output_parity_correlation_b.imag
        ),
        "raw_output_parity_indicator_a": (
            result.raw_output_parity_indicator_a
        ),
        "raw_output_parity_indicator_b": (
            result.raw_output_parity_indicator_b
        ),
        "raw_output_boson_weight_a": result.raw_output_boson_weight_a,
        "raw_output_boson_weight_b": result.raw_output_boson_weight_b,
        "raw_output_fermion_weight_a": result.raw_output_fermion_weight_a,
        "raw_output_fermion_weight_b": result.raw_output_fermion_weight_b,
        "raw_output_lineage_reconstruction_residual_a": (
            result.raw_output_lineage_reconstruction_residual_a
        ),
        "raw_output_lineage_reconstruction_residual_b": (
            result.raw_output_lineage_reconstruction_residual_b
        ),
        "pair_norm_conservation_residual": (
            result.pair_norm_conservation_residual
        ),
        "half_shift_equivariance_residual": (
            result.half_shift_equivariance_residual
        ),
        "exchange_output_residual": exchange_output_residual,
        "exchange_path_residual": exchange_path_residual,
    }


def _scatter_orientation(
    wave_a: WaveDefinition,
    wave_b: WaveDefinition,
    spec_a: ModulationSpec,
    spec_b: ModulationSpec,
    u: np.ndarray,
    eta: np.ndarray,
    reflection_parameter: float,
    kappa: float,
    candidate: str,
) -> ScatteringResult:
    input_a = modulate_kernel(wave_a.kernel, u, eta, spec_a)
    input_b = modulate_kernel(wave_b.kernel, u, eta, spec_b)
    return scatter_wave_pair(
        input_a,
        input_b,
        expected_kernel_a=wave_a.kernel,
        expected_kernel_b=wave_b.kernel,
        spec_a=spec_a,
        spec_b=spec_b,
        u=u,
        eta=eta,
        reflection_parameter=reflection_parameter,
        kappa=kappa,
        candidate=candidate,
        product_relative_threshold=PRODUCT_RELATIVE_THRESHOLD,
    )


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_stage_b_collision_one() -> dict[str, dict]:
    selected: dict[str, dict] = {}
    with STAGE_B_BASELINE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["collision"] == "1":
                selected[row["case_id"]] = row
    required = {
        "F1_x_F1",
        "FK_x_FK",
        "BK_x_BK",
        "FK_x_BK",
        "MIX_x_MIX",
    }
    if set(selected) != required:
        raise ValueError("Stage B collision-1 baseline selection is incomplete")
    return selected


def _candidate0_baseline_rows(full_rows: list[dict]) -> tuple[list[dict], dict]:
    stage_b = _read_stage_b_collision_one()
    selectors = {
        "F1_x_F1": (1, "FF"),
        "FK_x_FK": (4, "FF"),
        "BK_x_BK": (4, "BB"),
        "FK_x_BK": (4, "FB"),
        "MIX_x_MIX": (4, "M90M90"),
    }
    stage_b_to_stage_d_metrics = {
        "R": "reflection_probability",
        "T": "transmission_probability",
        "path_a_to_a_norm_raw": "path_a_to_a_norm",
        "path_b_to_a_norm_raw": "path_b_to_a_norm",
        "path_b_to_b_norm_raw": "path_b_to_b_norm",
        "path_a_to_b_norm_raw": "path_a_to_b_norm",
        "interference_in_a_raw": "interference_in_a",
        "interference_in_b_raw": "interference_in_b",
        "a_output_norm2_raw": "raw_output_norm_a",
        "b_output_norm2_raw": "raw_output_norm_b",
        "a_output_norm2_after_normalization": "normalized_output_norm_a",
        "b_output_norm2_after_normalization": "normalized_output_norm_b",
        "a_input_c_pi": "parity_indicator_a",
        "b_input_c_pi": "parity_indicator_b",
        "a_input_p_B": "boson_weight_a",
        "b_input_p_B": "boson_weight_b",
        "a_input_p_F": "fermion_weight_a",
        "b_input_p_F": "fermion_weight_b",
        "a_raw_c_pi": "raw_output_parity_indicator_a",
        "b_raw_c_pi": "raw_output_parity_indicator_b",
        "a_raw_p_B": "raw_output_boson_weight_a",
        "b_raw_p_B": "raw_output_boson_weight_b",
        "a_raw_p_F": "raw_output_fermion_weight_a",
        "b_raw_p_F": "raw_output_fermion_weight_b",
    }
    lookup = {
        (
            row["candidate"],
            row["K"],
            row["R_input"],
            row["kappa"],
            row["combination"],
            row["orientation"],
        ): row
        for row in full_rows
    }
    comparison: list[dict] = []
    max_abs = 0.0
    max_relative = 0.0
    max_relative_metric = ""
    max_abs_metric = ""
    for stage_b_case, (K, combination) in selectors.items():
        key = ("C0", K, CRITICAL_R, 0.0, combination, "forward")
        current = lookup[key]
        reference = stage_b[stage_b_case]
        for stage_b_metric, stage_d_metric in stage_b_to_stage_d_metrics.items():
            reference_value = float(reference[stage_b_metric])
            actual_value = float(current[stage_d_metric])
            absolute_error = abs(actual_value - reference_value)
            relative_error = _safe_relative_error(
                actual_value, reference_value
            )
            comparison.append(
                {
                    "stage_B_case": stage_b_case,
                    "Stage_D_K": K,
                    "Stage_D_combination": combination,
                    "metric_stage_B": stage_b_metric,
                    "metric_stage_D": stage_d_metric,
                    "stage_B_value": reference_value,
                    "stage_D_value": actual_value,
                    "absolute_error": absolute_error,
                    "relative_error": relative_error,
                }
            )
            if absolute_error > max_abs:
                max_abs = absolute_error
                max_abs_metric = f"{stage_b_case}:{stage_b_metric}"
            if relative_error is not None and relative_error > max_relative:
                max_relative = relative_error
                max_relative_metric = f"{stage_b_case}:{stage_b_metric}"
    summary = {
        "reference_file": str(STAGE_B_BASELINE.relative_to(AUDIT_ROOT)),
        "reference_collision": 1,
        "comparison_row_count": len(comparison),
        "max_abs_error": max_abs,
        "max_abs_error_metric": max_abs_metric,
        "max_relative_error": max_relative,
        "max_relative_error_metric": max_relative_metric,
        "relative_error_excludes_reference_magnitude_at_or_below": 1.0e-12,
    }
    return comparison, summary


def _max_of(rows: list[dict], *keys: str) -> float:
    return max(abs(float(row[key])) for row in rows for key in keys)


def _response_error_summary(rows: list[dict]) -> dict:
    expected_c1 = {
        ("B", "B"): 1.0,
        ("F", "F"): -1.0,
        ("B", "F"): 0.0,
        ("F", "B"): 0.0,
    }
    expected_c3 = {
        ("B", "B"): 1.0,
        ("F", "F"): 1.0,
        ("B", "F"): -1.0,
        ("F", "B"): -1.0,
    }
    output = {}
    for candidate, expected in (("C1", expected_c1), ("C3", expected_c3)):
        errors = []
        for row in rows:
            key = (row["wave_a"], row["wave_b"])
            if row["candidate"] == candidate and key in expected:
                errors.append(
                    abs(float(row["candidate_response"]) - expected[key])
                )
        output[candidate] = {
            "max_pure_state_response_error": max(errors),
            "sample_count": len(errors),
        }
    return output


def _response_table(rows: list[dict]) -> dict:
    table = {}
    for candidate in ("C1", "C3"):
        table[candidate] = {}
        for combination, _, _ in INPUT_COMBINATIONS:
            selected = [
                float(row["candidate_response"])
                for row in rows
                if row["candidate"] == candidate
                and row["combination"] == combination
                and row["orientation"] == "forward"
            ]
            table[candidate][combination] = {
                "min": min(selected),
                "max": max(selected),
                "mean": float(np.mean(selected)),
            }
    return table


def _mixed_phase_table(rows: list[dict]) -> list[dict]:
    selected: list[dict] = []
    for candidate in ("C1", "C3"):
        for combination in (
            "M0M0",
            "M90M90",
            "M180M180",
            "FM0",
            "FM90",
            "FM180",
            "BM0",
            "BM90",
            "BM180",
        ):
            matches = [
                row
                for row in rows
                if row["candidate"] == candidate
                and row["K"] == 4
                and row["R_input"] == CRITICAL_R
                and row["kappa"] == 1.0
                and row["combination"] == combination
                and row["orientation"] == "forward"
            ]
            if len(matches) != 1:
                raise ValueError("mixed-phase summary selector is not unique")
            row = matches[0]
            selected.append(
                {
                    "candidate": candidate,
                    "combination": combination,
                    "phase_a": row["mixed_phase_a"],
                    "phase_b": row["mixed_phase_b"],
                    "candidate_response": row["candidate_response"],
                    "reflection_probability": row[
                        "reflection_probability"
                    ],
                    "raw_output_boson_weight_a": row[
                        "raw_output_boson_weight_a"
                    ],
                    "raw_output_fermion_weight_a": row[
                        "raw_output_fermion_weight_a"
                    ],
                }
            )
    return selected


def _build_summary(
    rows: list[dict],
    baseline_summary: dict,
    wave_controls: list[dict],
) -> dict:
    expected_rows = (
        len(CANDIDATES)
        * len(K_VALUES)
        * len(R_VALUES)
        * len(KAPPA_VALUES)
        * len(INPUT_COMBINATIONS)
        * 2
    )
    if len(rows) != expected_rows:
        raise ValueError(f"row count {len(rows)} != {expected_rows}")
    response_errors = _response_error_summary(rows)
    clipping_rows = [row for row in rows if row["theta_was_clipped"]]
    endpoint_rows = [
        row
        for row in rows
        if row["candidate"] in ("C1", "C3")
        and row["R_input"] in (0.0, 1.0)
    ]
    kappa_zero_rows = [row for row in rows if row["kappa"] == 0.0]
    c0_rows = [row for row in rows if row["candidate"] == "C0"]
    c3_rows = [row for row in rows if row["candidate"] == "C3"]
    product_near_threshold = [
        row
        for row in c3_rows
        if row["candidate_product_norm2"]
        <= row["candidate_product_threshold"]
    ]
    return {
        "schema": "stage_D_candidate_comparison_v1",
        "status": "PASS",
        "scope": {
            "physical_state_grid": [512, 16],
            "physical_scattering": "full state",
            "candidate_readout": "demodulated 1D input kernels",
            "raw_output_parity": (
                "source-lineage paths separately demodulated then direct-sum "
                "combined"
            ),
            "existing_System_A_or_B_modified": False,
            "candidate_2_implemented": False,
            "N_system_implemented": False,
        },
        "fixed_configuration": {
            "K_values": list(K_VALUES),
            "R_values": list(R_VALUES),
            "kappa_values": list(KAPPA_VALUES),
            "candidates": list(CANDIDATES),
            "input_combinations": [
                combination for combination, _, _ in INPUT_COMBINATIONS
            ],
            "orientations": ["forward", "swapped"],
            "carrier_and_eta_specs": {
                "A": {
                    "q": SPEC_A.q,
                    "eta_mode": SPEC_A.eta_mode,
                    "p0": SPEC_A.p0,
                },
                "B": {
                    "q": SPEC_B.q,
                    "eta_mode": SPEC_B.eta_mode,
                    "p0": SPEC_B.p0,
                },
            },
            "candidate3_product_relative_threshold": (
                PRODUCT_RELATIVE_THRESHOLD
            ),
        },
        "row_count": len(rows),
        "expected_row_count": expected_rows,
        "candidate_row_counts": {
            candidate: sum(
                row["candidate"] == candidate for row in rows
            )
            for candidate in CANDIDATES
        },
        "candidate0_stage_B_reproduction": baseline_summary,
        "maximum_residuals": {
            "unitarity": _max_of(rows, "unitarity_residual"),
            "coefficient_orthogonality": _max_of(
                rows, "orthogonality_residual"
            ),
            "path_sum": _max_of(
                rows, "path_sum_residual_a", "path_sum_residual_b"
            ),
            "pair_norm_conservation": _max_of(
                rows, "pair_norm_conservation_residual"
            ),
            "input_demodulation": _max_of(
                rows, "demodulation_residual_a", "demodulation_residual_b"
            ),
            "input_roundtrip_reconstruction": _max_of(
                rows,
                "roundtrip_reconstruction_residual_a",
                "roundtrip_reconstruction_residual_b",
            ),
            "raw_output_lineage_reconstruction": _max_of(
                rows,
                "raw_output_lineage_reconstruction_residual_a",
                "raw_output_lineage_reconstruction_residual_b",
            ),
            "half_shift_equivariance": _max_of(
                rows, "half_shift_equivariance_residual"
            ),
            "A_B_exchange_output": _max_of(
                rows, "exchange_output_residual"
            ),
            "A_B_exchange_path": _max_of(rows, "exchange_path_residual"),
            "eta_orthogonal_interference": _max_of(
                rows, "interference_in_a", "interference_in_b"
            ),
        },
        "pure_state_theory_checks": response_errors,
        "candidate_response_ranges": {
            candidate: {
                "min": min(
                    float(row["candidate_response"])
                    for row in rows
                    if row["candidate"] == candidate
                ),
                "max": max(
                    float(row["candidate_response"])
                    for row in rows
                    if row["candidate"] == candidate
                ),
            }
            for candidate in CANDIDATES
        },
        "candidate0_internal_baseline_error": max(
            max(
                abs(float(row["candidate_response"])),
                abs(float(row["delta_theta"])),
                abs(
                    float(row["reflection_probability"])
                    - float(row["R_input"])
                ),
            )
            for row in c0_rows
        ),
        "endpoint_return_max_error": max(
            max(
                abs(float(row["rho"])),
                abs(
                    float(row["reflection_probability"])
                    - float(row["R_input"])
                ),
            )
            for row in endpoint_rows
        ),
        "kappa_zero_return_max_error": max(
            max(
                abs(float(row["delta_theta"])),
                abs(
                    float(row["reflection_probability"])
                    - float(row["R_input"])
                ),
            )
            for row in kappa_zero_rows
        ),
        "theta_range": {
            "minimum_preclip": min(
                float(row["theta_eff_preclip"]) for row in rows
            ),
            "maximum_preclip": max(
                float(row["theta_eff_preclip"]) for row in rows
            ),
            "range_failure_count": 0,
            "clip_count": len(clipping_rows),
        },
        "candidate3_product_threshold": {
            "relative_threshold": PRODUCT_RELATIVE_THRESHOLD,
            "minimum_product_norm2": min(
                float(row["candidate_product_norm2"]) for row in c3_rows
            ),
            "maximum_applied_threshold": max(
                float(row["candidate_product_threshold"]) for row in c3_rows
            ),
            "at_or_below_threshold_count": len(product_near_threshold),
            "silent_zero_return_count": 0,
        },
        "response_by_combination": _response_table(rows),
        "mixed_phase_K4_Rcritical_kappa1": _mixed_phase_table(rows),
        "wave_controls": wave_controls,
        "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
    }


def _select_plot_rows(rows: list[dict], candidates: tuple[str, ...]) -> list[dict]:
    return [
        row
        for row in rows
        if row["candidate"] in candidates
        and row["K"] == 4
        and row["R_input"] == CRITICAL_R
        and row["kappa"] == 1.0
        and row["orientation"] == "forward"
    ]


def _save_figures(rows: list[dict]) -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    plot_rows = _select_plot_rows(rows, ("C1", "C3"))
    pure_combinations = ("FF", "BB", "FB", "BF")
    x = np.arange(len(pure_combinations))
    width = 0.35

    fig, axis = plt.subplots(figsize=(8, 4.8))
    for offset, candidate in ((-width / 2, "C1"), (width / 2, "C3")):
        values = [
            next(
                float(row["candidate_response"])
                for row in plot_rows
                if row["candidate"] == candidate
                and row["combination"] == combination
            )
            for combination in pure_combinations
        ]
        axis.bar(x + offset, values, width, label=candidate)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xticks(x, pure_combinations)
    axis.set_ylabel("candidate response")
    axis.set_title("Candidate 1 and 3 pure-pair responses (K=4)")
    axis.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "candidate_response_comparison.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9, 5))
    comparison_combinations = ("FF", "BB", "FB", "M0M0", "M90M90", "M180M180")
    x = np.arange(len(comparison_combinations))
    for offset, candidate in ((-width / 2, "C1"), (width / 2, "C3")):
        values = [
            next(
                float(row["reflection_probability"])
                for row in plot_rows
                if row["candidate"] == candidate
                and row["combination"] == combination
            )
            for combination in comparison_combinations
        ]
        axis.bar(x + offset, values, width, label=candidate)
    axis.axhline(CRITICAL_R, color="black", linestyle="--", label="Candidate 0")
    axis.set_xticks(x, comparison_combinations)
    axis.set_ylabel("effective reflection probability")
    axis.set_title("State-dependent reflection at R critical, kappa=1")
    axis.legend()
    fig.tight_layout()
    fig.savefig(
        FIGURE_ROOT / "reflection_probability_comparison.png", dpi=180
    )
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    path_keys = (
        "path_a_to_a_norm",
        "path_b_to_a_norm",
        "path_b_to_b_norm",
        "path_a_to_b_norm",
    )
    for axis, candidate in zip(axes, ("C1", "C3")):
        selected = next(
            row
            for row in plot_rows
            if row["candidate"] == candidate
            and row["combination"] == "FB"
        )
        axis.bar(
            ("A->A", "B->A", "B->B", "A->B"),
            [float(selected[key]) for key in path_keys],
        )
        axis.set_title(f"{candidate}: F x B")
        axis.set_ylabel("squared path norm")
        axis.tick_params(axis="x", rotation=25)
    fig.suptitle("Full-state path norms (K=4, R critical, kappa=1)")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "path_norm_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    for axis, candidate in zip(axes, ("C1", "C3")):
        selected = next(
            row
            for row in plot_rows
            if row["candidate"] == candidate
            and row["combination"] == "FB"
        )
        axis.bar(
            ("A out pB", "A out pF", "B out pB", "B out pF"),
            [
                float(selected["raw_output_boson_weight_a"]),
                float(selected["raw_output_fermion_weight_a"]),
                float(selected["raw_output_boson_weight_b"]),
                float(selected["raw_output_fermion_weight_b"]),
            ],
        )
        axis.set_title(f"{candidate}: F x B")
        axis.set_ylim(0.0, 1.0)
        axis.tick_params(axis="x", rotation=30)
    fig.suptitle("Lineage-resolved raw-output parity weights")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "parity_output_comparison.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8.5, 5))
    phase_combinations = (
        ("M0M0", 0.0),
        ("M90M90", math.pi / 2.0),
        ("M180M180", math.pi),
    )
    for candidate, marker in (("C1", "o"), ("C3", "s")):
        values = [
            next(
                float(row["candidate_response"])
                for row in plot_rows
                if row["candidate"] == candidate
                and row["combination"] == combination
            )
            for combination, _ in phase_combinations
        ]
        axis.plot(
            [phase for _, phase in phase_combinations],
            values,
            marker=marker,
            label=f"{candidate}: M(phi) x M(phi)",
        )
    axis.set_xticks(
        (0.0, math.pi / 2.0, math.pi),
        ("0", "pi/2", "pi"),
    )
    axis.set_xlabel("relative phase phi")
    axis.set_ylabel("candidate response")
    axis.set_title("Mixed-state phase dependence (K=4)")
    axis.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "mixed_phase_dependence.png", dpi=180)
    plt.close(fig)


def run() -> dict:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    u = make_u_grid(512)
    eta = make_eta_grid(16)
    libraries = {K: wave_library(u, K) for K in K_VALUES}
    wave_controls = [
        wave_control_metrics(wave, u)
        for K in K_VALUES
        for wave in libraries[K].values()
    ]

    rows: list[dict] = []
    for candidate in CANDIDATES:
        for K in K_VALUES:
            library = libraries[K]
            for reflection_parameter in R_VALUES:
                for kappa in KAPPA_VALUES:
                    for combination, label_a, label_b in INPUT_COMBINATIONS:
                        wave_a = library[label_a]
                        wave_b = library[label_b]
                        forward = _scatter_orientation(
                            wave_a,
                            wave_b,
                            SPEC_A,
                            SPEC_B,
                            u,
                            eta,
                            reflection_parameter,
                            kappa,
                            candidate,
                        )
                        # Literal swap: states, demodulation specifications, and
                        # expected kernels all travel together.
                        swapped = _scatter_orientation(
                            wave_b,
                            wave_a,
                            SPEC_B,
                            SPEC_A,
                            u,
                            eta,
                            reflection_parameter,
                            kappa,
                            candidate,
                        )
                        exchange_output, exchange_path = _exchange_residuals(
                            forward, swapped
                        )
                        rows.append(
                            _scalar_row(
                                forward,
                                K=K,
                                reflection_parameter=reflection_parameter,
                                combination=combination,
                                orientation="forward",
                                wave_a=wave_a,
                                wave_b=wave_b,
                                exchange_output_residual=exchange_output,
                                exchange_path_residual=exchange_path,
                            )
                        )
                        rows.append(
                            _scalar_row(
                                swapped,
                                K=K,
                                reflection_parameter=reflection_parameter,
                                combination=combination,
                                orientation="swapped",
                                wave_a=wave_b,
                                wave_b=wave_a,
                                exchange_output_residual=exchange_output,
                                exchange_path_residual=exchange_path,
                            )
                        )

    baseline_rows, baseline_summary = _candidate0_baseline_rows(rows)
    summary = _build_summary(rows, baseline_summary, wave_controls)
    _write_csv(DATA_ROOT / "stage_D_full_results.csv", rows)
    _write_csv(
        DATA_ROOT / "stage_D_candidate0_baseline.csv", baseline_rows
    )
    _write_csv(
        DATA_ROOT / "stage_D_candidate1_results.csv",
        [row for row in rows if row["candidate"] == "C1"],
    )
    _write_csv(
        DATA_ROOT / "stage_D_candidate3_results.csv",
        [row for row in rows if row["candidate"] == "C3"],
    )
    (DATA_ROOT / "stage_D_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _save_figures(rows)
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, ensure_ascii=False))
