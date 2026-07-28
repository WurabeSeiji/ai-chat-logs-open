"""Common two-layer scattering API for Stage D Candidates 0, 1, and 3."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from candidate_responses import (
    CandidateName,
    candidate_response,
    effective_angle,
    scattering_coefficients,
)
from parity_metrics import (
    ModulationSpec,
    combined_lineage_parity,
    demodulate_state,
    lifted_half_shift,
    lineage_reconstruction_residual,
    modulate_kernel,
    norm2,
    normalize,
    parity_measurement,
)


@dataclass(frozen=True)
class ScatteringResult:
    candidate: str
    kappa: float

    input_a: np.ndarray
    input_b: np.ndarray
    input_norm_a: float
    input_norm_b: float

    demodulated_kernel_a: np.ndarray
    demodulated_kernel_b: np.ndarray
    demodulation_residual_a: float
    demodulation_residual_b: float
    roundtrip_reconstruction_residual_a: float
    roundtrip_reconstruction_residual_b: float
    modulation_q_a: float
    modulation_q_b: float
    eta_mode_a: int
    eta_mode_b: int

    theta_0: float
    rho: float
    candidate_response: float
    delta_theta: float
    theta_eff_preclip: float
    theta_eff: float
    theta_was_clipped: bool
    candidate_product_norm2: float | None
    candidate_product_threshold: float | None
    candidate_response_status: str

    r_eff: complex
    t_eff: complex
    reflection_probability: float
    transmission_probability: float
    unitarity_residual: float
    orthogonality_residual: float

    parity_correlation_raw_a: complex
    parity_correlation_raw_b: complex
    parity_indicator_a: float
    parity_indicator_b: float
    boson_weight_a: float
    boson_weight_b: float
    fermion_weight_a: float
    fermion_weight_b: float

    path_a_to_a_amplitude: np.ndarray
    path_b_to_a_amplitude: np.ndarray
    path_b_to_b_amplitude: np.ndarray
    path_a_to_b_amplitude: np.ndarray

    path_a_to_a_norm: float
    path_b_to_a_norm: float
    path_b_to_b_norm: float
    path_a_to_b_norm: float

    interference_in_a: float
    interference_in_b: float

    raw_output_a: np.ndarray
    raw_output_b: np.ndarray
    raw_output_norm_a: float
    raw_output_norm_b: float

    path_sum_residual_a: float
    path_sum_residual_b: float

    normalized_output_a: np.ndarray
    normalized_output_b: np.ndarray
    normalized_output_norm_a: float
    normalized_output_norm_b: float

    raw_output_parity_correlation_a: complex
    raw_output_parity_correlation_b: complex
    raw_output_parity_indicator_a: float
    raw_output_parity_indicator_b: float
    raw_output_boson_weight_a: float
    raw_output_boson_weight_b: float
    raw_output_fermion_weight_a: float
    raw_output_fermion_weight_b: float
    raw_output_lineage_reconstruction_residual_a: float
    raw_output_lineage_reconstruction_residual_b: float

    pair_norm_conservation_residual: float
    half_shift_equivariance_residual: float


def _core_scatter(
    input_a: np.ndarray,
    input_b: np.ndarray,
    *,
    expected_kernel_a: np.ndarray,
    expected_kernel_b: np.ndarray,
    spec_a: ModulationSpec,
    spec_b: ModulationSpec,
    u: np.ndarray,
    eta: np.ndarray,
    reflection_parameter: float,
    kappa: float,
    candidate: CandidateName,
    product_relative_threshold: float,
) -> dict:
    demod_a = demodulate_state(input_a, expected_kernel_a, u, eta, spec_a)
    demod_b = demodulate_state(input_b, expected_kernel_b, u, eta, spec_b)
    response = candidate_response(
        candidate,
        demod_a.kernel,
        demod_b.kernel,
        product_relative_threshold=product_relative_threshold,
    )
    angle = effective_angle(reflection_parameter, kappa, response.value)
    t_eff, r_eff = scattering_coefficients(angle.theta_eff)
    path_a_to_a = r_eff * input_a
    path_b_to_a = t_eff * input_b
    path_b_to_b = r_eff * input_b
    path_a_to_b = t_eff * input_a
    raw_a = path_a_to_a + path_b_to_a
    raw_b = path_b_to_b + path_a_to_b
    return {
        "demod_a": demod_a,
        "demod_b": demod_b,
        "response": response,
        "angle": angle,
        "t_eff": t_eff,
        "r_eff": r_eff,
        "path_a_to_a": path_a_to_a,
        "path_b_to_a": path_b_to_a,
        "path_b_to_b": path_b_to_b,
        "path_a_to_b": path_a_to_b,
        "raw_a": raw_a,
        "raw_b": raw_b,
    }


def _half_shift_equivariance_residual(
    core: dict,
    input_a: np.ndarray,
    input_b: np.ndarray,
    *,
    expected_kernel_a: np.ndarray,
    expected_kernel_b: np.ndarray,
    spec_a: ModulationSpec,
    spec_b: ModulationSpec,
    u: np.ndarray,
    eta: np.ndarray,
    reflection_parameter: float,
    kappa: float,
    candidate: CandidateName,
    product_relative_threshold: float,
) -> float:
    shifted_a, shifted_kernel_a = lifted_half_shift(
        input_a, expected_kernel_a, u, eta, spec_a
    )
    shifted_b, shifted_kernel_b = lifted_half_shift(
        input_b, expected_kernel_b, u, eta, spec_b
    )
    left = _core_scatter(
        shifted_a,
        shifted_b,
        expected_kernel_a=shifted_kernel_a,
        expected_kernel_b=shifted_kernel_b,
        spec_a=spec_a,
        spec_b=spec_b,
        u=u,
        eta=eta,
        reflection_parameter=reflection_parameter,
        kappa=kappa,
        candidate=candidate,
        product_relative_threshold=product_relative_threshold,
    )

    right_aa, _ = lifted_half_shift(
        core["path_a_to_a"],
        core["r_eff"] * core["demod_a"].kernel,
        u,
        eta,
        spec_a,
    )
    right_ba, _ = lifted_half_shift(
        core["path_b_to_a"],
        core["t_eff"] * core["demod_b"].kernel,
        u,
        eta,
        spec_b,
    )
    right_bb, _ = lifted_half_shift(
        core["path_b_to_b"],
        core["r_eff"] * core["demod_b"].kernel,
        u,
        eta,
        spec_b,
    )
    right_ab, _ = lifted_half_shift(
        core["path_a_to_b"],
        core["t_eff"] * core["demod_a"].kernel,
        u,
        eta,
        spec_a,
    )
    right_a = right_aa + right_ba
    right_b = right_bb + right_ab
    numerator = math.sqrt(
        norm2(left["raw_a"] - right_a) + norm2(left["raw_b"] - right_b)
    )
    denominator = max(
        math.sqrt(norm2(core["raw_a"]) + norm2(core["raw_b"])),
        1.0e-300,
    )
    return float(numerator / denominator)


def scatter_wave_pair(
    input_a: np.ndarray,
    input_b: np.ndarray,
    *,
    expected_kernel_a: np.ndarray,
    expected_kernel_b: np.ndarray,
    spec_a: ModulationSpec,
    spec_b: ModulationSpec,
    u: np.ndarray,
    eta: np.ndarray,
    reflection_parameter: float,
    kappa: float,
    candidate: CandidateName,
    product_relative_threshold: float = 1.0e-14,
) -> ScatteringResult:
    """Scatter full physical states; read candidate response from demodulated kernels."""
    core = _core_scatter(
        input_a,
        input_b,
        expected_kernel_a=expected_kernel_a,
        expected_kernel_b=expected_kernel_b,
        spec_a=spec_a,
        spec_b=spec_b,
        u=u,
        eta=eta,
        reflection_parameter=reflection_parameter,
        kappa=kappa,
        candidate=candidate,
        product_relative_threshold=product_relative_threshold,
    )
    r_eff = core["r_eff"]
    t_eff = core["t_eff"]
    angle = core["angle"]
    response = core["response"]
    demod_a = core["demod_a"]
    demod_b = core["demod_b"]
    input_parity_a = parity_measurement(demod_a.kernel)
    input_parity_b = parity_measurement(demod_b.kernel)

    path_a_to_a_norm = norm2(core["path_a_to_a"])
    path_b_to_a_norm = norm2(core["path_b_to_a"])
    path_b_to_b_norm = norm2(core["path_b_to_b"])
    path_a_to_b_norm = norm2(core["path_a_to_b"])
    interference_a = float(
        2.0 * np.vdot(core["path_a_to_a"], core["path_b_to_a"]).real
    )
    interference_b = float(
        2.0 * np.vdot(core["path_b_to_b"], core["path_a_to_b"]).real
    )
    raw_norm_a = norm2(core["raw_a"])
    raw_norm_b = norm2(core["raw_b"])
    normalized_a = normalize(core["raw_a"])
    normalized_b = normalize(core["raw_b"])

    kernel_path_aa = demodulate_state(
        core["path_a_to_a"],
        r_eff * demod_a.kernel,
        u,
        eta,
        spec_a,
    ).kernel
    kernel_path_ba = demodulate_state(
        core["path_b_to_a"],
        t_eff * demod_b.kernel,
        u,
        eta,
        spec_b,
    ).kernel
    kernel_path_bb = demodulate_state(
        core["path_b_to_b"],
        r_eff * demod_b.kernel,
        u,
        eta,
        spec_b,
    ).kernel
    kernel_path_ab = demodulate_state(
        core["path_a_to_b"],
        t_eff * demod_a.kernel,
        u,
        eta,
        spec_a,
    ).kernel
    raw_parity_a = combined_lineage_parity((kernel_path_aa, kernel_path_ba))
    raw_parity_b = combined_lineage_parity((kernel_path_bb, kernel_path_ab))
    lineage_residual_a = lineage_reconstruction_residual(
        (core["path_a_to_a"], core["path_b_to_a"]),
        (kernel_path_aa, kernel_path_ba),
        (spec_a, spec_b),
        u,
        eta,
    )
    lineage_residual_b = lineage_reconstruction_residual(
        (core["path_b_to_b"], core["path_a_to_b"]),
        (kernel_path_bb, kernel_path_ab),
        (spec_b, spec_a),
        u,
        eta,
    )

    reflection_probability = float(abs(r_eff) ** 2)
    transmission_probability = float(abs(t_eff) ** 2)
    unitarity_residual = abs(
        reflection_probability + transmission_probability - 1.0
    )
    orthogonality_residual = abs(
        np.conjugate(r_eff) * t_eff + np.conjugate(t_eff) * r_eff
    )
    path_sum_residual_a = abs(
        raw_norm_a
        - (path_a_to_a_norm + path_b_to_a_norm + interference_a)
    )
    path_sum_residual_b = abs(
        raw_norm_b
        - (path_b_to_b_norm + path_a_to_b_norm + interference_b)
    )
    pair_norm_residual = abs(
        raw_norm_a
        + raw_norm_b
        - norm2(input_a)
        - norm2(input_b)
    )
    half_shift_residual = _half_shift_equivariance_residual(
        core,
        input_a,
        input_b,
        expected_kernel_a=expected_kernel_a,
        expected_kernel_b=expected_kernel_b,
        spec_a=spec_a,
        spec_b=spec_b,
        u=u,
        eta=eta,
        reflection_parameter=reflection_parameter,
        kappa=kappa,
        candidate=candidate,
        product_relative_threshold=product_relative_threshold,
    )

    return ScatteringResult(
        candidate=candidate,
        kappa=float(kappa),
        input_a=input_a,
        input_b=input_b,
        input_norm_a=norm2(input_a),
        input_norm_b=norm2(input_b),
        demodulated_kernel_a=demod_a.kernel,
        demodulated_kernel_b=demod_b.kernel,
        demodulation_residual_a=demod_a.reference_residual,
        demodulation_residual_b=demod_b.reference_residual,
        roundtrip_reconstruction_residual_a=demod_a.roundtrip_residual,
        roundtrip_reconstruction_residual_b=demod_b.roundtrip_residual,
        modulation_q_a=spec_a.q,
        modulation_q_b=spec_b.q,
        eta_mode_a=spec_a.eta_mode,
        eta_mode_b=spec_b.eta_mode,
        theta_0=angle.theta_0,
        rho=angle.rho,
        candidate_response=angle.candidate_response,
        delta_theta=angle.delta_theta,
        theta_eff_preclip=angle.theta_eff_preclip,
        theta_eff=angle.theta_eff,
        theta_was_clipped=angle.was_clipped,
        candidate_product_norm2=response.product_norm2,
        candidate_product_threshold=response.product_threshold,
        candidate_response_status=response.status,
        r_eff=r_eff,
        t_eff=t_eff,
        reflection_probability=reflection_probability,
        transmission_probability=transmission_probability,
        unitarity_residual=float(unitarity_residual),
        orthogonality_residual=float(orthogonality_residual),
        parity_correlation_raw_a=input_parity_a.correlation_raw,
        parity_correlation_raw_b=input_parity_b.correlation_raw,
        parity_indicator_a=input_parity_a.indicator,
        parity_indicator_b=input_parity_b.indicator,
        boson_weight_a=input_parity_a.boson_weight,
        boson_weight_b=input_parity_b.boson_weight,
        fermion_weight_a=input_parity_a.fermion_weight,
        fermion_weight_b=input_parity_b.fermion_weight,
        path_a_to_a_amplitude=core["path_a_to_a"],
        path_b_to_a_amplitude=core["path_b_to_a"],
        path_b_to_b_amplitude=core["path_b_to_b"],
        path_a_to_b_amplitude=core["path_a_to_b"],
        path_a_to_a_norm=path_a_to_a_norm,
        path_b_to_a_norm=path_b_to_a_norm,
        path_b_to_b_norm=path_b_to_b_norm,
        path_a_to_b_norm=path_a_to_b_norm,
        interference_in_a=interference_a,
        interference_in_b=interference_b,
        raw_output_a=core["raw_a"],
        raw_output_b=core["raw_b"],
        raw_output_norm_a=raw_norm_a,
        raw_output_norm_b=raw_norm_b,
        path_sum_residual_a=float(path_sum_residual_a),
        path_sum_residual_b=float(path_sum_residual_b),
        normalized_output_a=normalized_a,
        normalized_output_b=normalized_b,
        normalized_output_norm_a=norm2(normalized_a),
        normalized_output_norm_b=norm2(normalized_b),
        raw_output_parity_correlation_a=raw_parity_a.correlation_raw,
        raw_output_parity_correlation_b=raw_parity_b.correlation_raw,
        raw_output_parity_indicator_a=raw_parity_a.indicator,
        raw_output_parity_indicator_b=raw_parity_b.indicator,
        raw_output_boson_weight_a=raw_parity_a.boson_weight,
        raw_output_boson_weight_b=raw_parity_b.boson_weight,
        raw_output_fermion_weight_a=raw_parity_a.fermion_weight,
        raw_output_fermion_weight_b=raw_parity_b.fermion_weight,
        raw_output_lineage_reconstruction_residual_a=lineage_residual_a,
        raw_output_lineage_reconstruction_residual_b=lineage_residual_b,
        pair_norm_conservation_residual=float(pair_norm_residual),
        half_shift_equivariance_residual=half_shift_residual,
    )
