"""Independent System A experimental copy with repeated raw-state scattering.

The current C0 kernel and the reversed Candidate 1 hypothesis share this
implementation.  Physical scattering acts on the complete 512 x 16 state.
Parity is recomputed at every collision after resolving both carrier/eta
lineages.  Raw outputs, not channel-normalized outputs, are propagated.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

import numpy as np


KernelName = Literal["C0", "C1_reversed"]


@dataclass(frozen=True)
class ExperimentConfig:
    u_grid_n: int = 512
    eta_grid_n: int = 16
    k_components: int = 4
    reflection_baseline: float = 0.6971778791282474
    collision_count: int = 32
    p0: float = 1.0


@dataclass(frozen=True)
class LineageSpec:
    label: str
    q: float
    eta_mode: int


@dataclass(frozen=True)
class ChannelParity:
    correlation_raw: complex
    indicator: float
    boson_weight: float
    fermion_weight: float
    norm2: float
    origin_a_weight: float
    origin_b_weight: float
    reconstruction_residual: float
    kernel_origin_a: np.ndarray
    kernel_origin_b: np.ndarray


@dataclass(frozen=True)
class AngleState:
    theta_0: float
    rho: float
    parity_mean: float
    delta_theta: float
    theta_eff: float
    reflection_probability: float
    transmission_probability: float
    r_eff: complex
    t_eff: complex
    unitarity_residual: float
    coefficient_orthogonality_residual: float


@dataclass(frozen=True)
class CollisionResult:
    input_a: np.ndarray
    input_b: np.ndarray
    input_parity_a: ChannelParity
    input_parity_b: ChannelParity
    angle: AngleState
    path_a_to_a: np.ndarray
    path_b_to_a: np.ndarray
    path_b_to_b: np.ndarray
    path_a_to_b: np.ndarray
    interference_a: float
    interference_b: float
    raw_output_a: np.ndarray
    raw_output_b: np.ndarray
    raw_parity_a: ChannelParity
    raw_parity_b: ChannelParity
    path_sum_residual_a: float
    path_sum_residual_b: float
    pair_norm_conservation_residual: float


SPEC_ORIGIN_A = LineageSpec("origin_A", q=1.0, eta_mode=1)
SPEC_ORIGIN_B = LineageSpec("origin_B", q=-1.0, eta_mode=2)


def norm2(vector: np.ndarray) -> float:
    return float(np.vdot(vector, vector).real)


def normalize(vector: np.ndarray) -> np.ndarray:
    value = math.sqrt(max(norm2(vector), 0.0))
    if value <= 0.0:
        raise ValueError("zero norm state")
    return vector / value


def make_grids(config: ExperimentConfig) -> tuple[np.ndarray, np.ndarray]:
    u = np.linspace(-math.pi, math.pi, config.u_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, config.eta_grid_n, endpoint=False)
    return u, eta


def fermion_kernel(u: np.ndarray, k_components: int) -> np.ndarray:
    wave = sum(
        (np.cos((2 * index + 1) * u) for index in range(k_components)),
        start=np.zeros_like(u),
    ) / math.sqrt(k_components)
    return normalize(np.asarray(wave, dtype=complex))


def boson_kernel(u: np.ndarray, k_components: int) -> np.ndarray:
    wave = sum(
        (
            np.cos(2 * (index + 1) * u)
            for index in range(k_components)
        ),
        start=np.zeros_like(u),
    ) / math.sqrt(k_components)
    return normalize(np.asarray(wave, dtype=complex))


def half_shift(kernel: np.ndarray) -> np.ndarray:
    if kernel.ndim != 1 or kernel.size % 2:
        raise ValueError("kernel must be a one-dimensional even-length array")
    return np.roll(kernel, -kernel.size // 2)


def modulate_kernel(
    kernel: np.ndarray,
    spec: LineageSpec,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> np.ndarray:
    carrier = np.exp(1j * spec.q * p0 * u)
    eta_embedding = (
        np.exp(1j * spec.eta_mode * eta) / math.sqrt(eta.size)
    )
    return ((kernel * carrier)[:, None] * eta_embedding[None, :]).reshape(-1)


def project_lineage_kernel(
    state: np.ndarray,
    spec: LineageSpec,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> np.ndarray:
    if state.size != u.size * eta.size:
        raise ValueError("state does not match the 512 x 16 physical grid")
    array = state.reshape(u.size, eta.size)
    eta_embedding = (
        np.exp(1j * spec.eta_mode * eta) / math.sqrt(eta.size)
    )
    carrier_state = np.sum(
        array * np.conjugate(eta_embedding)[None, :],
        axis=1,
    )
    return carrier_state * np.exp(-1j * spec.q * p0 * u)


def channel_parity(
    state: np.ndarray,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> ChannelParity:
    """Measure parity after resolving both persistent source lineages."""
    kernel_a = project_lineage_kernel(
        state, SPEC_ORIGIN_A, u, eta, p0
    )
    kernel_b = project_lineage_kernel(
        state, SPEC_ORIGIN_B, u, eta, p0
    )
    norm_a = norm2(kernel_a)
    norm_b = norm2(kernel_b)
    total_norm = norm_a + norm_b
    if total_norm <= 0.0:
        raise ValueError("zero norm channel")

    correlation = 0.0j
    boson_numerator = 0.0
    fermion_numerator = 0.0
    for kernel, kernel_norm in ((kernel_a, norm_a), (kernel_b, norm_b)):
        if kernel_norm == 0.0:
            continue
        shifted = half_shift(kernel)
        correlation += complex(np.vdot(kernel, shifted))
        boson_numerator += norm2(0.5 * (kernel + shifted))
        fermion_numerator += norm2(0.5 * (kernel - shifted))

    reconstructed = modulate_kernel(
        kernel_a, SPEC_ORIGIN_A, u, eta, p0
    ) + modulate_kernel(kernel_b, SPEC_ORIGIN_B, u, eta, p0)
    state_norm = max(math.sqrt(norm2(state)), 1.0e-300)
    reconstruction_residual = (
        math.sqrt(norm2(reconstructed - state)) / state_norm
    )
    return ChannelParity(
        correlation_raw=correlation,
        indicator=float(correlation.real / total_norm),
        boson_weight=float(boson_numerator / total_norm),
        fermion_weight=float(fermion_numerator / total_norm),
        norm2=float(total_norm),
        origin_a_weight=float(norm_a / total_norm),
        origin_b_weight=float(norm_b / total_norm),
        reconstruction_residual=float(reconstruction_residual),
        kernel_origin_a=kernel_a,
        kernel_origin_b=kernel_b,
    )


def effective_angle(
    kernel_name: KernelName,
    reflection_baseline: float,
    kappa: float,
    c_a: float,
    c_b: float,
) -> AngleState:
    if not 0.0 <= reflection_baseline <= 1.0:
        raise ValueError("reflection_baseline must be in [0,1]")
    if not 0.0 <= kappa <= 1.0:
        raise ValueError("kappa must be in [0,1]")
    theta_0 = math.asin(math.sqrt(reflection_baseline))
    rho = (2.0 / math.pi) * theta_0 * (math.pi / 2.0 - theta_0)
    parity_mean = 0.5 * (c_a + c_b)
    if kernel_name == "C0":
        delta_theta = 0.0
    elif kernel_name == "C1_reversed":
        delta_theta = -kappa * rho * parity_mean
    else:
        raise ValueError(f"unknown kernel: {kernel_name}")
    theta_eff = theta_0 + delta_theta
    tolerance = 2.0e-15
    if theta_eff < -tolerance or theta_eff > math.pi / 2.0 + tolerance:
        raise ValueError(f"theta_eff outside [0,pi/2]: {theta_eff}")
    theta_eff = min(max(theta_eff, 0.0), math.pi / 2.0)
    t_eff = complex(np.exp(1j * theta_eff) * math.cos(theta_eff))
    r_eff = complex(-1j * np.exp(1j * theta_eff) * math.sin(theta_eff))
    reflection = float(abs(r_eff) ** 2)
    transmission = float(abs(t_eff) ** 2)
    return AngleState(
        theta_0=float(theta_0),
        rho=float(rho),
        parity_mean=float(parity_mean),
        delta_theta=float(delta_theta),
        theta_eff=float(theta_eff),
        reflection_probability=reflection,
        transmission_probability=transmission,
        r_eff=r_eff,
        t_eff=t_eff,
        unitarity_residual=float(abs(reflection + transmission - 1.0)),
        coefficient_orthogonality_residual=float(
            abs(
                np.conjugate(r_eff) * t_eff
                + np.conjugate(t_eff) * r_eff
            )
        ),
    )


def scatter_once(
    input_a: np.ndarray,
    input_b: np.ndarray,
    *,
    kernel_name: KernelName,
    reflection_baseline: float,
    kappa: float,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> CollisionResult:
    parity_a = channel_parity(input_a, u, eta, p0)
    parity_b = channel_parity(input_b, u, eta, p0)
    angle = effective_angle(
        kernel_name,
        reflection_baseline,
        kappa,
        parity_a.indicator,
        parity_b.indicator,
    )
    path_a_to_a = angle.r_eff * input_a
    path_b_to_a = angle.t_eff * input_b
    path_b_to_b = angle.r_eff * input_b
    path_a_to_b = angle.t_eff * input_a
    interference_a = float(
        2.0 * np.vdot(path_a_to_a, path_b_to_a).real
    )
    interference_b = float(
        2.0 * np.vdot(path_b_to_b, path_a_to_b).real
    )
    raw_a = path_a_to_a + path_b_to_a
    raw_b = path_b_to_b + path_a_to_b
    raw_parity_a = channel_parity(raw_a, u, eta, p0)
    raw_parity_b = channel_parity(raw_b, u, eta, p0)
    raw_norm_a = norm2(raw_a)
    raw_norm_b = norm2(raw_b)
    path_sum_a = abs(
        raw_norm_a
        - (
            norm2(path_a_to_a)
            + norm2(path_b_to_a)
            + interference_a
        )
    )
    path_sum_b = abs(
        raw_norm_b
        - (
            norm2(path_b_to_b)
            + norm2(path_a_to_b)
            + interference_b
        )
    )
    pair_norm_residual = abs(
        raw_norm_a
        + raw_norm_b
        - norm2(input_a)
        - norm2(input_b)
    )
    return CollisionResult(
        input_a=input_a,
        input_b=input_b,
        input_parity_a=parity_a,
        input_parity_b=parity_b,
        angle=angle,
        path_a_to_a=path_a_to_a,
        path_b_to_a=path_b_to_a,
        path_b_to_b=path_b_to_b,
        path_a_to_b=path_a_to_b,
        interference_a=interference_a,
        interference_b=interference_b,
        raw_output_a=raw_a,
        raw_output_b=raw_b,
        raw_parity_a=raw_parity_a,
        raw_parity_b=raw_parity_b,
        path_sum_residual_a=float(path_sum_a),
        path_sum_residual_b=float(path_sum_b),
        pair_norm_conservation_residual=float(pair_norm_residual),
    )


def initial_state_pair(
    case_id: str,
    config: ExperimentConfig,
    u: np.ndarray,
    eta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    fermion = fermion_kernel(u, config.k_components)
    boson = boson_kernel(u, config.k_components)
    cases = {
        "F_x_F": (fermion, fermion),
        "B_x_B": (boson, boson),
        "F_x_B": (fermion, boson),
    }
    if case_id not in cases:
        raise ValueError(f"unknown case: {case_id}")
    kernel_a, kernel_b = cases[case_id]
    return (
        modulate_kernel(kernel_a, SPEC_ORIGIN_A, u, eta, config.p0),
        modulate_kernel(kernel_b, SPEC_ORIGIN_B, u, eta, config.p0),
    )


def collision_row(
    run_id: str,
    case_id: str,
    kernel_name: KernelName,
    kappa: float,
    collision: int,
    result: CollisionResult,
    config: ExperimentConfig,
) -> dict:
    angle = result.angle
    input_a = result.input_parity_a
    input_b = result.input_parity_b
    raw_a = result.raw_parity_a
    raw_b = result.raw_parity_b
    return {
        "run_id": run_id,
        "case_id": case_id,
        "kernel": kernel_name,
        "kappa_requested": kappa,
        "kappa_applied": 0.0 if kernel_name == "C0" else kappa,
        "collision": collision,
        "collision_count": config.collision_count,
        "K": config.k_components,
        "u_grid_n": config.u_grid_n,
        "eta_grid_n": config.eta_grid_n,
        "R_0": config.reflection_baseline,
        "raw_output_is_next_physical_state": True,
        "channel_normalization_applied": False,
        "c_A": input_a.indicator,
        "c_B": input_b.indicator,
        "c_mean": angle.parity_mean,
        "theta_0": angle.theta_0,
        "rho": angle.rho,
        "delta_theta": angle.delta_theta,
        "theta_eff": angle.theta_eff,
        "R_eff": angle.reflection_probability,
        "T_eff": angle.transmission_probability,
        "r_eff_real": angle.r_eff.real,
        "r_eff_imag": angle.r_eff.imag,
        "t_eff_real": angle.t_eff.real,
        "t_eff_imag": angle.t_eff.imag,
        "unitarity_residual": angle.unitarity_residual,
        "coefficient_orthogonality_residual": (
            angle.coefficient_orthogonality_residual
        ),
        "p_B_A_input": input_a.boson_weight,
        "p_F_A_input": input_a.fermion_weight,
        "p_B_B_input": input_b.boson_weight,
        "p_F_B_input": input_b.fermion_weight,
        "origin_A_weight_in_A_input": input_a.origin_a_weight,
        "origin_B_weight_in_A_input": input_a.origin_b_weight,
        "origin_A_weight_in_B_input": input_b.origin_a_weight,
        "origin_B_weight_in_B_input": input_b.origin_b_weight,
        "input_norm2_A": norm2(result.input_a),
        "input_norm2_B": norm2(result.input_b),
        "input_reconstruction_residual_A": (
            input_a.reconstruction_residual
        ),
        "input_reconstruction_residual_B": (
            input_b.reconstruction_residual
        ),
        "path_a_to_a_norm2": norm2(result.path_a_to_a),
        "path_b_to_a_norm2": norm2(result.path_b_to_a),
        "path_b_to_b_norm2": norm2(result.path_b_to_b),
        "path_a_to_b_norm2": norm2(result.path_a_to_b),
        "I_A": result.interference_a,
        "I_B": result.interference_b,
        "a_raw_norm2": norm2(result.raw_output_a),
        "b_raw_norm2": norm2(result.raw_output_b),
        "p_B_A_raw": raw_a.boson_weight,
        "p_F_A_raw": raw_a.fermion_weight,
        "p_B_B_raw": raw_b.boson_weight,
        "p_F_B_raw": raw_b.fermion_weight,
        "c_A_raw": raw_a.indicator,
        "c_B_raw": raw_b.indicator,
        "origin_A_weight_in_A_raw": raw_a.origin_a_weight,
        "origin_B_weight_in_A_raw": raw_a.origin_b_weight,
        "origin_A_weight_in_B_raw": raw_b.origin_a_weight,
        "origin_B_weight_in_B_raw": raw_b.origin_b_weight,
        "raw_reconstruction_residual_A": raw_a.reconstruction_residual,
        "raw_reconstruction_residual_B": raw_b.reconstruction_residual,
        "path_sum_residual_A": result.path_sum_residual_a,
        "path_sum_residual_B": result.path_sum_residual_b,
        "pair_norm_conservation_residual": (
            result.pair_norm_conservation_residual
        ),
    }


def run_repeated(
    case_id: str,
    kernel_name: KernelName,
    kappa: float,
    config: ExperimentConfig,
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    u, eta = make_grids(config)
    state_a, state_b = initial_state_pair(case_id, config, u, eta)
    run_id = f"{kernel_name}__{case_id}__kappa_{kappa:g}"
    rows = []
    for collision in range(1, config.collision_count + 1):
        result = scatter_once(
            state_a,
            state_b,
            kernel_name=kernel_name,
            reflection_baseline=config.reflection_baseline,
            kappa=kappa,
            u=u,
            eta=eta,
            p0=config.p0,
        )
        rows.append(
            collision_row(
                run_id,
                case_id,
                kernel_name,
                kappa,
                collision,
                result,
                config,
            )
        )
        # This is the defining Stage E update: do not normalize the channels.
        state_a = result.raw_output_a
        state_b = result.raw_output_b
    return rows, state_a, state_b
