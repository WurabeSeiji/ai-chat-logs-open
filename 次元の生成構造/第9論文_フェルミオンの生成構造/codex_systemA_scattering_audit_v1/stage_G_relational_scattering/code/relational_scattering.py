"""C0, reversed_C1, and the single Stage G relational_C1 candidate."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

import numpy as np

from parity_demodulation import ParityReadout, channel_parity, norm2


ScatteringMode = Literal["C0", "reversed_C1", "relational_C1"]
RELATION_NORM2_THRESHOLD = 1.0e-24
GAMMA_RANGE_TOLERANCE = 5.0e-12


class RelationUndefinedError(ValueError):
    def __init__(self, norm2_a: float, norm2_b: float):
        super().__init__(
            "demodulated coherent relation wave below threshold: "
            f"norm2_a={norm2_a:.17g}, norm2_b={norm2_b:.17g}, "
            f"threshold={RELATION_NORM2_THRESHOLD:.17g}"
        )
        self.norm2_a = norm2_a
        self.norm2_b = norm2_b


@dataclass(frozen=True)
class RelationReadout:
    parity_a: ParityReadout
    parity_b: ParityReadout
    coherent_kernel_a: np.ndarray
    coherent_kernel_b: np.ndarray
    coherent_norm2_a: float
    coherent_norm2_b: float
    overlap_complex: complex
    overlap_abs: float
    gamma_ab: float


@dataclass(frozen=True)
class AngleState:
    theta0: float
    rho: float
    c_mean: float
    gamma_ab: float
    candidate_response: float
    delta_theta: float
    theta_eff: float
    r_eff: complex
    t_eff: complex
    R_eff: float
    T_eff: float
    unitarity_residual: float
    orthogonality_residual: float


@dataclass(frozen=True)
class Collision:
    relation: RelationReadout
    angle: AngleState
    path_a_to_a: np.ndarray
    path_b_to_a: np.ndarray
    path_b_to_b: np.ndarray
    path_a_to_b: np.ndarray
    interference_a: float
    interference_b: float
    raw_a: np.ndarray
    raw_b: np.ndarray
    raw_parity_a: ParityReadout
    raw_parity_b: ParityReadout
    path_sum_residual_a: float
    path_sum_residual_b: float
    total_norm_residual: float


def relation_readout(
    state_a: np.ndarray,
    state_b: np.ndarray,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> RelationReadout:
    """Read overlap after correct lineage demodulation and coherent summing."""
    parity_a = channel_parity(state_a, u, eta, p0)
    parity_b = channel_parity(state_b, u, eta, p0)
    kernel_a = parity_a.kernel_origin_a + parity_a.kernel_origin_b
    kernel_b = parity_b.kernel_origin_a + parity_b.kernel_origin_b
    norm_a = norm2(kernel_a)
    norm_b = norm2(kernel_b)
    if (
        norm_a < RELATION_NORM2_THRESHOLD
        or norm_b < RELATION_NORM2_THRESHOLD
    ):
        raise RelationUndefinedError(norm_a, norm_b)
    overlap = complex(np.vdot(kernel_a, kernel_b)) / math.sqrt(
        norm_a * norm_b
    )
    overlap_abs = float(abs(overlap))
    gamma = float(overlap_abs**2)
    if (
        gamma < -GAMMA_RANGE_TOLERANCE
        or gamma > 1.0 + GAMMA_RANGE_TOLERANCE
    ):
        raise ValueError(f"Gamma_AB outside [0,1]: {gamma:.17g}")
    return RelationReadout(
        parity_a=parity_a,
        parity_b=parity_b,
        coherent_kernel_a=kernel_a,
        coherent_kernel_b=kernel_b,
        coherent_norm2_a=float(norm_a),
        coherent_norm2_b=float(norm_b),
        overlap_complex=overlap,
        overlap_abs=overlap_abs,
        gamma_ab=gamma,
    )


def effective_angle(
    scattering_mode: ScatteringMode,
    reflection_baseline: float,
    kappa: float,
    c_a: float,
    c_b: float,
    gamma_ab: float,
) -> AngleState:
    if not 0.0 <= reflection_baseline <= 1.0:
        raise ValueError("reflection_baseline must be in [0,1]")
    if kappa not in (0.01, 0.1, 1.0):
        raise ValueError("Stage G kappa must be one of 0.01, 0.1, 1")
    if (
        gamma_ab < -GAMMA_RANGE_TOLERANCE
        or gamma_ab > 1.0 + GAMMA_RANGE_TOLERANCE
    ):
        raise ValueError("Gamma_AB outside numerical tolerance")
    theta0 = math.asin(math.sqrt(reflection_baseline))
    rho = (2.0 / math.pi) * theta0 * (math.pi / 2.0 - theta0)
    c_mean = 0.5 * (c_a + c_b)
    if scattering_mode == "C0":
        response = 0.0
    elif scattering_mode == "reversed_C1":
        response = -c_mean
    elif scattering_mode == "relational_C1":
        response = -c_mean * gamma_ab
    else:
        raise ValueError(f"unknown scattering mode: {scattering_mode}")
    delta_theta = kappa * rho * response
    theta_eff = theta0 + delta_theta
    endpoint_tolerance = 2.0e-15
    if (
        theta_eff < -endpoint_tolerance
        or theta_eff > math.pi / 2.0 + endpoint_tolerance
    ):
        raise ValueError(
            "theta_eff outside [0,pi/2]; Stage G forbids clipping: "
            f"{theta_eff:.17g}"
        )
    if abs(theta_eff) <= endpoint_tolerance:
        theta_eff = 0.0
    if abs(theta_eff - math.pi / 2.0) <= endpoint_tolerance:
        theta_eff = math.pi / 2.0
    t_eff = complex(np.exp(1j * theta_eff) * math.cos(theta_eff))
    r_eff = complex(
        -1j * np.exp(1j * theta_eff) * math.sin(theta_eff)
    )
    R_eff = float(abs(r_eff) ** 2)
    T_eff = float(abs(t_eff) ** 2)
    return AngleState(
        theta0=float(theta0),
        rho=float(rho),
        c_mean=float(c_mean),
        gamma_ab=float(gamma_ab),
        candidate_response=float(response),
        delta_theta=float(delta_theta),
        theta_eff=float(theta_eff),
        r_eff=r_eff,
        t_eff=t_eff,
        R_eff=R_eff,
        T_eff=T_eff,
        unitarity_residual=float(abs(R_eff + T_eff - 1.0)),
        orthogonality_residual=float(
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
    scattering_mode: ScatteringMode,
    reflection_baseline: float,
    kappa: float,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> Collision:
    relation = relation_readout(input_a, input_b, u, eta, p0)
    angle = effective_angle(
        scattering_mode,
        reflection_baseline,
        kappa,
        relation.parity_a.c_pi,
        relation.parity_b.c_pi,
        relation.gamma_ab,
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
    return Collision(
        relation=relation,
        angle=angle,
        path_a_to_a=path_a_to_a,
        path_b_to_a=path_b_to_a,
        path_b_to_b=path_b_to_b,
        path_a_to_b=path_a_to_b,
        interference_a=interference_a,
        interference_b=interference_b,
        raw_a=raw_a,
        raw_b=raw_b,
        raw_parity_a=raw_parity_a,
        raw_parity_b=raw_parity_b,
        path_sum_residual_a=float(
            abs(
                raw_norm_a
                - norm2(path_a_to_a)
                - norm2(path_b_to_a)
                - interference_a
            )
        ),
        path_sum_residual_b=float(
            abs(
                raw_norm_b
                - norm2(path_b_to_b)
                - norm2(path_a_to_b)
                - interference_b
            )
        ),
        total_norm_residual=float(
            abs(
                raw_norm_a
                + raw_norm_b
                - norm2(input_a)
                - norm2(input_b)
            )
        ),
    )
