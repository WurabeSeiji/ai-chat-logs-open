"""C0 and reversed Candidate 1 scattering coefficients for Stage F."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

import numpy as np

from parity_demodulation import ParityReadout, channel_parity, norm2


ScatteringMode = Literal["C0", "reversed_C1"]


@dataclass(frozen=True)
class AngleState:
    theta0: float
    rho: float
    c_mean: float
    delta_theta: float
    theta_eff: float
    r_eff: complex
    t_eff: complex
    R_eff: float
    T_eff: float
    unitarity_residual: float
    coefficient_orthogonality_residual: float


@dataclass(frozen=True)
class Collision:
    input_parity_a: ParityReadout
    input_parity_b: ParityReadout
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


def effective_angle(
    scattering_mode: ScatteringMode,
    reflection_baseline: float,
    kappa: float,
    c_a: float,
    c_b: float,
) -> AngleState:
    if not 0.0 <= reflection_baseline <= 1.0:
        raise ValueError("reflection_baseline must be in [0,1]")
    if kappa not in (0.01, 0.1, 1.0):
        raise ValueError("Stage F kappa must be one of 0.01, 0.1, 1")
    theta0 = math.asin(math.sqrt(reflection_baseline))
    rho = (2.0 / math.pi) * theta0 * (math.pi / 2.0 - theta0)
    c_mean = 0.5 * (c_a + c_b)
    if scattering_mode == "C0":
        delta_theta = 0.0
    elif scattering_mode == "reversed_C1":
        delta_theta = -kappa * rho * c_mean
    else:
        raise ValueError(f"unknown scattering mode: {scattering_mode}")
    theta_eff = theta0 + delta_theta
    tolerance = 2.0e-15
    if theta_eff < -tolerance or theta_eff > math.pi / 2.0 + tolerance:
        raise ValueError(
            "theta_eff outside [0,pi/2]; Stage F forbids clipping: "
            f"{theta_eff}"
        )
    # Round-off at the two exact endpoints is not a physical clip.
    if abs(theta_eff) <= tolerance:
        theta_eff = 0.0
    if abs(theta_eff - math.pi / 2.0) <= tolerance:
        theta_eff = math.pi / 2.0
    t_eff = complex(np.exp(1j * theta_eff) * math.cos(theta_eff))
    r_eff = complex(-1j * np.exp(1j * theta_eff) * math.sin(theta_eff))
    R_eff = float(abs(r_eff) ** 2)
    T_eff = float(abs(t_eff) ** 2)
    return AngleState(
        theta0=float(theta0),
        rho=float(rho),
        c_mean=float(c_mean),
        delta_theta=float(delta_theta),
        theta_eff=float(theta_eff),
        r_eff=r_eff,
        t_eff=t_eff,
        R_eff=R_eff,
        T_eff=T_eff,
        unitarity_residual=float(abs(R_eff + T_eff - 1.0)),
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
    scattering_mode: ScatteringMode,
    reflection_baseline: float,
    kappa: float,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> Collision:
    parity_a = channel_parity(input_a, u, eta, p0)
    parity_b = channel_parity(input_b, u, eta, p0)
    angle = effective_angle(
        scattering_mode,
        reflection_baseline,
        kappa,
        parity_a.c_pi,
        parity_b.c_pi,
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
    path_sum_residual_a = abs(
        raw_norm_a
        - norm2(path_a_to_a)
        - norm2(path_b_to_a)
        - interference_a
    )
    path_sum_residual_b = abs(
        raw_norm_b
        - norm2(path_b_to_b)
        - norm2(path_a_to_b)
        - interference_b
    )
    total_norm_residual = abs(
        raw_norm_a
        + raw_norm_b
        - norm2(input_a)
        - norm2(input_b)
    )
    return Collision(
        input_parity_a=parity_a,
        input_parity_b=parity_b,
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
        path_sum_residual_a=float(path_sum_residual_a),
        path_sum_residual_b=float(path_sum_residual_b),
        total_norm_residual=float(total_norm_residual),
    )
