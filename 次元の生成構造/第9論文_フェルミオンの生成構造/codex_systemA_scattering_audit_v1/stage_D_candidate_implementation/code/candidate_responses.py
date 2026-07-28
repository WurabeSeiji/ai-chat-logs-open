"""Candidate 0, 1, and 3 response functions and shared angle map."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

import numpy as np

from parity_metrics import norm2, parity_measurement


CandidateName = Literal["C0", "C1", "C3"]


class ProductNormTooSmallError(ValueError):
    pass


@dataclass(frozen=True)
class CandidateResponse:
    candidate: CandidateName
    value: float
    product_norm2: float | None
    product_threshold: float | None
    status: str


@dataclass(frozen=True)
class EffectiveAngle:
    theta_0: float
    rho: float
    candidate_response: float
    delta_theta: float
    theta_eff_preclip: float
    theta_eff: float
    was_clipped: bool


def candidate_response(
    candidate: CandidateName,
    kernel_a: np.ndarray,
    kernel_b: np.ndarray,
    *,
    product_relative_threshold: float = 1.0e-14,
) -> CandidateResponse:
    if candidate == "C0":
        return CandidateResponse(candidate, 0.0, None, None, "baseline")
    if candidate == "C1":
        c_a = parity_measurement(kernel_a).indicator
        c_b = parity_measurement(kernel_b).indicator
        value = 0.5 * (c_a + c_b)
        return CandidateResponse(candidate, float(value), None, None, "ok")
    if candidate == "C3":
        product = kernel_a * kernel_b
        product_norm = norm2(product)
        threshold = (
            product_relative_threshold * norm2(kernel_a) * norm2(kernel_b)
        )
        if product_norm <= threshold:
            raise ProductNormTooSmallError(
                "Candidate 3 product norm is at or below the preregistered "
                f"threshold: norm2={product_norm}, threshold={threshold}"
            )
        value = parity_measurement(product).indicator
        return CandidateResponse(
            candidate,
            float(value),
            float(product_norm),
            float(threshold),
            "ok",
        )
    raise ValueError(f"unknown candidate: {candidate}")


def effective_angle(
    reflection_parameter: float,
    kappa: float,
    response: float,
) -> EffectiveAngle:
    if not 0.0 <= reflection_parameter <= 1.0:
        raise ValueError("reflection_parameter must be in [0,1]")
    if not 0.0 <= kappa <= 1.0:
        raise ValueError("Stage D kappa must be in [0,1]")
    if abs(response) > 1.0 + 1.0e-12:
        raise ValueError("candidate response is outside [-1,1]")
    theta_0 = math.asin(math.sqrt(reflection_parameter))
    rho = (2.0 / math.pi) * theta_0 * (math.pi / 2.0 - theta_0)
    delta_theta = kappa * rho * response
    preclip = theta_0 + delta_theta
    tolerance = 2.0e-15
    if preclip < -tolerance or preclip > math.pi / 2.0 + tolerance:
        raise ValueError(f"theta_eff range failure: {preclip}")
    theta_eff = min(max(preclip, 0.0), math.pi / 2.0)
    was_clipped = theta_eff != preclip
    return EffectiveAngle(
        theta_0=float(theta_0),
        rho=float(rho),
        candidate_response=float(response),
        delta_theta=float(delta_theta),
        theta_eff_preclip=float(preclip),
        theta_eff=float(theta_eff),
        was_clipped=was_clipped,
    )


def scattering_coefficients(theta_eff: float) -> tuple[complex, complex]:
    t_eff = np.exp(1j * theta_eff) * math.cos(theta_eff)
    r_eff = -1j * np.exp(1j * theta_eff) * math.sin(theta_eff)
    return complex(t_eff), complex(r_eff)
