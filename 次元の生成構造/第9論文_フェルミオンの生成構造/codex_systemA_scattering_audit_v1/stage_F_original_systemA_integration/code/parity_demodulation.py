"""Lineage-resolved demodulation used as the Stage F parity definition."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class LineageSpec:
    label: str
    q: float
    eta_mode: int


@dataclass(frozen=True)
class ParityReadout:
    correlation_raw: complex
    c_pi: float
    boson_weight: float
    fermion_weight: float
    projected_norm2: float
    state_norm2: float
    origin_a_weight: float
    origin_b_weight: float
    reconstruction_residual: float
    projection_sum_residual: float
    kernel_origin_a: np.ndarray
    kernel_origin_b: np.ndarray


ORIGIN_A = LineageSpec("origin_A", q=1.0, eta_mode=1)
ORIGIN_B = LineageSpec("origin_B", q=-1.0, eta_mode=2)


def norm2(vector: np.ndarray) -> float:
    return float(np.vdot(vector, vector).real)


def half_shift(kernel: np.ndarray) -> np.ndarray:
    if kernel.ndim != 1 or kernel.size % 2:
        raise ValueError("kernel must be a one-dimensional even-length array")
    return np.roll(kernel, -kernel.size // 2)


def eta_basis(spec: LineageSpec, eta: np.ndarray) -> np.ndarray:
    return np.exp(1j * spec.eta_mode * eta) / math.sqrt(eta.size)


def modulate_kernel(
    kernel: np.ndarray,
    spec: LineageSpec,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> np.ndarray:
    carrier = np.exp(1j * spec.q * p0 * u)
    return (
        (kernel * carrier)[:, None] * eta_basis(spec, eta)[None, :]
    ).reshape(-1)


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
    carrier_state = np.sum(
        array * np.conjugate(eta_basis(spec, eta))[None, :],
        axis=1,
    )
    return carrier_state * np.exp(-1j * spec.q * p0 * u)


def channel_parity(
    state: np.ndarray,
    u: np.ndarray,
    eta: np.ndarray,
    p0: float,
) -> ParityReadout:
    """Demodulate both persistent eta/carrier lineages before parity readout."""
    kernel_a = project_lineage_kernel(state, ORIGIN_A, u, eta, p0)
    kernel_b = project_lineage_kernel(state, ORIGIN_B, u, eta, p0)
    norm_a = norm2(kernel_a)
    norm_b = norm2(kernel_b)
    projected_norm = norm_a + norm_b
    state_norm = norm2(state)
    if projected_norm <= 0.0 or state_norm <= 0.0:
        raise ValueError("zero norm channel")

    correlation = 0.0j
    boson_numerator = 0.0
    fermion_numerator = 0.0
    for kernel, kernel_norm in ((kernel_a, norm_a), (kernel_b, norm_b)):
        if kernel_norm <= 0.0:
            continue
        shifted = half_shift(kernel)
        correlation += complex(np.vdot(kernel, shifted))
        boson_numerator += norm2(0.5 * (kernel + shifted))
        fermion_numerator += norm2(0.5 * (kernel - shifted))

    reconstructed = modulate_kernel(
        kernel_a, ORIGIN_A, u, eta, p0
    ) + modulate_kernel(kernel_b, ORIGIN_B, u, eta, p0)
    reconstruction_residual = math.sqrt(norm2(reconstructed - state)) / max(
        math.sqrt(state_norm), 1.0e-300
    )
    boson_weight = boson_numerator / projected_norm
    fermion_weight = fermion_numerator / projected_norm
    return ParityReadout(
        correlation_raw=correlation,
        c_pi=float(correlation.real / projected_norm),
        boson_weight=float(boson_weight),
        fermion_weight=float(fermion_weight),
        projected_norm2=float(projected_norm),
        state_norm2=float(state_norm),
        origin_a_weight=float(norm_a / projected_norm),
        origin_b_weight=float(norm_b / projected_norm),
        reconstruction_residual=float(reconstruction_residual),
        projection_sum_residual=float(
            max(
                abs(projected_norm - state_norm) / max(state_norm, 1.0e-300),
                abs(boson_weight + fermion_weight - 1.0),
            )
        ),
        kernel_origin_a=kernel_a,
        kernel_origin_b=kernel_b,
    )
