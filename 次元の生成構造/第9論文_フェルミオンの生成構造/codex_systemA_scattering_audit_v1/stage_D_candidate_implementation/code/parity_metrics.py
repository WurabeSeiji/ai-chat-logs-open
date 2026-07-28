"""Two-layer modulation, demodulation, and parity measurements for Stage D."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class ModulationSpec:
    label: str
    q: float
    eta_mode: int
    p0: float = 1.0


@dataclass(frozen=True)
class ParityMeasurement:
    correlation_raw: complex
    indicator: float
    boson_weight: float
    fermion_weight: float
    norm2: float


@dataclass(frozen=True)
class DemodulationResult:
    kernel: np.ndarray
    reference_residual: float
    roundtrip_residual: float


def norm2(vector: np.ndarray) -> float:
    return float(np.vdot(vector, vector).real)


def normalize(vector: np.ndarray) -> np.ndarray:
    value = math.sqrt(max(norm2(vector), 0.0))
    if value <= 0.0:
        raise ValueError("zero norm state")
    return vector / value


def half_shift_kernel(kernel: np.ndarray) -> np.ndarray:
    if kernel.ndim != 1 or kernel.size % 2:
        raise ValueError("kernel must be a one-dimensional even-length array")
    return np.roll(kernel, -kernel.size // 2)


def parity_measurement(kernel: np.ndarray) -> ParityMeasurement:
    denominator = norm2(kernel)
    if denominator <= 0.0:
        raise ValueError("zero norm parity input")
    shifted = half_shift_kernel(kernel)
    correlation = complex(np.vdot(kernel, shifted))
    boson = 0.5 * (kernel + shifted)
    fermion = 0.5 * (kernel - shifted)
    return ParityMeasurement(
        correlation_raw=correlation,
        indicator=float(correlation.real / denominator),
        boson_weight=float(norm2(boson) / denominator),
        fermion_weight=float(norm2(fermion) / denominator),
        norm2=denominator,
    )


def combined_lineage_parity(
    kernels: Iterable[np.ndarray],
) -> ParityMeasurement:
    """Combine mutually orthogonal source lineages without cross interference."""
    kernels = tuple(kernels)
    if not kernels:
        raise ValueError("at least one lineage kernel is required")
    # At R=0 or R=1 one lineage coefficient is exactly zero.  Such a zero path
    # contributes neither norm nor correlation and is omitted from the direct
    # sum; the combined physical output itself remains nonzero.
    measurements = tuple(
        parity_measurement(kernel) for kernel in kernels if norm2(kernel) > 0.0
    )
    if not measurements:
        raise ValueError("all lineage kernels have zero norm")
    total_norm = sum(item.norm2 for item in measurements)
    if total_norm <= 0.0:
        raise ValueError("zero norm combined lineage state")
    correlation = sum(
        (item.correlation_raw for item in measurements), start=0.0j
    )
    boson_numerator = sum(
        item.boson_weight * item.norm2 for item in measurements
    )
    fermion_numerator = sum(
        item.fermion_weight * item.norm2 for item in measurements
    )
    return ParityMeasurement(
        correlation_raw=complex(correlation),
        indicator=float(correlation.real / total_norm),
        boson_weight=float(boson_numerator / total_norm),
        fermion_weight=float(fermion_numerator / total_norm),
        norm2=float(total_norm),
    )


def modulate_kernel(
    kernel: np.ndarray,
    u: np.ndarray,
    eta: np.ndarray,
    spec: ModulationSpec,
) -> np.ndarray:
    if kernel.ndim != 1 or kernel.shape != u.shape:
        raise ValueError("kernel and u must be matching one-dimensional arrays")
    carrier = np.exp(1j * spec.q * spec.p0 * u)
    eta_embedding = np.exp(1j * spec.eta_mode * eta) / math.sqrt(eta.size)
    return ((kernel * carrier)[:, None] * eta_embedding[None, :]).reshape(-1)


def demodulate_state(
    state: np.ndarray,
    expected_kernel: np.ndarray,
    u: np.ndarray,
    eta: np.ndarray,
    spec: ModulationSpec,
) -> DemodulationResult:
    expected_shape = (u.size, eta.size)
    if state.size != u.size * eta.size:
        raise ValueError("physical state size does not match u x eta grid")
    array = state.reshape(expected_shape)
    eta_embedding = np.exp(1j * spec.eta_mode * eta) / math.sqrt(eta.size)
    carrier_state = np.sum(
        array * np.conjugate(eta_embedding)[None, :], axis=1
    )
    kernel = carrier_state * np.exp(-1j * spec.q * spec.p0 * u)
    expected_norm = max(math.sqrt(norm2(expected_kernel)), 1.0e-300)
    state_norm = max(math.sqrt(norm2(state)), 1.0e-300)
    reference_residual = math.sqrt(norm2(kernel - expected_kernel)) / expected_norm
    reconstructed = modulate_kernel(kernel, u, eta, spec)
    roundtrip_residual = math.sqrt(norm2(reconstructed - state)) / state_norm
    return DemodulationResult(
        kernel=kernel,
        reference_residual=float(reference_residual),
        roundtrip_residual=float(roundtrip_residual),
    )


def lifted_half_shift(
    state: np.ndarray,
    expected_kernel: np.ndarray,
    u: np.ndarray,
    eta: np.ndarray,
    spec: ModulationSpec,
) -> tuple[np.ndarray, np.ndarray]:
    demodulated = demodulate_state(state, expected_kernel, u, eta, spec)
    shifted_kernel = half_shift_kernel(demodulated.kernel)
    shifted_state = modulate_kernel(shifted_kernel, u, eta, spec)
    return shifted_state, shifted_kernel


def lineage_reconstruction_residual(
    physical_paths: Iterable[np.ndarray],
    kernels: Iterable[np.ndarray],
    specs: Iterable[ModulationSpec],
    u: np.ndarray,
    eta: np.ndarray,
) -> float:
    physical_paths = tuple(physical_paths)
    kernels = tuple(kernels)
    specs = tuple(specs)
    if not (len(physical_paths) == len(kernels) == len(specs)):
        raise ValueError("lineage path, kernel, and spec counts differ")
    reconstructed = sum(
        (
            modulate_kernel(kernel, u, eta, spec)
            for kernel, spec in zip(kernels, specs)
        ),
        start=np.zeros_like(physical_paths[0], dtype=complex),
    )
    physical = sum(
        physical_paths,
        start=np.zeros_like(physical_paths[0], dtype=complex),
    )
    denominator = max(math.sqrt(norm2(physical)), 1.0e-300)
    return float(math.sqrt(norm2(reconstructed - physical)) / denominator)
