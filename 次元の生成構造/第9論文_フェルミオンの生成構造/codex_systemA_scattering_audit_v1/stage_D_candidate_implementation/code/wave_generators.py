"""Fixed Stage D kernel-wave generators."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from parity_metrics import norm2, normalize


@dataclass(frozen=True)
class WaveDefinition:
    label: str
    family: str
    k_components: int
    phase: float | None
    kernel: np.ndarray


def make_u_grid(sample_count: int = 512) -> np.ndarray:
    if sample_count <= 0 or sample_count % 2:
        raise ValueError("sample_count must be a positive even integer")
    return np.linspace(-math.pi, math.pi, sample_count, endpoint=False)


def make_eta_grid(sample_count: int = 16) -> np.ndarray:
    if sample_count <= 0:
        raise ValueError("eta sample_count must be positive")
    return np.linspace(-math.pi, math.pi, sample_count, endpoint=False)


def fermion_kernel(u: np.ndarray, k_components: int) -> np.ndarray:
    if k_components <= 0:
        raise ValueError("k_components must be positive")
    kernel = sum(
        (np.cos((2 * index + 1) * u) for index in range(k_components)),
        start=np.zeros_like(u),
    ) / math.sqrt(k_components)
    return normalize(np.asarray(kernel, dtype=complex))


def boson_kernel(u: np.ndarray, k_components: int) -> np.ndarray:
    if k_components <= 0:
        raise ValueError("k_components must be positive")
    kernel = sum(
        (np.cos(2 * index * u) for index in range(1, k_components + 1)),
        start=np.zeros_like(u),
    ) / math.sqrt(k_components)
    return normalize(np.asarray(kernel, dtype=complex))


def mixed_kernel(
    boson: np.ndarray,
    fermion: np.ndarray,
    phase: float,
) -> np.ndarray:
    return normalize(boson + np.exp(1j * phase) * fermion)


def wave_library(u: np.ndarray, k_components: int) -> dict[str, WaveDefinition]:
    fermion = fermion_kernel(u, k_components)
    boson = boson_kernel(u, k_components)
    definitions = {
        "F": WaveDefinition("F", "fermion", k_components, None, fermion),
        "B": WaveDefinition("B", "boson", k_components, None, boson),
    }
    for label, phase in (
        ("M0", 0.0),
        ("M90", math.pi / 2.0),
        ("M180", math.pi),
    ):
        definitions[label] = WaveDefinition(
            label,
            "mixed",
            k_components,
            phase,
            mixed_kernel(boson, fermion, phase),
        )
    return definitions


def wave_control_metrics(wave: WaveDefinition, u: np.ndarray) -> dict:
    kernel = wave.kernel
    probability = np.abs(kernel) ** 2 / norm2(kernel)
    du = 2.0 * math.pi / u.size
    transformed = np.fft.fft(kernel, norm="ortho")
    power = np.abs(transformed) ** 2
    frequencies = np.fft.fftfreq(u.size, d=1.0 / u.size)
    total_power = float(np.sum(power))
    active = np.flatnonzero(power > max(1.0e-15, total_power * 1.0e-12))
    highest = float(np.max(np.abs(frequencies[active]))) if active.size else 0.0
    return {
        "label": wave.label,
        "family": wave.family,
        "K": wave.k_components,
        "phase": wave.phase,
        "sample_count": u.size,
        "domain": "[-pi,pi)",
        "norm2": norm2(kernel),
        "peak_amplitude": float(np.max(np.abs(kernel))),
        "RMS_amplitude": float(np.sqrt(np.mean(np.abs(kernel) ** 2))),
        "localization_width_ipr": float(du / np.sum(probability**2)),
        "spectral_RMS_wavenumber": float(
            math.sqrt(np.sum((frequencies**2) * power) / total_power)
        ),
        "highest_wavenumber": highest,
    }
