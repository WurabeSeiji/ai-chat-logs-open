#!/usr/bin/env python3
"""Independent Stage B reproduction of the current System A scattering kernel.

The audited System A/System B source files are never imported or executed.
Their already-audited formulas are independently transcribed here.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np


SEMANTIC_NOTICE = (
    "B_to_A_transfer is spectral cosine similarity of the A-channel state "
    "to the initial B spectrum; it is NOT a path-exchange norm."
)

AUDIT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = AUDIT_ROOT / "data" / "stage_B"
PATH_ARRAY_ROOT = DATA_ROOT / "path_arrays"


@dataclass(frozen=True)
class Config:
    chi_grid_n: int = 512
    eta_grid_n: int = 16
    high_n: int = 63
    k_components: int = 4
    mix_p: float = 0.5
    mix_phi: float = math.pi / 2.0
    reflection_rate: float = 0.6971778791282474
    max_collision: int = 32
    q_A: float = 1.0
    q_B: float = -1.0
    m_A: int = 1
    m_B: int = 2
    p0: float = 1.0


@dataclass(frozen=True)
class Case:
    case_id: str
    wave_A: str
    wave_B: str
    type_A: str
    type_B: str


CASES = (
    Case("F1_x_F1", "F1", "F1", "odd", "odd"),
    Case("FK_x_FK", "FK", "FK", "odd", "odd"),
    Case("BK_x_BK", "BK", "BK", "even", "even"),
    Case("FK_x_BK", "FK", "BK", "odd", "even"),
    Case("MIX_x_MIX", "MIX", "MIX", "mixed", "mixed"),
)


def norm2(vector: np.ndarray) -> float:
    return float(np.vdot(vector, vector).real)


def normalize(vector: np.ndarray) -> np.ndarray:
    value = math.sqrt(max(norm2(vector), 0.0))
    if value <= 0.0:
        raise ValueError("zero norm state")
    return vector / value


def grids(config: Config) -> tuple[np.ndarray, np.ndarray]:
    u = np.linspace(-math.pi, math.pi, config.chi_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, config.eta_grid_n, endpoint=False)
    return u, eta


def wave_kernels(config: Config, u: np.ndarray) -> dict[str, np.ndarray]:
    odd = tuple(2 * index + 1 for index in range(config.k_components))
    even = tuple(2 * index for index in range(1, config.k_components + 1))
    f1 = np.cos(u).astype(complex)
    fk = sum((np.cos(n * u) for n in odd), start=np.zeros_like(u)) / config.k_components
    bk = sum((np.cos(n * u) for n in even), start=np.zeros_like(u)) / config.k_components
    mixed = (
        math.sqrt(config.mix_p) * fk
        + np.exp(1j * config.mix_phi) * math.sqrt(1.0 - config.mix_p) * bk
    )
    return {
        "F1": np.asarray(f1, dtype=complex),
        "FK": np.asarray(fk, dtype=complex),
        "BK": np.asarray(bk, dtype=complex),
        "MIX": np.asarray(mixed, dtype=complex),
    }


def make_state(
    kernel: np.ndarray,
    q: float,
    m: int,
    config: Config,
    u: np.ndarray,
    eta: np.ndarray,
) -> np.ndarray:
    carrier = np.exp(1j * q * config.p0 * u)
    eta_phase = np.exp(1j * m * eta)
    state = (kernel * carrier)[:, None] * eta_phase[None, :]
    return normalize(state.reshape(-1))


def scattering_coefficients(reflection_rate: float) -> tuple[float, complex, complex, float, float]:
    if not 0.0 <= reflection_rate <= 1.0:
        raise ValueError("reflection_rate must be in [0, 1]")
    delta = 2.0 * math.asin(math.sqrt(reflection_rate))
    t = np.exp(0.5j * delta) * math.cos(0.5 * delta)
    r = -1j * np.exp(0.5j * delta) * math.sin(0.5 * delta)
    return delta, complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


def reshape(vector: np.ndarray, config: Config) -> np.ndarray:
    return vector.reshape(config.chi_grid_n, config.eta_grid_n)


def inverse_carrier_by_eta_lineage(
    vector: np.ndarray,
    config: Config,
    u: np.ndarray,
    eta: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Remove q by using the orthogonal eta mode as the source-lineage label."""
    array = reshape(vector, config)
    reconstructed = np.zeros_like(array, dtype=complex)
    decarried = np.zeros_like(array, dtype=complex)
    for m, q in ((config.m_A, config.q_A), (config.m_B, config.q_B)):
        eta_phase = np.exp(1j * m * eta)
        coefficient = np.sum(array * np.conjugate(eta_phase)[None, :], axis=1) / config.eta_grid_n
        component = coefficient[:, None] * eta_phase[None, :]
        reconstructed += component
        decarried += component * np.exp(-1j * q * config.p0 * u)[:, None]
    denominator = max(math.sqrt(norm2(vector)), 1.0e-300)
    residual = math.sqrt(norm2(array - reconstructed)) / denominator
    return decarried.reshape(-1), residual


def half_shift(vector: np.ndarray, config: Config) -> np.ndarray:
    array = reshape(vector, config)
    return np.roll(array, -config.chi_grid_n // 2, axis=0).reshape(-1)


def parity_metrics(vector: np.ndarray, config: Config) -> dict[str, float]:
    denominator = norm2(vector)
    if denominator <= 0.0:
        raise ValueError("zero norm parity input")
    shifted = half_shift(vector, config)
    correlation = np.vdot(vector, shifted)
    even = 0.5 * (vector + shifted)
    odd = 0.5 * (vector - shifted)
    return {
        "C_pi_raw_real": float(correlation.real),
        "C_pi_raw_imag": float(correlation.imag),
        "c_pi": float(correlation.real / denominator),
        "p_B": float(norm2(even) / denominator),
        "p_F": float(norm2(odd) / denominator),
    }


def chi_density(vector: np.ndarray, config: Config) -> np.ndarray:
    array = reshape(vector, config)
    density = np.sum(np.abs(array) ** 2, axis=1)
    total = float(np.sum(density))
    if total <= 0.0:
        raise ValueError("zero density")
    return density / total


def localization_metrics(vector: np.ndarray, config: Config) -> dict[str, float]:
    probability = np.abs(vector) ** 2 / norm2(vector)
    density = chi_density(vector, config)
    du = 2.0 * math.pi / config.chi_grid_n
    return {
        "localization_L": float(np.sum(probability**2)),
        "localization_width_ipr": float(du / np.sum(density**2)),
        "chi_peak_density": float(np.max(density)),
    }


def harmonic_distribution(vector: np.ndarray, config: Config) -> dict[int, float]:
    denominator = norm2(vector)
    array = reshape(vector, config)
    transformed = np.fft.fft(array, axis=0, norm="ortho")
    power = np.sum(np.abs(transformed) ** 2, axis=1)
    frequencies = np.fft.fftfreq(config.chi_grid_n, d=1.0 / config.chi_grid_n)
    rounded = np.rint(frequencies).astype(int)
    max_n = min(config.chi_grid_n // 2, config.high_n + 2)
    raw: dict[int, float] = {}
    total = 0.0
    for n_abs in range(max_n + 1):
        indices = np.flatnonzero(np.abs(rounded) == n_abs)
        amount = float(np.sum(power[indices]).real) / denominator if indices.size else 0.0
        if amount > 1.0e-14:
            raw[n_abs] = max(amount, 0.0)
            total += raw[n_abs]
    if total <= 0.0:
        return {0: 1.0}
    return {key: value / total for key, value in raw.items()}


def distribution_similarity(a: dict[int, float], b: dict[int, float]) -> float:
    keys = sorted(set(a) | set(b))
    va = np.asarray([a.get(key, 0.0) for key in keys], dtype=float)
    vb = np.asarray([b.get(key, 0.0) for key in keys], dtype=float)
    denominator = float(np.linalg.norm(va) * np.linalg.norm(vb))
    if denominator <= 0.0:
        return float("nan")
    return float(max(0.0, min(1.0, np.dot(va, vb) / denominator)))


def state_metrics(
    vector: np.ndarray,
    config: Config,
    u: np.ndarray,
    eta: np.ndarray,
) -> dict[str, float]:
    decarried, residual = inverse_carrier_by_eta_lineage(vector, config, u, eta)
    values = {
        "norm2": norm2(vector),
        "inverse_carrier_lineage_residual": residual,
    }
    values.update(localization_metrics(vector, config))
    values.update(parity_metrics(decarried, config))
    return values


def prefixed(prefix: str, values: dict[str, float]) -> dict[str, float]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def overlap_phase(reference: np.ndarray, vector: np.ndarray) -> float:
    value = np.vdot(reference, vector)
    if abs(value) <= 1.0e-14:
        return float("nan")
    return float(np.angle(value))


def kernel_control_metrics(
    name: str,
    kernel: np.ndarray,
    component_count: int,
    type_label: str,
    config: Config,
) -> dict[str, float | int | str]:
    normalized = normalize(kernel)
    power = np.abs(np.fft.fft(normalized, norm="ortho")) ** 2
    frequencies = np.fft.fftfreq(config.chi_grid_n, d=1.0 / config.chi_grid_n)
    total = float(np.sum(power))
    spectral_rms = math.sqrt(float(np.sum((frequencies**2) * power) / total))
    active = np.flatnonzero(power > max(total * 1.0e-12, 1.0e-15))
    highest = float(np.max(np.abs(frequencies[active]))) if active.size else 0.0
    density = np.abs(normalized) ** 2
    du = 2.0 * math.pi / config.chi_grid_n
    return {
        "wave": name,
        "type": type_label,
        "sample_count": config.chi_grid_n,
        "spatial_domain": "[-pi, pi)",
        "component_count": component_count,
        "kernel_norm2": norm2(normalized),
        "peak_amplitude": float(np.max(np.abs(normalized))),
        "RMS_amplitude": float(np.sqrt(np.mean(np.abs(normalized) ** 2))),
        "localization_width_ipr": float(du / np.sum(density**2)),
        "spectral_RMS_wavenumber": spectral_rms,
        "highest_wavenumber": highest,
        "phase_origin": 0.0,
        "carrier_treatment": (
            "A:q=+1,m=1; B:q=-1,m=2; inverse carrier resolved by eta lineage before parity"
        ),
    }


def scatter_pair(
    a: np.ndarray,
    b: np.ndarray,
    r: complex,
    t: complex,
) -> tuple[np.ndarray, np.ndarray]:
    return normalize(r * a + t * b), normalize(t * a + r * b)


def commutator_response(
    a: np.ndarray,
    b: np.ndarray,
    r: complex,
    t: complex,
    config: Config,
    u: np.ndarray,
    eta: np.ndarray,
) -> float:
    a_kernel, _ = inverse_carrier_by_eta_lineage(a, config, u, eta)
    b_kernel, _ = inverse_carrier_by_eta_lineage(b, config, u, eta)
    p_a = half_shift(a_kernel, config)
    p_b = half_shift(b_kernel, config)
    left_a, left_b = scatter_pair(p_a, p_b, r, t)
    out_a, out_b = scatter_pair(a_kernel, b_kernel, r, t)
    right_a = half_shift(out_a, config)
    right_b = half_shift(out_b, config)
    numerator = math.sqrt(norm2(left_a - right_a) + norm2(left_b - right_b))
    denominator = math.sqrt(norm2(a_kernel) + norm2(b_kernel))
    return float(numerator / denominator)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def range_across_cases(rows: Iterable[dict], key: str) -> float:
    by_collision: dict[int, list[float]] = {}
    for row in rows:
        by_collision.setdefault(int(row["collision"]), []).append(float(row[key]))
    return max(max(values) - min(values) for values in by_collision.values())


def main() -> int:
    config = Config()
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    PATH_ARRAY_ROOT.mkdir(parents=True, exist_ok=True)
    u, eta = grids(config)
    kernels = wave_kernels(config, u)
    delta, t, r, T, R = scattering_coefficients(config.reflection_rate)

    kernel_meta = {
        "F1": (1, "odd"),
        "FK": (config.k_components, "odd"),
        "BK": (config.k_components, "even"),
        "MIX": (2 * config.k_components, "mixed"),
    }
    wave_controls = {
        name: kernel_control_metrics(
            name,
            kernel,
            kernel_meta[name][0],
            kernel_meta[name][1],
            config,
        )
        for name, kernel in kernels.items()
    }
    control_rows: list[dict] = []
    transition_rows: list[dict] = []
    state_rows: list[dict] = []
    case_summaries: list[dict] = []
    npz_outputs: list[Path] = []

    for case in CASES:
        control_a = wave_controls[case.wave_A]
        control_b = wave_controls[case.wave_B]
        control_rows.append(
            {
                "case_id": case.case_id,
                "wave_A": case.wave_A,
                "wave_B": case.wave_B,
                "type_A": case.type_A,
                "type_B": case.type_B,
                "sample_count_A": control_a["sample_count"],
                "sample_count_B": control_b["sample_count"],
                "spatial_domain_A": control_a["spatial_domain"],
                "spatial_domain_B": control_b["spatial_domain"],
                "component_count_A": control_a["component_count"],
                "component_count_B": control_b["component_count"],
                "combined_initial_norm2": 2.0,
                "peak_amplitude_A": control_a["peak_amplitude"],
                "peak_amplitude_B": control_b["peak_amplitude"],
                "peak_amplitude_difference": abs(
                    float(control_a["peak_amplitude"]) - float(control_b["peak_amplitude"])
                ),
                "RMS_amplitude_A": control_a["RMS_amplitude"],
                "RMS_amplitude_B": control_b["RMS_amplitude"],
                "RMS_amplitude_difference": abs(
                    float(control_a["RMS_amplitude"]) - float(control_b["RMS_amplitude"])
                ),
                "localization_width_ipr_A": control_a["localization_width_ipr"],
                "localization_width_ipr_B": control_b["localization_width_ipr"],
                "localization_width_difference": abs(
                    float(control_a["localization_width_ipr"])
                    - float(control_b["localization_width_ipr"])
                ),
                "spectral_RMS_wavenumber_A": control_a["spectral_RMS_wavenumber"],
                "spectral_RMS_wavenumber_B": control_b["spectral_RMS_wavenumber"],
                "spectral_RMS_wavenumber_difference": abs(
                    float(control_a["spectral_RMS_wavenumber"])
                    - float(control_b["spectral_RMS_wavenumber"])
                ),
                "highest_wavenumber_A": control_a["highest_wavenumber"],
                "highest_wavenumber_B": control_b["highest_wavenumber"],
                "highest_wavenumber_difference": abs(
                    float(control_a["highest_wavenumber"])
                    - float(control_b["highest_wavenumber"])
                ),
                "carrier_treatment_A": "q=+1,m=1; inverse by eta lineage",
                "carrier_treatment_B": "q=-1,m=2; inverse by eta lineage",
                "phase_origin_A": 0.0,
                "phase_origin_B": 0.0,
                "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
            }
        )

        a = make_state(kernels[case.wave_A], config.q_A, config.m_A, config, u, eta)
        b = make_state(kernels[case.wave_B], config.q_B, config.m_B, config, u, eta)
        initial_a = a.copy()
        initial_b = b.copy()
        initial_h_a = harmonic_distribution(initial_a, config)
        initial_h_b = harmonic_distribution(initial_b, config)

        def append_state(collision: int, channel: str, vector: np.ndarray) -> None:
            metrics = state_metrics(vector, config, u, eta)
            distribution = harmonic_distribution(vector, config)
            row = {
                "case_id": case.case_id,
                "collision": collision,
                "channel": channel,
                "wave_A": case.wave_A,
                "wave_B": case.wave_B,
                "type_A": case.type_A,
                "type_B": case.type_B,
                **metrics,
                "sim_to_A0": distribution_similarity(distribution, initial_h_a),
                "sim_to_B0": distribution_similarity(distribution, initial_h_b),
                "B_to_A_transfer": (
                    distribution_similarity(distribution, initial_h_b)
                    if channel == "A_channel"
                    else ""
                ),
                "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
            }
            state_rows.append(row)

        append_state(0, "A_channel", a)
        append_state(0, "B_channel", b)
        path_buffers: dict[str, list[np.ndarray]] = {
            "path_a_to_a": [],
            "path_b_to_a": [],
            "path_b_to_b": [],
            "path_a_to_b": [],
            "interference_in_a_density": [],
            "interference_in_b_density": [],
        }
        case_transition_start = len(transition_rows)

        for collision in range(1, config.max_collision + 1):
            a_before = a
            b_before = b
            path_a_to_a = r * a_before
            path_b_to_a = t * b_before
            path_b_to_b = r * b_before
            path_a_to_b = t * a_before
            interference_a_density = 2.0 * np.real(
                np.conjugate(path_a_to_a) * path_b_to_a
            )
            interference_b_density = 2.0 * np.real(
                np.conjugate(path_b_to_b) * path_a_to_b
            )
            interference_a = float(np.sum(interference_a_density))
            interference_b = float(np.sum(interference_b_density))
            a_raw = path_a_to_a + path_b_to_a
            b_raw = path_b_to_b + path_a_to_b
            a_raw_norm2 = norm2(a_raw)
            b_raw_norm2 = norm2(b_raw)
            a = normalize(a_raw)
            b = normalize(b_raw)

            a_before_values = state_metrics(a_before, config, u, eta)
            b_before_values = state_metrics(b_before, config, u, eta)
            a_raw_values = state_metrics(a_raw, config, u, eta)
            b_raw_values = state_metrics(b_raw, config, u, eta)
            a_after_values = state_metrics(a, config, u, eta)
            b_after_values = state_metrics(b, config, u, eta)
            a_after_h = harmonic_distribution(a, config)
            b_after_h = harmonic_distribution(b, config)

            path_a_to_a_norm = norm2(path_a_to_a)
            path_b_to_a_norm = norm2(path_b_to_a)
            path_b_to_b_norm = norm2(path_b_to_b)
            path_a_to_b_norm = norm2(path_a_to_b)
            row = {
                "case_id": case.case_id,
                "collision": collision,
                "input_state_collision": collision - 1,
                "output_state_collision": collision,
                "wave_A": case.wave_A,
                "wave_B": case.wave_B,
                "type_A": case.type_A,
                "type_B": case.type_B,
                "R_input": config.reflection_rate,
                "R": R,
                "T": T,
                "Delta_F": delta,
                "r_real": r.real,
                "r_imag": r.imag,
                "t_real": t.real,
                "t_imag": t.imag,
                "phase_r": float(np.angle(r)),
                "phase_t": float(np.angle(t)),
                "phase_t_minus_r": float(np.angle(t / r)),
                "path_a_to_a_norm_raw": path_a_to_a_norm,
                "path_b_to_a_norm_raw": path_b_to_a_norm,
                "path_b_to_b_norm_raw": path_b_to_b_norm,
                "path_a_to_b_norm_raw": path_a_to_b_norm,
                "interference_in_a_raw": interference_a,
                "interference_in_b_raw": interference_b,
                "a_output_norm2_from_decomposition": (
                    path_a_to_a_norm + path_b_to_a_norm + interference_a
                ),
                "b_output_norm2_from_decomposition": (
                    path_b_to_b_norm + path_a_to_b_norm + interference_b
                ),
                "a_output_norm2_raw": a_raw_norm2,
                "b_output_norm2_raw": b_raw_norm2,
                "a_output_norm2_after_normalization": norm2(a),
                "b_output_norm2_after_normalization": norm2(b),
                "path_a_to_a_norm_after_channel_normalization": (
                    path_a_to_a_norm / a_raw_norm2
                ),
                "path_b_to_a_norm_after_channel_normalization": (
                    path_b_to_a_norm / a_raw_norm2
                ),
                "path_b_to_b_norm_after_channel_normalization": (
                    path_b_to_b_norm / b_raw_norm2
                ),
                "path_a_to_b_norm_after_channel_normalization": (
                    path_a_to_b_norm / b_raw_norm2
                ),
                "interference_in_a_after_channel_normalization": (
                    interference_a / a_raw_norm2
                ),
                "interference_in_b_after_channel_normalization": (
                    interference_b / b_raw_norm2
                ),
                "combined_input_norm2": norm2(a_before) + norm2(b_before),
                "combined_raw_output_norm2": a_raw_norm2 + b_raw_norm2,
                "combined_normalized_output_norm2": norm2(a) + norm2(b),
                **prefixed("a_input", a_before_values),
                **prefixed("b_input", b_before_values),
                **prefixed("a_raw", a_raw_values),
                **prefixed("b_raw", b_raw_values),
                **prefixed("a_normalized", a_after_values),
                **prefixed("b_normalized", b_after_values),
                "delta_a_p_B": a_after_values["p_B"] - a_before_values["p_B"],
                "delta_a_p_F": a_after_values["p_F"] - a_before_values["p_F"],
                "delta_b_p_B": b_after_values["p_B"] - b_before_values["p_B"],
                "delta_b_p_F": b_after_values["p_F"] - b_before_values["p_F"],
                "eta_commutator_current_kernel": commutator_response(
                    a_before, b_before, r, t, config, u, eta
                ),
                "a_phase_to_initial_A": overlap_phase(initial_a, a),
                "a_phase_to_initial_B": overlap_phase(initial_b, a),
                "b_phase_to_initial_A": overlap_phase(initial_a, b),
                "b_phase_to_initial_B": overlap_phase(initial_b, b),
                "B_to_A_transfer": distribution_similarity(a_after_h, initial_h_b),
                "A_to_B_spectral_similarity": distribution_similarity(
                    b_after_h, initial_h_a
                ),
                "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
            }
            transition_rows.append(row)
            append_state(collision, "A_channel", a)
            append_state(collision, "B_channel", b)
            for name, value in (
                ("path_a_to_a", path_a_to_a),
                ("path_b_to_a", path_b_to_a),
                ("path_b_to_b", path_b_to_b),
                ("path_a_to_b", path_a_to_b),
                ("interference_in_a_density", interference_a_density),
                ("interference_in_b_density", interference_b_density),
            ):
                path_buffers[name].append(reshape(value, config))

        case_rows = transition_rows[case_transition_start:]
        output_path = PATH_ARRAY_ROOT / f"{case.case_id}_collision_paths.npz"
        np.savez_compressed(
            output_path,
            collisions=np.arange(1, config.max_collision + 1, dtype=np.int64),
            path_a_to_a=np.stack(path_buffers["path_a_to_a"]),
            path_b_to_a=np.stack(path_buffers["path_b_to_a"]),
            path_b_to_b=np.stack(path_buffers["path_b_to_b"]),
            path_a_to_b=np.stack(path_buffers["path_a_to_b"]),
            interference_in_a_density=np.stack(path_buffers["interference_in_a_density"]),
            interference_in_b_density=np.stack(path_buffers["interference_in_b_density"]),
            case_id=np.asarray(case.case_id),
            R=np.asarray(R),
            T=np.asarray(T),
            B_to_A_transfer_semantics=np.asarray(SEMANTIC_NOTICE),
        )
        npz_outputs.append(output_path)
        case_summaries.append(
            {
                "case_id": case.case_id,
                "wave_A": case.wave_A,
                "wave_B": case.wave_B,
                "type_A": case.type_A,
                "type_B": case.type_B,
                "max_abs_interference_in_a": max(
                    abs(float(row["interference_in_a_raw"])) for row in case_rows
                ),
                "max_abs_interference_in_b": max(
                    abs(float(row["interference_in_b_raw"])) for row in case_rows
                ),
                "max_path_decomposition_error": max(
                    max(
                        abs(
                            float(row["a_output_norm2_from_decomposition"])
                            - float(row["a_output_norm2_raw"])
                        ),
                        abs(
                            float(row["b_output_norm2_from_decomposition"])
                            - float(row["b_output_norm2_raw"])
                        ),
                    )
                    for row in case_rows
                ),
                "max_raw_channel_norm2_error_from_one": max(
                    max(
                        abs(float(row["a_output_norm2_raw"]) - 1.0),
                        abs(float(row["b_output_norm2_raw"]) - 1.0),
                    )
                    for row in case_rows
                ),
                "max_normalized_channel_norm2_error_from_one": max(
                    max(
                        abs(float(row["a_output_norm2_after_normalization"]) - 1.0),
                        abs(float(row["b_output_norm2_after_normalization"]) - 1.0),
                    )
                    for row in case_rows
                ),
                "max_eta_commutator_current_kernel": max(
                    float(row["eta_commutator_current_kernel"]) for row in case_rows
                ),
                "B_to_A_transfer_min": min(
                    float(row["B_to_A_transfer"]) for row in case_rows
                ),
                "B_to_A_transfer_max": max(
                    float(row["B_to_A_transfer"]) for row in case_rows
                ),
                "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
            }
        )

    baseline_path = DATA_ROOT / "current_behavior_baseline.csv"
    state_path = DATA_ROOT / "state_parity_metrics.csv"
    controls_path = DATA_ROOT / "input_control_metrics.csv"
    write_csv(baseline_path, transition_rows)
    write_csv(state_path, state_rows)
    write_csv(controls_path, control_rows)

    maximums = {
        "path_a_to_a_norm_cross_case_spread": range_across_cases(
            transition_rows, "path_a_to_a_norm_raw"
        ),
        "path_b_to_a_norm_cross_case_spread": range_across_cases(
            transition_rows, "path_b_to_a_norm_raw"
        ),
        "path_b_to_b_norm_cross_case_spread": range_across_cases(
            transition_rows, "path_b_to_b_norm_raw"
        ),
        "path_a_to_b_norm_cross_case_spread": range_across_cases(
            transition_rows, "path_a_to_b_norm_raw"
        ),
        "raw_channel_norm2_error_from_one": max(
            max(
                abs(float(row["a_output_norm2_raw"]) - 1.0),
                abs(float(row["b_output_norm2_raw"]) - 1.0),
            )
            for row in transition_rows
        ),
        "normalized_channel_norm2_error_from_one": max(
            max(
                abs(float(row["a_output_norm2_after_normalization"]) - 1.0),
                abs(float(row["b_output_norm2_after_normalization"]) - 1.0),
            )
            for row in transition_rows
        ),
        "path_decomposition_error": max(
            max(
                abs(
                    float(row["a_output_norm2_from_decomposition"])
                    - float(row["a_output_norm2_raw"])
                ),
                abs(
                    float(row["b_output_norm2_from_decomposition"])
                    - float(row["b_output_norm2_raw"])
                ),
            )
            for row in transition_rows
        ),
        "absolute_interference": max(
            max(
                abs(float(row["interference_in_a_raw"])),
                abs(float(row["interference_in_b_raw"])),
            )
            for row in transition_rows
        ),
        "commutator_response": max(
            float(row["eta_commutator_current_kernel"]) for row in transition_rows
        ),
        "inverse_carrier_lineage_residual": max(
            max(
                float(row["a_normalized_inverse_carrier_lineage_residual"]),
                float(row["b_normalized_inverse_carrier_lineage_residual"]),
            )
            for row in transition_rows
        ),
    }

    state_by_case = {
        case.case_id: [row for row in state_rows if row["case_id"] == case.case_id]
        for case in CASES
    }
    parity_checks = {
        "F1_x_F1_max_even_leakage": max(
            float(row["p_B"]) for row in state_by_case["F1_x_F1"]
        ),
        "FK_x_FK_max_even_leakage": max(
            float(row["p_B"]) for row in state_by_case["FK_x_FK"]
        ),
        "BK_x_BK_max_odd_leakage": max(
            float(row["p_F"]) for row in state_by_case["BK_x_BK"]
        ),
        "MIX_x_MIX_max_half_weight_error": max(
            max(abs(float(row["p_B"]) - 0.5), abs(float(row["p_F"]) - 0.5))
            for row in state_by_case["MIX_x_MIX"]
        ),
    }
    fk_bk_by_collision: dict[int, list[dict]] = {}
    for row in state_by_case["FK_x_BK"]:
        fk_bk_by_collision.setdefault(int(row["collision"]), []).append(row)
    parity_checks["FK_x_BK_pair_total_sector_error"] = max(
        max(
            abs(sum(float(row["p_B"]) for row in rows) - 1.0),
            abs(sum(float(row["p_F"]) for row in rows) - 1.0),
        )
        for rows in fk_bk_by_collision.values()
    )
    parity_checks["c_pi_equals_p_B_minus_p_F_error"] = max(
        abs(float(row["c_pi"]) - (float(row["p_B"]) - float(row["p_F"])))
        for row in state_rows
    )
    tolerance = 2.0e-12
    assertions = {
        "all_path_norms_case_independent": max(
            maximums["path_a_to_a_norm_cross_case_spread"],
            maximums["path_b_to_a_norm_cross_case_spread"],
            maximums["path_b_to_b_norm_cross_case_spread"],
            maximums["path_a_to_b_norm_cross_case_spread"],
        )
        <= tolerance,
        "raw_channel_norms_conserved": maximums["raw_channel_norm2_error_from_one"]
        <= tolerance,
        "normalized_channel_norms_equal_one": maximums[
            "normalized_channel_norm2_error_from_one"
        ]
        <= tolerance,
        "path_decomposition_identity_holds": maximums["path_decomposition_error"]
        <= tolerance,
        "eta_orthogonality_removes_path_interference": maximums[
            "absolute_interference"
        ]
        <= tolerance,
        "current_kernel_commutes_with_half_shift": maximums["commutator_response"]
        <= tolerance,
        "inverse_carrier_lineage_projection_complete": maximums[
            "inverse_carrier_lineage_residual"
        ]
        <= tolerance,
        "pure_and_mixed_parity_weights_preserved": max(parity_checks.values())
        <= tolerance,
    }
    if not all(assertions.values()):
        failed = [name for name, passed in assertions.items() if not passed]
        raise AssertionError(f"Stage B reproduction assertions failed: {failed}")

    output_files = [baseline_path, state_path, controls_path, *npz_outputs]
    diagnostics = {
        "schema": "current_System_A_stage_B_reproduction_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_policy": (
            "Independent reproduction only. Audited System A/System B originals "
            "were not imported, executed, or modified."
        ),
        "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
        "config": asdict(config),
        "cases": [asdict(case) for case in CASES],
        "formula_provenance": {
            "System_A_SHA256": (
                "91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec"
            ),
            "scattering_source_SHA256": (
                "f815320f5632ae1b23ccade3a53b01e9110ad770da8407e2b10fa6065ef1695c"
            ),
            "formula": {
                "delta": "2*asin(sqrt(R_input))",
                "t": "exp(i*delta/2)*cos(delta/2)",
                "r": "-i*exp(i*delta/2)*sin(delta/2)",
                "a_raw": "r*a + t*b",
                "b_raw": "t*a + r*b",
                "normalization": "each output channel independently normalized",
            },
        },
        "coefficients": {
            "delta": delta,
            "r": {"real": r.real, "imag": r.imag, "phase": float(np.angle(r))},
            "t": {"real": t.real, "imag": t.imag, "phase": float(np.angle(t))},
            "R": R,
            "T": T,
        },
        "input_wave_controls": wave_controls,
        "case_summaries": case_summaries,
        "maximum_errors_and_spreads": maximums,
        "parity_checks": parity_checks,
        "assertions": assertions,
        "result_interpretation": {
            "type_preservation": (
                "The current kernel commutes with half-period shift and preserves "
                "the total parity sectors represented in the two-channel state."
            ),
            "type_dependent_response": (
                "No type-dependent path-norm response was generated; the same r,t "
                "acted on all cases."
            ),
            "mixed_channel_caution": (
                "For odd x even input, each named output channel becomes a mixture "
                "of preserved odd/even source sectors; this is channel mixing, not "
                "parity conversion."
            ),
        },
        "outputs": [
            {
                "path": path.relative_to(AUDIT_ROOT).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
            }
            for path in output_files
        ],
    }
    diagnostics_path = DATA_ROOT / "current_behavior_diagnostics.json"
    diagnostics_path.write_text(
        json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"B_to_A_transfer semantics: {SEMANTIC_NOTICE}")
    print(diagnostics_path)
    print(json.dumps(assertions, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
