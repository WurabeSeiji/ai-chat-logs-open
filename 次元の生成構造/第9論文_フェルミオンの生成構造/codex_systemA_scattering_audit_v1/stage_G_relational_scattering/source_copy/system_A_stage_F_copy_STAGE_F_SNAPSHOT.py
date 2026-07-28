"""Stage F independent integration copy of the existing System A experiment.

Only the N_A=1, N_B=63 recursive System A path is reproduced here.  The
physical state remains the original 512 x 16 complex array.  C0 and reversed
Candidate 1 differ only through the effective angle.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Literal

import numpy as np

from parity_demodulation import channel_parity, norm2
from state_dependent_scattering import ScatteringMode, scatter_once


NormalizationMode = Literal["existing_normalization", "raw_update"]
KAPPA_VALUES = (0.01, 0.1, 1.0)
ORIGINAL_R_VALUES = (0.00, 0.51, 0.55, 0.60, 0.70, 0.90, 1.00)


@dataclass(frozen=True)
class Params:
    chi_grid_n: int = 512
    eta_grid_n: int = 16
    high_n: int = 63
    chi_center: float = 0.0
    p0: float = 1.0
    q_A: float = 1.0
    q_B: float = -1.0
    A_A: float = 1.0
    A_B: float = 1.0
    m_A: int = 1
    m_B: int = 2
    recursive_collision_count: int = 128
    r_sweep_values: tuple[float, ...] = ORIGINAL_R_VALUES


@dataclass(frozen=True)
class RunResult:
    rows: list[dict[str, Any]]
    initial_a: np.ndarray
    initial_b: np.ndarray
    final_a: np.ndarray
    final_b: np.ndarray
    snapshots: dict[int, tuple[np.ndarray, np.ndarray]]


def normalize(vector: np.ndarray) -> np.ndarray:
    value = np.linalg.norm(vector)
    if value <= 0.0:
        raise ValueError("zero norm state")
    return vector / value


def make_grids(params: Params) -> tuple[np.ndarray, np.ndarray]:
    chi = np.linspace(-math.pi, math.pi, params.chi_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, params.eta_grid_n, endpoint=False)
    return chi, eta


def odd_harmonic_kernel(u: np.ndarray, nh: int) -> np.ndarray:
    numerator = np.sin((nh + 1) * u)
    denominator = (nh + 1) * np.sin(u)
    out = np.empty_like(u, dtype=float)
    regular = np.abs(np.sin(u)) > 1.0e-12
    out[regular] = numerator[regular] / denominator[regular]
    if np.any(~regular):
        k = np.rint(u[~regular] / math.pi).astype(int)
        out[~regular] = np.where(k % 2 == 0, 1.0, -1.0)
    return out


def make_state(
    params: Params,
    n_chi: int,
    q: float,
    m: int,
    amplitude: float,
) -> np.ndarray:
    chi, eta = make_grids(params)
    kernel = odd_harmonic_kernel(chi - params.chi_center, n_chi)
    carrier = np.exp(1j * q * params.p0 * (chi - params.chi_center))
    eta_phase = np.exp(1j * m * eta)
    state = (kernel * carrier)[:, None] * eta_phase[None, :]
    return amplitude * normalize(state.reshape(-1))


def initial_state_pair(params: Params) -> tuple[np.ndarray, np.ndarray]:
    return (
        make_state(params, 1, params.q_A, params.m_A, params.A_A),
        make_state(
            params,
            params.high_n,
            params.q_B,
            params.m_B,
            params.A_B,
        ),
    )


def make_explicit_packet_state(
    params: Params,
    harmonics: tuple[int, ...],
    *,
    q: float,
    m: int,
) -> np.ndarray:
    chi, eta = make_grids(params)
    u = chi - params.chi_center
    kernel = sum(
        (np.cos(float(harmonic) * u) for harmonic in harmonics),
        start=np.zeros_like(u),
    ) / math.sqrt(len(harmonics))
    carrier = np.exp(1j * q * params.p0 * u)
    eta_phase = np.exp(1j * m * eta)
    state = (kernel * carrier)[:, None] * eta_phase[None, :]
    return normalize(state.reshape(-1))


def custom_31_initial_state_pair(
    params: Params,
) -> tuple[np.ndarray, np.ndarray]:
    """Existing 31-series condition: custom packet A=1, B=1+2."""
    return (
        make_explicit_packet_state(
            params, (1,), q=params.q_A, m=params.m_A
        ),
        make_explicit_packet_state(
            params, (1, 2), q=params.q_B, m=params.m_B
        ),
    )


def harmonic_distribution(
    params: Params, vector: np.ndarray
) -> dict[int, float]:
    denominator = norm2(vector)
    if denominator <= 0.0:
        return {0: 1.0}
    array = vector.reshape(params.chi_grid_n, params.eta_grid_n)
    transformed = np.fft.fft(array, axis=0, norm="ortho")
    frequencies = np.fft.fftfreq(
        params.chi_grid_n, d=1.0 / params.chi_grid_n
    )
    rounded = np.rint(frequencies).astype(int)
    max_n = min(params.chi_grid_n // 2, params.high_n + 2)
    raw: dict[int, float] = {}
    total = 0.0
    for n_abs in range(max_n + 1):
        indices = np.flatnonzero(np.abs(rounded) == n_abs)
        amount = (
            float(np.sum(np.abs(transformed[indices, :]) ** 2).real)
            / denominator
            if indices.size
            else 0.0
        )
        if amount > 1.0e-14:
            raw[n_abs] = max(amount, 0.0)
            total += raw[n_abs]
    if total <= 0.0:
        return {0: 1.0}
    return {key: value / total for key, value in raw.items()}


def effective_n(distribution: dict[int, float]) -> float:
    return float(
        sum(float(n) * float(weight) for n, weight in distribution.items())
    )


def localization(vector: np.ndarray) -> float:
    denominator = norm2(vector)
    if denominator <= 0.0:
        return float("nan")
    probability = np.abs(vector) ** 2 / denominator
    return float(np.sum(probability**2))


def distribution_similarity(
    left: dict[int, float], right: dict[int, float]
) -> float:
    keys = sorted(set(left) | set(right))
    a = np.asarray([float(left.get(key, 0.0)) for key in keys])
    b = np.asarray([float(right.get(key, 0.0)) for key in keys])
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denominator <= 0.0:
        return float("nan")
    return float(max(0.0, min(1.0, np.dot(a, b) / denominator)))


def state_metrics(
    params: Params,
    vector: np.ndarray,
    initial_h_a: dict[int, float],
    initial_h_b: dict[int, float],
) -> dict[str, float]:
    distribution = harmonic_distribution(params, vector)
    return {
        "L": localization(vector),
        "N_eff": effective_n(distribution),
        "similarity_to_initial_A": distribution_similarity(
            distribution, initial_h_a
        ),
        "similarity_to_initial_B": distribution_similarity(
            distribution, initial_h_b
        ),
        "norm2": norm2(vector),
    }


def state_series(
    params: Params,
    reflection_baseline: float,
    *,
    collision_count: int | None = None,
) -> list[dict[str, float | int]]:
    """Exact C0 existing-normalization series used by the local integration."""
    count = (
        params.recursive_collision_count
        if collision_count is None
        else int(collision_count)
    )
    a, b = initial_state_pair(params)
    initial_h_a = harmonic_distribution(params, a)
    initial_h_b = harmonic_distribution(params, b)
    theta0 = math.asin(math.sqrt(reflection_baseline))
    t = complex(np.exp(1j * theta0) * math.cos(theta0))
    r = complex(-1j * np.exp(1j * theta0) * math.sin(theta0))
    rows: list[dict[str, float | int]] = []
    for collision in range(count + 1):
        metrics_a = state_metrics(params, a, initial_h_a, initial_h_b)
        metrics_b = state_metrics(params, b, initial_h_a, initial_h_b)
        rows.append(
            {
                "collision_index": collision,
                "L_A": metrics_a["L"],
                "L_B": metrics_b["L"],
                "N_eff_A": metrics_a["N_eff"],
                "N_eff_B": metrics_b["N_eff"],
                "spectral_similarity_A_to_initial_A": metrics_a[
                    "similarity_to_initial_A"
                ],
                "spectral_similarity_A_to_initial_B": metrics_a[
                    "similarity_to_initial_B"
                ],
                "spectral_similarity_B_to_initial_A": metrics_b[
                    "similarity_to_initial_A"
                ],
                "spectral_similarity_B_to_initial_B": metrics_b[
                    "similarity_to_initial_B"
                ],
                "channel_norm_A": metrics_a["norm2"],
                "channel_norm_B": metrics_b["norm2"],
            }
        )
        if collision < count:
            a, b = normalize(r * a + t * b), normalize(t * a + r * b)
    return rows


def _finite_count(values: Iterable[Any]) -> int:
    count = 0
    for value in values:
        if isinstance(value, (int, float, np.integer, np.floating)):
            if not np.isfinite(float(value)):
                count += 1
    return count


def run_series(
    params: Params,
    *,
    reflection_baseline: float,
    scattering_mode: ScatteringMode,
    normalization_mode: NormalizationMode,
    kappa: float,
    collision_count: int | None = None,
    initial_state_override: tuple[np.ndarray, np.ndarray] | None = None,
    run_label: str | None = None,
    snapshot_indices: Iterable[int] = (),
) -> RunResult:
    if normalization_mode not in ("existing_normalization", "raw_update"):
        raise ValueError(f"unknown normalization mode: {normalization_mode}")
    count = (
        params.recursive_collision_count
        if collision_count is None
        else int(collision_count)
    )
    u, eta = make_grids(params)
    if initial_state_override is None:
        a, b = initial_state_pair(params)
    else:
        a, b = (
            np.asarray(initial_state_override[0], dtype=complex).copy(),
            np.asarray(initial_state_override[1], dtype=complex).copy(),
        )
    initial_a, initial_b = a.copy(), b.copy()
    initial_h_a = harmonic_distribution(params, initial_a)
    initial_h_b = harmonic_distribution(params, initial_b)
    rows: list[dict[str, Any]] = []
    run_id = (
        f"{scattering_mode}__{normalization_mode}__"
        f"kappa_{kappa:g}__R0_{reflection_baseline:g}"
    )
    if run_label:
        run_id = f"{run_label}__{run_id}"
    targets = {int(value) for value in snapshot_indices}
    snapshots: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    if 0 in targets:
        snapshots[0] = (a.copy(), b.copy())
    for collision_index in range(1, count + 1):
        collision = scatter_once(
            a,
            b,
            scattering_mode=scattering_mode,
            reflection_baseline=reflection_baseline,
            kappa=kappa,
            u=u,
            eta=eta,
            p0=params.p0,
        )
        raw_norm_a = norm2(collision.raw_a)
        raw_norm_b = norm2(collision.raw_b)
        if normalization_mode == "existing_normalization":
            next_a = normalize(collision.raw_a)
            next_b = normalize(collision.raw_b)
            scale_a = 1.0 / math.sqrt(raw_norm_a)
            scale_b = 1.0 / math.sqrt(raw_norm_b)
        else:
            next_a = collision.raw_a
            next_b = collision.raw_b
            scale_a = 1.0
            scale_b = 1.0

        next_parity_a = channel_parity(next_a, u, eta, params.p0)
        next_parity_b = channel_parity(next_b, u, eta, params.p0)
        metrics_a = state_metrics(
            params, next_a, initial_h_a, initial_h_b
        )
        metrics_b = state_metrics(
            params, next_b, initial_h_a, initial_h_b
        )
        angle = collision.angle
        row: dict[str, Any] = {
            "run_id": run_id,
            "collision_index": collision_index,
            "scattering_mode": scattering_mode,
            "normalization_mode": normalization_mode,
            "kappa": kappa,
            "R0": reflection_baseline,
            "theta0": angle.theta0,
            "rho": angle.rho,
            "c_A": collision.input_parity_a.c_pi,
            "c_B": collision.input_parity_b.c_pi,
            "c_mean": angle.c_mean,
            "delta_theta": angle.delta_theta,
            "theta_eff": angle.theta_eff,
            "R_eff": angle.R_eff,
            "T_eff": angle.T_eff,
            "L_A": metrics_a["L"],
            "L_B": metrics_b["L"],
            "N_eff_A": metrics_a["N_eff"],
            "N_eff_B": metrics_b["N_eff"],
            "spectral_similarity_A_to_initial_A": metrics_a[
                "similarity_to_initial_A"
            ],
            "spectral_similarity_A_to_initial_B": metrics_a[
                "similarity_to_initial_B"
            ],
            "spectral_similarity_B_to_initial_A": metrics_b[
                "similarity_to_initial_A"
            ],
            "spectral_similarity_B_to_initial_B": metrics_b[
                "similarity_to_initial_B"
            ],
            "B_to_A_transfer": metrics_a["similarity_to_initial_B"],
            "B_to_A_transfer_definition": (
                "spectral_similarity_to_initial_B; not path flux"
            ),
            "path_A_to_A_norm": norm2(collision.path_a_to_a),
            "path_B_to_A_norm": norm2(collision.path_b_to_a),
            "path_B_to_B_norm": norm2(collision.path_b_to_b),
            "path_A_to_B_norm": norm2(collision.path_a_to_b),
            "interference_A": collision.interference_a,
            "interference_B": collision.interference_b,
            "raw_norm_A": raw_norm_a,
            "raw_norm_B": raw_norm_b,
            "next_state_norm_A": norm2(next_a),
            "next_state_norm_B": norm2(next_b),
            "boson_weight_A": collision.input_parity_a.boson_weight,
            "fermion_weight_A": collision.input_parity_a.fermion_weight,
            "boson_weight_B": collision.input_parity_b.boson_weight,
            "fermion_weight_B": collision.input_parity_b.fermion_weight,
            "next_c_A": next_parity_a.c_pi,
            "next_c_B": next_parity_b.c_pi,
            "unitarity_residual": angle.unitarity_residual,
            "coefficient_orthogonality_residual": (
                angle.coefficient_orthogonality_residual
            ),
            "path_sum_residual_A": collision.path_sum_residual_a,
            "path_sum_residual_B": collision.path_sum_residual_b,
            "total_norm_conservation_residual": (
                collision.total_norm_residual
            ),
            "demodulation_reconstruction_residual": max(
                collision.input_parity_a.reconstruction_residual,
                collision.input_parity_b.reconstruction_residual,
                collision.raw_parity_a.reconstruction_residual,
                collision.raw_parity_b.reconstruction_residual,
                next_parity_a.reconstruction_residual,
                next_parity_b.reconstruction_residual,
            ),
            "parity_projection_sum_residual": max(
                collision.input_parity_a.projection_sum_residual,
                collision.input_parity_b.projection_sum_residual,
                collision.raw_parity_a.projection_sum_residual,
                collision.raw_parity_b.projection_sum_residual,
                next_parity_a.projection_sum_residual,
                next_parity_b.projection_sum_residual,
            ),
            "normalization_scale_A": scale_a,
            "normalization_scale_B": scale_b,
            "theta_range_violation": not (
                0.0 <= angle.theta_eff <= math.pi / 2.0
            ),
        }
        row["nan_inf_count"] = _finite_count(row.values())
        rows.append(row)
        a, b = next_a, next_b
        if collision_index in targets:
            snapshots[collision_index] = (a.copy(), b.copy())
    return RunResult(
        rows=rows,
        initial_a=initial_a,
        initial_b=initial_b,
        final_a=a,
        final_b=b,
        snapshots=snapshots,
    )
