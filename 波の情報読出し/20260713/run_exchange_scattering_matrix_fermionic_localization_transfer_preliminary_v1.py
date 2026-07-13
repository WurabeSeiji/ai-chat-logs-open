from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)
ACCELERATION_BASE_PATH = (
    BASE_DIR.parent
    / "20260711"
    / "ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v2"
    / "ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v2.json"
)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
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
    delta_f_fermion: float = math.pi
    delta_f_boson: float = 0.0
    delta_f_half: float = math.pi / 2.0
    p_tol: float = 1.0e-2
    d_q_tol: float = 1.0e-2
    name_tol: float = 1.0e-8
    norm_tol: float = 1.0e-10
    recursive_collision_count: int = 128
    r_sweep_values: Tuple[float, ...] = (0.00, 0.51, 0.55, 0.60, 0.70, 0.90, 1.00)


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


def normalize(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    if norm <= 0.0:
        raise ValueError("zero norm state")
    return v / norm


def wrap_angle(value: np.ndarray | float) -> np.ndarray | float:
    return (value + math.pi) % (2.0 * math.pi) - math.pi


def make_grids(params: Params) -> Tuple[np.ndarray, np.ndarray]:
    chi = np.linspace(-math.pi, math.pi, params.chi_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, params.eta_grid_n, endpoint=False)
    return chi, eta


def make_state(
    params: Params,
    n_chi: int,
    q: float,
    m: int,
    hair_enabled: bool,
    amplitude: float,
) -> np.ndarray:
    chi, eta = make_grids(params)
    chi_part = odd_harmonic_kernel(chi - params.chi_center, n_chi)
    phase_chi = np.exp(1j * q * params.p0 * (chi - params.chi_center))
    eta_phase = np.exp(1j * m * eta) if hair_enabled else np.ones_like(eta, dtype=complex)
    psi = (chi_part * phase_chi)[:, None] * eta_phase[None, :]
    return amplitude * normalize(psi.reshape(-1))


def reshape_state(params: Params, v: np.ndarray) -> np.ndarray:
    return v.reshape(params.chi_grid_n, params.eta_grid_n)


def chi_density(params: Params, v: np.ndarray) -> np.ndarray:
    psi = reshape_state(params, v)
    rho = np.sum(np.abs(psi) ** 2, axis=1)
    total = float(np.sum(rho))
    if total <= 0.0:
        raise ValueError("zero chi density")
    return rho / total


def inner(v: np.ndarray, w: np.ndarray) -> complex:
    return complex(np.vdot(v, w))


def norm2(v: np.ndarray) -> float:
    return float(np.vdot(v, v).real)


def normalized_distance(v: np.ndarray, w: np.ndarray) -> float:
    nv = math.sqrt(max(norm2(v), 0.0))
    nw = math.sqrt(max(norm2(w), 0.0))
    if nv <= 0.0 or nw <= 0.0:
        return float("inf")
    overlap = abs(inner(v / nv, w / nw))
    return math.sqrt(max(2.0 - 2.0 * overlap, 0.0))


def scattering_coefficients(delta_f: float) -> Tuple[complex, complex, float, float]:
    t = np.exp(0.5j * delta_f) * math.cos(0.5 * delta_f)
    r = -1j * np.exp(0.5j * delta_f) * math.sin(0.5 * delta_f)
    return complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


def delta_from_reflection_rate(r_value: float) -> float:
    if r_value < 0.0 or r_value > 1.0:
        raise ValueError(f"reflection rate must be in [0, 1]: {r_value}")
    return 2.0 * math.asin(math.sqrt(r_value))


def scattering_outputs(
    params: Params,
    n_a: int,
    n_b: int,
    hair_enabled: bool,
    delta_f: float,
) -> Dict[str, Any]:
    a_trans = make_state(params, n_a, params.q_A, params.m_A, hair_enabled, params.A_A)
    a_ref = make_state(params, n_a, -params.q_A, params.m_A, hair_enabled, params.A_A)
    b_trans = make_state(params, n_b, params.q_B, params.m_B, hair_enabled, params.A_B)
    b_ref = make_state(params, n_b, -params.q_B, params.m_B, hair_enabled, params.A_B)
    t, r, T, R = scattering_coefficients(delta_f)

    # Outgoing channels retain the two-channel structure.
    # minus: A reflected to the negative direction + B transmitted.
    # plus:  A transmitted + B reflected to the positive direction.
    out_minus = r * a_ref + t * b_trans
    out_plus = t * a_trans + r * b_ref

    return {
        "a_trans": a_trans,
        "a_ref": a_ref,
        "b_trans": b_trans,
        "b_ref": b_ref,
        "out_minus": out_minus,
        "out_plus": out_plus,
        "t": t,
        "r": r,
        "T": T,
        "R": R,
    }


def spectral_p_chi(params: Params, w: np.ndarray) -> np.ndarray:
    arr = reshape_state(params, w)
    freqs = np.fft.fftfreq(params.chi_grid_n, d=1.0 / params.chi_grid_n)
    transformed = np.fft.fft(arr, axis=0, norm="ortho")
    p_arr = np.fft.ifft(freqs[:, None] * transformed, axis=0, norm="ortho")
    return p_arr.reshape(-1)


def pure_expect_p(params: Params, vector: np.ndarray) -> float:
    denom = norm2(vector)
    if denom <= 0.0:
        return float("nan")
    return float(inner(vector, spectral_p_chi(params, vector)).real / denom)


def chi_frequency_components(params: Params, v: np.ndarray) -> np.ndarray:
    arr = reshape_state(params, v)
    return np.fft.fft(arr, axis=0, norm="ortho")


def eta_frequency_components(params: Params, v: np.ndarray) -> np.ndarray:
    arr = reshape_state(params, v)
    return np.fft.fft(arr, axis=1, norm="ortho")


def frequency_index(freqs: np.ndarray, target: int) -> List[int]:
    return [int(i) for i, freq in enumerate(freqs) if int(round(freq)) == target]


def harmonic_distribution(params: Params, vector: np.ndarray) -> Dict[int, float]:
    denom = norm2(vector)
    if denom <= 0.0:
        return {0: 1.0}
    fv = chi_frequency_components(params, vector)
    freqs = np.fft.fftfreq(params.chi_grid_n, d=1.0 / params.chi_grid_n)
    max_n = min(params.chi_grid_n // 2, params.high_n + 2)
    raw: Dict[int, float] = {}
    total = 0.0
    for n_abs in range(max_n + 1):
        indices = frequency_index(freqs, n_abs)
        if n_abs != 0:
            indices += frequency_index(freqs, -n_abs)
        amount = 0.0
        for idx in indices:
            amount += float(np.vdot(fv[idx, :], fv[idx, :]).real)
        amount = max(amount / denom, 0.0)
        if amount > 1.0e-14:
            raw[n_abs] = amount
            total += amount
    if total <= 0.0:
        return {0: 1.0}
    return {k: v / total for k, v in raw.items()}


def eta_distribution(params: Params, vector: np.ndarray, modes: Iterable[int]) -> Dict[int, float]:
    denom = norm2(vector)
    if denom <= 0.0:
        return {mode: 0.0 for mode in modes}
    fv = eta_frequency_components(params, vector)
    freqs = np.fft.fftfreq(params.eta_grid_n, d=1.0 / params.eta_grid_n)
    raw: Dict[int, float] = {}
    total = 0.0
    for mode in modes:
        amount = 0.0
        for idx in frequency_index(freqs, mode):
            amount += float(np.vdot(fv[:, idx], fv[:, idx]).real)
        amount = max(amount / denom, 0.0)
        raw[mode] = amount
        total += amount
    if total <= 0.0:
        return {mode: 0.0 for mode in modes}
    return {k: v / total for k, v in raw.items()}


def effective_n(distribution: Dict[int, float]) -> Tuple[float, float]:
    n_eff = sum(float(n) * float(weight) for n, weight in distribution.items())
    n_eff_2 = math.sqrt(sum(float(n * n) * float(weight) for n, weight in distribution.items()))
    return n_eff, n_eff_2


def l1_distance(a: Dict[int, float], b: Dict[int, float]) -> float:
    keys = set(a) | set(b)
    return float(sum(abs(float(a.get(k, 0.0)) - float(b.get(k, 0.0))) for k in keys))


def localization(vector: np.ndarray) -> float:
    denom = norm2(vector)
    if denom <= 0.0:
        return float("nan")
    prob = np.abs(vector) ** 2 / denom
    return float(np.sum(prob**2))


def chi_position_metrics(params: Params, vector: np.ndarray, n_chi: int) -> Dict[str, float]:
    denom = norm2(vector)
    if denom <= 0.0:
        return {
            "chi_peak": float("nan"),
            "chi_peak_abs_error": float("nan"),
            "chi_peak_contrast": float("nan"),
            "chi_center_cell_mass": float("nan"),
            "chi_effective_cell_count": float("nan"),
        }
    chi, _ = make_grids(params)
    arr = reshape_state(params, vector)
    prob_chi = np.sum(np.abs(arr) ** 2, axis=1) / denom
    peak_index = int(np.argmax(prob_chi))
    chi_peak = float(chi[peak_index])
    mean_cell = 1.0 / float(params.chi_grid_n)
    peak_contrast = float(prob_chi[peak_index] / mean_cell)
    delta = np.asarray(wrap_angle(chi - params.chi_center), dtype=float)
    cell_half_width = math.pi / float(n_chi + 1)
    center_cell_mass = float(np.sum(prob_chi[np.abs(delta) <= cell_half_width]))
    entropy = -float(np.sum([p * math.log(max(float(p), 1.0e-300)) for p in prob_chi]))
    return {
        "chi_peak": chi_peak,
        "chi_peak_abs_error": abs(float(wrap_angle(chi_peak - params.chi_center))),
        "chi_peak_contrast": peak_contrast,
        "chi_center_cell_mass": center_cell_mass,
        "chi_effective_cell_count": float(math.exp(entropy)),
    }


def name_readout_metrics(
    params: Params,
    p_m_a: float,
    p_m_b: float,
    expected_origin_A: float,
    expected_origin_B: float,
    hair_enabled: bool,
) -> Dict[str, Any]:
    if not hair_enabled:
        return {
            "name_hair_total": float(p_m_a + p_m_b),
            "name_readout_l1_error": float("nan"),
            "name_readout_ok": False,
        }
    error = abs(float(p_m_a) - float(expected_origin_A)) + abs(float(p_m_b) - float(expected_origin_B))
    return {
        "name_hair_total": float(p_m_a + p_m_b),
        "name_readout_l1_error": error,
        "name_readout_ok": bool(error <= params.name_tol),
    }


def channel_metrics(
    params: Params,
    stage: str,
    model: str,
    delta_f: float,
    n_a: int,
    n_b: int,
    hair_enabled: bool,
    channel: str,
    vector: np.ndarray,
    p_target: float | None,
    expected_origin_A: float,
    expected_origin_B: float,
    T: float,
    R: float,
    readout_enabled: bool = True,
    control_family: str = "scattering_matrix",
) -> Dict[str, Any]:
    h = harmonic_distribution(params, vector)
    n_eff, n_eff_2 = effective_n(h)
    eta = eta_distribution(params, vector, [params.m_A, params.m_B])
    p_m_a = eta.get(params.m_A, 0.0)
    p_m_b = eta.get(params.m_B, 0.0)
    chi_metrics = chi_position_metrics(params, vector, max(n_a, n_b))
    name_metrics = name_readout_metrics(params, p_m_a, p_m_b, expected_origin_A, expected_origin_B, hair_enabled)
    p = pure_expect_p(params, vector)
    p_abs_error = abs(p - p_target) if p_target is not None else float("nan")
    compressed_q_A = params.q_A * (T - R)
    compressed_q_B = params.q_B * (T - R)
    compressed_q_channel = compressed_q_A if channel == "minus_out" else compressed_q_B
    d_q = abs(p - compressed_q_channel) / (abs(p) + 1.0e-12)
    return {
        "stage": stage,
        "model": model,
        "delta_f": delta_f,
        "T": T,
        "R": R,
        "channel": channel,
        "N_A": n_a,
        "N_B": n_b,
        "hair_enabled": hair_enabled,
        "norm": norm2(vector),
        "p_chi": p,
        "p_target": p_target if p_target is not None else float("nan"),
        "p_abs_error": p_abs_error,
        "d_q": d_q,
        "L": localization(vector),
        **chi_metrics,
        "N_eff": n_eff,
        "N_eff_2": n_eff_2,
        "P_m_A": p_m_a,
        "P_m_B": p_m_b,
        **name_metrics,
        "expected_origin_A": expected_origin_A,
        "expected_origin_B": expected_origin_B,
        "readout_enabled": readout_enabled,
        "control_family": control_family,
        "compressed_q_A": compressed_q_A,
        "compressed_q_B": compressed_q_B,
        "compressed_q_channel": compressed_q_channel,
    }


def projection_weight(vector: np.ndarray, basis: np.ndarray) -> float:
    denom = norm2(vector)
    base_denom = norm2(basis)
    if denom <= 0.0 or base_denom <= 0.0:
        return float("nan")
    v = vector / math.sqrt(denom)
    b = basis / math.sqrt(base_denom)
    return float(abs(inner(b, v)) ** 2)


def recursive_state_metrics(
    params: Params,
    stage: str,
    model: str,
    delta_f: float,
    T: float,
    R: float,
    n_a_initial: int,
    n_b_initial: int,
    hair_enabled: bool,
    collision_index: int,
    channel: str,
    vector: np.ndarray,
    initial_a: np.ndarray,
    initial_b: np.ndarray,
) -> Dict[str, Any]:
    h = harmonic_distribution(params, vector)
    n_eff, n_eff_2 = effective_n(h)
    eta = eta_distribution(params, vector, [params.m_A, params.m_B])
    p_m_a = eta.get(params.m_A, 0.0)
    p_m_b = eta.get(params.m_B, 0.0)
    chi_metrics = chi_position_metrics(params, vector, max(n_a_initial, n_b_initial))
    return {
        "stage": stage,
        "model": model,
        "delta_f": delta_f,
        "T": T,
        "R": R,
        "channel": channel,
        "N_A": n_a_initial,
        "N_B": n_b_initial,
        "hair_enabled": hair_enabled,
        "collision_index": collision_index,
        "norm": norm2(vector),
        "p_chi": pure_expect_p(params, vector),
        "p_target": float("nan"),
        "p_abs_error": float("nan"),
        "d_q": float("nan"),
        "L": localization(vector),
        **chi_metrics,
        "N_eff": n_eff,
        "N_eff_2": n_eff_2,
        "P_m_A": p_m_a,
        "P_m_B": p_m_b,
        "name_hair_total": float(p_m_a + p_m_b),
        "name_readout_l1_error": float("nan"),
        "name_readout_ok": bool(hair_enabled and (p_m_a + p_m_b) > 0.0),
        "expected_origin_A": projection_weight(vector, initial_a),
        "expected_origin_B": projection_weight(vector, initial_b),
        "readout_enabled": True,
        "control_family": "recursive_scattering_matrix",
        "compressed_q_A": float("nan"),
        "compressed_q_B": float("nan"),
        "compressed_q_channel": float("nan"),
    }


def model_delta(name: str, delta_f: float) -> Tuple[str, float]:
    return name, delta_f


def rows_for_case(
    params: Params,
    stage: str,
    model: str,
    delta_f: float,
    n_a: int,
    n_b: int,
    hair_enabled: bool,
    compressed_targets: bool,
    readout_enabled: bool = True,
) -> List[Dict[str, Any]]:
    out = scattering_outputs(params, n_a, n_b, hair_enabled, delta_f)
    T = out["T"]
    R = out["R"]
    rows: List[Dict[str, Any]] = []
    rows.append(
        channel_metrics(
            params,
            stage,
            model,
            delta_f,
            n_a,
            n_b,
            hair_enabled,
            "minus_out",
            out["out_minus"],
            -params.q_A if compressed_targets else None,
            expected_origin_A=R,
            expected_origin_B=T,
            T=T,
            R=R,
            readout_enabled=readout_enabled,
            control_family="scattering_matrix",
        )
    )
    rows.append(
        channel_metrics(
            params,
            stage,
            model,
            delta_f,
            n_a,
            n_b,
            hair_enabled,
            "plus_out",
            out["out_plus"],
            -params.q_B if compressed_targets else None,
            expected_origin_A=T,
            expected_origin_B=R,
            T=T,
            R=R,
            readout_enabled=readout_enabled,
            control_family="scattering_matrix",
        )
    )
    return rows


def rows_for_compressed_control(
    params: Params,
    stage: str,
    model: str,
    n_a: int,
    n_b: int,
    hair_enabled: bool,
    mode: str,
) -> List[Dict[str, Any]]:
    a_trans = make_state(params, n_a, params.q_A, params.m_A, hair_enabled, params.A_A)
    a_ref = make_state(params, n_a, -params.q_A, params.m_A, hair_enabled, params.A_A)
    b_trans = make_state(params, n_b, params.q_B, params.m_B, hair_enabled, params.A_B)
    b_ref = make_state(params, n_b, -params.q_B, params.m_B, hair_enabled, params.A_B)

    if mode in {"compressed_reflection", "simple_reflection"}:
        minus = a_ref
        plus = b_ref
        expected_minus = (1.0, 0.0)
        expected_plus = (0.0, 1.0)
        T = 0.0
        R = 1.0
    elif mode == "compressed_transmission":
        minus = b_trans
        plus = a_trans
        expected_minus = (0.0, 1.0)
        expected_plus = (1.0, 0.0)
        T = 1.0
        R = 0.0
    else:
        raise ValueError(f"unknown compressed control mode: {mode}")

    return [
        channel_metrics(
            params,
            stage,
            model,
            float("nan"),
            n_a,
            n_b,
            hair_enabled,
            "minus_out",
            minus,
            None,
            expected_origin_A=expected_minus[0],
            expected_origin_B=expected_minus[1],
            T=T,
            R=R,
            readout_enabled=True,
            control_family="compressed_control",
        ),
        channel_metrics(
            params,
            stage,
            model,
            float("nan"),
            n_a,
            n_b,
            hair_enabled,
            "plus_out",
            plus,
            None,
            expected_origin_A=expected_plus[0],
            expected_origin_B=expected_plus[1],
            T=T,
            R=R,
            readout_enabled=True,
            control_family="compressed_control",
        ),
    ]


def stage0_rows(params: Params) -> List[Dict[str, Any]]:
    n = params.high_n
    rows: List[Dict[str, Any]] = []
    cases = [
        model_delta("fermionic_scattering_complete_reflection", params.delta_f_fermion),
        model_delta("bosonic_scattering_transmission", params.delta_f_boson),
        model_delta("partial_scattering_half", params.delta_f_half),
        model_delta("partial_scattering_R055", 2.0 * math.asin(math.sqrt(0.55))),
    ]
    for model, delta_f in cases:
        rows.extend(
            rows_for_case(
                params,
                "stage0_full_reflection_base",
                model,
                delta_f,
                n,
                n,
                True,
                compressed_targets=(model == "fermionic_scattering_complete_reflection"),
            )
        )
    return rows


def low_n_rows(params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for hair_enabled in [True, False]:
        for n in [99, 63, 31, 15, 7, 3, 1]:
            rows.extend(
                rows_for_case(
                    params,
                    "stage1_odd_harmonic_bottom",
                    "fermionic_scattering_complete_reflection",
                    params.delta_f_fermion,
                    n,
                    n,
                    hair_enabled,
                    compressed_targets=True,
                )
            )
    return rows


def observation_stop_rows(params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for hair_enabled in [True, False]:
        for readout_enabled in [True, False]:
            for n in [63, 15, 3, 1]:
                rows.extend(
                    rows_for_case(
                        params,
                        "stage2_observation_stop_control",
                        "fermionic_scattering_complete_reflection",
                        params.delta_f_fermion,
                        n,
                        n,
                        hair_enabled,
                        compressed_targets=True,
                        readout_enabled=readout_enabled,
                    )
                )
    return rows


def one_side_high_harmonic_rows(params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    pairs = [(3, 63), (7, 63), (15, 63), (3, 31)]
    cases = [
        model_delta("bosonic_scattering_transmission", params.delta_f_boson),
        model_delta("partial_scattering_half", params.delta_f_half),
        model_delta("partial_scattering_R055", 2.0 * math.asin(math.sqrt(0.55))),
        model_delta("fermionic_scattering_complete_reflection", params.delta_f_fermion),
    ]
    for hair_enabled in [True, False]:
        for n_a, n_b in pairs:
            for model, delta_f in cases:
                rows.extend(
                    rows_for_case(
                        params,
                        "stage4_one_side_high_harmonic",
                        model,
                        delta_f,
                        n_a,
                        n_b,
                        hair_enabled,
                        compressed_targets=False,
                    )
                )
    return rows


def control_comparison_rows(params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    pairs = [(3, 63), (7, 63), (3, 31), (1, 31)]
    scattering_cases = [
        model_delta("fermionic_scattering_complete_reflection", params.delta_f_fermion),
        model_delta("bosonic_scattering_transmission", params.delta_f_boson),
        model_delta("partial_scattering_R055", 2.0 * math.asin(math.sqrt(0.55))),
    ]
    compressed_cases = ["compressed_reflection", "simple_reflection", "compressed_transmission"]
    for hair_enabled in [True, False]:
        for n_a, n_b in pairs:
            for model, delta_f in scattering_cases:
                rows.extend(
                    rows_for_case(
                        params,
                        "stage5_control_group_comparison",
                        model,
                        delta_f,
                        n_a,
                        n_b,
                        hair_enabled,
                        compressed_targets=False,
                    )
                )
            for model in compressed_cases:
                rows.extend(
                    rows_for_compressed_control(
                        params,
                        "stage5_control_group_comparison",
                        model,
                        n_a,
                        n_b,
                        hair_enabled,
                        model,
                    )
                )
    return rows


def recursive_transfer_rows(params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    pairs = [(1, 63), (63, 1), (1, 31), (31, 1)]
    cases = [
        model_delta("fermionic_scattering_complete_reflection", params.delta_f_fermion),
        model_delta("bosonic_scattering_transmission", params.delta_f_boson),
        model_delta("partial_scattering_half", params.delta_f_half),
        model_delta("partial_scattering_R055", 2.0 * math.asin(math.sqrt(0.55))),
    ]
    hair_enabled = True
    for n_a, n_b in pairs:
        for model, delta_f in cases:
            t, r, T, R = scattering_coefficients(delta_f)
            a = make_state(params, n_a, params.q_A, params.m_A, hair_enabled, params.A_A)
            b = make_state(params, n_b, params.q_B, params.m_B, hair_enabled, params.A_B)
            initial_a = a.copy()
            initial_b = b.copy()
            for collision_index in range(params.recursive_collision_count + 1):
                rows.append(
                    recursive_state_metrics(
                        params,
                        "stage6_recursive_one_side_high_harmonic",
                        model,
                        delta_f,
                        T,
                        R,
                        n_a,
                        n_b,
                        hair_enabled,
                        collision_index,
                        "A_channel",
                        a,
                        initial_a,
                        initial_b,
                    )
                )
                rows.append(
                    recursive_state_metrics(
                        params,
                        "stage6_recursive_one_side_high_harmonic",
                        model,
                        delta_f,
                        T,
                        R,
                        n_a,
                        n_b,
                        hair_enabled,
                        collision_index,
                        "B_channel",
                        b,
                        initial_a,
                        initial_b,
                    )
                )
                if collision_index >= params.recursive_collision_count:
                    break
                a_next = normalize(r * a + t * b)
                b_next = normalize(t * a + r * b)
                a = a_next
                b = b_next
    return rows


def recursive_r_sweep_rows(params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    n_a = 1
    n_b = params.high_n
    hair_enabled = True
    for r_value in params.r_sweep_values:
        delta_f = delta_from_reflection_rate(r_value)
        model = f"recursive_scattering_R{int(round(r_value * 100)):03d}"
        t, r, T, R = scattering_coefficients(delta_f)
        a = make_state(params, n_a, params.q_A, params.m_A, hair_enabled, params.A_A)
        b = make_state(params, n_b, params.q_B, params.m_B, hair_enabled, params.A_B)
        initial_a = a.copy()
        initial_b = b.copy()
        for collision_index in range(params.recursive_collision_count + 1):
            rows.append(
                recursive_state_metrics(
                    params,
                    "stage7_recursive_R_sweep",
                    model,
                    delta_f,
                    T,
                    R,
                    n_a,
                    n_b,
                    hair_enabled,
                    collision_index,
                    "A_channel",
                    a,
                    initial_a,
                    initial_b,
                )
            )
            rows.append(
                recursive_state_metrics(
                    params,
                    "stage7_recursive_R_sweep",
                    model,
                    delta_f,
                    T,
                    R,
                    n_a,
                    n_b,
                    hair_enabled,
                    collision_index,
                    "B_channel",
                    b,
                    initial_a,
                    initial_b,
                )
            )
            if collision_index >= params.recursive_collision_count:
                break
            a_next = normalize(r * a + t * b)
            b_next = normalize(t * a + r * b)
            a = a_next
            b = b_next
    return rows


def rows_by_collision(rows: List[Dict[str, Any]]) -> Dict[int, Dict[str, Dict[str, Any]]]:
    grouped: Dict[int, Dict[str, Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(int(row["collision_index"]), {})[str(row["channel"])] = row
    return grouped


def recursive_pair_summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    grouped = rows_by_collision(rows)
    first = grouped[min(grouped)]
    last = grouped[max(grouped)]
    gaps = []
    for collision_index, by_channel in grouped.items():
        if "A_channel" in by_channel and "B_channel" in by_channel:
            l_gap = abs(float(by_channel["A_channel"]["L"]) - float(by_channel["B_channel"]["L"]))
            n_gap = abs(float(by_channel["A_channel"]["N_eff"]) - float(by_channel["B_channel"]["N_eff"]))
            gaps.append((l_gap, n_gap, collision_index))
    min_l_gap, min_n_gap, min_gap_collision = min(gaps, key=lambda item: item[0]) if gaps else (float("nan"), float("nan"), None)
    max_collision = max(grouped)
    tail_from = max(0, max_collision - 24)
    tail_gaps = [(l_gap, n_gap, collision_index) for l_gap, n_gap, collision_index in gaps if collision_index >= tail_from]
    tail_l_values = [item[0] for item in tail_gaps]
    tail_n_values = [item[1] for item in tail_gaps]
    return {
        "L_A_initial": float(first["A_channel"]["L"]),
        "L_B_initial": float(first["B_channel"]["L"]),
        "L_A_final": float(last["A_channel"]["L"]),
        "L_B_final": float(last["B_channel"]["L"]),
        "N_eff_A_initial": float(first["A_channel"]["N_eff"]),
        "N_eff_B_initial": float(first["B_channel"]["N_eff"]),
        "N_eff_A_final": float(last["A_channel"]["N_eff"]),
        "N_eff_B_final": float(last["B_channel"]["N_eff"]),
        "L_gap_initial": abs(float(first["A_channel"]["L"]) - float(first["B_channel"]["L"])),
        "L_gap_final": abs(float(last["A_channel"]["L"]) - float(last["B_channel"]["L"])),
        "L_gap_min": float(min_l_gap),
        "N_eff_gap_initial": abs(float(first["A_channel"]["N_eff"]) - float(first["B_channel"]["N_eff"])),
        "N_eff_gap_final": abs(float(last["A_channel"]["N_eff"]) - float(last["B_channel"]["N_eff"])),
        "N_eff_gap_at_L_gap_min": float(min_n_gap),
        "L_gap_min_collision": min_gap_collision,
        "tail_from_collision": tail_from,
        "tail_L_gap_min": min(tail_l_values) if tail_l_values else float("nan"),
        "tail_L_gap_max": max(tail_l_values) if tail_l_values else float("nan"),
        "tail_N_eff_gap_min": min(tail_n_values) if tail_n_values else float("nan"),
        "tail_N_eff_gap_max": max(tail_n_values) if tail_n_values else float("nan"),
        "A_localization_increased": bool(float(last["A_channel"]["L"]) > float(first["A_channel"]["L"])),
        "B_localization_decreased": bool(float(last["B_channel"]["L"]) < float(first["B_channel"]["L"])),
    }


def recursive_sweep_summary(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    sweep_rows = [row for row in rows if row["stage"] == "stage7_recursive_R_sweep"]
    summaries: List[Dict[str, Any]] = []
    for model in sorted({str(row["model"]) for row in sweep_rows}):
        model_rows = [row for row in sweep_rows if row["model"] == model]
        if not model_rows:
            continue
        summary = recursive_pair_summary(model_rows)
        first = model_rows[0]
        summaries.append(
            {
                "model": model,
                "R": float(first["R"]),
                "T": float(first["T"]),
                "L_gap_initial": summary["L_gap_initial"],
                "L_gap_final": summary["L_gap_final"],
                "L_gap_min": summary["L_gap_min"],
                "L_gap_min_collision": summary["L_gap_min_collision"],
                "N_eff_gap_initial": summary["N_eff_gap_initial"],
                "N_eff_gap_final": summary["N_eff_gap_final"],
                "N_eff_gap_at_L_gap_min": summary["N_eff_gap_at_L_gap_min"],
                "tail_from_collision": summary["tail_from_collision"],
                "tail_L_gap_min": summary["tail_L_gap_min"],
                "tail_L_gap_max": summary["tail_L_gap_max"],
                "tail_N_eff_gap_min": summary["tail_N_eff_gap_min"],
                "tail_N_eff_gap_max": summary["tail_N_eff_gap_max"],
            }
        )
    return sorted(summaries, key=lambda row: float(row["R"]))


def recursive_snapshot_states(
    params: Params,
    r_value: float,
    collisions: Iterable[int],
    n_a: int = 1,
    n_b: int = 63,
) -> Dict[int, Dict[str, Any]]:
    targets = set(int(collision) for collision in collisions)
    delta_f = delta_from_reflection_rate(r_value)
    t, r, T, R = scattering_coefficients(delta_f)
    hair_enabled = True
    a = make_state(params, n_a, params.q_A, params.m_A, hair_enabled, params.A_A)
    b = make_state(params, n_b, params.q_B, params.m_B, hair_enabled, params.A_B)
    initial_a = a.copy()
    initial_b = b.copy()
    out: Dict[int, Dict[str, Any]] = {}
    for collision_index in range(max(targets) + 1):
        if collision_index in targets:
            a_metrics = recursive_state_metrics(
                params,
                "snapshot",
                f"R{int(round(r_value * 100)):03d}",
                delta_f,
                T,
                R,
                n_a,
                n_b,
                hair_enabled,
                collision_index,
                "A_channel",
                a,
                initial_a,
                initial_b,
            )
            b_metrics = recursive_state_metrics(
                params,
                "snapshot",
                f"R{int(round(r_value * 100)):03d}",
                delta_f,
                T,
                R,
                n_a,
                n_b,
                hair_enabled,
                collision_index,
                "B_channel",
                b,
                initial_a,
                initial_b,
            )
            out[collision_index] = {
                "A_channel": a_metrics,
                "B_channel": b_metrics,
                "rho_A": chi_density(params, a),
                "rho_B": chi_density(params, b),
            }
        if collision_index >= max(targets):
            break
        a_next = normalize(r * a + t * b)
        b_next = normalize(t * a + r * b)
        a = a_next
        b = b_next
    return out


def compute_verdict(params: Params, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    stage0 = [
        row
        for row in rows
        if row["stage"] == "stage0_full_reflection_base"
        and row["model"] == "fermionic_scattering_complete_reflection"
    ]
    stage0_p_ok = all(float(row["p_abs_error"]) <= params.p_tol for row in stage0)
    stage0_d_q_ok = all(float(row["d_q"]) <= params.d_q_tol for row in stage0)
    stage0_norm_ok = all(abs(float(row["norm"]) - 1.0) <= params.norm_tol for row in stage0)
    stage0_name_ok = all(bool(row["name_readout_ok"]) for row in stage0 if bool(row["hair_enabled"]))

    low = [row for row in rows if row["stage"] == "stage1_odd_harmonic_bottom"]
    good_ns = sorted(
        {
            int(row["N_A"])
            for row in low
            if bool(row["hair_enabled"])
            and float(row["p_abs_error"]) <= params.p_tol
            and float(row["d_q"]) <= params.d_q_tol
            and bool(row["name_readout_ok"])
        }
    )
    good_ns_no_hair = sorted(
        {
            int(row["N_A"])
            for row in low
            if not bool(row["hair_enabled"])
            and float(row["p_abs_error"]) <= params.p_tol
            and float(row["d_q"]) <= params.d_q_tol
            and float(row["name_hair_total"]) <= params.name_tol
        }
    )
    bottom_with_hair = min(good_ns) if good_ns else None
    bottom_without_hair = min(good_ns_no_hair) if good_ns_no_hair else None
    low_with_hair_bottom_rows = [
        row for row in low if bottom_with_hair is not None and int(row["N_A"]) == bottom_with_hair and bool(row["hair_enabled"])
    ]
    low_without_hair_bottom_rows = [
        row for row in low if bottom_without_hair is not None and int(row["N_A"]) == bottom_without_hair and not bool(row["hair_enabled"])
    ]
    name_hair_removed_total_max = max(
        [float(row["name_hair_total"]) for row in low if not bool(row["hair_enabled"])] or [float("nan")]
    )

    one_side = [row for row in rows if row["stage"] == "stage4_one_side_high_harmonic"]
    partial = [row for row in one_side if row["model"] == "partial_scattering_R055"]
    one_side_recorded = any(
        0.0 < float(row["expected_origin_A"]) < 1.0 and 0.0 < float(row["expected_origin_B"]) < 1.0
        for row in partial
    )

    obs = [row for row in rows if row["stage"] == "stage2_observation_stop_control"]
    obs_executed = bool(obs)
    obs_deltas: List[float] = []
    obs_keys = sorted({(row["N_A"], row["N_B"], row["hair_enabled"], row["channel"]) for row in obs})
    for key in obs_keys:
        on = [
            row
            for row in obs
            if (row["N_A"], row["N_B"], row["hair_enabled"], row["channel"]) == key
            and bool(row["readout_enabled"])
        ]
        off = [
            row
            for row in obs
            if (row["N_A"], row["N_B"], row["hair_enabled"], row["channel"]) == key
            and not bool(row["readout_enabled"])
        ]
        if on and off:
            obs_deltas.append(abs(float(on[0]["L"]) - float(off[0]["L"])))

    controls = [row for row in rows if row["stage"] == "stage5_control_group_comparison"]
    control_models = sorted({row["model"] for row in controls})

    recursive = [row for row in rows if row["stage"] == "stage6_recursive_one_side_high_harmonic"]
    recursive_main = [
        row
        for row in recursive
        if row["model"] == "partial_scattering_R055"
        and int(row["N_A"]) == 1
        and int(row["N_B"]) == 63
    ]
    recursive_summary = recursive_pair_summary(recursive_main) if recursive_main else {}
    r_sweep_summary = recursive_sweep_summary(rows)

    return {
        "stage0_reproduced": bool(stage0_p_ok and stage0_d_q_ok and stage0_norm_ok),
        "stage0_p_reflection_ok": bool(stage0_p_ok),
        "stage0_d_q_small": bool(stage0_d_q_ok),
        "stage0_norm_ok": bool(stage0_norm_ok),
        "stage0_name_readout_ok": bool(stage0_name_ok),
        "odd_harmonic_bottom_with_hair": bottom_with_hair,
        "odd_harmonic_bottom_without_hair": bottom_without_hair,
        "bottom_with_hair_name_readout_ok": bool(
            low_with_hair_bottom_rows and all(bool(row["name_readout_ok"]) for row in low_with_hair_bottom_rows)
        ),
        "bottom_without_hair_name_removed": bool(
            low_without_hair_bottom_rows and all(float(row["name_hair_total"]) <= params.name_tol for row in low_without_hair_bottom_rows)
        ),
        "name_hair_removed_total_max": name_hair_removed_total_max,
        "bottom_with_hair_chi_center_cell_mass_min": min(
            [float(row["chi_center_cell_mass"]) for row in low_with_hair_bottom_rows] or [float("nan")]
        ),
        "bottom_without_hair_chi_center_cell_mass_min": min(
            [float(row["chi_center_cell_mass"]) for row in low_without_hair_bottom_rows] or [float("nan")]
        ),
        "one_side_high_harmonic_recorded": bool(one_side_recorded),
        "observation_stop_executed": obs_executed,
        "observation_stop_max_L_delta": max(obs_deltas) if obs_deltas else None,
        "stage5_control_comparison_executed": bool(controls),
        "stage5_control_models": control_models,
        "recursive_one_side_high_harmonic_executed": bool(recursive),
        "recursive_partial_R055_L_A_initial": recursive_summary.get("L_A_initial"),
        "recursive_partial_R055_L_A_final": recursive_summary.get("L_A_final"),
        "recursive_partial_R055_L_B_initial": recursive_summary.get("L_B_initial"),
        "recursive_partial_R055_L_B_final": recursive_summary.get("L_B_final"),
        "recursive_partial_R055_L_gap_initial": recursive_summary.get("L_gap_initial"),
        "recursive_partial_R055_L_gap_final": recursive_summary.get("L_gap_final"),
        "recursive_partial_R055_L_gap_min": recursive_summary.get("L_gap_min"),
        "recursive_partial_R055_L_gap_min_collision": recursive_summary.get("L_gap_min_collision"),
        "recursive_partial_R055_N_eff_A_initial": recursive_summary.get("N_eff_A_initial"),
        "recursive_partial_R055_N_eff_A_final": recursive_summary.get("N_eff_A_final"),
        "recursive_partial_R055_N_eff_B_initial": recursive_summary.get("N_eff_B_initial"),
        "recursive_partial_R055_N_eff_B_final": recursive_summary.get("N_eff_B_final"),
        "recursive_partial_R055_N_eff_gap_initial": recursive_summary.get("N_eff_gap_initial"),
        "recursive_partial_R055_N_eff_gap_final": recursive_summary.get("N_eff_gap_final"),
        "recursive_partial_R055_tail_from_collision": recursive_summary.get("tail_from_collision"),
        "recursive_partial_R055_tail_L_gap_min": recursive_summary.get("tail_L_gap_min"),
        "recursive_partial_R055_tail_L_gap_max": recursive_summary.get("tail_L_gap_max"),
        "recursive_partial_R055_tail_N_eff_gap_min": recursive_summary.get("tail_N_eff_gap_min"),
        "recursive_partial_R055_tail_N_eff_gap_max": recursive_summary.get("tail_N_eff_gap_max"),
        "recursive_partial_R055_A_localization_increased": recursive_summary.get("A_localization_increased"),
        "recursive_partial_R055_B_localization_decreased": recursive_summary.get("B_localization_decreased"),
        "recursive_R_sweep_executed": bool(r_sweep_summary),
        "recursive_R_sweep": r_sweep_summary,
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def clean_value(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return float(value)
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def serialise_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [{key: clean_value(value) for key, value in row.items()} for row in rows]


def load_acceleration_base() -> Dict[str, Any]:
    if not ACCELERATION_BASE_PATH.exists():
        return {
            "path": str(ACCELERATION_BASE_PATH),
            "loaded": False,
            "acceleration_base_ok": False,
            "verdict": {},
        }
    data = json.loads(ACCELERATION_BASE_PATH.read_text(encoding="utf-8"))
    verdict = data.get("aggregate_verdict", {})
    ok = bool(
        abs(float(verdict.get("fermionic_reflection_rate", float("nan"))) - 1.0) <= 1.0e-12
        and abs(float(verdict.get("fermionic_transmission_rate", float("nan")))) <= 1.0e-20
        and abs(float(verdict.get("fermionic_q_out_factor", float("nan"))) + 1.0) <= 1.0e-12
        and abs(float(verdict.get("max_Q_closed_abs", float("nan")))) <= 1.0e-12
        and bool(verdict.get("label_free_pass_vs_fermionic_match_all_cases", False))
        and bool(verdict.get("fermionic_regular_cell_harmonic_consistent_nonstrong_modes", False))
        and bool(verdict.get("fermionic_c1_area_sweep_detected_all_cases", False))
    )
    keys = [
        "fermionic_reflection_rate",
        "fermionic_transmission_rate",
        "fermionic_q_out_factor",
        "max_Q_closed_abs",
        "label_free_pass_vs_fermionic_match_all_cases",
        "fermionic_regular_cell_harmonic_consistent_nonstrong_modes",
        "fermionic_c1_area_sweep_detected_all_cases",
        "fermionic_c1_readout_off_max_epsilon_c_abs",
    ]
    return {
        "path": str(ACCELERATION_BASE_PATH),
        "loaded": True,
        "acceleration_base_ok": ok,
        "verdict": {key: verdict.get(key) for key in keys},
    }


def attach_acceleration_base(rows: List[Dict[str, Any]], acceleration_base_ok: bool) -> None:
    for row in rows:
        row["acceleration_base_ok"] = bool(acceleration_base_ok)


def make_plots(rows: List[Dict[str, Any]]) -> Dict[str, str]:
    outputs: Dict[str, str] = {}

    stage0 = [row for row in rows if row["stage"] == "stage0_full_reflection_base"]
    labels = [f"{row['model']}\n{row['channel']}" for row in stage0]
    p_values = [float(row["p_chi"]) for row in stage0]
    d_q_values = [0.0 if math.isnan(float(row["d_q"])) else float(row["d_q"]) for row in stage0]
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    axes[0].bar(labels, p_values)
    axes[0].axhline(-1.0, color="black", linestyle="--", linewidth=1)
    axes[0].axhline(1.0, color="black", linestyle=":", linewidth=1)
    axes[0].set_ylabel("p_chi")
    axes[0].set_title("Stage 0 scattering-matrix direction readout")
    axes[1].bar(labels, d_q_values)
    axes[1].set_ylabel("d_q")
    axes[1].tick_params(axis="x", rotation=80)
    path = OUT_DIR / "exchange_scattering_matrix_stage0_diagnostics_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["stage0_plot"] = path.name

    low = [row for row in rows if row["stage"] == "stage1_odd_harmonic_bottom"]
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for hair_enabled, marker in [(True, "o"), (False, "s")]:
        subset = [row for row in low if bool(row["hair_enabled"]) == hair_enabled and row["channel"] == "minus_out"]
        subset = sorted(subset, key=lambda row: int(row["N_A"]))
        ax.plot([int(row["N_A"]) for row in subset], [float(row["d_q"]) for row in subset], marker=marker, label=f"hair={hair_enabled}")
    ax.set_xscale("log", base=2)
    ax.invert_xaxis()
    ax.set_xlabel("N")
    ax.set_ylabel("d_q")
    ax.set_title("Odd-harmonic bottom, minus channel")
    ax.legend()
    path = OUT_DIR / "exchange_scattering_matrix_low_n_bottom_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["low_n_plot"] = path.name

    asym = [
        row
        for row in rows
        if row["stage"] == "stage4_one_side_high_harmonic"
        and row["hair_enabled"] is True
        and row["channel"] == "minus_out"
    ]
    labels = [f"{row['model']}\nN{row['N_A']}-{row['N_B']}" for row in asym]
    n_values = [float(row["N_eff"]) for row in asym]
    fig, ax = plt.subplots(figsize=(13, 5), constrained_layout=True)
    ax.bar(labels, n_values)
    ax.set_ylabel("N_eff")
    ax.set_title("Asymmetric harmonic transfer, hair enabled, minus channel")
    ax.tick_params(axis="x", rotation=80)
    path = OUT_DIR / "exchange_scattering_matrix_one_side_high_harmonic_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["one_side_high_harmonic_plot"] = path.name

    recursive = [
        row
        for row in rows
        if row["stage"] == "stage6_recursive_one_side_high_harmonic"
        and row["model"] == "partial_scattering_R055"
        and int(row["N_A"]) == 1
        and int(row["N_B"]) == 63
    ]
    if recursive:
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
        for channel, marker in [("A_channel", "o"), ("B_channel", "s")]:
            subset = [row for row in recursive if row["channel"] == channel]
            subset = sorted(subset, key=lambda row: int(row["collision_index"]))
            x = [int(row["collision_index"]) for row in subset]
            axes[0].plot(x, [float(row["L"]) for row in subset], marker=marker, label=channel)
            axes[1].plot(x, [float(row["N_eff"]) for row in subset], marker=marker, label=channel)
        axes[0].set_ylabel("L")
        axes[0].set_title("Recursive one-side high harmonic: localization")
        axes[0].legend()
        axes[1].set_xlabel("collision index")
        axes[1].set_ylabel("N_eff")
        axes[1].set_title("Recursive one-side high harmonic: N_eff")
        axes[1].legend()
        path = OUT_DIR / "exchange_scattering_matrix_recursive_localization_transfer_v1.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        outputs["recursive_localization_transfer_plot"] = path.name

    r_sweep_summary = recursive_sweep_summary(rows)
    if r_sweep_summary:
        fig, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
        x = [float(row["R"]) for row in r_sweep_summary]
        axes[0].plot(x, [float(row["L_gap_min"]) for row in r_sweep_summary], marker="o", label="min L gap")
        axes[0].plot(x, [float(row["tail_L_gap_max"]) for row in r_sweep_summary], marker="s", label="tail max L gap")
        axes[0].set_ylabel("L gap")
        axes[0].set_title("Recursive R sweep: localization gap")
        axes[0].legend()
        axes[1].plot(x, [float(row["N_eff_gap_at_L_gap_min"]) for row in r_sweep_summary], marker="o", label="N_eff gap at min L gap")
        axes[1].plot(x, [float(row["tail_N_eff_gap_max"]) for row in r_sweep_summary], marker="s", label="tail max N_eff gap")
        axes[1].set_xlabel("R")
        axes[1].set_ylabel("N_eff gap")
        axes[1].set_title("Recursive R sweep: effective harmonic gap")
        axes[1].legend()
        path = OUT_DIR / "exchange_scattering_matrix_recursive_R_sweep_v1.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        outputs["recursive_R_sweep_plot"] = path.name

    chi, _ = make_grids(Params())
    x = chi / math.pi
    snapshot_specs = [
        ("R=0 transmission endpoint", 0.0, 1),
        ("R=0.70 intermediate scattering", 0.70, 42),
        ("R=1 reflection endpoint", 1.0, 1),
    ]
    fig, axes = plt.subplots(len(snapshot_specs), 2, figsize=(12, 9), sharex=True, constrained_layout=True)
    for row_index, (label, r_value, event_collision) in enumerate(snapshot_specs):
        snapshots = recursive_snapshot_states(Params(), r_value, [0, event_collision])
        for col_index, collision_index in enumerate([0, event_collision]):
            ax = axes[row_index][col_index]
            snap = snapshots[collision_index]
            a_metrics = snap["A_channel"]
            b_metrics = snap["B_channel"]
            rho_a = snap["rho_A"] / np.max(snap["rho_A"])
            rho_b = snap["rho_B"] / np.max(snap["rho_B"])
            ax.plot(x, rho_a, label=f"A L={a_metrics['L']:.3g}, N={a_metrics['N_eff']:.3g}")
            ax.plot(x, rho_b, label=f"B L={b_metrics['L']:.3g}, N={b_metrics['N_eff']:.3g}")
            ax.set_title(f"{label}, collision={collision_index}")
            ax.set_ylabel("rho_chi / max")
            ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")
    path = OUT_DIR / "exchange_scattering_matrix_waveform_localization_snapshots_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["waveform_localization_snapshots_plot"] = path.name

    evolution_collisions = [0, 1, 2, 3, 5, 10, 20, 42]
    snapshots = recursive_snapshot_states(Params(), 0.70, evolution_collisions)
    fig, axes = plt.subplots(4, 2, figsize=(12, 12), sharex=True, sharey=True, constrained_layout=True)
    for ax, collision_index in zip(axes.flatten(), evolution_collisions):
        snap = snapshots[collision_index]
        a_metrics = snap["A_channel"]
        b_metrics = snap["B_channel"]
        rho_a = snap["rho_A"] / np.max(snap["rho_A"])
        rho_b = snap["rho_B"] / np.max(snap["rho_B"])
        ax.plot(x, rho_a, label=f"A L={a_metrics['L']:.3g}, N={a_metrics['N_eff']:.3g}")
        ax.plot(x, rho_b, label=f"B L={b_metrics['L']:.3g}, N={b_metrics['N_eff']:.3g}")
        ax.set_title(f"R=0.70, collision={collision_index}")
        ax.set_ylabel("rho_chi / max")
        ax.legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")
    path = OUT_DIR / "exchange_scattering_matrix_R070_waveform_evolution_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["R070_waveform_evolution_plot"] = path.name

    return outputs


def build_report(result: Dict[str, Any]) -> str:
    verdict = result["verdict"]
    outputs = result["outputs"]
    rows = result["rows"]
    acceleration_base = result["acceleration_base"]

    stage0_rows = [row for row in rows if row["stage"] == "stage0_full_reflection_base"]
    stage0_table = "\n".join(
        "| {model} | {channel} | {delta_f:.6g} | {R:.6g} | {T:.6g} | {p_chi:.6g} | {p_target} | {d_q:.6g} | {chi_center_cell_mass:.6g} | {P_m_A:.6g} | {P_m_B:.6g} | {name_readout_l1_error:.6g} | {name_readout_ok} |".format(
            **row
        )
        for row in stage0_rows
    )

    low_rows = [
        row
        for row in rows
        if row["stage"] == "stage1_odd_harmonic_bottom"
        and row["channel"] == "minus_out"
        and int(row["N_A"]) in {63, 15, 7, 3, 1}
    ]
    low_rows = sorted(low_rows, key=lambda row: (not bool(row["hair_enabled"]), -int(row["N_A"])))
    low_table = "\n".join(
        "| {N_A} | {hair_enabled} | {p_chi:.6g} | {d_q:.6g} | {chi_center_cell_mass:.6g} | {chi_peak_contrast:.6g} | {L:.6g} | {P_m_A:.6g} | {P_m_B:.6g} | {name_hair_total:.6g} | {name_readout_l1_error:.6g} | {name_readout_ok} |".format(
            **row
        )
        for row in low_rows
    )

    partial_rows = [
        row
        for row in rows
        if row["stage"] == "stage4_one_side_high_harmonic"
        and row["model"] == "partial_scattering_R055"
        and row["hair_enabled"] is True
    ][:8]
    partial_table = "\n".join(
        "| {channel} | {N_A} | {N_B} | {R:.6g} | {T:.6g} | {N_eff:.6g} | {L:.6g} | {expected_origin_A:.6g} | {expected_origin_B:.6g} |".format(
            **row
        )
        for row in partial_rows
    )

    obs_rows = [
        row
        for row in rows
        if row["stage"] == "stage2_observation_stop_control"
        and row["channel"] == "minus_out"
        and int(row["N_A"]) in {63, 3, 1}
    ][:12]
    obs_table = "\n".join(
        "| {N_A} | {hair_enabled} | {readout_enabled} | {p_chi:.6g} | {L:.6g} | {d_q:.6g} | {name_hair_total:.6g} | {name_readout_ok} |".format(
            **row
        )
        for row in obs_rows
    )

    control_rows = [
        row
        for row in rows
        if row["stage"] == "stage5_control_group_comparison"
        and row["hair_enabled"] is True
        and row["N_A"] in [3, 7]
        and row["N_B"] == 63
        and row["channel"] == "minus_out"
    ]
    control_table = "\n".join(
        "| {model} | {N_A} | {N_B} | {channel} | {R:.6g} | {T:.6g} | {N_eff:.6g} | {L:.6g} | {expected_origin_A:.6g} | {expected_origin_B:.6g} |".format(
            **row
        )
        for row in control_rows
    )

    recursive_rows = [
        row
        for row in rows
        if row["stage"] == "stage6_recursive_one_side_high_harmonic"
        and row["model"] == "partial_scattering_R055"
        and int(row["N_A"]) == 1
        and int(row["N_B"]) == 63
        and int(row["collision_index"]) in {0, 1, 2, 4, 8, 12, 16, 20, 24, 32, 48, 64, 80, 96, 112, 128}
    ]
    recursive_rows = sorted(recursive_rows, key=lambda row: (int(row["collision_index"]), str(row["channel"])))
    recursive_table = "\n".join(
        "| {collision_index} | {channel} | {L:.6g} | {N_eff:.6g} | {chi_center_cell_mass:.6g} | {expected_origin_A:.6g} | {expected_origin_B:.6g} | {P_m_A:.6g} | {P_m_B:.6g} |".format(
            **row
        )
        for row in recursive_rows
    )
    r_sweep_table = "\n".join(
        "| {R:.6g} | {T:.6g} | {L_gap_min:.6g} | {L_gap_min_collision} | {N_eff_gap_at_L_gap_min:.6g} | {tail_L_gap_min:.6g} | {tail_L_gap_max:.6g} | {tail_N_eff_gap_min:.6g} | {tail_N_eff_gap_max:.6g} |".format(
            **row
        )
        for row in verdict["recursive_R_sweep"]
    )

    return f"""# 交換干渉散乱行列フェルミオン的衝突 予備実験検証メモ v1

## 実行条件

V2 散乱行列基準で、加速度基底、低奇数倍音底、片側高次倍音条件を実行した。

## 判定

| 項目 | 結果 |
|---|---:|
| Stage 0 完全反射基底確認 | `{str(verdict['stage0_reproduced']).lower()}` |
| Stage 0 p 反転 | `{str(verdict['stage0_p_reflection_ok']).lower()}` |
| Stage 0 d_q | `{str(verdict['stage0_d_q_small']).lower()}` |
| Stage 0 ノルム | `{str(verdict['stage0_norm_ok']).lower()}` |
| Stage 0 名前毛読出し | `{str(verdict['stage0_name_readout_ok']).lower()}` |
| 加速度V2基底読込 | `{str(verdict['acceleration_base_loaded']).lower()}` |
| 加速度V2基底 | `{str(verdict['acceleration_base_ok']).lower()}` |
| 奇数倍音底 hairあり | `{verdict['odd_harmonic_bottom_with_hair']}` |
| 奇数倍音底 hairなし | `{verdict['odd_harmonic_bottom_without_hair']}` |
| 底 hairあり 名前毛 | `{str(verdict['bottom_with_hair_name_readout_ok']).lower()}` |
| 底 hairなし 名前毛除去 | `{str(verdict['bottom_without_hair_name_removed']).lower()}` |
| hairなし 名前毛総量最大 | `{verdict['name_hair_removed_total_max']}` |
| 底 hairあり chi中心セル質量最小 | `{verdict['bottom_with_hair_chi_center_cell_mass_min']}` |
| 底 hairなし chi中心セル質量最小 | `{verdict['bottom_without_hair_chi_center_cell_mass_min']}` |
| 片側高次倍音条件記録 | `{str(verdict['one_side_high_harmonic_recorded']).lower()}` |
| 観測停止対照 | `{str(verdict['observation_stop_executed']).lower()}` |
| 観測停止 L 最大差分 | `{verdict['observation_stop_max_L_delta']}` |
| Stage5 対照群比較 | `{str(verdict['stage5_control_comparison_executed']).lower()}` |
| 再帰片側高次倍音 | `{str(verdict['recursive_one_side_high_harmonic_executed']).lower()}` |
| 再帰 R055 L_A 初期 | `{verdict['recursive_partial_R055_L_A_initial']}` |
| 再帰 R055 L_A 最終 | `{verdict['recursive_partial_R055_L_A_final']}` |
| 再帰 R055 L_B 初期 | `{verdict['recursive_partial_R055_L_B_initial']}` |
| 再帰 R055 L_B 最終 | `{verdict['recursive_partial_R055_L_B_final']}` |
| 再帰 R055 L差 初期 | `{verdict['recursive_partial_R055_L_gap_initial']}` |
| 再帰 R055 L差 最終 | `{verdict['recursive_partial_R055_L_gap_final']}` |
| 再帰 R055 L差 最小 | `{verdict['recursive_partial_R055_L_gap_min']}` |
| 再帰 R055 L差 最小衝突回 | `{verdict['recursive_partial_R055_L_gap_min_collision']}` |
| 再帰 R055 N_eff差 初期 | `{verdict['recursive_partial_R055_N_eff_gap_initial']}` |
| 再帰 R055 N_eff差 最終 | `{verdict['recursive_partial_R055_N_eff_gap_final']}` |
| 再帰 R055 末尾区間開始 | `{verdict['recursive_partial_R055_tail_from_collision']}` |
| 再帰 R055 末尾 L差 最小 | `{verdict['recursive_partial_R055_tail_L_gap_min']}` |
| 再帰 R055 末尾 L差 最大 | `{verdict['recursive_partial_R055_tail_L_gap_max']}` |
| 再帰 R055 末尾 N_eff差 最小 | `{verdict['recursive_partial_R055_tail_N_eff_gap_min']}` |
| 再帰 R055 末尾 N_eff差 最大 | `{verdict['recursive_partial_R055_tail_N_eff_gap_max']}` |
| 再帰 R/T スイープ | `{str(verdict['recursive_R_sweep_executed']).lower()}` |

## 加速度V2基底

| 量 | 値 |
|---|---:|
| loaded | `{str(acceleration_base['loaded']).lower()}` |
| acceleration_base_ok | `{str(acceleration_base['acceleration_base_ok']).lower()}` |
| fermionic_reflection_rate | `{acceleration_base['verdict'].get('fermionic_reflection_rate')}` |
| fermionic_transmission_rate | `{acceleration_base['verdict'].get('fermionic_transmission_rate')}` |
| fermionic_q_out_factor | `{acceleration_base['verdict'].get('fermionic_q_out_factor')}` |
| max_Q_closed_abs | `{acceleration_base['verdict'].get('max_Q_closed_abs')}` |
| label_free_pass_vs_fermionic_match_all_cases | `{acceleration_base['verdict'].get('label_free_pass_vs_fermionic_match_all_cases')}` |
| fermionic_regular_cell_harmonic_consistent_nonstrong_modes | `{acceleration_base['verdict'].get('fermionic_regular_cell_harmonic_consistent_nonstrong_modes')}` |
| fermionic_c1_area_sweep_detected_all_cases | `{acceleration_base['verdict'].get('fermionic_c1_area_sweep_detected_all_cases')}` |

## Stage 0

| model | channel | delta_f | R | T | p_chi | p_target | d_q | chi_center_cell_mass | P_m_A | P_m_B | name_l1_error | name_ok |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
{stage0_table}

![stage0]({outputs['stage0_plot']})

## 低奇数倍音底と内在読出し

| N | hair_enabled | p_chi | d_q | chi_center_cell_mass | chi_peak_contrast | L | P_m_A | P_m_B | name_hair_total | name_l1_error | name_ok |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
{low_table}

## 片側高次倍音条件

| channel | N_A | N_B | R | T | N_eff | L | expected_origin_A | expected_origin_B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{partial_table}

![one side high harmonic]({outputs['one_side_high_harmonic_plot']})

## 観測停止対照

観測あり/なしを診断対照として記録した。

| N | hair_enabled | readout_enabled | p_chi | L | d_q | name_hair_total | name_ok |
|---:|---|---|---:|---:|---:|---:|---|
{obs_table}

## Stage5 対照群比較

片側高次倍音条件で、散乱行列版と圧縮表示対照を比較した。

| model | N_A | N_B | channel | R | T | N_eff | L | expected_origin_A | expected_origin_B |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
{control_table}

## 再帰片側高次倍音

出射チャネルを次回入力へ渡し、複数回の散乱で `L` と `N_eff` を読む。

| collision | channel | L | N_eff | chi_center_cell_mass | origin_A | origin_B | P_m_A | P_m_B |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
{recursive_table}

![recursive localization transfer]({outputs['recursive_localization_transfer_plot']})

## 再帰 R/T スイープ

`N_A=1`, `N_B=63` の片側高次倍音条件で、反射率 `R` を変えて再帰散乱を実行した。

| R | T | L差最小 | L差最小衝突回 | N_eff差 at L差最小 | 末尾L差最小 | 末尾L差最大 | 末尾N_eff差最小 | 末尾N_eff差最大 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{r_sweep_table}

![recursive R sweep]({outputs['recursive_R_sweep_plot']})

## 波形局在スナップショット

`chi` 方向の縮約密度 `rho_chi` を、各線の最大値で振幅正規化し、完全透過端点、中間散乱、完全反射端点で比較した。

完全透過端点 `R=0` と完全反射端点 `R=1` では、低次数波形と高次数波形の差が保存される。

中間散乱 `R=0.70` では、衝突回 `42` で両チャネルの `L` と `N_eff` が近接する。

![waveform localization snapshots]({outputs['waveform_localization_snapshots_plot']})

## R=0.70 波形局在化の再帰発展

中間散乱 `R=0.70` について、衝突回 `0,1,2,3,5,10,20,42` の `rho_chi / max` を図化した。

![R070 waveform evolution]({outputs['R070_waveform_evolution_plot']})

## 奇数倍音底

完全反射 `Delta_F=pi` で `N` を下げ、`N_min` を読む。

![low N]({outputs['low_n_plot']})

## 二チャネル出力

```text
minus_out = r A_ref + t B_trans
plus_out  = t A_trans + r B_ref
```

上記の二チャネル出力から `L` と `H(n)` を読む。

## 出力

| 種別 | ファイル |
|---|---|
| JSON | `{outputs['json']}` |
| CSV | `{outputs['csv']}` |
| Stage0 図 | `{outputs['stage0_plot']}` |
| 低N図 | `{outputs['low_n_plot']}` |
| 片側高次倍音図 | `{outputs['one_side_high_harmonic_plot']}` |
| 再帰局在性図 | `{outputs['recursive_localization_transfer_plot']}` |
| 再帰Rスイープ図 | `{outputs['recursive_R_sweep_plot']}` |
| 波形局在スナップショット図 | `{outputs['waveform_localization_snapshots_plot']}` |
| R=0.70 波形発展図 | `{outputs['R070_waveform_evolution_plot']}` |
| report | `{outputs['report']}` |
"""


def run() -> Dict[str, Any]:
    params = Params()
    acceleration_base = load_acceleration_base()
    rows: List[Dict[str, Any]] = []
    rows.extend(stage0_rows(params))
    rows.extend(low_n_rows(params))
    rows.extend(observation_stop_rows(params))
    rows.extend(one_side_high_harmonic_rows(params))
    rows.extend(control_comparison_rows(params))
    rows.extend(recursive_transfer_rows(params))
    rows.extend(recursive_r_sweep_rows(params))
    attach_acceleration_base(rows, bool(acceleration_base["acceleration_base_ok"]))
    verdict = compute_verdict(params, rows)
    verdict["acceleration_base_loaded"] = bool(acceleration_base["loaded"])
    verdict["acceleration_base_ok"] = bool(acceleration_base["acceleration_base_ok"])

    outputs = {
        "json": "exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1.json",
        "csv": "exchange_scattering_matrix_fermionic_localization_transfer_rows_v1.csv",
        "report": "exchange_scattering_matrix_fermionic_localization_transfer_report_v1.md",
    }
    outputs.update(make_plots(rows))

    result: Dict[str, Any] = {
        "experiment": "exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1",
        "params": asdict(params),
        "acceleration_base": acceleration_base,
        "verdict": verdict,
        "rows": serialise_rows(rows),
        "outputs": outputs,
    }

    write_csv(OUT_DIR / outputs["csv"], rows)
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    report = build_report({"verdict": verdict, "outputs": outputs, "rows": rows, "acceleration_base": acceleration_base})
    (OUT_DIR / outputs["report"]).write_text(report, encoding="utf-8")
    (BASE_DIR / "交換干渉散乱行列フェルミオン的衝突における加速度基底・低奇数倍音底・片側高次倍音予備実験検証メモ_v1.md").write_text(
        report,
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))
