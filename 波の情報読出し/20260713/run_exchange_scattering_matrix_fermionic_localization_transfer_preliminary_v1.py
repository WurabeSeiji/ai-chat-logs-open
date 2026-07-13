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
    copy_distance_tol: float = 1.0e-2
    norm_tol: float = 1.0e-10


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
    copy_target: np.ndarray | None,
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
    p = pure_expect_p(params, vector)
    copy_d = normalized_distance(vector, copy_target) if copy_target is not None else float("nan")
    p_abs_error = abs(p - p_target) if p_target is not None else float("nan")
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
        "copy_distance_d": copy_d,
        "L": localization(vector),
        "N_eff": n_eff,
        "N_eff_2": n_eff_2,
        "P_m_A": eta.get(params.m_A, 0.0),
        "P_m_B": eta.get(params.m_B, 0.0),
        "expected_origin_A": expected_origin_A,
        "expected_origin_B": expected_origin_B,
        "readout_enabled": readout_enabled,
        "control_family": control_family,
        "compressed_q_A": params.q_A * (T - R),
        "compressed_q_B": params.q_B * (T - R),
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
    copy_targets: bool,
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
            -params.q_A if copy_targets else None,
            out["a_ref"] if copy_targets else None,
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
            -params.q_B if copy_targets else None,
            out["b_ref"] if copy_targets else None,
            expected_origin_A=T,
            expected_origin_B=R,
            T=T,
            R=R,
            readout_enabled=readout_enabled,
            control_family="scattering_matrix",
        )
    )
    return rows


def rows_for_copy_control(
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

    if mode in {"copy_reflection", "simple_reflection"}:
        minus = a_ref
        plus = b_ref
        expected_minus = (1.0, 0.0)
        expected_plus = (0.0, 1.0)
        T = 0.0
        R = 1.0
    elif mode == "copy_transmission":
        minus = b_trans
        plus = a_trans
        expected_minus = (0.0, 1.0)
        expected_plus = (1.0, 0.0)
        T = 1.0
        R = 0.0
    else:
        raise ValueError(f"unknown copy control mode: {mode}")

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
            None,
            expected_origin_A=expected_minus[0],
            expected_origin_B=expected_minus[1],
            T=T,
            R=R,
            readout_enabled=True,
            control_family="copy_control",
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
            None,
            expected_origin_A=expected_plus[0],
            expected_origin_B=expected_plus[1],
            T=T,
            R=R,
            readout_enabled=True,
            control_family="copy_control",
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
                "stage0_old_condition_reproduction",
                model,
                delta_f,
                n,
                n,
                True,
                copy_targets=(model == "fermionic_scattering_complete_reflection"),
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
                    "stage1_low_localization_bottom",
                    "fermionic_scattering_complete_reflection",
                    params.delta_f_fermion,
                    n,
                    n,
                    hair_enabled,
                    copy_targets=True,
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
                        copy_targets=True,
                        readout_enabled=readout_enabled,
                    )
                )
    return rows


def asymmetric_transfer_rows(params: Params) -> List[Dict[str, Any]]:
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
                        "stage2_asymmetric_harmonic_transfer",
                        model,
                        delta_f,
                        n_a,
                        n_b,
                        hair_enabled,
                        copy_targets=False,
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
    copy_cases = ["copy_reflection", "simple_reflection", "copy_transmission"]
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
                        copy_targets=False,
                    )
                )
            for model in copy_cases:
                rows.extend(
                    rows_for_copy_control(
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


def compute_verdict(params: Params, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    stage0 = [
        row
        for row in rows
        if row["stage"] == "stage0_old_condition_reproduction"
        and row["model"] == "fermionic_scattering_complete_reflection"
    ]
    stage0_p_ok = all(float(row["p_abs_error"]) <= params.p_tol for row in stage0)
    stage0_copy_ok = all(float(row["copy_distance_d"]) <= params.copy_distance_tol for row in stage0)
    stage0_norm_ok = all(abs(float(row["norm"]) - 1.0) <= params.norm_tol for row in stage0)

    low = [row for row in rows if row["stage"] == "stage1_low_localization_bottom"]
    good_ns = sorted(
        {
            int(row["N_A"])
            for row in low
            if bool(row["hair_enabled"])
            and float(row["p_abs_error"]) <= params.p_tol
            and float(row["copy_distance_d"]) <= params.copy_distance_tol
        }
    )
    good_ns_no_hair = sorted(
        {
            int(row["N_A"])
            for row in low
            if not bool(row["hair_enabled"])
            and float(row["p_abs_error"]) <= params.p_tol
            and float(row["copy_distance_d"]) <= params.copy_distance_tol
        }
    )

    asym = [row for row in rows if row["stage"] == "stage2_asymmetric_harmonic_transfer"]
    partial = [row for row in asym if row["model"] == "partial_scattering_R055"]
    transfer_seen = any(
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

    return {
        "stage0_reproduced": bool(stage0_p_ok and stage0_copy_ok and stage0_norm_ok),
        "stage0_p_reflection_ok": bool(stage0_p_ok),
        "stage0_copy_distance_small": bool(stage0_copy_ok),
        "stage0_norm_ok": bool(stage0_norm_ok),
        "low_localization_bottom_with_hair": min(good_ns) if good_ns else None,
        "low_localization_bottom_without_hair": min(good_ns_no_hair) if good_ns_no_hair else None,
        "asymmetric_partial_transfer_recorded": bool(transfer_seen),
        "observation_stop_executed": obs_executed,
        "observation_stop_max_L_delta": max(obs_deltas) if obs_deltas else None,
        "observation_stop_note": "No damping model is included; on/off rows verify that diagnostic readout does not alter the scattering-matrix output.",
        "stage5_control_comparison_executed": bool(controls),
        "stage5_control_models": control_models,
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def clean_value(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return float(value)
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def serialise_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [{key: clean_value(value) for key, value in row.items()} for row in rows]


def make_plots(rows: List[Dict[str, Any]]) -> Dict[str, str]:
    outputs: Dict[str, str] = {}

    stage0 = [row for row in rows if row["stage"] == "stage0_old_condition_reproduction"]
    labels = [f"{row['model']}\n{row['channel']}" for row in stage0]
    p_values = [float(row["p_chi"]) for row in stage0]
    copy_values = [0.0 if math.isnan(float(row["copy_distance_d"])) else float(row["copy_distance_d"]) for row in stage0]
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    axes[0].bar(labels, p_values)
    axes[0].axhline(-1.0, color="black", linestyle="--", linewidth=1)
    axes[0].axhline(1.0, color="black", linestyle=":", linewidth=1)
    axes[0].set_ylabel("p_chi")
    axes[0].set_title("Stage 0 scattering-matrix direction readout")
    axes[1].bar(labels, copy_values)
    axes[1].set_ylabel("copy distance to full reflection")
    axes[1].tick_params(axis="x", rotation=80)
    path = OUT_DIR / "exchange_scattering_matrix_stage0_diagnostics_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["stage0_plot"] = path.name

    low = [row for row in rows if row["stage"] == "stage1_low_localization_bottom"]
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for hair_enabled, marker in [(True, "o"), (False, "s")]:
        subset = [row for row in low if bool(row["hair_enabled"]) == hair_enabled and row["channel"] == "minus_out"]
        subset = sorted(subset, key=lambda row: int(row["N_A"]))
        ax.plot([int(row["N_A"]) for row in subset], [float(row["copy_distance_d"]) for row in subset], marker=marker, label=f"hair={hair_enabled}")
    ax.set_xscale("log", base=2)
    ax.invert_xaxis()
    ax.set_xlabel("N")
    ax.set_ylabel("copy distance d")
    ax.set_title("Low-localization bottom, minus channel")
    ax.legend()
    path = OUT_DIR / "exchange_scattering_matrix_low_n_bottom_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["low_n_plot"] = path.name

    asym = [
        row
        for row in rows
        if row["stage"] == "stage2_asymmetric_harmonic_transfer"
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
    path = OUT_DIR / "exchange_scattering_matrix_asymmetric_transfer_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    outputs["asymmetric_transfer_plot"] = path.name

    return outputs


def build_report(result: Dict[str, Any]) -> str:
    verdict = result["verdict"]
    outputs = result["outputs"]
    rows = result["rows"]

    stage0_rows = [row for row in rows if row["stage"] == "stage0_old_condition_reproduction"]
    stage0_table = "\n".join(
        "| {model} | {channel} | {delta_f:.6g} | {R:.6g} | {T:.6g} | {p_chi:.6g} | {p_target} | {copy_distance_d} | {P_m_A:.6g} | {P_m_B:.6g} |".format(
            **row
        )
        for row in stage0_rows
    )

    partial_rows = [
        row
        for row in rows
        if row["stage"] == "stage2_asymmetric_harmonic_transfer"
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
        "| {N_A} | {hair_enabled} | {readout_enabled} | {p_chi:.6g} | {L:.6g} | {copy_distance_d:.6g} |".format(
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

    return f"""# 交換干渉散乱行列フェルミオン的衝突 予備実験検証メモ v1

## 目的

20260713 の準備論文を散乱行列版へ修正したため、予備実験を同じ方針で実行し直した。

本実験では、線形重ね合わせ `A+B` から直接 `A',B'` を復元するのではなく、交換干渉位相から反射振幅 `r` と透過振幅 `t` を計算し、分離された入射チャネルへ二チャネル散乱行列として作用させた。

## 判定

| 項目 | 結果 |
|---|---:|
| Stage 0 旧完全反射条件再現 | `{str(verdict['stage0_reproduced']).lower()}` |
| Stage 0 p 反転 | `{str(verdict['stage0_p_reflection_ok']).lower()}` |
| Stage 0 保存コピー距離 | `{str(verdict['stage0_copy_distance_small']).lower()}` |
| Stage 0 ノルム | `{str(verdict['stage0_norm_ok']).lower()}` |
| 低局在性底 hairあり | `{verdict['low_localization_bottom_with_hair']}` |
| 低局在性底 hairなし | `{verdict['low_localization_bottom_without_hair']}` |
| 非対称次数の部分移乗記録 | `{str(verdict['asymmetric_partial_transfer_recorded']).lower()}` |
| 観測停止対照 | `{str(verdict['observation_stop_executed']).lower()}` |
| 観測停止 L 最大差分 | `{verdict['observation_stop_max_L_delta']}` |
| Stage5 対照群比較 | `{str(verdict['stage5_control_comparison_executed']).lower()}` |

## Stage 0

| model | channel | delta_f | R | T | p_chi | p_target | copy_distance_d | P_m_A | P_m_B |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
{stage0_table}

![stage0]({outputs['stage0_plot']})

## 非対称次数の部分移乗例

| channel | N_A | N_B | R | T | N_eff | L | expected_origin_A | expected_origin_B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{partial_table}

![asymmetric transfer]({outputs['asymmetric_transfer_plot']})

## 観測停止対照

本スクリプトでは読出し波による減衰モデルを入れていない。

そのため、観測あり/なしは散乱行列出力を変えない診断対照として記録した。

| N | hair_enabled | readout_enabled | p_chi | L | copy_distance_d |
|---:|---|---|---:|---:|---:|
{obs_table}

## Stage5 対照群比較

非対称次数条件で、散乱行列版と保存コピー型対照を比較した。

| model | N_A | N_B | channel | R | T | N_eff | L | expected_origin_A | expected_origin_B |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
{control_table}

## 低局在性底

散乱行列版では、完全反射 `Delta_F=pi` の場合、低い `N` でも反射コピー条件は保たれた。

![low N]({outputs['low_n_plot']})

## 解釈

今回の結果は、前回の静的な縮約密度版とは異なる。

前回は `A+B` 型の交換合成から縮約密度を作ったため、二つの出射チャネルを復元できなかった。

今回の散乱行列版では、反射振幅 `r` と透過振幅 `t` を用いて、

```text
minus_out = r A_ref + t B_trans
plus_out  = t A_trans + r B_ref
```

を明示的に保持した。

このため、完全反射では旧保存コピー反射条件を再現し、部分反射では A 起因成分と B 起因成分の混合を出射チャネル上に残せた。

## 注意

本実験では、観測停止による減衰差はまだ扱っていない。

散乱行列写像は衝突セル内の局所写像であり、読出し波による包絡減衰モデルを別に入れていないためである。

## 出力

| 種別 | ファイル |
|---|---|
| JSON | `{outputs['json']}` |
| CSV | `{outputs['csv']}` |
| Stage0 図 | `{outputs['stage0_plot']}` |
| 低N図 | `{outputs['low_n_plot']}` |
| 非対称次数図 | `{outputs['asymmetric_transfer_plot']}` |
| report | `{outputs['report']}` |
"""


def run() -> Dict[str, Any]:
    params = Params()
    rows: List[Dict[str, Any]] = []
    rows.extend(stage0_rows(params))
    rows.extend(low_n_rows(params))
    rows.extend(observation_stop_rows(params))
    rows.extend(asymmetric_transfer_rows(params))
    rows.extend(control_comparison_rows(params))
    verdict = compute_verdict(params, rows)

    outputs = {
        "json": "exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1.json",
        "csv": "exchange_scattering_matrix_fermionic_localization_transfer_rows_v1.csv",
        "report": "exchange_scattering_matrix_fermionic_localization_transfer_report_v1.md",
    }
    outputs.update(make_plots(rows))

    result: Dict[str, Any] = {
        "experiment": "exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1",
        "params": asdict(params),
        "verdict": verdict,
        "rows": serialise_rows(rows),
        "outputs": outputs,
    }

    write_csv(OUT_DIR / outputs["csv"], rows)
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    report = build_report({"verdict": verdict, "outputs": outputs, "rows": rows})
    (OUT_DIR / outputs["report"]).write_text(report, encoding="utf-8")
    (BASE_DIR / "交換干渉散乱行列フェルミオン的衝突における低局在性・倍音移乗予備実験検証メモ_v1.md").write_text(
        report,
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))
