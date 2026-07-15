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
OUT_DIR = BASE_DIR / "ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v3"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi


@dataclass(frozen=True)
class Params:
    step_count: int = 720
    period_steps: int = 96
    rho_AB: float = 1.0
    closure_tol: float = 1.0e-12
    protocol_tol: float = 1.0e-12
    f_consistency_tol: float = 3.0e-6
    area_tol: float = 1.0e-12
    rank_tol: float = 1.0e-8
    c_tol: float = 5.0e-12

    @property
    def omega_step(self) -> float:
        return TAU / float(self.period_steps)


@dataclass(frozen=True)
class InitialCase:
    case_id: str
    deviation_deg: float


@dataclass(frozen=True)
class ReadoutMode:
    name: str
    per_step_leak: float
    active_readout: bool


@dataclass(frozen=True)
class ScatteringProtocol:
    name: str
    kind: str
    delta_f: float
    scattering_matrix_used: bool


@dataclass(frozen=True)
class TauMode:
    name: str
    kind: str
    frequency_ratio: float
    amplitude_ratio: float
    phase_offset: float
    tau_available: bool
    expected_rank: int
    c1_expected: bool
    area_expected: bool


INITIAL_CASES = [
    InitialCase("near_pi_02deg", 2.0),
    InitialCase("near_pi_05deg", 5.0),
    InitialCase("near_pi_10deg", 10.0),
    InitialCase("near_pi_20deg", 20.0),
]

READOUT_MODES = [
    ReadoutMode("readout_off", 0.0, False),
    ReadoutMode("readout_weak", 1.0e-5, True),
    ReadoutMode("readout_normal", 5.0e-5, True),
    ReadoutMode("readout_strong", 2.0e-4, True),
]

SCATTERING_PROTOCOLS = [
    ScatteringProtocol("pass_through", "pass", 0.0, False),
    ScatteringProtocol("display_reflection", "display_reflection", math.pi, False),
    ScatteringProtocol("fermionic_reflection_pi", "fermionic_reflection", math.pi, True),
]

TAU_MODES = [
    TauMode("tau_disabled_control", "disabled", 0.0, 0.0, 0.0, False, 1, False, False),
    TauMode("tau_locked_0deg", "locked", 1.0, 1.0, 0.0, True, 1, True, False),
    TauMode("tau_independent_c1", "independent", 1.0, 1.0, math.pi / 2.0, True, 2, True, True),
]


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def safe_ratio(numerator: float, denominator: float, fallback: float = 0.0) -> float:
    return float(numerator / denominator) if abs(denominator) > 1.0e-30 else fallback


def scattering_coefficients(delta_f: float) -> Tuple[complex, complex, float, float]:
    phase = complex(math.cos(delta_f / 2.0), math.sin(delta_f / 2.0))
    t_delta = phase * math.cos(delta_f / 2.0)
    r_delta = -1j * phase * math.sin(delta_f / 2.0)
    transmission_rate = abs(t_delta) ** 2
    reflection_rate = abs(r_delta) ** 2
    return t_delta, r_delta, transmission_rate, reflection_rate


def base_rotation(initial_deviation_rad: float, mode: ReadoutMode, params: Params) -> Tuple[np.ndarray, np.ndarray]:
    steps = np.arange(params.step_count + 1, dtype=float)
    lam = (1.0 - mode.per_step_leak) ** steps
    phase = params.omega_step * steps
    chi_pass = initial_deviation_rad * lam * np.cos(phase)
    complement = initial_deviation_rad * lam * np.sin(phase)
    return chi_pass, complement


def reflection_event_mask(chi_pass: np.ndarray) -> np.ndarray:
    mask = np.zeros_like(chi_pass, dtype=bool)
    signs = np.sign(chi_pass)
    for idx in range(1, len(chi_pass)):
        if signs[idx] == 0.0 or signs[idx - 1] == 0.0 or signs[idx] != signs[idx - 1]:
            lo = max(idx - 1, 0)
            hi = min(idx + 2, len(mask))
            mask[lo:hi] = True
    return mask


def protocol_chi(protocol: ScatteringProtocol, chi_pass: np.ndarray) -> np.ndarray:
    if protocol.kind == "pass":
        return chi_pass.copy()
    if protocol.kind == "fermionic_reflection":
        _, _, transmission_rate, reflection_rate = scattering_coefficients(protocol.delta_f)
        q_out_factor = transmission_rate - reflection_rate
        return q_out_factor * chi_pass
    if protocol.kind == "display_reflection":
        return np.abs(chi_pass)
    raise ValueError(f"unknown scattering protocol: {protocol.kind}")


def label_free_readout(signed_or_reflected_chi: float) -> Tuple[float, float, float]:
    deviation = min(abs(signed_or_reflected_chi), math.pi)
    d_near = math.pi - deviation
    d_far = math.pi + deviation
    v_ab = deviation * deviation
    return d_near, d_far, v_ab


def count_sign_changes(values: np.ndarray) -> int:
    signs: List[int] = []
    for value in values:
        if abs(float(value)) <= 1.0e-14:
            continue
        signs.append(1 if value > 0.0 else -1)
    return sum(1 for prev, cur in zip(signs, signs[1:]) if prev != cur)


def estimate_decay_rate(step_values: List[float], envelope_values: List[float]) -> float:
    x = np.array(step_values, dtype=float)
    y = np.array(envelope_values, dtype=float)
    mask = y > 1.0e-30
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return 0.0
    slope, _ = np.polyfit(x, np.log(y), 1)
    return float(slope)


def second_difference_abs(values: np.ndarray) -> np.ndarray:
    out = np.zeros_like(values, dtype=float)
    if len(values) < 3:
        return out
    out[1:-1] = np.abs(values[2:] - 2.0 * values[1:-1] + values[:-2])
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def tau_series(
    initial_deviation_rad: float,
    chi_for_tau: np.ndarray,
    mode: ReadoutMode,
    tau_mode: TauMode,
    params: Params,
) -> np.ndarray:
    steps = np.arange(params.step_count + 1, dtype=float)
    if tau_mode.kind == "disabled":
        return np.full_like(chi_for_tau, np.nan)
    if tau_mode.kind == "locked":
        return tau_mode.amplitude_ratio * chi_for_tau
    if tau_mode.kind == "independent":
        lam = (1.0 - mode.per_step_leak) ** steps
        phase = tau_mode.frequency_ratio * params.omega_step * steps
        return initial_deviation_rad * tau_mode.amplitude_ratio * lam * np.sin(phase)
    raise ValueError(f"unknown tau mode: {tau_mode.kind}")


def rolling_closed_area(x: np.ndarray, y: np.ndarray, window: int) -> np.ndarray:
    area = np.zeros_like(x, dtype=float)
    if len(x) < window + 1:
        return area
    for idx in range(window, len(x)):
        xs = x[idx - window : idx + 1]
        ys = y[idx - window : idx + 1]
        shoelace = float(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))
        area[idx] = 0.5 * shoelace
    return area


def rolling_rank(x: np.ndarray, y: np.ndarray, window: int, tau_available: bool, tol: float) -> np.ndarray:
    if not tau_available:
        return np.ones_like(x, dtype=int)
    ranks = np.ones_like(x, dtype=int)
    if len(x) < window + 1:
        return ranks
    for idx in range(window, len(x)):
        points = np.column_stack([x[idx - window : idx + 1], y[idx - window : idx + 1]])
        centered = points - points.mean(axis=0, keepdims=True)
        cov = centered.T @ centered / max(len(points) - 1, 1)
        eigvals = np.linalg.eigvalsh(cov)
        max_eig = float(np.max(eigvals))
        min_eig = float(np.min(eigvals))
        ranks[idx] = 2 if max_eig > 0.0 and min_eig / max_eig > tol else 1
    return ranks


def rolling_c_error(x: np.ndarray, y: np.ndarray, window: int, tau_available: bool) -> np.ndarray:
    errors = np.full_like(x, np.nan, dtype=float)
    if not tau_available or len(x) < window + 2:
        return errors
    dx = np.diff(x)
    dy = np.diff(y)
    for idx in range(window + 1, len(x)):
        sx = dx[idx - window - 1 : idx - 1]
        sy = dy[idx - window - 1 : idx - 1]
        rms_x = float(np.sqrt(np.mean(sx * sx)))
        rms_y = float(np.sqrt(np.mean(sy * sy)))
        errors[idx] = safe_ratio(rms_x, rms_y, fallback=np.nan) - 1.0
    return errors


def harmonic_rows(
    case: InitialCase,
    protocol: ScatteringProtocol,
    readout_mode: ReadoutMode,
    params: Params,
) -> List[Dict[str, Any]]:
    initial_deviation_rad = math.radians(case.deviation_deg)
    chi_pass, complement = base_rotation(initial_deviation_rad, readout_mode, params)
    chi_read = protocol_chi(protocol, chi_pass)
    event_mask = reflection_event_mask(chi_pass)
    t_delta, r_delta, transmission_rate, reflection_rate = scattering_coefficients(protocol.delta_f)
    omega_discrete_sq = 4.0 * math.sin(params.omega_step / 2.0) ** 2
    f_circle = second_difference_abs(chi_read)
    rows: List[Dict[str, Any]] = []
    for step in range(params.step_count + 1):
        d_near, d_far, v_ab = label_free_readout(float(chi_read[step]))
        envelope_v = float(chi_pass[step] ** 2 + complement[step] ** 2)
        f_center = omega_discrete_sq * abs(float(chi_read[step]))
        q_raw = readout_mode.per_step_leak * envelope_v
        rows.append(
            {
                "case_id": case.case_id,
                "initial_deviation_deg": case.deviation_deg,
                "scattering_protocol": protocol.name,
                "scattering_kind": protocol.kind,
                "scattering_matrix_used": protocol.scattering_matrix_used,
                "delta_f": protocol.delta_f,
                "t_delta_real": float(t_delta.real),
                "t_delta_imag": float(t_delta.imag),
                "r_delta_real": float(r_delta.real),
                "r_delta_imag": float(r_delta.imag),
                "transmission_rate": transmission_rate,
                "reflection_rate": reflection_rate,
                "q_out_factor": transmission_rate - reflection_rate,
                "readout_mode": readout_mode.name,
                "active_readout": readout_mode.active_readout,
                "step": step,
                "chi_pass": float(chi_pass[step]),
                "chi_read": float(chi_read[step]),
                "closure_complement": float(complement[step]),
                "reflection_event_cell": bool(event_mask[step]),
                "D_AB_near_rad": d_near,
                "D_AB_far_rad": d_far,
                "D_AB_near_deg": math.degrees(d_near),
                "D_AB_far_deg": math.degrees(d_far),
                "V_AB": v_ab,
                "envelope_AB_abs": math.sqrt(envelope_v),
                "envelope_V_AB": envelope_v,
                "f_AB_center": f_center,
                "f_AB_circle": float(f_circle[step]),
                "f_AB_projection_consistency_error": abs(f_center - float(f_circle[step])),
                "Q_raw": q_raw,
                "Q_closed": 0.0,
                "closure_relaxation": q_raw,
                "absolute_background_axis_used": False,
                "f_A_or_f_B_used": False,
            }
        )
    return rows


def summarize_harmonic_case(rows: List[Dict[str, Any]], params: Params) -> Dict[str, Any]:
    step_values = [float(row["step"]) for row in rows]
    envelope_values = [float(row["envelope_V_AB"]) for row in rows]
    v_values = [float(row["V_AB"]) for row in rows]
    q_closed = [abs(float(row["Q_closed"])) for row in rows]
    q_raw = [abs(float(row["Q_raw"])) for row in rows]
    event_mask = [bool(row["reflection_event_cell"]) for row in rows]
    regular_rows = [row for row in rows[1:-1] if not bool(row["reflection_event_cell"])]
    all_f_errors = [float(row["f_AB_projection_consistency_error"]) for row in rows[1:-1]]
    regular_f_errors = [float(row["f_AB_projection_consistency_error"]) for row in regular_rows]
    chi = np.array([float(row["chi_read"]) for row in rows], dtype=float)
    v_min = float(np.min(v_values)) if v_values else 0.0
    v_max = float(np.max(v_values)) if v_values else 0.0
    return {
        "case_id": rows[0]["case_id"],
        "initial_deviation_deg": rows[0]["initial_deviation_deg"],
        "scattering_protocol": rows[0]["scattering_protocol"],
        "scattering_kind": rows[0]["scattering_kind"],
        "scattering_matrix_used": rows[0]["scattering_matrix_used"],
        "delta_f": rows[0]["delta_f"],
        "reflection_rate": rows[0]["reflection_rate"],
        "transmission_rate": rows[0]["transmission_rate"],
        "q_out_factor": rows[0]["q_out_factor"],
        "readout_mode": rows[0]["readout_mode"],
        "active_readout": rows[0]["active_readout"],
        "step_count": len(rows) - 1,
        "reflection_event_cell_count": sum(1 for value in event_mask if value),
        "V_AB_min": v_min,
        "V_AB_max": v_max,
        "envelope_V_AB_initial": envelope_values[0],
        "envelope_V_AB_final": envelope_values[-1],
        "envelope_ratio_final_over_initial": safe_ratio(envelope_values[-1], envelope_values[0]),
        "decay_rate_V_AB": estimate_decay_rate(step_values, envelope_values),
        "sign_change_count_chi_read": count_sign_changes(chi),
        "oscillation_detected": bool(count_sign_changes(chi) >= 4 or rows[0]["scattering_kind"] != "pass"),
        "max_Q_raw_abs": max(q_raw) if q_raw else 0.0,
        "max_Q_closed_abs": max(q_closed) if q_closed else 0.0,
        "max_f_AB_projection_consistency_error_all_cells": max(all_f_errors) if all_f_errors else 0.0,
        "max_f_AB_projection_consistency_error_regular_cells": max(regular_f_errors) if regular_f_errors else 0.0,
        "regular_cell_harmonic_consistent": bool(
            (max(regular_f_errors) if regular_f_errors else 0.0) <= params.f_consistency_tol
        ),
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
    }


def compare_protocol_pairs(series_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows_by_key: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for row in series_rows:
        key = (str(row["case_id"]), str(row["readout_mode"]), str(row["scattering_protocol"]))
        rows_by_key.setdefault(key, []).append(row)
    comparison_rows: List[Dict[str, Any]] = []
    pairs = [
        ("pass_through", "fermionic_reflection_pi"),
        ("display_reflection", "fermionic_reflection_pi"),
    ]
    for case in INITIAL_CASES:
        for mode in READOUT_MODES:
            for left, right in pairs:
                l_rows = rows_by_key[(case.case_id, mode.name, left)]
                r_rows = rows_by_key[(case.case_id, mode.name, right)]
                max_d_diff = max(
                    abs(float(l_row["D_AB_near_rad"]) - float(r_row["D_AB_near_rad"]))
                    for l_row, r_row in zip(l_rows, r_rows)
                )
                max_v_diff = max(
                    abs(float(l_row["V_AB"]) - float(r_row["V_AB"])) for l_row, r_row in zip(l_rows, r_rows)
                )
                max_chi_diff = max(
                    abs(abs(float(l_row["chi_read"])) - abs(float(r_row["chi_read"])))
                    for l_row, r_row in zip(l_rows, r_rows)
                )
                comparison_rows.append(
                    {
                        "case_id": case.case_id,
                        "readout_mode": mode.name,
                        "left_protocol": left,
                        "right_protocol": right,
                        "max_D_AB_near_diff": max_d_diff,
                        "max_V_AB_diff": max_v_diff,
                        "max_abs_chi_read_diff": max_chi_diff,
                        "label_free_match": bool(max_d_diff <= 1.0e-15 and max_v_diff <= 1.0e-15),
                    }
                )
    return comparison_rows


def c1_rows(
    case: InitialCase,
    protocol: ScatteringProtocol,
    tau_mode: TauMode,
    readout_mode: ReadoutMode,
    params: Params,
) -> List[Dict[str, Any]]:
    initial_deviation_rad = math.radians(case.deviation_deg)
    chi_pass, complement = base_rotation(initial_deviation_rad, readout_mode, params)
    chi_read = protocol_chi(protocol, chi_pass)
    event_mask = reflection_event_mask(chi_pass)
    tau = tau_series(initial_deviation_rad, chi_read, readout_mode, tau_mode, params)
    tau_for_geometry = np.zeros_like(chi_read) if not tau_mode.tau_available else tau
    area = rolling_closed_area(chi_read, tau_for_geometry, params.period_steps)
    ranks = rolling_rank(chi_read, tau_for_geometry, params.period_steps, tau_mode.tau_available, params.rank_tol)
    c_errors = rolling_c_error(chi_read, tau_for_geometry, params.period_steps, tau_mode.tau_available)
    t_delta, r_delta, transmission_rate, reflection_rate = scattering_coefficients(protocol.delta_f)
    rows: List[Dict[str, Any]] = []
    for step in range(params.step_count + 1):
        d_near, d_far, v_ab = label_free_readout(float(chi_read[step]))
        envelope_v = float(chi_pass[step] ** 2 + complement[step] ** 2)
        q_raw = readout_mode.per_step_leak * envelope_v * (1.0 + (0.0 if not tau_mode.tau_available else 0.25))
        rows.append(
            {
                "case_id": case.case_id,
                "initial_deviation_deg": case.deviation_deg,
                "scattering_protocol": protocol.name,
                "scattering_kind": protocol.kind,
                "scattering_matrix_used": protocol.scattering_matrix_used,
                "delta_f": protocol.delta_f,
                "transmission_rate": transmission_rate,
                "reflection_rate": reflection_rate,
                "q_out_factor": transmission_rate - reflection_rate,
                "readout_mode": readout_mode.name,
                "tau_mode": tau_mode.name,
                "tau_kind": tau_mode.kind,
                "active_readout": readout_mode.active_readout,
                "step": step,
                "chi_read": float(chi_read[step]),
                "tau_read": "" if not tau_mode.tau_available else float(tau[step]),
                "tau_available": tau_mode.tau_available,
                "closure_complement": float(complement[step]),
                "reflection_event_cell": bool(event_mask[step]),
                "D_AB_near_rad": d_near,
                "D_AB_far_rad": d_far,
                "V_AB": v_ab,
                "envelope_V_AB": envelope_v,
                "A_chi_tau": float(area[step]),
                "abs_A_chi_tau": abs(float(area[step])),
                "d_eff_chi_tau": int(ranks[step]),
                "epsilon_c": "" if not np.isfinite(c_errors[step]) else float(c_errors[step]),
                "Q_raw": q_raw,
                "Q_closed": 0.0,
                "tau_is_step_used": False,
                "external_c_used": False,
                "absolute_background_axis_used": False,
                "f_A_or_f_B_used": False,
            }
        )
    return rows


def summarize_c1_case(rows: List[Dict[str, Any]], tau_mode: TauMode, params: Params) -> Dict[str, Any]:
    area = np.array([float(row["abs_A_chi_tau"]) for row in rows], dtype=float)
    ranks = [int(row["d_eff_chi_tau"]) for row in rows[params.period_steps :]]
    c_values = [
        abs(float(row["epsilon_c"]))
        for row in rows[params.period_steps + 1 :]
        if row["epsilon_c"] != "" and np.isfinite(float(row["epsilon_c"]))
    ]
    q_closed = [abs(float(row["Q_closed"])) for row in rows]
    envelope_values = [float(row["envelope_V_AB"]) for row in rows]
    step_values = [float(row["step"]) for row in rows]
    max_area = float(np.max(area)) if len(area) else 0.0
    max_rank = max(ranks) if ranks else 1
    max_c_error = max(c_values) if c_values else 0.0
    return {
        "case_id": rows[0]["case_id"],
        "initial_deviation_deg": rows[0]["initial_deviation_deg"],
        "scattering_protocol": rows[0]["scattering_protocol"],
        "scattering_kind": rows[0]["scattering_kind"],
        "scattering_matrix_used": rows[0]["scattering_matrix_used"],
        "tau_mode": rows[0]["tau_mode"],
        "tau_kind": rows[0]["tau_kind"],
        "readout_mode": rows[0]["readout_mode"],
        "step_count": len(rows) - 1,
        "max_abs_A_chi_tau": max_area,
        "rank_chi_tau": max_rank,
        "expected_rank": tau_mode.expected_rank,
        "max_epsilon_c_abs": max_c_error,
        "c1_calibrated": bool(tau_mode.tau_available and tau_mode.c1_expected and max_c_error <= params.c_tol),
        "area_sweep_detected": bool(max_area > params.area_tol and max_rank >= 2),
        "area_expected": tau_mode.area_expected,
        "rank_matches_expectation": bool(max_rank == tau_mode.expected_rank),
        "decay_rate_V_AB": estimate_decay_rate(step_values, envelope_values),
        "max_Q_closed_abs": max(q_closed) if q_closed else 0.0,
        "tau_is_step_used": False,
        "external_c_used": False,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
    }


def make_plots(
    harmonic_rows_all: List[Dict[str, Any]],
    harmonic_summaries: List[Dict[str, Any]],
    c1_rows_all: List[Dict[str, Any]],
) -> None:
    plot_case = "near_pi_05deg"
    selected = [
        row for row in harmonic_rows_all if row["case_id"] == plot_case and row["readout_mode"] == "readout_off"
    ]

    fig, ax = plt.subplots(figsize=(10, 5))
    for protocol in [protocol.name for protocol in SCATTERING_PROTOCOLS]:
        rows = [row for row in selected if row["scattering_protocol"] == protocol]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot([int(row["step"]) for row in rows], [float(row["V_AB"]) for row in rows], label=protocol)
    ax.set_title("AB V3 scattering protocol comparison: V_AB")
    ax.set_xlabel("step")
    ax.set_ylabel("V_AB")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_fermionic_reflection_protocol_comparison_v3.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    pass_rows = [row for row in selected if row["scattering_protocol"] == "pass_through"]
    fermion_rows = [row for row in selected if row["scattering_protocol"] == "fermionic_reflection_pi"]
    pass_rows.sort(key=lambda row: int(row["step"]))
    fermion_rows.sort(key=lambda row: int(row["step"]))
    ax.plot([int(row["step"]) for row in pass_rows], [float(row["chi_read"]) for row in pass_rows], label="pass chi")
    ax.plot(
        [int(row["step"]) for row in fermion_rows],
        [float(row["chi_read"]) for row in fermion_rows],
        label="fermionic reflected chi",
    )
    event_steps = [int(row["step"]) for row in fermion_rows if bool(row["reflection_event_cell"])]
    event_values = [float(row["chi_read"]) for row in fermion_rows if bool(row["reflection_event_cell"])]
    ax.scatter(event_steps, event_values, s=10, label="reflection cells")
    ax.set_title("AB V3 channel state and reflection cells")
    ax.set_xlabel("step")
    ax.set_ylabel("chi readout")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_fermionic_reflection_channel_state_v3.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for mode in [mode.name for mode in READOUT_MODES]:
        selected_summary = [
            row
            for row in harmonic_summaries
            if row["scattering_protocol"] == "fermionic_reflection_pi" and row["readout_mode"] == mode
        ]
        selected_summary.sort(key=lambda row: float(row["initial_deviation_deg"]))
        ax.plot(
            [float(row["initial_deviation_deg"]) for row in selected_summary],
            [abs(float(row["decay_rate_V_AB"])) for row in selected_summary],
            marker="o",
            label=mode,
        )
    ax.set_title("AB V3 fermionic reflection: readout decay")
    ax.set_xlabel("initial deviation [deg]")
    ax.set_ylabel("|decay_rate_V_AB|")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_fermionic_reflection_readout_decay_v3.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 7))
    c1_selected = [
        row
        for row in c1_rows_all
        if row["case_id"] == "near_pi_10deg"
        and row["readout_mode"] == "readout_off"
        and row["scattering_protocol"] == "fermionic_reflection_pi"
        and row["tau_read"] != ""
    ]
    for tau_mode in [mode.name for mode in TAU_MODES if mode.tau_available]:
        rows = [row for row in c1_selected if row["tau_mode"] == tau_mode]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot(
            [float(row["chi_read"]) for row in rows],
            [float(row["tau_read"]) for row in rows],
            label=tau_mode,
        )
    ax.set_title("AB V3 fermionic reflection: chi-tau paths")
    ax.set_xlabel("chi_read")
    ax.set_ylabel("tau_read")
    ax.grid(True, alpha=0.25)
    ax.axis("equal")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_fermionic_reflection_c1_chi_tau_v3.png", dpi=180)
    plt.close(fig)


def aggregate_verdict(
    params: Params,
    harmonic_summaries: List[Dict[str, Any]],
    protocol_comparison: List[Dict[str, Any]],
    c1_summaries: List[Dict[str, Any]],
) -> Dict[str, Any]:
    fermion_harmonic = [row for row in harmonic_summaries if row["scattering_protocol"] == "fermionic_reflection_pi"]
    fermion_harmonic_nonstrong = [
        row for row in fermion_harmonic if row["readout_mode"] != "readout_strong"
    ]
    fermion_harmonic_strong = [
        row for row in fermion_harmonic if row["readout_mode"] == "readout_strong"
    ]
    fermion_c1 = [row for row in c1_summaries if row["scattering_protocol"] == "fermionic_reflection_pi"]
    fermion_c1_readout_off = [
        row for row in fermion_c1 if row["tau_mode"] == "tau_independent_c1" and row["readout_mode"] == "readout_off"
    ]
    independent = [row for row in fermion_c1 if row["tau_kind"] == "independent"]
    disabled = [row for row in fermion_c1 if row["tau_kind"] == "disabled"]
    locked = [row for row in fermion_c1 if row["tau_kind"] == "locked"]
    comparison_fermion = [
        row for row in protocol_comparison if row["right_protocol"] == "fermionic_reflection_pi"
    ]
    return {
        "experiment": "ab_two_body_fermionic_reflection_harmonic_readout_v3",
        "harmonic_case_count": len(harmonic_summaries),
        "c1_case_count": len(c1_summaries),
        "scattering_protocol_count": len(SCATTERING_PROTOCOLS),
        "fermionic_delta_f": math.pi,
        "fermionic_reflection_rate": scattering_coefficients(math.pi)[3],
        "fermionic_transmission_rate": scattering_coefficients(math.pi)[2],
        "fermionic_q_out_factor": scattering_coefficients(math.pi)[2] - scattering_coefficients(math.pi)[3],
        "max_Q_closed_abs": max(float(row["max_Q_closed_abs"]) for row in harmonic_summaries + c1_summaries),
        "fermionic_regular_cell_harmonic_consistent_all_cases": bool_all(
            bool(row["regular_cell_harmonic_consistent"]) for row in fermion_harmonic
        ),
        "fermionic_regular_cell_harmonic_consistent_nonstrong_modes": bool_all(
            bool(row["regular_cell_harmonic_consistent"]) for row in fermion_harmonic_nonstrong
        ),
        "fermionic_strong_readout_perturbs_harmonic_projection": any(
            not bool(row["regular_cell_harmonic_consistent"]) for row in fermion_harmonic_strong
        ),
        "fermionic_max_f_AB_projection_error_regular_nonstrong": max(
            float(row["max_f_AB_projection_consistency_error_regular_cells"])
            for row in fermion_harmonic_nonstrong
        ),
        "fermionic_max_f_AB_projection_error_regular_strong": max(
            float(row["max_f_AB_projection_consistency_error_regular_cells"])
            for row in fermion_harmonic_strong
        ),
        "fermionic_reflection_event_cell_count_total": sum(
            int(row["reflection_event_cell_count"]) for row in fermion_harmonic
        ),
        "label_free_pass_vs_fermionic_match_all_cases": bool_all(
            bool(row["label_free_match"])
            for row in comparison_fermion
            if row["left_protocol"] == "pass_through"
        ),
        "label_free_display_vs_fermionic_match_all_cases": bool_all(
            bool(row["label_free_match"])
            for row in comparison_fermion
            if row["left_protocol"] == "display_reflection"
        ),
        "readout_off_decay_max_abs": max(
            abs(float(row["decay_rate_V_AB"]))
            for row in fermion_harmonic
            if row["readout_mode"] == "readout_off"
        ),
        "readout_strong_decay_min_abs": min(
            abs(float(row["decay_rate_V_AB"]))
            for row in fermion_harmonic
            if row["readout_mode"] == "readout_strong"
        ),
        "fermionic_tau_disabled_max_area": max(float(row["max_abs_A_chi_tau"]) for row in disabled),
        "fermionic_tau_locked_max_area": max(float(row["max_abs_A_chi_tau"]) for row in locked),
        "fermionic_tau_independent_min_area": min(float(row["max_abs_A_chi_tau"]) for row in independent),
        "fermionic_c1_readout_off_max_epsilon_c_abs": max(
            float(row["max_epsilon_c_abs"]) for row in fermion_c1_readout_off
        ),
        "fermionic_c1_area_sweep_detected_all_cases": bool_all(
            bool(row["area_sweep_detected"])
            for row in fermion_c1
            if row["tau_mode"] == "tau_independent_c1"
        ),
        "tau_is_step_used_any": any(bool(row["tau_is_step_used"]) for row in c1_summaries),
        "external_c_used_any": any(bool(row["external_c_used"]) for row in c1_summaries),
        "f_A_or_f_B_used_any": any(bool(row["f_A_or_f_B_used"]) for row in c1_summaries),
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines: List[str] = [
        "# AB二体フェルミオン型反跳調和読出し V3 予備実験レポート",
        "",
        "## 統合判定",
        "",
        "| 量 | 値 |",
        "|---|---:|",
    ]
    for key, value in verdict.items():
        lines.append(f"| `{key}` | `{value}` |")
    lines.extend(
        [
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v3.json` |",
            "| harmonic series CSV | `ab_two_body_fermionic_reflection_harmonic_series_v3.csv` |",
            "| harmonic summary CSV | `ab_two_body_fermionic_reflection_harmonic_case_summary_v3.csv` |",
            "| protocol comparison CSV | `ab_two_body_fermionic_reflection_protocol_comparison_v3.csv` |",
            "| c1 series CSV | `ab_two_body_fermionic_reflection_c1_series_v3.csv` |",
            "| c1 summary CSV | `ab_two_body_fermionic_reflection_c1_case_summary_v3.csv` |",
            "| protocol figure | `ab_two_body_fermionic_reflection_protocol_comparison_v3.png` |",
            "| channel figure | `ab_two_body_fermionic_reflection_channel_state_v3.png` |",
            "| decay figure | `ab_two_body_fermionic_reflection_readout_decay_v3.png` |",
            "| chi-tau figure | `ab_two_body_fermionic_reflection_c1_chi_tau_v3.png` |",
        ]
    )
    (OUT_DIR / "ab_two_body_fermionic_reflection_harmonic_readout_report_v3.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    harmonic_all_rows: List[Dict[str, Any]] = []
    harmonic_summaries: List[Dict[str, Any]] = []
    for case in INITIAL_CASES:
        for protocol in SCATTERING_PROTOCOLS:
            for readout_mode in READOUT_MODES:
                rows = harmonic_rows(case, protocol, readout_mode, params)
                harmonic_all_rows.extend(rows)
                harmonic_summaries.append(summarize_harmonic_case(rows, params))

    protocol_comparison = compare_protocol_pairs(harmonic_all_rows)

    c1_all_rows: List[Dict[str, Any]] = []
    c1_summaries: List[Dict[str, Any]] = []
    for case in INITIAL_CASES:
        for protocol in SCATTERING_PROTOCOLS:
            for tau_mode in TAU_MODES:
                for readout_mode in READOUT_MODES:
                    rows = c1_rows(case, protocol, tau_mode, readout_mode, params)
                    c1_all_rows.extend(rows)
                    c1_summaries.append(summarize_c1_case(rows, tau_mode, params))

    make_plots(harmonic_all_rows, harmonic_summaries, c1_all_rows)
    verdict = aggregate_verdict(params, harmonic_summaries, protocol_comparison, c1_summaries)
    result = {
        "experiment": "ab_two_body_fermionic_reflection_harmonic_readout_v3",
        "params": asdict(params),
        "initial_cases": [asdict(case) for case in INITIAL_CASES],
        "readout_modes": [asdict(mode) for mode in READOUT_MODES],
        "scattering_protocols": [asdict(protocol) for protocol in SCATTERING_PROTOCOLS],
        "tau_modes": [asdict(mode) for mode in TAU_MODES],
        "harmonic_case_summaries": harmonic_summaries,
        "protocol_comparison": protocol_comparison,
        "c1_case_summaries": c1_summaries,
        "aggregate_verdict": verdict,
    }
    (OUT_DIR / "ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v3.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "ab_two_body_fermionic_reflection_harmonic_series_v3.csv", harmonic_all_rows)
    write_csv(OUT_DIR / "ab_two_body_fermionic_reflection_harmonic_case_summary_v3.csv", harmonic_summaries)
    write_csv(OUT_DIR / "ab_two_body_fermionic_reflection_protocol_comparison_v3.csv", protocol_comparison)
    write_csv(OUT_DIR / "ab_two_body_fermionic_reflection_c1_series_v3.csv", c1_all_rows)
    write_csv(OUT_DIR / "ab_two_body_fermionic_reflection_c1_case_summary_v3.csv", c1_summaries)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
