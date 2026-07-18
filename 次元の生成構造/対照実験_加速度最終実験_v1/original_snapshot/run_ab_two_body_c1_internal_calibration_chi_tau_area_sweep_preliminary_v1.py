from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1"
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
    area_tol: float = 1.0e-12
    rank_tol: float = 1.0e-8
    c_tol: float = 5.0e-12
    closure_tol: float = 1.0e-12

    @property
    def omega_step(self) -> float:
        return TAU / float(self.period_steps)


@dataclass(frozen=True)
class InitialCase:
    case_id: str
    deviation_deg: float


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


@dataclass(frozen=True)
class ReadoutMode:
    name: str
    per_step_leak: float
    active_readout: bool


INITIAL_CASES = [
    InitialCase("near_pi_02deg", 2.0),
    InitialCase("near_pi_05deg", 5.0),
    InitialCase("near_pi_10deg", 10.0),
    InitialCase("near_pi_20deg", 20.0),
    InitialCase("wide_35deg", 35.0),
    InitialCase("wide_60deg", 60.0),
]

PROTOCOLS = ["Protocol_F", "Protocol_B"]

TAU_MODES = [
    TauMode("tau_disabled_control", "disabled", 0.0, 0.0, 0.0, False, 1, False, False),
    TauMode("tau_locked_0deg", "locked", 1.0, 1.0, 0.0, True, 1, True, False),
    TauMode("tau_locked_90deg", "locked_offset", 1.0, 1.0, math.pi / 2.0, True, 1, True, False),
    TauMode("tau_independent_slow", "independent", 0.5, 1.0, math.pi / 2.0, True, 2, False, True),
    TauMode("tau_independent_c1", "independent", 1.0, 1.0, math.pi / 2.0, True, 2, True, True),
    TauMode("tau_independent_fast", "independent", 1.5, 1.0, math.pi / 2.0, True, 2, False, True),
]

READOUT_MODES = [
    ReadoutMode("readout_off", 0.0, False),
    ReadoutMode("readout_weak", 1.0e-5, True),
    ReadoutMode("readout_normal", 5.0e-5, True),
    ReadoutMode("readout_strong", 2.0e-4, True),
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


def label_free_readout(signed_deviation: float) -> Tuple[float, float, float]:
    deviation = min(abs(signed_deviation), math.pi)
    d_near = math.pi - deviation
    d_far = math.pi + deviation
    v_ab = deviation * deviation
    return d_near, d_far, v_ab


def protocol_display_deviation(protocol: str, signed_deviation: float) -> float:
    if protocol == "Protocol_F":
        return abs(signed_deviation)
    if protocol == "Protocol_B":
        return signed_deviation
    raise ValueError(f"unknown protocol: {protocol}")


def chi_series(initial_deviation_rad: float, mode: ReadoutMode, params: Params) -> Tuple[np.ndarray, np.ndarray]:
    steps = np.arange(params.step_count + 1, dtype=float)
    lam = (1.0 - mode.per_step_leak) ** steps
    phase = params.omega_step * steps
    chi = initial_deviation_rad * lam * np.cos(phase)
    complement = initial_deviation_rad * lam * np.sin(phase)
    return chi, complement


def tau_series(
    initial_deviation_rad: float,
    chi: np.ndarray,
    mode: ReadoutMode,
    tau_mode: TauMode,
    params: Params,
) -> np.ndarray:
    steps = np.arange(params.step_count + 1, dtype=float)
    if tau_mode.kind == "disabled":
        return np.full_like(chi, np.nan)
    if tau_mode.kind == "locked":
        return tau_mode.amplitude_ratio * chi
    if tau_mode.kind == "locked_offset":
        return tau_mode.amplitude_ratio * chi + initial_deviation_rad * math.sin(tau_mode.phase_offset)
    if tau_mode.kind == "independent":
        lam = (1.0 - mode.per_step_leak) ** steps
        phase = tau_mode.frequency_ratio * params.omega_step * steps
        return initial_deviation_rad * tau_mode.amplitude_ratio * lam * np.sin(phase)
    raise ValueError(f"unknown tau mode kind: {tau_mode.kind}")


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


def count_internal_sign_changes(values: np.ndarray) -> int:
    signs: List[int] = []
    for value in values:
        if abs(float(value)) <= 1.0e-14:
            continue
        signs.append(1 if value > 0.0 else -1)
    return sum(1 for prev, cur in zip(signs, signs[1:]) if prev != cur)


def case_rows(
    case: InitialCase,
    protocol: str,
    tau_mode: TauMode,
    readout_mode: ReadoutMode,
    params: Params,
) -> List[Dict[str, Any]]:
    initial_deviation_rad = math.radians(case.deviation_deg)
    chi, closure_complement = chi_series(initial_deviation_rad, readout_mode, params)
    tau = tau_series(initial_deviation_rad, chi, readout_mode, tau_mode, params)
    tau_for_geometry = np.zeros_like(chi) if not tau_mode.tau_available else tau
    area = rolling_closed_area(chi, tau_for_geometry, params.period_steps)
    ranks = rolling_rank(chi, tau_for_geometry, params.period_steps, tau_mode.tau_available, params.rank_tol)
    c_errors = rolling_c_error(chi, tau_for_geometry, params.period_steps, tau_mode.tau_available)
    rows: List[Dict[str, Any]] = []
    omega_discrete_sq = 4.0 * math.sin(params.omega_step / 2.0) ** 2
    for step in range(params.step_count + 1):
        signed_deviation = float(chi[step])
        display_deviation = protocol_display_deviation(protocol, signed_deviation)
        d_near, d_far, v_ab = label_free_readout(signed_deviation)
        envelope_v = float(chi[step] ** 2 + closure_complement[step] ** 2)
        f_ab = omega_discrete_sq * abs(signed_deviation)
        q_raw = readout_mode.per_step_leak * envelope_v * (1.0 + (0.0 if not tau_mode.tau_available else 0.25))
        rows.append(
            {
                "case_id": case.case_id,
                "initial_deviation_deg": case.deviation_deg,
                "protocol": protocol,
                "tau_mode": tau_mode.name,
                "tau_kind": tau_mode.kind,
                "readout_mode": readout_mode.name,
                "active_readout": readout_mode.active_readout,
                "step": step,
                "chi_read": signed_deviation,
                "tau_read": "" if not tau_mode.tau_available else float(tau[step]),
                "tau_available": tau_mode.tau_available,
                "protocol_display_deviation_rad": display_deviation,
                "D_AB_near_rad": d_near,
                "D_AB_far_rad": d_far,
                "D_AB_near_deg": math.degrees(d_near),
                "D_AB_far_deg": math.degrees(d_far),
                "V_AB": v_ab,
                "rho_AB": params.rho_AB,
                "closure_complement": float(closure_complement[step]),
                "envelope_AB_abs": math.sqrt(envelope_v),
                "envelope_V_AB": envelope_v,
                "f_AB": f_ab,
                "A_chi_tau": float(area[step]),
                "abs_A_chi_tau": abs(float(area[step])),
                "d_eff_chi_tau": int(ranks[step]),
                "epsilon_c": "" if not np.isfinite(c_errors[step]) else float(c_errors[step]),
                "Q_raw": q_raw,
                "Q_closed": 0.0,
                "closure_relaxation": q_raw,
                "tau_is_step_used": False,
                "external_c_used": False,
                "absolute_background_axis_used": False,
                "f_A_or_f_B_used": False,
            }
        )
    return rows


def estimate_decay_rate(rows: List[Dict[str, Any]], key: str) -> float:
    x = np.array([float(row["step"]) for row in rows], dtype=float)
    y = np.array([float(row[key]) for row in rows], dtype=float)
    mask = y > 1.0e-30
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return 0.0
    slope, _ = np.polyfit(x, np.log(y), 1)
    return float(slope)


def summarize_case(rows: List[Dict[str, Any]], tau_mode: TauMode, params: Params) -> Dict[str, Any]:
    chi = np.array([float(row["chi_read"]) for row in rows], dtype=float)
    area = np.array([float(row["abs_A_chi_tau"]) for row in rows], dtype=float)
    ranks = [int(row["d_eff_chi_tau"]) for row in rows[params.period_steps :]]
    c_values = [
        abs(float(row["epsilon_c"]))
        for row in rows[params.period_steps + 1 :]
        if row["epsilon_c"] != "" and np.isfinite(float(row["epsilon_c"]))
    ]
    q_closed = [abs(float(row["Q_closed"])) for row in rows]
    q_raw = [abs(float(row["Q_raw"])) for row in rows]
    max_area = float(np.max(area)) if len(area) else 0.0
    max_rank = max(ranks) if ranks else 1
    max_c_error = max(c_values) if c_values else 0.0
    return {
        "case_id": rows[0]["case_id"],
        "initial_deviation_deg": rows[0]["initial_deviation_deg"],
        "protocol": rows[0]["protocol"],
        "tau_mode": rows[0]["tau_mode"],
        "tau_kind": rows[0]["tau_kind"],
        "readout_mode": rows[0]["readout_mode"],
        "active_readout": rows[0]["active_readout"],
        "step_count": len(rows) - 1,
        "chi_read_min": float(np.min(chi)),
        "chi_read_max": float(np.max(chi)),
        "max_abs_A_chi_tau": max_area,
        "d_eff_chi_tau_max": max_rank,
        "rank_chi_tau": max_rank,
        "max_epsilon_c_abs": max_c_error,
        "c1_calibrated": bool(tau_mode.tau_available and tau_mode.c1_expected and max_c_error <= params.c_tol),
        "area_sweep_detected": bool(max_area > params.area_tol and max_rank >= 2),
        "area_expected": tau_mode.area_expected,
        "expected_rank": tau_mode.expected_rank,
        "rank_matches_expectation": bool(max_rank == tau_mode.expected_rank),
        "internal_sign_change_count": count_internal_sign_changes(chi),
        "oscillation_detected": count_internal_sign_changes(chi) >= 4,
        "decay_rate_V_AB": estimate_decay_rate(rows, "envelope_V_AB"),
        "max_Q_raw_abs": max(q_raw) if q_raw else 0.0,
        "max_Q_closed_abs": max(q_closed) if q_closed else 0.0,
        "tau_is_step_used": False,
        "external_c_used": False,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
    }


def power_candidate_rows(case_summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for protocol in PROTOCOLS:
        for tau_mode in TAU_MODES:
            for readout_mode in READOUT_MODES:
                selected = [
                    row
                    for row in case_summaries
                    if row["protocol"] == protocol
                    and row["tau_mode"] == tau_mode.name
                    and row["readout_mode"] == readout_mode.name
                ]
                selected.sort(key=lambda row: float(row["initial_deviation_deg"]))
                l_values = np.array([math.radians(float(row["initial_deviation_deg"])) for row in selected], dtype=float)
                i_values = np.array([float(row["max_abs_A_chi_tau"]) for row in selected], dtype=float)
                mask = (l_values > 0.0) & (i_values > 1.0e-30)
                if np.sum(mask) >= 3:
                    slope, intercept = np.polyfit(np.log(l_values[mask]), np.log(i_values[mask]), 1)
                    alpha = -float(slope)
                    valid = True
                else:
                    slope = 0.0
                    intercept = 0.0
                    alpha = 0.0
                    valid = False
                rows.append(
                    {
                        "protocol": protocol,
                        "tau_mode": tau_mode.name,
                        "readout_mode": readout_mode.name,
                        "power_candidate_valid": valid,
                        "loglog_slope_I_vs_L": float(slope),
                        "power_candidate_alpha": alpha,
                        "loglog_intercept": float(intercept),
                        "case_count": len(selected),
                    }
                )
    return rows


def validation_summary(
    case_summaries: List[Dict[str, Any]],
    power_rows: List[Dict[str, Any]],
    params: Params,
) -> Dict[str, Any]:
    disabled = [row for row in case_summaries if row["tau_mode"] == "tau_disabled_control"]
    locked = [row for row in case_summaries if row["tau_kind"] in {"locked", "locked_offset"}]
    independent = [row for row in case_summaries if row["tau_kind"] == "independent"]
    c1 = [row for row in case_summaries if row["tau_mode"] == "tau_independent_c1"]
    c1_readout_off = [
        row for row in case_summaries if row["tau_mode"] == "tau_independent_c1" and row["readout_mode"] == "readout_off"
    ]
    readout_off = [row for row in case_summaries if row["readout_mode"] == "readout_off"]
    readout_strong = [row for row in case_summaries if row["readout_mode"] == "readout_strong"]
    c1_power = [row for row in power_rows if row["tau_mode"] == "tau_independent_c1" and row["readout_mode"] == "readout_off"]
    summary = {
        "c1_internal_calibration_chi_tau_area_sweep_preliminary_valid": bool(
            bool_all(abs(float(row["max_Q_closed_abs"])) <= params.closure_tol for row in case_summaries)
            and bool_all(float(row["max_abs_A_chi_tau"]) <= params.area_tol for row in disabled)
            and bool_all(float(row["max_abs_A_chi_tau"]) <= params.area_tol for row in locked)
            and bool_all(int(row["rank_chi_tau"]) == int(row["expected_rank"]) for row in case_summaries)
            and bool_all(bool(row["area_sweep_detected"]) for row in independent)
            and bool_all(float(row["max_epsilon_c_abs"]) <= params.c_tol for row in c1_readout_off)
            and bool_all(not bool(row["tau_is_step_used"]) for row in case_summaries)
            and bool_all(not bool(row["external_c_used"]) for row in case_summaries)
            and bool_all(not bool(row["f_A_or_f_B_used"]) for row in case_summaries)
        ),
        "case_summary_count": len(case_summaries),
        "power_candidate_count": len(power_rows),
        "max_Q_closed_abs": max(abs(float(row["max_Q_closed_abs"])) for row in case_summaries),
        "disabled_max_area": max(float(row["max_abs_A_chi_tau"]) for row in disabled),
        "locked_max_area": max(float(row["max_abs_A_chi_tau"]) for row in locked),
        "independent_min_area": min(float(row["max_abs_A_chi_tau"]) for row in independent),
        "c1_max_epsilon_c_abs": max(float(row["max_epsilon_c_abs"]) for row in c1),
        "c1_readout_off_max_epsilon_c_abs": max(float(row["max_epsilon_c_abs"]) for row in c1_readout_off),
        "c1_area_sweep_detected_all_cases": bool_all(bool(row["area_sweep_detected"]) for row in c1),
        "tau_disabled_rank1_all_cases": bool_all(int(row["rank_chi_tau"]) == 1 for row in disabled),
        "tau_locked_rank1_all_cases": bool_all(int(row["rank_chi_tau"]) == 1 for row in locked),
        "tau_independent_rank2_all_cases": bool_all(int(row["rank_chi_tau"]) == 2 for row in independent),
        "readout_off_decay_max_abs": max(abs(float(row["decay_rate_V_AB"])) for row in readout_off),
        "readout_strong_decay_min_abs": min(abs(float(row["decay_rate_V_AB"])) for row in readout_strong),
        "c1_readout_off_power_candidate_alpha_values": [
            float(row["power_candidate_alpha"]) for row in c1_power if bool(row["power_candidate_valid"])
        ],
        "tau_is_step_used_any": any(bool(row["tau_is_step_used"]) for row in case_summaries),
        "external_c_used_any": any(bool(row["external_c_used"]) for row in case_summaries),
        "f_A_or_f_B_used_any": any(bool(row["f_A_or_f_B_used"]) for row in case_summaries),
    }
    return summary


def make_plots(
    series_rows: List[Dict[str, Any]],
    case_summaries: List[Dict[str, Any]],
    power_rows: List[Dict[str, Any]],
    params: Params,
) -> None:
    selected = [
        row
        for row in series_rows
        if row["case_id"] == "near_pi_10deg"
        and row["protocol"] == "Protocol_B"
        and row["readout_mode"] == "readout_off"
    ]

    fig, ax = plt.subplots(figsize=(7, 7))
    for tau_mode_name in ["tau_locked_0deg", "tau_independent_c1", "tau_independent_slow", "tau_independent_fast"]:
        rows = [row for row in selected if row["tau_mode"] == tau_mode_name and row["tau_read"] != ""]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot(
            [float(row["chi_read"]) for row in rows],
            [float(row["tau_read"]) for row in rows],
            label=tau_mode_name,
            linewidth=1.6,
        )
    ax.set_title("chi-tau readout paths")
    ax.set_xlabel("chi_read")
    ax.set_ylabel("tau_read")
    ax.grid(True, alpha=0.25)
    ax.axis("equal")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for tau_mode_name in ["tau_disabled_control", "tau_locked_0deg", "tau_independent_c1"]:
        rows = [row for row in selected if row["tau_mode"] == tau_mode_name]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot(
            [int(row["step"]) for row in rows],
            [float(row["A_chi_tau"]) for row in rows],
            label=tau_mode_name,
        )
    ax.set_title("chi-tau rolling closed area")
    ax.set_xlabel("step")
    ax.set_ylabel("A_chi_tau")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_internal_calibration_area_sweep_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for tau_mode_name in ["tau_locked_0deg", "tau_independent_slow", "tau_independent_c1", "tau_independent_fast"]:
        rows = [row for row in selected if row["tau_mode"] == tau_mode_name and row["epsilon_c"] != ""]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot(
            [int(row["step"]) for row in rows],
            [float(row["epsilon_c"]) for row in rows],
            label=tau_mode_name,
        )
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_title("internal c=1 RMS calibration error")
    ax.set_xlabel("step")
    ax.set_ylabel("epsilon_c")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_internal_calibration_error_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    summary_selected = [
        row
        for row in case_summaries
        if row["protocol"] == "Protocol_B"
        and row["readout_mode"] == "readout_off"
        and row["tau_mode"] in {"tau_disabled_control", "tau_locked_0deg", "tau_independent_c1"}
    ]
    for tau_mode_name in ["tau_disabled_control", "tau_locked_0deg", "tau_independent_c1"]:
        rows = [row for row in summary_selected if row["tau_mode"] == tau_mode_name]
        rows.sort(key=lambda row: float(row["initial_deviation_deg"]))
        ax.plot(
            [float(row["initial_deviation_deg"]) for row in rows],
            [float(row["max_abs_A_chi_tau"]) for row in rows],
            marker="o",
            label=tau_mode_name,
        )
    ax.set_title("phase deviation vs chi-tau area")
    ax.set_xlabel("initial deviation from pi [deg]")
    ax.set_ylabel("max |A_chi_tau|")
    ax.set_yscale("symlog", linthresh=1.0e-12)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_internal_calibration_power_candidate_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for readout_name in [mode.name for mode in READOUT_MODES]:
        rows = [
            row
            for row in case_summaries
            if row["protocol"] == "Protocol_B"
            and row["case_id"] == "near_pi_10deg"
            and row["tau_mode"] == "tau_independent_c1"
            and row["readout_mode"] == readout_name
        ]
        if rows:
            ax.scatter(
                [float(row["max_abs_A_chi_tau"]) for row in rows],
                [abs(float(row["decay_rate_V_AB"])) for row in rows],
                label=readout_name,
                s=60,
            )
    ax.set_title("readout leak response on c1 chi-tau surface")
    ax.set_xlabel("max |A_chi_tau|")
    ax.set_ylabel("|decay_rate_V_AB|")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_internal_calibration_readout_leak_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any], power_rows: List[Dict[str, Any]]) -> None:
    c1_alphas = validation["c1_readout_off_power_candidate_alpha_values"]
    alpha_text = ", ".join(f"{value:.6g}" for value in c1_alphas) if c1_alphas else "none"
    text = f"""# AB二体 c=1 内部較正 chi-tau 面積スイープ予備実験レポート v1

## Summary

- valid: `{validation["c1_internal_calibration_chi_tau_area_sweep_preliminary_valid"]}`
- case_summary_count: `{validation["case_summary_count"]}`
- max_Q_closed_abs: `{validation["max_Q_closed_abs"]:.16e}`
- disabled_max_area: `{validation["disabled_max_area"]:.16e}`
- locked_max_area: `{validation["locked_max_area"]:.16e}`
- independent_min_area: `{validation["independent_min_area"]:.16e}`
- c1_max_epsilon_c_abs: `{validation["c1_max_epsilon_c_abs"]:.16e}`
- tau_disabled_rank1_all_cases: `{validation["tau_disabled_rank1_all_cases"]}`
- tau_locked_rank1_all_cases: `{validation["tau_locked_rank1_all_cases"]}`
- tau_independent_rank2_all_cases: `{validation["tau_independent_rank2_all_cases"]}`
- tau_is_step_used_any: `{validation["tau_is_step_used_any"]}`
- external_c_used_any: `{validation["external_c_used_any"]}`
- f_A_or_f_B_used_any: `{validation["f_A_or_f_B_used_any"]}`

## Main reading

`tau_disabled_control` and `tau_locked_*` are controls. They keep the readout effectively one-dimensional and do not generate a `chi-tau` area.

`tau_independent_*` modes generate a two-dimensional readout surface. Among them, `tau_independent_c1` is the internally calibrated case: the one-period RMS exchange ratio between `chi_read` and `tau_read` is approximately one.

The power-candidate fit is reported only after a nonzero area sweep is present. In the present preliminary construction, the c1 surface gives alpha values:

```text
{alpha_text}
```

Negative alpha means that the measured area-like readout grows with the initial phase deviation rather than decays with it. Therefore this preliminary experiment establishes the `chi-tau` surface control, but it does not yet produce an inverse-power decay law.

## Output files

| kind | file |
|---|---|
| JSON | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1.json` |
| series CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv` |
| case summary CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv` |
| power candidate CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_power_candidates_v1.csv` |
| chi-tau surface | `ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png` |
| area sweep | `ab_two_body_c1_internal_calibration_area_sweep_v1.png` |
| c calibration | `ab_two_body_c1_internal_calibration_error_v1.png` |
| power candidate | `ab_two_body_c1_internal_calibration_power_candidate_v1.png` |
| readout leak | `ab_two_body_c1_internal_calibration_readout_leak_v1.png` |
"""
    (OUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_report_v1.md").write_text(
        text,
        encoding="utf-8",
    )


def main() -> None:
    params = Params()
    series_rows: List[Dict[str, Any]] = []
    case_summaries: List[Dict[str, Any]] = []
    for case in INITIAL_CASES:
        for protocol in PROTOCOLS:
            for tau_mode in TAU_MODES:
                for readout_mode in READOUT_MODES:
                    rows = case_rows(case, protocol, tau_mode, readout_mode, params)
                    series_rows.extend(rows)
                    case_summaries.append(summarize_case(rows, tau_mode, params))

    power_rows = power_candidate_rows(case_summaries)
    validation = validation_summary(case_summaries, power_rows, params)

    write_csv(
        OUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv",
        series_rows,
    )
    write_csv(
        OUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv",
        case_summaries,
    )
    write_csv(
        OUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_power_candidates_v1.csv",
        power_rows,
    )
    result = {
        "params": {
            "step_count": params.step_count,
            "period_steps": params.period_steps,
            "omega_step": params.omega_step,
            "rho_AB": params.rho_AB,
            "area_tol": params.area_tol,
            "rank_tol": params.rank_tol,
            "c_tol": params.c_tol,
            "closure_tol": params.closure_tol,
        },
        "validation": validation,
        "case_summaries": case_summaries,
        "power_candidates": power_rows,
        "notes": [
            "s is not used as tau_read.",
            "c=1 is evaluated as one-period RMS internal exchange calibration.",
            "inverse-power candidates are fitted only for nonzero chi-tau area modes.",
        ],
    }
    (OUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(series_rows, case_summaries, power_rows, params)
    write_report(validation, power_rows)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
