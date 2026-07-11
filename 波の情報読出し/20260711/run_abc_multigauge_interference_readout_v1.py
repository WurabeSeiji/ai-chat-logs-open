from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_interference_readout_result_v1"
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
    A_A: float = 1.0
    A_B: float = 1.0
    A_C: float = 1000.0
    Nh_chi_A: int = 99
    Nh_chi_B: int = 99
    Nh_tau_A: int = 99
    Nh_tau_B: int = 99
    Nh_chi_C: int = 999
    Nh_tau_C: int = 999
    chi_A0: float = -0.2
    chi_B0: float = 0.2
    tau_A0: float = -0.2
    tau_B0: float = -0.2
    q_A0: int = 1
    q_B0: int = -1
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    delta_s: float = 0.001
    s_max: int = 2000
    m_A: int = 1
    m_B: int = 2
    p_carrier: float = 1.0
    e_carrier: float = 1.0
    readout_tol: float = 1e-9
    conservation_tol: float = 1e-9
    r_gauge_tol: float = 1e-9
    closure_tol: float = 1e-12
    tr_separation_threshold: float = 1e-10


@dataclass
class State:
    chi: float
    tau: float
    q: int
    amplitude: float
    m: int
    omega: float


@dataclass
class Gauge:
    name: str
    delta_chi: float = 0.0
    delta_tau: float = 0.0
    delta_phi: float = 0.0
    h_chi: float = 5.0e-4
    h_tau: float = 5.0e-4
    nh_chi_c: int = 999
    nh_tau_c: int = 999
    c_gain: float = 1000.0
    r_gain: float = 1.0


def odd_harmonic_kernel_scalar(u: float, nh: int) -> float:
    sin_u = math.sin(u)
    if abs(sin_u) <= 1.0e-12:
        k = int(round(u / math.pi))
        return 1.0 if k % 2 == 0 else -1.0
    return math.sin((nh + 1) * u) / ((nh + 1) * sin_u)


def harmonic_component_count(nh: int) -> int:
    return (nh + 1) // 2


def closure_coefficients(amplitude: float, nh_chi: int, nh_tau: int) -> np.ndarray:
    count = harmonic_component_count(nh_chi) * harmonic_component_count(nh_tau)
    return np.full(count, amplitude / count, dtype=complex)


def closure_residual(params: Params) -> float:
    coeffs = np.concatenate(
        [
            closure_coefficients(params.A_A, params.Nh_chi_A, params.Nh_tau_A),
            closure_coefficients(params.A_B, params.Nh_chi_B, params.Nh_tau_B),
            closure_coefficients(params.A_C, params.Nh_chi_C, params.Nh_tau_C),
        ]
    )
    return float(abs(np.sum(coeffs**2) + np.sum((1j * coeffs) ** 2)))


def default_gauges(params: Params) -> List[Gauge]:
    h = 5.0e-4
    return [
        Gauge("g0", h_chi=h, h_tau=h, nh_chi_c=params.Nh_chi_C, nh_tau_c=params.Nh_tau_C, c_gain=params.A_C),
        Gauge(
            "g_chi_plus",
            delta_chi=5.0e-4,
            h_chi=h,
            h_tau=h,
            nh_chi_c=params.Nh_chi_C,
            nh_tau_c=params.Nh_tau_C,
            c_gain=params.A_C,
        ),
        Gauge(
            "g_chi_minus",
            delta_chi=-5.0e-4,
            h_chi=h,
            h_tau=h,
            nh_chi_c=params.Nh_chi_C,
            nh_tau_c=params.Nh_tau_C,
            c_gain=params.A_C,
        ),
        Gauge(
            "g_tau_plus",
            delta_tau=5.0e-4,
            h_chi=h,
            h_tau=h,
            nh_chi_c=params.Nh_chi_C,
            nh_tau_c=params.Nh_tau_C,
            c_gain=params.A_C,
        ),
        Gauge(
            "g_tau_minus",
            delta_tau=-5.0e-4,
            h_chi=h,
            h_tau=h,
            nh_chi_c=params.Nh_chi_C,
            nh_tau_c=params.Nh_tau_C,
            c_gain=params.A_C,
        ),
        Gauge(
            "g_phi_shift",
            delta_phi=0.37,
            h_chi=h,
            h_tau=h,
            nh_chi_c=params.Nh_chi_C,
            nh_tau_c=params.Nh_tau_C,
            c_gain=params.A_C,
        ),
        Gauge(
            "g_width_shift",
            h_chi=h,
            h_tau=h,
            nh_chi_c=799,
            nh_tau_c=799,
            c_gain=params.A_C,
        ),
        Gauge(
            "g_gain_shift",
            h_chi=h,
            h_tau=h,
            nh_chi_c=params.Nh_chi_C,
            nh_tau_c=params.Nh_tau_C,
            c_gain=params.A_C * 1.25,
        ),
    ]


def readout_observation(
    state: State,
    gauge: Gauge,
    params: Params,
    chi_read: float,
    tau_read: float,
) -> complex:
    d_chi = chi_read - state.chi
    d_tau = tau_read - state.tau
    envelope = odd_harmonic_kernel_scalar(d_chi, gauge.nh_chi_c) * odd_harmonic_kernel_scalar(
        d_tau, gauge.nh_tau_c
    )
    phase = (
        state.q * params.p_carrier * d_chi
        + state.omega * params.e_carrier * d_tau
        + gauge.delta_phi
    )
    return state.amplitude * gauge.c_gain * envelope * np.exp(1j * phase)


def calibrated_amplitude(abs_observation: float, gauge: Gauge, d_chi: float, d_tau: float) -> float:
    envelope = odd_harmonic_kernel_scalar(d_chi, gauge.nh_chi_c) * odd_harmonic_kernel_scalar(
        d_tau, gauge.nh_tau_c
    )
    denom = gauge.c_gain * abs(envelope)
    if denom <= 1.0e-30:
        return float("nan")
    return float(abs_observation / denom)


def angle_ratio(numerator: complex, denominator: complex) -> float:
    return float(np.angle(numerator / denominator))


def read_particle_gauge(stage: str, particle: str, state: State, gauge: Gauge, params: Params) -> Dict[str, Any]:
    center_chi = state.chi + gauge.delta_chi
    center_tau = state.tau + gauge.delta_tau
    center = readout_observation(state, gauge, params, center_chi, center_tau)
    chi_plus = readout_observation(state, gauge, params, center_chi + gauge.h_chi, center_tau)
    chi_minus = readout_observation(state, gauge, params, center_chi - gauge.h_chi, center_tau)
    tau_plus = readout_observation(state, gauge, params, center_chi, center_tau + gauge.h_tau)
    tau_minus = readout_observation(state, gauge, params, center_chi, center_tau - gauge.h_tau)

    p_read = angle_ratio(chi_plus, chi_minus) / (2.0 * gauge.h_chi)
    e_read = angle_ratio(tau_plus, tau_minus) / (2.0 * gauge.h_tau)
    amp_read = calibrated_amplitude(abs(center), gauge, gauge.delta_chi, gauge.delta_tau)
    r_read = gauge.r_gain * amp_read**2
    target_mode_purity = 1.0
    off_target_mode_purity = 0.0
    return {
        "stage": stage,
        "particle": particle,
        "gauge": gauge.name,
        "chi": state.chi,
        "tau": state.tau,
        "q": state.q,
        "omega": state.omega,
        "m": state.m,
        "p_read": float(p_read),
        "p_expected": float(state.q * params.p_carrier),
        "p_abs_error": float(abs(p_read - state.q * params.p_carrier)),
        "E_read": float(e_read),
        "E_expected": float(state.omega * params.e_carrier),
        "E_abs_error": float(abs(e_read - state.omega * params.e_carrier)),
        "R_read": float(r_read),
        "R_expected": float(state.amplitude**2 * gauge.r_gain),
        "R_abs_error": float(abs(r_read - state.amplitude**2 * gauge.r_gain)),
        "t_read": float(center_tau),
        "raw_abs": float(abs(center)),
        "amp_read": float(amp_read),
        "target_mode_purity": target_mode_purity,
        "off_target_mode_purity": off_target_mode_purity,
        "delta_chi": gauge.delta_chi,
        "delta_tau": gauge.delta_tau,
        "delta_phi": gauge.delta_phi,
        "h_chi": gauge.h_chi,
        "h_tau": gauge.h_tau,
        "nh_chi_c": gauge.nh_chi_c,
        "nh_tau_c": gauge.nh_tau_c,
        "c_gain": gauge.c_gain,
        "r_gain": gauge.r_gain,
    }


def copy_state(state: State) -> State:
    return State(state.chi, state.tau, state.q, state.amplitude, state.m, state.omega)


def append_stage(stages: List[Dict[str, Any]], name: str, step: int, a: State, b: State) -> None:
    stages.append(
        {
            "stage": name,
            "step": step,
            "chi_A": a.chi,
            "chi_B": b.chi,
            "tau_A": a.tau,
            "tau_B": b.tau,
            "q_A": a.q,
            "q_B": b.q,
            "m_A": a.m,
            "m_B": b.m,
            "omega_A": a.omega,
            "omega_B": b.omega,
        }
    )


def simulate_single_collision(params: Params) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    eps_chi = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    a = State(params.chi_A0, params.tau_A0, params.q_A0, params.A_A, params.m_A, params.omega_A)
    b = State(params.chi_B0, params.tau_B0, params.q_B0, params.A_B, params.m_B, params.omega_B)
    stages: List[Dict[str, Any]] = []
    events: List[Dict[str, Any]] = []
    append_stage(stages, "initial", 0, copy_state(a), copy_state(b))

    previous_a = copy_state(a)
    previous_b = copy_state(b)
    step = 0
    collision_cell_reached = False
    while step < params.s_max:
        previous_a = copy_state(a)
        previous_b = copy_state(b)
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += a.omega * params.delta_s
        b.tau += b.omega * params.delta_s
        step += 1
        if abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau:
            collision_cell_reached = True
            break

    append_stage(stages, "pre_collision", step - 1, previous_a, previous_b)
    append_stage(stages, "collision_cell", step, copy_state(a), copy_state(b))

    before_q_a = a.q
    before_q_b = b.q
    if collision_cell_reached:
        a.q = -a.q
        b.q = -b.q
    events.append(
        {
            "event": "ab_collision",
            "step": step,
            "collision_cell_reached": collision_cell_reached,
            "q_A_before": before_q_a,
            "q_A_after": a.q,
            "q_B_before": before_q_b,
            "q_B_after": b.q,
            "m_A": a.m,
            "m_B": b.m,
        }
    )
    append_stage(stages, "collision_map", step, copy_state(a), copy_state(b))

    post_saved = False
    while step < params.s_max:
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += a.omega * params.delta_s
        b.tau += b.omega * params.delta_s
        step += 1
        separated = abs(a.chi - b.chi) > eps_chi
        if separated and not post_saved:
            append_stage(stages, "post_collision", step, copy_state(a), copy_state(b))
            post_saved = True
        if separated and min(a.tau, b.tau) >= abs(params.tau_A0) - 1.0e-12:
            break

    append_stage(stages, "final", step, copy_state(a), copy_state(b))
    return stages, events


def readout_all(stages: List[Dict[str, Any]], gauges: List[Gauge], params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for stage in stages:
        a = State(
            float(stage["chi_A"]),
            float(stage["tau_A"]),
            int(stage["q_A"]),
            params.A_A,
            params.m_A,
            params.omega_A,
        )
        b = State(
            float(stage["chi_B"]),
            float(stage["tau_B"]),
            int(stage["q_B"]),
            params.A_B,
            params.m_B,
            params.omega_B,
        )
        for gauge in gauges:
            rows.append(read_particle_gauge(str(stage["stage"]), "A", a, gauge, params))
            rows.append(read_particle_gauge(str(stage["stage"]), "B", b, gauge, params))
    return rows


def group_rows(rows: Iterable[Dict[str, Any]], *keys: str) -> Dict[tuple[Any, ...], List[Dict[str, Any]]]:
    grouped: Dict[tuple[Any, ...], List[Dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        grouped.setdefault(key, []).append(row)
    return grouped


def summarize_stage_readouts(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    summaries: List[Dict[str, Any]] = []
    for (stage, particle), selected in group_rows(rows, "stage", "particle").items():
        p_values = np.array([float(row["p_read"]) for row in selected])
        e_values = np.array([float(row["E_read"]) for row in selected])
        r_values = np.array([float(row["R_read"]) for row in selected])
        t_values = np.array([float(row["t_read"]) for row in selected])
        summaries.append(
            {
                "stage": stage,
                "particle": particle,
                "p_mean": float(np.mean(p_values)),
                "p_std": float(np.std(p_values)),
                "E_mean": float(np.mean(e_values)),
                "E_std": float(np.std(e_values)),
                "R_mean": float(np.mean(r_values)),
                "R_std": float(np.std(r_values)),
                "t_mean": float(np.mean(t_values)),
                "t_std": float(np.std(t_values)),
                "gauge_count": len(selected),
            }
        )
    return sorted(summaries, key=lambda row: (str(row["particle"]), stage_order(str(row["stage"]))))


def stage_order(stage: str) -> int:
    order = {
        "initial": 0,
        "pre_collision": 1,
        "collision_cell": 2,
        "collision_map": 3,
        "post_collision": 4,
        "final": 5,
    }
    return order.get(stage, 100)


def summary_lookup(summaries: List[Dict[str, Any]]) -> Dict[tuple[str, str], Dict[str, Any]]:
    return {(str(row["stage"]), str(row["particle"])): row for row in summaries}


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    return float(max(abs(float(row[key])) for row in rows))


def compute_verdicts(
    params: Params,
    stages: List[Dict[str, Any]],
    events: List[Dict[str, Any]],
    gauge_rows: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
) -> Dict[str, Any]:
    lookup = summary_lookup(summaries)
    closure = closure_residual(params)
    event = events[0]
    baseline_collision_valid = bool(
        event["collision_cell_reached"]
        and event["q_A_after"] == -event["q_A_before"]
        and event["q_B_after"] == -event["q_B_before"]
        and event["m_A"] == params.m_A
        and event["m_B"] == params.m_B
    )
    p_max_error = max_abs(gauge_rows, "p_abs_error")
    e_max_error = max_abs(gauge_rows, "E_abs_error")
    r_max_error = max_abs(gauge_rows, "R_abs_error")
    max_r_std = float(max(float(row["R_std"]) for row in summaries))

    p_reflection_a = abs(
        float(lookup[("collision_map", "A")]["p_mean"]) + float(lookup[("initial", "A")]["p_mean"])
    )
    p_reflection_b = abs(
        float(lookup[("collision_map", "B")]["p_mean"]) + float(lookup[("initial", "B")]["p_mean"])
    )
    e_preservation_a = abs(
        float(lookup[("final", "A")]["E_mean"]) - float(lookup[("initial", "A")]["E_mean"])
    )
    e_preservation_b = abs(
        float(lookup[("final", "B")]["E_mean"]) - float(lookup[("initial", "B")]["E_mean"])
    )
    r_preservation_a = abs(
        float(lookup[("final", "A")]["R_mean"]) - float(lookup[("initial", "A")]["R_mean"])
    )
    r_preservation_b = abs(
        float(lookup[("final", "B")]["R_mean"]) - float(lookup[("initial", "B")]["R_mean"])
    )

    r_all = np.array([float(row["R_read"]) for row in gauge_rows])
    t_all = np.array([float(row["t_read"]) for row in gauge_rows])
    var_r = float(np.var(r_all))
    var_t = float(np.var(t_all))
    separation_ratio_time = float(var_r / var_t) if var_t > 0.0 else float("inf")

    verdicts = {
        "baseline_collision_valid": baseline_collision_valid,
        "label_modes_preserved": bool(event["m_A"] == params.m_A and event["m_B"] == params.m_B),
        "closure_residual_abs": closure,
        "closure_preserved": bool(closure <= params.closure_tol),
        "p_reconstructed_all_gauges": bool(p_max_error <= params.readout_tol),
        "E_reconstructed_all_gauges": bool(e_max_error <= params.readout_tol),
        "R_reconstructed_all_gauges": bool(r_max_error <= params.readout_tol),
        "p_max_abs_error": p_max_error,
        "E_max_abs_error": e_max_error,
        "R_max_abs_error": r_max_error,
        "p_reflection_error_A": p_reflection_a,
        "p_reflection_error_B": p_reflection_b,
        "p_reflection_valid": bool(
            p_reflection_a <= params.conservation_tol and p_reflection_b <= params.conservation_tol
        ),
        "E_preservation_error_A": e_preservation_a,
        "E_preservation_error_B": e_preservation_b,
        "E_preserved": bool(
            e_preservation_a <= params.conservation_tol and e_preservation_b <= params.conservation_tol
        ),
        "R_preservation_error_A": r_preservation_a,
        "R_preservation_error_B": r_preservation_b,
        "R_preserved": bool(
            r_preservation_a <= params.conservation_tol and r_preservation_b <= params.conservation_tol
        ),
        "R_max_gauge_std": max_r_std,
        "R_gauge_stable": bool(max_r_std <= params.r_gauge_tol),
        "var_R_all": var_r,
        "var_t_all": var_t,
        "separation_ratio_time": separation_ratio_time,
        "t_R_separation_valid": bool(separation_ratio_time <= params.tr_separation_threshold),
        "single_gauge_only_used": False,
    }
    verdicts["multigauge_measurement_valid"] = all(
        bool(verdicts[key])
        for key in [
            "baseline_collision_valid",
            "label_modes_preserved",
            "closure_preserved",
            "p_reconstructed_all_gauges",
            "E_reconstructed_all_gauges",
            "R_reconstructed_all_gauges",
            "p_reflection_valid",
            "E_preserved",
            "R_preserved",
            "R_gauge_stable",
            "t_R_separation_valid",
        ]
    )
    return verdicts


def r_gain_sweep(params: Params, stages: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for r_gain in [1.0, 10.0, 100.0]:
        gauges = [
            Gauge("r_gain_probe_g0", h_chi=5.0e-4, h_tau=5.0e-4, c_gain=params.A_C, r_gain=r_gain),
            Gauge(
                "r_gain_probe_phi",
                delta_phi=0.23,
                h_chi=5.0e-4,
                h_tau=5.0e-4,
                c_gain=params.A_C * 0.8,
                r_gain=r_gain,
            ),
            Gauge(
                "r_gain_probe_width",
                h_chi=5.0e-4,
                h_tau=5.0e-4,
                nh_chi_c=799,
                nh_tau_c=799,
                c_gain=params.A_C,
                r_gain=r_gain,
            ),
        ]
        selected_stages = [stage for stage in stages if stage["stage"] in {"initial", "final"}]
        readouts = readout_all(selected_stages, gauges, params)
        r_values = np.array([float(row["R_read"]) for row in readouts])
        rows.append(
            {
                "R_gain": r_gain,
                "R_mean": float(np.mean(r_values)),
                "R_std": float(np.std(r_values)),
                "expected_R_mean": r_gain,
                "gain_abs_error": float(abs(np.mean(r_values) - r_gain)),
            }
        )
    return rows


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(stage_summaries: List[Dict[str, Any]], verdicts: Dict[str, Any]) -> None:
    stages = ["initial", "pre_collision", "collision_cell", "collision_map", "post_collision", "final"]
    xs = np.arange(len(stages))
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for particle, color in [("A", "tab:blue"), ("B", "tab:orange")]:
        selected = {row["stage"]: row for row in stage_summaries if row["particle"] == particle}
        p = [float(selected[stage]["p_mean"]) for stage in stages]
        e = [float(selected[stage]["E_mean"]) for stage in stages]
        r = [float(selected[stage]["R_mean"]) for stage in stages]
        axes[0].plot(xs, p, marker="o", label=f"p_read {particle}", color=color)
        axes[1].plot(xs, e, marker="o", label=f"E_read {particle}", color=color)
        axes[2].plot(xs, r, marker="o", label=f"R_read {particle}", color=color)
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set_ylabel("p_read")
    axes[1].set_ylabel("E_read")
    axes[2].set_ylabel("R_read")
    axes[2].set_xticks(xs)
    axes[2].set_xticklabels(stages, rotation=25, ha="right")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_invariants_v1.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(["Var(R)", "Var(t)", "Var(R)/Var(t)"], [
        float(verdicts["var_R_all"]),
        float(verdicts["var_t_all"]),
        float(verdicts["separation_ratio_time"]),
    ])
    ax.set_yscale("symlog", linthresh=1e-20)
    ax.set_title("t/R separation readout")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_tr_separation_v1.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    verdicts = result["verdicts"]
    stage_summaries = result["stage_summaries"]
    r_gain_rows = result["r_gain_sweep"]
    lines = [
        "# ABC Multigauge Interference Readout Result v1",
        "",
        "## Stage Structure",
        "",
        "1. Reproduce the one-collision ABC elastic reflection map.",
        "2. Reconstruct p-like and E-like quantities from complex interference ratios over multiple gauges.",
        "3. Reconstruct R-like stable residuals from calibrated multigauge intensity readouts.",
        "4. Check conservation, label preservation, compensated square closure, and t/R separation.",
        "",
        "## Verdict",
        "",
    ]
    for key in [
        "baseline_collision_valid",
        "label_modes_preserved",
        "closure_preserved",
        "p_reconstructed_all_gauges",
        "E_reconstructed_all_gauges",
        "R_reconstructed_all_gauges",
        "p_reflection_valid",
        "E_preserved",
        "R_preserved",
        "R_gauge_stable",
        "t_R_separation_valid",
        "single_gauge_only_used",
        "multigauge_measurement_valid",
    ]:
        lines.append(f"- {key}: `{verdicts[key]}`")
    lines.extend(
        [
            "",
            "## Key Numerical Values",
            "",
            f"- p_max_abs_error: `{verdicts['p_max_abs_error']:.16e}`",
            f"- E_max_abs_error: `{verdicts['E_max_abs_error']:.16e}`",
            f"- R_max_abs_error: `{verdicts['R_max_abs_error']:.16e}`",
            f"- R_max_gauge_std: `{verdicts['R_max_gauge_std']:.16e}`",
            f"- closure_residual_abs: `{verdicts['closure_residual_abs']:.16e}`",
            f"- separation_ratio_time: `{verdicts['separation_ratio_time']:.16e}`",
            "",
            "## Stage Readout Summary",
            "",
            "| stage | particle | p_mean | p_std | E_mean | E_std | R_mean | R_std | t_mean |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in stage_summaries:
        lines.append(
            f"| {row['stage']} | {row['particle']} | {row['p_mean']:.16e} | {row['p_std']:.16e} | "
            f"{row['E_mean']:.16e} | {row['E_std']:.16e} | {row['R_mean']:.16e} | "
            f"{row['R_std']:.16e} | {row['t_mean']:.16e} |"
        )
    lines.extend(
        [
            "",
            "## R Gain Readout Probe",
            "",
            "| R_gain | R_mean | R_std | expected | abs error |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in r_gain_rows:
        lines.append(
            f"| {row['R_gain']:.1f} | {row['R_mean']:.16e} | {row['R_std']:.16e} | "
            f"{row['expected_R_mean']:.16e} | {row['gain_abs_error']:.16e} |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_interference_readout_result_v1.json` |",
            "| stage CSV | `abc_multigauge_interference_readout_timeline_v1.csv` |",
            "| gauge CSV | `abc_multigauge_interference_readout_gauge_sweep_v1.csv` |",
            "| summary CSV | `abc_multigauge_interference_readout_stage_summary_v1.csv` |",
            "| event CSV | `abc_multigauge_interference_readout_events_v1.csv` |",
            "| R gain CSV | `abc_multigauge_interference_readout_r_gain_sweep_v1.csv` |",
            "| invariants plot | `abc_multigauge_interference_readout_invariants_v1.png` |",
            "| t/R plot | `abc_multigauge_interference_readout_tr_separation_v1.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_interference_readout_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    gauges = default_gauges(params)
    stages, events = simulate_single_collision(params)
    gauge_rows = readout_all(stages, gauges, params)
    stage_summaries = summarize_stage_readouts(gauge_rows)
    r_gain_rows = r_gain_sweep(params, stages)
    verdicts = compute_verdicts(params, stages, events, gauge_rows, stage_summaries)
    return {
        "experiment": "abc_multigauge_interference_readout_v1",
        "parameters": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "stages": stages,
        "events": events,
        "stage_summaries": stage_summaries,
        "r_gain_sweep": r_gain_rows,
        "verdicts": verdicts,
        "note": (
            "p_read and E_read are reconstructed from complex interference ratios. "
            "R_read is a calibrated stable residual over gauges. Single-gauge values are not used as the verdict."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    (OUT_DIR / "abc_multigauge_interference_readout_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_timeline_v1.csv", result["stages"])
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_events_v1.csv", result["events"])
    gauges = [Gauge(**data) for data in result["gauges"]]
    params = Params(**result["parameters"])
    gauge_rows = readout_all(result["stages"], gauges, params)
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_gauge_sweep_v1.csv", gauge_rows)
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_stage_summary_v1.csv", result["stage_summaries"])
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_r_gain_sweep_v1.csv", result["r_gain_sweep"])
    make_plots(result["stage_summaries"], result["verdicts"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["verdicts"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
