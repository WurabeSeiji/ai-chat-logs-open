#!/usr/bin/env python3
"""Minimal direct System B gray-state check, two-phase v4.

This is a two-complex-amplitude bug-check bed. It intentionally removes labels,
harmonic cases, plots, samples, packet readout, and the 36-condition exploratory
grid. The only direct conditions retained are

    s0 = 0.01
    stability_gain = 0
    phi = 0 or pi

The only dynamics is

    (a, b) -> normalize_pair(r a + t b, t a + r b)

and the only classification is the gray-state score used by System B.

V4 keeps R as an explicit direct input. The input text is preserved in CSV, JSON,
and graph titles so that control runs can be checked against the exact R label
used in the previous experiments.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/codex-cache")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_R_TEXT = "0.697177879231003"


@dataclass
class Params:
    steps: int = 4096
    tail_fraction: float = 0.5
    s_gray_limit: float = 0.05
    eigen_amp_tol: float = 1.0e-4
    mean_tol: float = 5.0e-3
    drift_tol: float = 1.0e-2
    selection_limit: float = 0.95
    q_tol: float = 1.0e-10
    metastable_amp_target: float = 0.02
    phi_values: Tuple[float, ...] = (0.0, math.pi)
    s0_values: Tuple[float, ...] = (1.0e-2,)
    stability_gain_values: Tuple[float, ...] = (0.0,)


def normalize_pair(a: complex, b: complex) -> Tuple[complex, complex]:
    q = abs(a) ** 2 + abs(b) ** 2
    if q <= 0.0:
        raise ValueError("zero AB norm")
    scale = 1.0 / math.sqrt(q)
    return a * scale, b * scale


def state_from_s_phi(s: float, phi: float) -> Tuple[complex, complex]:
    if abs(s) >= 0.5:
        raise ValueError(f"|s| must be < 0.5: {s}")
    a = math.sqrt(0.5 + s)
    b = math.sqrt(0.5 - s) * complex(math.cos(phi), math.sin(phi))
    return normalize_pair(a, b)


def scattering_coefficients(reflection_rate: float) -> Tuple[complex, complex, float, float]:
    if reflection_rate < 0.0 or reflection_rate > 1.0:
        raise ValueError(f"reflection rate must be in [0, 1]: {reflection_rate}")
    delta_f = 2.0 * math.asin(math.sqrt(reflection_rate))
    t = np.exp(0.5j * delta_f) * math.cos(0.5 * delta_f)
    r = -1j * np.exp(0.5j * delta_f) * math.sin(0.5 * delta_f)
    return complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


def apply_scattering(a: complex, b: complex, reflection_rate: float) -> Tuple[complex, complex]:
    t, r, _t_power, _r_power = scattering_coefficients(reflection_rate)
    return normalize_pair(r * a + t * b, t * a + r * b)


def apply_stability_gain(a: complex, b: complex, gain: float) -> Tuple[complex, complex]:
    if gain == 0.0:
        return a, b
    p_a = abs(a) ** 2
    p_b = abs(b) ** 2
    q = p_a + p_b
    if q <= 0.0:
        raise ValueError("zero AB norm")
    s = (p_a - p_b) / q
    s_next = s + gain * s * (1.0 - s * s)
    s_next = max(min(s_next, 1.0 - 1.0e-12), -1.0 + 1.0e-12)
    phase_a = a / abs(a) if abs(a) > 0.0 else 1.0 + 0.0j
    phase_b = b / abs(b) if abs(b) > 0.0 else 1.0 + 0.0j
    a_next = math.sqrt(0.5 * (1.0 + s_next)) * phase_a
    b_next = math.sqrt(0.5 * (1.0 - s_next)) * phase_b
    return normalize_pair(a_next, b_next)


def classify_series(params: Params, s_values: np.ndarray, q_values: np.ndarray) -> Dict[str, Any]:
    tail_start = int(round(params.steps * (1.0 - params.tail_fraction)))
    tail = s_values[tail_start:]
    q_tail = q_values[tail_start:]
    first_half_tail = tail[: len(tail) // 2]
    second_half_tail = tail[len(tail) // 2 :]
    s_mean = float(np.mean(tail))
    s_amp = float(0.5 * (np.max(tail) - np.min(tail)))
    s_max_abs = float(np.max(np.abs(tail)))
    s_final = float(s_values[-1])
    s_drift = float(abs(np.mean(second_half_tail) - np.mean(first_half_tail))) if len(first_half_tail) else 0.0
    q_max_error = float(np.max(np.abs(q_values - 1.0)))
    q_tail_max_error = float(np.max(np.abs(q_tail - 1.0)))

    if q_max_error > params.q_tol:
        phase = "norm_error"
    elif abs(s_mean) >= params.selection_limit and s_amp < params.s_gray_limit and abs(s_final) >= params.selection_limit:
        phase = "natural_selection"
    elif abs(s_mean) <= params.mean_tol and s_amp <= params.eigen_amp_tol and s_drift <= params.drift_tol:
        phase = "gray_eigen"
    elif abs(s_mean) <= params.s_gray_limit and params.eigen_amp_tol < s_amp < params.s_gray_limit and s_drift <= params.drift_tol:
        phase = "gray_metastable"
    elif s_amp >= params.s_gray_limit and s_drift <= params.drift_tol:
        phase = "large_oscillation"
    else:
        phase = "unstable_or_drifting"

    if phase == "gray_metastable":
        phase_penalty = 0.0
    elif phase == "gray_eigen":
        phase_penalty = 0.25
    elif phase == "large_oscillation":
        phase_penalty = 0.5
    else:
        phase_penalty = 1.0
    gray_error = abs(s_mean) + abs(s_amp - params.metastable_amp_target) + s_drift + phase_penalty
    gray_depth = -math.log10(max(gray_error, 1.0e-300))
    return {
        "phase": phase,
        "S_mean": s_mean,
        "S_amp": s_amp,
        "S_max_abs": s_max_abs,
        "S_final": s_final,
        "S_drift": s_drift,
        "Q_max_error": q_max_error,
        "Q_tail_max_error": q_tail_max_error,
        "gray_error": gray_error,
        "gray_depth": gray_depth,
    }


def gray_error_no_phase_from_values(params: Params, s_window: np.ndarray) -> Tuple[float, float, float, float]:
    if len(s_window) == 0:
        raise ValueError("empty S window")
    first_half = s_window[: len(s_window) // 2]
    second_half = s_window[len(s_window) // 2 :]
    s_mean = float(np.mean(s_window))
    s_amp = float(0.5 * (np.max(s_window) - np.min(s_window)))
    s_drift = float(abs(np.mean(second_half) - np.mean(first_half))) if len(first_half) else 0.0
    gray_error_no_phase = abs(s_mean) + abs(s_amp - params.metastable_amp_target) + s_drift
    return s_mean, s_amp, s_drift, gray_error_no_phase


def wrapped_phase(value: float) -> float:
    return (value + math.pi) % (2.0 * math.pi) - math.pi


def make_window_metric_row(
    params: Params,
    reflection_rate: float,
    reflection_rate_text: str,
    condition_id: str,
    phi: float,
    s0: float,
    stability_gain: float,
    window_name: str,
    step_start: int,
    step_end: int,
    s_values: np.ndarray,
    q_values: np.ndarray,
) -> Dict[str, Any]:
    s_window = s_values[step_start : step_end + 1]
    q_window = q_values[step_start : step_end + 1]
    if len(s_window) == 0:
        raise ValueError(f"empty window: {window_name}")
    s_mean, s_amp, s_drift, gray_error_no_phase = gray_error_no_phase_from_values(params, s_window)
    return {
        "R": reflection_rate,
        "R_input_text": reflection_rate_text,
        "condition_id": condition_id,
        "window": window_name,
        "step_start": step_start,
        "step_end": step_end,
        "n": len(s_window),
        "phi": phi,
        "phi_over_pi": phi / math.pi,
        "s0": s0,
        "stability_gain": stability_gain,
        "S_mean": s_mean,
        "S_amp": s_amp,
        "S_min": float(np.min(s_window)),
        "S_max": float(np.max(s_window)),
        "S_drift": s_drift,
        "Q_max_error": float(np.max(np.abs(q_window - 1.0))),
        "gray_error_no_phase": gray_error_no_phase,
        "gray_depth_no_phase": -math.log10(max(gray_error_no_phase, 1.0e-300)),
    }


def make_window_metric_rows(
    params: Params,
    reflection_rate: float,
    reflection_rate_text: str,
    condition_id: str,
    phi: float,
    s0: float,
    stability_gain: float,
    s_values: np.ndarray,
    q_values: np.ndarray,
) -> List[Dict[str, Any]]:
    last_step = len(s_values) - 1
    tail_start = int(round(params.steps * (1.0 - params.tail_fraction)))
    windows: List[Tuple[str, int, int]] = [
        ("all", 0, last_step),
        ("first16", 0, min(15, last_step)),
        ("first64", 0, min(63, last_step)),
        ("first256", 0, min(255, last_step)),
        ("first_half", 0, max(0, tail_start - 1)),
        ("tail", tail_start, last_step),
    ]
    rows = [
        make_window_metric_row(
            params,
            reflection_rate,
            reflection_rate_text,
            condition_id,
            phi,
            s0,
            stability_gain,
            window_name,
            start,
            end,
            s_values,
            q_values,
        )
        for window_name, start, end in windows
        if start <= end
    ]
    block_size = 256
    for start in range(0, last_step + 1, block_size):
        end = min(start + block_size - 1, last_step)
        rows.append(
            make_window_metric_row(
                params,
                reflection_rate,
                reflection_rate_text,
                condition_id,
                phi,
                s0,
                stability_gain,
                f"block{block_size}_{start:04d}_{end:04d}",
                start,
                end,
                s_values,
                q_values,
            )
        )
    return rows


def run_condition(
    params: Params,
    reflection_rate: float,
    reflection_rate_text: str,
    phi: float,
    s0: float,
    stability_gain: float,
) -> Dict[str, Any]:
    a, b = state_from_s_phi(s0, phi)
    condition_id = f"phi{phi / math.pi:.6g}_s{s0:.6g}_g{stability_gain:.6g}"
    s_values = np.empty(params.steps + 1)
    q_values = np.empty(params.steps + 1)
    p_a_values = np.empty(params.steps + 1)
    p_b_values = np.empty(params.steps + 1)
    a_values: List[complex] = []
    b_values: List[complex] = []
    time_series_rows: List[Dict[str, Any]] = []
    for step in range(params.steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        p_a_values[step] = p_a
        p_b_values[step] = p_b
        s_values[step] = (p_a - p_b) / q
        q_values[step] = q
        a_values.append(a)
        b_values.append(b)
        if step >= params.steps:
            break
        a, b = apply_scattering(a, b, reflection_rate)
        a, b = apply_stability_gain(a, b, stability_gain)
    cumulative_abs_dS = 0.0
    cumulative_abs_dp = 0.0
    cumulative_s_zero_crossings = 0
    cumulative_dS_turns = 0
    best_prefix_step = 0
    best_prefix_error = float("inf")
    raw_relative_phases = [
        wrapped_phase(math.atan2(b_values[step].imag, b_values[step].real) - math.atan2(a_values[step].imag, a_values[step].real))
        for step in range(params.steps + 1)
    ]
    unwrapped_relative_phases = np.unwrap(np.array(raw_relative_phases))
    relative_phase_cycles = (unwrapped_relative_phases - unwrapped_relative_phases[0]) / (2.0 * math.pi)
    for step in range(params.steps + 1):
        dS = 0.0 if step == 0 else float(s_values[step] - s_values[step - 1])
        prev_dS = 0.0 if step <= 1 else float(s_values[step - 1] - s_values[step - 2])
        dp_a = 0.0 if step == 0 else float(p_a_values[step] - p_a_values[step - 1])
        dp_b = 0.0 if step == 0 else float(p_b_values[step] - p_b_values[step - 1])
        cumulative_abs_dS += abs(dS)
        cumulative_abs_dp += abs(dp_a) + abs(dp_b)
        s_zero_crossing_event = 0
        if step > 0 and s_values[step] != 0.0 and s_values[step - 1] != 0.0:
            s_zero_crossing_event = int(math.copysign(1.0, s_values[step]) != math.copysign(1.0, s_values[step - 1]))
            cumulative_s_zero_crossings += s_zero_crossing_event
        dS_turn_event = 0
        if step > 1 and dS != 0.0 and prev_dS != 0.0:
            dS_turn_event = int(math.copysign(1.0, dS) != math.copysign(1.0, prev_dS))
            cumulative_dS_turns += dS_turn_event
        prefix_mean, prefix_amp, prefix_drift, prefix_error = gray_error_no_phase_from_values(
            params, s_values[: step + 1]
        )
        if prefix_error < best_prefix_error:
            best_prefix_error = prefix_error
            best_prefix_step = step
        phase_a = math.atan2(a_values[step].imag, a_values[step].real)
        phase_b = math.atan2(b_values[step].imag, b_values[step].real)
        relative_phase = raw_relative_phases[step]
        time_series_rows.append(
            {
                "R": reflection_rate,
                "R_input_text": reflection_rate_text,
                "condition_id": condition_id,
                "step": step,
                "phi": phi,
                "phi_over_pi": phi / math.pi,
                "s0": s0,
                "stability_gain": stability_gain,
                "p_a": p_a_values[step],
                "p_b": p_b_values[step],
                "q": q_values[step],
                "S": s_values[step],
                "dS_from_prev": dS,
                "abs_dS_from_prev": abs(dS),
                "cumulative_abs_dS": cumulative_abs_dS,
                "S_zero_crossing_event": s_zero_crossing_event,
                "cumulative_S_zero_crossings": cumulative_s_zero_crossings,
                "S_zero_crossing_cycle_estimate": 0.5 * cumulative_s_zero_crossings,
                "dS_turn_event": dS_turn_event,
                "cumulative_dS_turns": cumulative_dS_turns,
                "dS_turn_cycle_estimate": 0.5 * cumulative_dS_turns,
                "dp_a_from_prev": dp_a,
                "dp_b_from_prev": dp_b,
                "abs_dp_total_from_prev": abs(dp_a) + abs(dp_b),
                "cumulative_abs_dp_total": cumulative_abs_dp,
                "phase_a": phase_a,
                "phase_b": phase_b,
                "relative_phase": relative_phase,
                "relative_phase_over_pi": relative_phase / math.pi,
                "unwrapped_relative_phase": unwrapped_relative_phases[step],
                "unwrapped_relative_phase_over_pi": unwrapped_relative_phases[step] / math.pi,
                "relative_phase_cycles": relative_phase_cycles[step],
                "prefix_S_mean": prefix_mean,
                "prefix_S_amp": prefix_amp,
                "prefix_S_drift": prefix_drift,
                "prefix_gray_error_no_phase": prefix_error,
                "prefix_gray_depth_no_phase": -math.log10(max(prefix_error, 1.0e-300)),
                "best_prefix_step_so_far": best_prefix_step,
                "best_prefix_gray_error_no_phase_so_far": best_prefix_error,
                "a_real": a_values[step].real,
                "a_imag": a_values[step].imag,
                "b_real": b_values[step].real,
                "b_imag": b_values[step].imag,
            }
        )
    metrics = classify_series(params, s_values, q_values)
    window_metric_rows = make_window_metric_rows(
        params,
        reflection_rate,
        reflection_rate_text,
        condition_id,
        phi,
        s0,
        stability_gain,
        s_values,
        q_values,
    )
    return {
        "R": reflection_rate,
        "R_input_text": reflection_rate_text,
        "phi": phi,
        "phi_over_pi": phi / math.pi,
        "s0": s0,
        "stability_gain": stability_gain,
        "condition_id": condition_id,
        "time_series_rows": time_series_rows,
        "window_metric_rows": window_metric_rows,
        **metrics,
    }


def aggregate_conditions(params: Params, reflection_rate: float, reflection_rate_text: str) -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    time_series_rows: List[Dict[str, Any]] = []
    window_metric_rows: List[Dict[str, Any]] = []
    for phi in params.phi_values:
        for s0 in params.s0_values:
            for stability_gain in params.stability_gain_values:
                row = run_condition(params, reflection_rate, reflection_rate_text, phi, s0, stability_gain)
                time_series_rows.extend(row["time_series_rows"])
                window_metric_rows.extend(row["window_metric_rows"])
                rows.append(
                    {
                        key: value
                        for key, value in row.items()
                        if key not in {"time_series_rows", "window_metric_rows"}
                    }
                )
    phase_counts: Dict[str, int] = {}
    for row in rows:
        phase_counts[str(row["phase"])] = phase_counts.get(str(row["phase"]), 0) + 1
    best = min(rows, key=lambda item: float(item["gray_error"]))
    gray_count = phase_counts.get("gray_eigen", 0) + phase_counts.get("gray_metastable", 0)
    selection_count = phase_counts.get("natural_selection", 0)
    joint_gray_score = float(best["gray_depth"]) + 0.25 * gray_count - 0.25 * selection_count
    _t, _r, t_power, r_power = scattering_coefficients(reflection_rate)
    return {
        "R": reflection_rate,
        "R_input_text": reflection_rate_text,
        "T": t_power,
        "reflection_power": r_power,
        "gray_eigen_count": phase_counts.get("gray_eigen", 0),
        "gray_metastable_count": phase_counts.get("gray_metastable", 0),
        "large_oscillation_count": phase_counts.get("large_oscillation", 0),
        "natural_selection_count": phase_counts.get("natural_selection", 0),
        "unstable_or_drifting_count": phase_counts.get("unstable_or_drifting", 0),
        "norm_error_count": phase_counts.get("norm_error", 0),
        "best_condition_id": str(best["condition_id"]),
        "best_phase": str(best["phase"]),
        "best_gray_depth": float(best["gray_depth"]),
        "best_gray_error": float(best["gray_error"]),
        "best_S_mean": float(best["S_mean"]),
        "best_S_amp": float(best["S_amp"]),
        "best_S_drift": float(best["S_drift"]),
        "joint_gray_score": joint_gray_score,
        "condition_rows": rows,
        "time_series_rows": time_series_rows,
        "window_metric_rows": window_metric_rows,
    }


def alpha_inv_to_r(alpha_inv: float) -> float:
    return 1.0 - math.sqrt(4.0 * math.pi / alpha_inv)


def compact_checked(checked: Dict[str, Any]) -> Dict[str, Any]:
    return {
        key: value
        for key, value in checked.items()
        if key not in {"condition_rows", "time_series_rows", "window_metric_rows"}
    }


def plot_condition_diagnostics(
    params: Params,
    out_dir: Path,
    reflection_rate_text: str,
    condition_id: str,
    rows: List[Dict[str, Any]],
) -> Path:
    steps = np.array([int(row["step"]) for row in rows])
    abs_dS = np.array([float(row["abs_dS_from_prev"]) for row in rows])
    cumulative_abs_dS = np.array([float(row["cumulative_abs_dS"]) for row in rows])
    prefix_error = np.array([float(row["prefix_gray_error_no_phase"]) for row in rows])
    zero_cross_cycles = np.array([float(row["S_zero_crossing_cycle_estimate"]) for row in rows])
    turn_cycles = np.array([float(row["dS_turn_cycle_estimate"]) for row in rows])
    best_index = int(np.argmin(prefix_error))
    best_step = int(steps[best_index])
    tail_start = int(round(params.steps * (1.0 - params.tail_fraction)))

    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    fig.suptitle(f"{condition_id}: direct time-evolution diagnostics\nR={reflection_rate_text}")

    axes[0].plot(steps, abs_dS, color="tab:blue", linewidth=0.9)
    axes[0].set_ylabel("|dS|")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(steps, zero_cross_cycles, color="tab:purple", linewidth=1.0, label="S zero-cross cycles")
    axes[1].plot(steps, turn_cycles, color="tab:brown", linewidth=0.9, alpha=0.75, label="dS turn cycles")
    axes[1].set_ylabel("cycle count")
    axes[1].legend(loc="upper left", fontsize=8)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(steps, cumulative_abs_dS, color="tab:orange", linewidth=1.0)
    axes[2].set_ylabel("cum |dS|")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(steps, prefix_error, color="tab:green", linewidth=0.9)
    axes[3].scatter([best_step], [prefix_error[best_index]], color="red", s=28, zorder=5)
    axes[3].set_yscale("log")
    axes[3].set_ylabel("prefix gray error")
    axes[3].set_xlabel("scattering-map iteration")
    axes[3].grid(True, alpha=0.3, which="both")

    for ax in axes:
        ax.axvline(best_step, color="red", linestyle="-", linewidth=1.4, alpha=0.85)
        ax.axvline(tail_start, color="black", linestyle="--", linewidth=1.0, alpha=0.45)
    axes[3].annotate(
        f"best step={best_step}\nerror={prefix_error[best_index]:.3e}",
        xy=(best_step, prefix_error[best_index]),
        xytext=(10, 18),
        textcoords="offset points",
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "red", "lw": 0.8},
    )

    fig.tight_layout()
    safe_id = condition_id.replace(".", "p").replace("-", "m")
    plot_path = out_dir / f"minimal_system_B_gray_direct_diagnostics_{safe_id}_v4.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    return plot_path


def plot_condition_zoom_diagnostics(
    out_dir: Path,
    reflection_rate_text: str,
    condition_id: str,
    rows: List[Dict[str, Any]],
    half_width: int = 50,
) -> Path:
    steps = np.array([int(row["step"]) for row in rows])
    abs_dS = np.array([float(row["abs_dS_from_prev"]) for row in rows])
    prefix_error = np.array([float(row["prefix_gray_error_no_phase"]) for row in rows])
    best_index = int(np.argmin(prefix_error))
    best_step = int(steps[best_index])
    lo = max(0, best_step - half_width)
    hi = min(int(steps[-1]), best_step + half_width)
    mask = (steps >= lo) & (steps <= hi)

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(steps[mask], abs_dS[mask], color="tab:blue", linewidth=1.2)
    ax.axvline(best_step, color="red", linestyle="-", linewidth=1.4, alpha=0.85)
    ax.scatter([best_step], [abs_dS[best_index]], color="red", s=28, zorder=5)
    ax.set_title(f"{condition_id}: |dS| zoom around best prefix step\nR={reflection_rate_text}")
    ax.set_xlabel("scattering-map iteration")
    ax.set_ylabel("|dS|")
    ax.grid(True, alpha=0.3)
    ax.annotate(
        f"best step={best_step}",
        xy=(best_step, abs_dS[best_index]),
        xytext=(10, 18),
        textcoords="offset points",
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "red", "lw": 0.8},
    )

    fig.tight_layout()
    safe_id = condition_id.replace(".", "p").replace("-", "m")
    plot_path = out_dir / f"minimal_system_B_gray_direct_dS_zoom_{safe_id}_v4.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    return plot_path


def plot_all_diagnostics(
    params: Params,
    out_dir: Path,
    reflection_rate_text: str,
    time_series_rows: List[Dict[str, Any]],
) -> List[Path]:
    plot_paths: List[Path] = []
    condition_ids = sorted({str(row["condition_id"]) for row in time_series_rows})
    for condition_id in condition_ids:
        rows = [row for row in time_series_rows if str(row["condition_id"]) == condition_id]
        plot_paths.append(plot_condition_diagnostics(params, out_dir, reflection_rate_text, condition_id, rows))
        plot_paths.append(plot_condition_zoom_diagnostics(out_dir, reflection_rate_text, condition_id, rows))
    return plot_paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--R",
        type=str,
        default=DEFAULT_R_TEXT,
        help=f"R value to check directly; default is previous control R={DEFAULT_R_TEXT}",
    )
    parser.add_argument("--alpha-inv", type=float, default=137.035999177)
    parser.add_argument("--use-alpha-r", action="store_true", help="override --R with R = 1 - sqrt(4pi / alpha_inv)")
    parser.add_argument("--steps", type=int, default=4096)
    parser.add_argument("--out-dir", type=Path, help="write direct summary JSON and two-condition CSV")
    args = parser.parse_args()
    params = Params(steps=args.steps)
    if args.use_alpha_r:
        r_value = alpha_inv_to_r(args.alpha_inv)
        r_input_text = format(r_value, ".17g")
        r_source = "alpha_inv"
    else:
        r_input_text = str(args.R)
        r_value = float(r_input_text)
        r_source = "direct_R"
    checked = aggregate_conditions(params, r_value, r_input_text)
    payload = {
        "model": "minimal_system_B_gray_direct_v4_direct_R_control",
        "steps": params.steps,
        "alpha_inv": args.alpha_inv,
        "R_from_alpha_inv": alpha_inv_to_r(args.alpha_inv),
        "R_source": r_source,
        "R_input_text": r_input_text,
        "R_input_float_17g": format(r_value, ".17g"),
        "checked": compact_checked(checked),
    }
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        condition_rows = checked["condition_rows"]
        time_series_rows = checked["time_series_rows"]
        window_metric_rows = checked["window_metric_rows"]
        condition_csv = args.out_dir / "minimal_system_B_gray_direct_condition_rows_v4.csv"
        with condition_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(condition_rows[0].keys()))
            writer.writeheader()
            writer.writerows(condition_rows)
        time_series_csv = args.out_dir / "minimal_system_B_gray_direct_time_series_v4.csv"
        with time_series_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(time_series_rows[0].keys()))
            writer.writeheader()
            writer.writerows(time_series_rows)
        window_metrics_csv = args.out_dir / "minimal_system_B_gray_direct_window_metrics_v4.csv"
        with window_metrics_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(window_metric_rows[0].keys()))
            writer.writeheader()
            writer.writerows(window_metric_rows)
        plot_paths = plot_all_diagnostics(params, args.out_dir, r_input_text, time_series_rows)
        summary_payload = {
            **payload,
            "condition_rows_csv": str(condition_csv),
            "time_series_csv": str(time_series_csv),
            "window_metrics_csv": str(window_metrics_csv),
            "diagnostic_plots": [str(path) for path in plot_paths],
        }
        summary_json = args.out_dir / "minimal_system_B_gray_direct_summary_v4.json"
        summary_json.write_text(json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
