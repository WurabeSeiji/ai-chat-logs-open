#!/usr/bin/env python3
"""Minimal direct System B gray-state check.

This is a two-complex-amplitude bug-check bed. It intentionally removes labels,
harmonic cases, plots, samples, and packet readout. The only dynamics is

    (a, b) -> normalize_pair(r a + t b, t a + r b)

and the only classification is the gray-state score used by System B.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


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
    phi_values: Tuple[float, ...] = (0.0, math.pi / 12.0, math.pi / 2.0, math.pi)
    s0_values: Tuple[float, ...] = (0.0, 1.0e-3, 1.0e-2)
    stability_gain_values: Tuple[float, ...] = (-2.0e-3, 0.0, 2.0e-3)


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


def run_condition(params: Params, reflection_rate: float, phi: float, s0: float, stability_gain: float) -> Dict[str, Any]:
    a, b = state_from_s_phi(s0, phi)
    s_values = np.empty(params.steps + 1)
    q_values = np.empty(params.steps + 1)
    for step in range(params.steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s_values[step] = (p_a - p_b) / q
        q_values[step] = q
        if step >= params.steps:
            break
        a, b = apply_scattering(a, b, reflection_rate)
        a, b = apply_stability_gain(a, b, stability_gain)
    metrics = classify_series(params, s_values, q_values)
    return {
        "R": reflection_rate,
        "phi": phi,
        "phi_over_pi": phi / math.pi,
        "s0": s0,
        "stability_gain": stability_gain,
        "condition_id": f"phi{phi / math.pi:.6g}_s{s0:.6g}_g{stability_gain:.6g}",
        **metrics,
    }


def aggregate_conditions(params: Params, reflection_rate: float) -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    for phi in params.phi_values:
        for s0 in params.s0_values:
            for stability_gain in params.stability_gain_values:
                rows.append(run_condition(params, reflection_rate, phi, s0, stability_gain))
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
    }


def alpha_inv_to_r(alpha_inv: float) -> float:
    return 1.0 - math.sqrt(4.0 * math.pi / alpha_inv)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--R", type=float, help="R value to check directly")
    parser.add_argument("--alpha-inv", type=float, default=137.035999177)
    parser.add_argument("--use-alpha-r", action="store_true", help="use R = 1 - sqrt(4pi / alpha_inv)")
    parser.add_argument("--steps", type=int, default=4096)
    parser.add_argument("--out-dir", type=Path, help="write direct summary JSON and 36-condition CSV")
    args = parser.parse_args()
    params = Params(steps=args.steps)
    r_value = alpha_inv_to_r(args.alpha_inv) if args.use_alpha_r or args.R is None else float(args.R)
    checked = aggregate_conditions(params, r_value)
    payload = {
        "model": "minimal_system_B_gray_direct",
        "steps": params.steps,
        "alpha_inv": args.alpha_inv,
        "R_from_alpha_inv": alpha_inv_to_r(args.alpha_inv),
        "checked": checked,
    }
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        condition_rows = checked["condition_rows"]
        condition_csv = args.out_dir / "minimal_system_B_gray_direct_condition_rows_v1.csv"
        with condition_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(condition_rows[0].keys()))
            writer.writeheader()
            writer.writerows(condition_rows)
        summary_payload = {
            **payload,
            "checked": {
                key: value for key, value in checked.items() if key != "condition_rows"
            },
            "condition_rows_csv": str(condition_csv),
        }
        summary_json = args.out_dir / "minimal_system_B_gray_direct_summary_v1.json"
        summary_json.write_text(json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
