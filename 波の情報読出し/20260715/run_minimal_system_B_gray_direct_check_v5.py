#!/usr/bin/env python3
"""Minimal System B gray-depth probe for one direct R value.

This V5 script is intentionally narrower than the V4 direct checker.
It keeps only the calculation needed for fixed-point scouting:

  R -> scattering coefficients
  (a, b) -> normalize_pair(r a + t b, t a + r b)
  prefix gray error over selected phi condition(s)

It does not generate plots, CSV files, window metrics, or time-series records.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, List, Tuple


DEFAULT_R_TEXT = "0.697177879231003"
DEFAULT_STEPS = 4096
METASTABLE_AMP_TARGET = 0.02


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
    half_delta = 0.5 * delta_f
    phase = complex(math.cos(half_delta), math.sin(half_delta))
    t = phase * math.cos(half_delta)
    r = -1j * phase * math.sin(half_delta)
    return complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


def condition_id(phi: float, s0: float, stability_gain: float) -> str:
    return f"phi{phi / math.pi:.6g}_s{s0:.6g}_g{stability_gain:.6g}"


def prefix_error(prefix_sums: List[float], s_min: float, s_max: float) -> Tuple[float, float, float, float]:
    n = len(prefix_sums) - 1
    if n <= 0:
        raise ValueError("empty prefix")
    total = prefix_sums[n]
    s_mean = total / n
    s_amp = 0.5 * (s_max - s_min)
    mid = n // 2
    if mid == 0:
        s_drift = 0.0
    else:
        first_mean = prefix_sums[mid] / mid
        second_mean = (total - prefix_sums[mid]) / (n - mid)
        s_drift = abs(second_mean - first_mean)
    gray_error = abs(s_mean) + abs(s_amp - METASTABLE_AMP_TARGET) + s_drift
    return gray_error, s_mean, s_amp, s_drift


def run_condition(
    reflection_rate: float,
    phi: float,
    steps: int,
    min_steps: int,
    early_stop_patience: int,
    s0: float = 1.0e-2,
) -> Dict[str, Any]:
    a, b = state_from_s_phi(s0, phi)
    t, r, _t_power, _r_power = scattering_coefficients(reflection_rate)

    prefix_sums = [0.0]
    s_min = float("inf")
    s_max = float("-inf")
    best = {
        "best_prefix_gray_error_no_phase": float("inf"),
        "best_prefix_gray_depth_no_phase": 0.0,
        "best_step": -1,
        "S_mean": 0.0,
        "S_amp": 0.0,
        "S_drift": 0.0,
    }

    stopped_at_step = steps
    stop_reason = "max_steps"
    for step in range(steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s_value = (p_a - p_b) / q
        s_min = min(s_min, s_value)
        s_max = max(s_max, s_value)
        prefix_sums.append(prefix_sums[-1] + s_value)

        gray_error, s_mean, s_amp, s_drift = prefix_error(prefix_sums, s_min, s_max)
        if gray_error < best["best_prefix_gray_error_no_phase"]:
            best = {
                "best_prefix_gray_error_no_phase": gray_error,
                "best_prefix_gray_depth_no_phase": -math.log10(max(gray_error, 1.0e-300)),
                "best_step": step,
                "S_mean": s_mean,
                "S_amp": s_amp,
                "S_drift": s_drift,
            }

        if step >= steps:
            break
        if (
            early_stop_patience > 0
            and step >= min_steps
            and step - int(best["best_step"]) >= early_stop_patience
        ):
            stopped_at_step = step
            stop_reason = "early_stop_no_best_update"
            break
        a, b = normalize_pair(r * a + t * b, t * a + r * b)

    return {
        "condition_id": condition_id(phi, s0, 0.0),
        "phi": phi,
        "phi_over_pi": phi / math.pi,
        "s0": s0,
        "stopped_at_step": stopped_at_step,
        "stop_reason": stop_reason,
        **best,
    }


def phi_values_from_mode(phi_mode: str) -> Tuple[float, ...]:
    if phi_mode == "zero":
        return (0.0,)
    if phi_mode == "pi":
        return (math.pi,)
    if phi_mode == "both":
        return (0.0, math.pi)
    raise ValueError(f"unsupported phi mode: {phi_mode}")


def probe_r(
    reflection_rate: float,
    reflection_rate_text: str,
    steps: int,
    phi_mode: str,
    min_steps: int,
    early_stop_patience: int,
) -> Dict[str, Any]:
    _t, _r, t_power, r_power = scattering_coefficients(reflection_rate)
    rows = [
        run_condition(reflection_rate, phi, steps, min_steps, early_stop_patience)
        for phi in phi_values_from_mode(phi_mode)
    ]
    best = min(rows, key=lambda row: float(row["best_prefix_gray_error_no_phase"]))
    return {
        "model": "minimal_system_B_gray_direct_depth_probe_v5",
        "R": reflection_rate,
        "R_input_text": reflection_rate_text,
        "steps": steps,
        "phi_mode": phi_mode,
        "min_steps": min_steps,
        "early_stop_patience": early_stop_patience,
        "T": t_power,
        "reflection_power": r_power,
        "best_condition_id": best["condition_id"],
        "best_step": best["best_step"],
        "best_prefix_gray_error_no_phase": best["best_prefix_gray_error_no_phase"],
        "best_prefix_gray_depth_no_phase": best["best_prefix_gray_depth_no_phase"],
        "stopped_at_step": best["stopped_at_step"],
        "stop_reason": best["stop_reason"],
        "best_S_mean": best["S_mean"],
        "best_S_amp": best["S_amp"],
        "best_S_drift": best["S_drift"],
        "condition_rows": rows,
    }


def stable_fixed_point_candidate(payload: Dict[str, Any]) -> bool:
    patience = int(payload["early_stop_patience"])
    min_steps = int(payload["min_steps"])
    if patience <= 0 or min_steps <= 0:
        return False
    best_step = int(payload["best_step"])
    return (
        payload["stop_reason"] == "early_stop_no_best_update"
        and min_steps - patience <= best_step <= min_steps
    )


def candidate_row(payload: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "R": payload["R"],
        "R_input_text": payload["R_input_text"],
        "best_step": payload["best_step"],
        "best_prefix_gray_error_no_phase": payload["best_prefix_gray_error_no_phase"],
        "best_prefix_gray_depth_no_phase": payload["best_prefix_gray_depth_no_phase"],
        "best_condition_id": payload["best_condition_id"],
        "best_S_mean": payload["best_S_mean"],
        "best_S_amp": payload["best_S_amp"],
        "best_S_drift": payload["best_S_drift"],
        "stopped_at_step": payload["stopped_at_step"],
        "stop_reason": payload["stop_reason"],
        "phi_mode": payload["phi_mode"],
        "steps": payload["steps"],
        "min_steps": payload["min_steps"],
        "early_stop_patience": payload["early_stop_patience"],
    }


def decimal_range(min_r: str, max_r: str, delta_r: str) -> List[str]:
    start = Decimal(min_r)
    stop = Decimal(max_r)
    step = Decimal(delta_r)
    if step <= 0:
        raise ValueError("--delta-R must be positive")
    if stop < start:
        raise ValueError("--max-R must be greater than or equal to --min-R")

    values: List[str] = []
    current = start
    guard = 0
    while current <= stop:
        values.append(format(current, "f"))
        current += step
        guard += 1
        if guard > 10_000_000:
            raise RuntimeError("too many R sweep points")
    return values


def run_sweep(args: argparse.Namespace) -> Dict[str, Any]:
    if args.min_R is None or args.max_R is None or args.delta_R is None:
        raise ValueError("--min-R, --max-R, and --delta-R must be specified together")
    if args.out_csv is None:
        raise ValueError("--out-csv is required for R sweep output")

    r_values = decimal_range(args.min_R, args.max_R, args.delta_R)
    candidate_rows: List[Dict[str, Any]] = []
    for index, r_text in enumerate(r_values, start=1):
        payload = probe_r(
            float(r_text),
            r_text,
            args.steps,
            args.phi_mode,
            args.min_steps,
            args.early_stop_patience,
        )
        if stable_fixed_point_candidate(payload):
            candidate_rows.append(candidate_row(payload))
            if args.progress_every > 0:
                print(
                    "candidate "
                    f"[{index}/{len(r_values)}] "
                    f"R={payload['R_input_text']} "
                    f"depth={payload['best_prefix_gray_depth_no_phase']:.12f} "
                    f"step={payload['best_step']} "
                    f"candidates={len(candidate_rows)}",
                    file=sys.stderr,
                    flush=True,
                )
        if args.progress_every > 0 and index % args.progress_every == 0:
            print(
                f"[{index}/{len(r_values)}] "
                f"R={r_text} "
                f"candidates={len(candidate_rows)}",
                file=sys.stderr,
                flush=True,
            )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "R",
        "R_input_text",
        "best_step",
        "best_prefix_gray_error_no_phase",
        "best_prefix_gray_depth_no_phase",
        "best_condition_id",
        "best_S_mean",
        "best_S_amp",
        "best_S_drift",
        "stopped_at_step",
        "stop_reason",
        "phi_mode",
        "steps",
        "min_steps",
        "early_stop_patience",
    ]
    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(candidate_rows)

    return {
        "model": "minimal_system_B_gray_direct_depth_probe_v5_R_sweep",
        "min_R": args.min_R,
        "max_R": args.max_R,
        "delta_R": args.delta_R,
        "n_R": len(r_values),
        "n_candidates": len(candidate_rows),
        "out_csv": str(args.out_csv),
        "candidate_rows": candidate_rows,
        "candidate_rule": "stop_reason == early_stop_no_best_update and min_steps - early_stop_patience <= best_step <= min_steps",
        "phi_mode": args.phi_mode,
        "steps": args.steps,
        "min_steps": args.min_steps,
        "early_stop_patience": args.early_stop_patience,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--R", type=str, default=DEFAULT_R_TEXT)
    parser.add_argument("--min-R", dest="min_R", type=str)
    parser.add_argument("--max-R", dest="max_R", type=str)
    parser.add_argument("--delta-R", dest="delta_R", type=str)
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    parser.add_argument("--phi-mode", choices=("both", "zero", "pi"), default="both")
    parser.add_argument("--min-steps", type=int, default=0)
    parser.add_argument("--early-stop-patience", type=int, default=0)
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-csv", type=Path)
    parser.add_argument("--progress-every", type=int, default=0)
    args = parser.parse_args()

    if args.min_steps < 0:
        raise ValueError("--min-steps must be non-negative")
    if args.early_stop_patience < 0:
        raise ValueError("--early-stop-patience must be non-negative")
    if args.progress_every < 0:
        raise ValueError("--progress-every must be non-negative")

    if args.min_R is not None or args.max_R is not None or args.delta_R is not None:
        payload = run_sweep(args)
    else:
        r_text = str(args.R)
        payload = probe_r(
            float(r_text),
            r_text,
            args.steps,
            args.phi_mode,
            args.min_steps,
            args.early_stop_patience,
        )
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
