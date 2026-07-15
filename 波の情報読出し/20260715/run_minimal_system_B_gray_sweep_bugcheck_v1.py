#!/usr/bin/env python3
"""Minimal System B gray-state R sweep.

This script sweeps R for the two-complex-amplitude System B debug bed. It keeps
the original gray-state scoring but removes case labels, harmonic labels,
plotting, samples, and report generation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, List

from run_minimal_system_B_gray_direct_check_v1 import Params, aggregate_conditions, alpha_inv_to_r


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "minimal_system_B_gray_bugcheck_result_v1"


def default_r_values(alpha_inv: float) -> List[float]:
    r137 = alpha_inv_to_r(alpha_inv)
    r128 = alpha_inv_to_r(128.0)
    values = [0.0, 0.5, 1.0, r137, r128]
    values.extend(0.600 + 0.010 * i for i in range(int(round((0.900 - 0.600) / 0.010)) + 1))
    values.extend(0.680 + 0.001 * i for i in range(int(round((0.710 - 0.680) / 0.001)) + 1))
    return sorted({round(float(v), 12) for v in values})


def uniform_r_values(start: float, stop: float, step: float) -> List[float]:
    if step <= 0.0:
        raise ValueError("R step must be positive")
    if stop < start:
        raise ValueError("R stop must be greater than or equal to R start")
    count = int(math.floor((stop - start) / step + 1.0e-12))
    values = [round(start + i * step, 15) for i in range(count + 1)]
    if not values or abs(values[-1] - stop) > 0.5 * step:
        values.append(round(stop, 15))
    return sorted(set(values))


def sweep_region(r_value: float) -> str:
    if 0.680 <= r_value <= 0.710:
        return "fine"
    if 0.600 <= r_value <= 0.900:
        return "coarse"
    return "control"


def local_maxima(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    ordered = sorted(rows, key=lambda row: float(row["R"]))
    peaks: List[Dict[str, Any]] = []
    for i, row in enumerate(ordered):
        score = float(row["joint_gray_score"])
        left = float(ordered[i - 1]["joint_gray_score"]) if i > 0 else -float("inf")
        right = float(ordered[i + 1]["joint_gray_score"]) if i < len(ordered) - 1 else -float("inf")
        if score >= left and score >= right:
            peaks.append(row)
    return sorted(peaks, key=lambda row: -float(row["joint_gray_score"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha-inv", type=float, default=137.035999177)
    parser.add_argument("--r-min", type=float)
    parser.add_argument("--r-max", type=float)
    parser.add_argument("--r-step", type=float)
    parser.add_argument("--r-values", help="comma-separated explicit R values")
    parser.add_argument("--steps", type=int, default=4096)
    parser.add_argument("--output-dir", default=str(OUT_DIR))
    args = parser.parse_args()

    if args.r_values:
        values = sorted({float(part.strip()) for part in args.r_values.split(",") if part.strip()})
    elif args.r_min is not None or args.r_max is not None or args.r_step is not None:
        if args.r_min is None or args.r_max is None or args.r_step is None:
            raise ValueError("--r-min, --r-max, and --r-step must be supplied together")
        values = uniform_r_values(args.r_min, args.r_max, args.r_step)
    else:
        values = default_r_values(args.alpha_inv)

    params = Params(steps=args.steps)
    rows = []
    for r_value in values:
        row = aggregate_conditions(params, r_value)
        row["sweep_region"] = sweep_region(float(r_value))
        row["distance_to_R_137"] = abs(float(r_value) - alpha_inv_to_r(args.alpha_inv))
        row["distance_to_R_128"] = abs(float(r_value) - alpha_inv_to_r(128.0))
        rows.append(row)

    sweep_rows = [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}] or rows
    best = max(sweep_rows, key=lambda row: float(row["joint_gray_score"]))
    peaks = local_maxima(sweep_rows)[:12]

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "minimal_system_B_gray_sweep_summary_v1.csv"
    peaks_path = out_dir / "minimal_system_B_gray_sweep_peaks_v1.csv"
    best_path = out_dir / "minimal_system_B_gray_sweep_best_v1.json"

    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with peaks_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(peaks[0].keys()) if peaks else list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(peaks)

    payload = {
        "model": "minimal_system_B_gray_sweep",
        "steps": params.steps,
        "alpha_inv": args.alpha_inv,
        "R_from_alpha_inv": alpha_inv_to_r(args.alpha_inv),
        "row_count": len(rows),
        "best": best,
        "top_peaks": peaks,
        "summary_csv": str(summary_path),
        "peaks_csv": str(peaks_path),
    }
    best_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
