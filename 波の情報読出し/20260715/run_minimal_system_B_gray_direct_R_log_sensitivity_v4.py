#!/usr/bin/env python3
"""Log-offset R sensitivity check for minimal System B gray-state direct runs.

This script uses the V4 direct checker as the calculation kernel. It probes both
sides of center R with logarithmic offsets and plots the normalized depth loss.
When baseline R is provided, the fixed step and normalization error are kept
from that baseline for all center-R comparison runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/codex-cache")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_minimal_system_B_gray_direct_check_v4 as direct


def format_r(value: float) -> str:
    return format(value, ".17g")


def format_delta(value: float) -> str:
    return f"{value:.0e}".replace("+", "").replace("-", "m").replace("e", "e")


def fixed_step_error(checked: Dict[str, Any], fixed_step: int) -> Tuple[float, str]:
    candidates = [
        row
        for row in checked["time_series_rows"]
        if int(row["step"]) == fixed_step
    ]
    if not candidates:
        raise ValueError(f"fixed step not found: {fixed_step}")
    best = min(candidates, key=lambda row: float(row["prefix_gray_error_no_phase"]))
    return float(best["prefix_gray_error_no_phase"]), str(best["condition_id"])


def global_prefix_min(checked: Dict[str, Any]) -> Tuple[float, int, str]:
    best = min(
        checked["time_series_rows"],
        key=lambda row: float(row["prefix_gray_error_no_phase"]),
    )
    return (
        float(best["prefix_gray_error_no_phase"]),
        int(best["step"]),
        str(best["condition_id"]),
    )


def write_case_outputs(
    params: direct.Params,
    out_dir: Path,
    case_dir: Path,
    checked: Dict[str, Any],
    r_text: str,
) -> List[str]:
    case_dir.mkdir(parents=True, exist_ok=True)
    condition_rows = checked["condition_rows"]
    time_series_rows = checked["time_series_rows"]
    window_metric_rows = checked["window_metric_rows"]

    condition_csv = case_dir / "minimal_system_B_gray_direct_condition_rows_v4.csv"
    with condition_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(condition_rows[0].keys()))
        writer.writeheader()
        writer.writerows(condition_rows)

    time_series_csv = case_dir / "minimal_system_B_gray_direct_time_series_v4.csv"
    with time_series_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(time_series_rows[0].keys()))
        writer.writeheader()
        writer.writerows(time_series_rows)

    window_csv = case_dir / "minimal_system_B_gray_direct_window_metrics_v4.csv"
    with window_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(window_metric_rows[0].keys()))
        writer.writeheader()
        writer.writerows(window_metric_rows)

    plot_paths = direct.plot_all_diagnostics(params, case_dir, r_text, time_series_rows)
    summary_json = case_dir / "minimal_system_B_gray_direct_summary_v4.json"
    summary_json.write_text(
        json.dumps(
            {
                "model": "minimal_system_B_gray_direct_v4_log_sensitivity_case",
                "R_input_text": r_text,
                "checked": direct.compact_checked(checked),
                "condition_rows_csv": str(condition_csv),
                "time_series_csv": str(time_series_csv),
                "window_metrics_csv": str(window_csv),
                "diagnostic_plots": [str(path) for path in plot_paths],
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return [str(condition_csv), str(time_series_csv), str(window_csv), str(summary_json), *map(str, plot_paths)]


def choose_edge_delta(
    params: direct.Params,
    center_r: float,
    center_error: float,
    fixed_step: int,
    sign: float,
    start_delta: float,
    max_delta: float,
    target_ratio: float,
) -> float:
    candidates: List[Tuple[float, float]] = []
    delta = start_delta
    while delta <= max_delta * (1.0 + 1.0e-12):
        r_value = center_r + sign * delta
        checked = direct.aggregate_conditions(params, r_value, format_r(r_value))
        err, _condition_id = fixed_step_error(checked, fixed_step)
        ratio = center_error / err if err > 0.0 else float("inf")
        candidates.append((delta, ratio))
        delta *= 10.0
    return min(candidates, key=lambda item: abs(item[1] - target_ratio))[0]


def build_offsets(edge_delta: float, points_per_side: int) -> np.ndarray:
    if points_per_side <= 1:
        return np.array([edge_delta], dtype=float)
    start = edge_delta / (10.0 ** (points_per_side - 1))
    return np.geomspace(start, edge_delta, points_per_side)


def plot_depth_summary(
    out_dir: Path,
    rows: List[Dict[str, Any]],
    center_r_text: str,
    baseline_r_text: str,
    fixed_step: int,
    baseline_error: float,
    target_ratio: float,
) -> Path:
    fig, ax = plt.subplots(figsize=(10, 5.8))
    plot_rows = [row for row in rows if row["case"] != "center"]
    x = np.array([float(row["delta_R"]) for row in plot_rows])
    y = np.array([float(row["normalized_inverse_depth"]) for row in plot_rows])
    colors = ["tab:blue" if value < 0 else "tab:orange" for value in x]

    ax.scatter(x, y, c=colors, s=46, zorder=3)
    ax.plot(x[x < 0], y[x < 0], color="tab:blue", linewidth=1.0, alpha=0.7)
    ax.plot(x[x > 0], y[x > 0], color="tab:orange", linewidth=1.0, alpha=0.7)
    ax.scatter([0.0], [1.0], color="black", s=60, zorder=4, label="center")
    ax.axhline(target_ratio, color="red", linestyle="--", linewidth=1.0, alpha=0.75, label=f"target {target_ratio:g}")
    ax.axvline(0.0, color="black", linestyle="-", linewidth=0.8, alpha=0.45)

    min_abs_delta = min(abs(value) for value in x if value != 0.0)
    ax.set_xscale("symlog", linthresh=min_abs_delta * 0.5)
    ax.set_ylim(bottom=0.0, top=max(1.05, float(np.max(y)) * 1.08))
    ax.set_xlabel("signed R offset from center (symlog)")
    ax.set_ylabel("normalized depth = baseline error / shifted error")
    ax.set_title(
        "System B direct R log-sensitivity\n"
        f"center R={center_r_text}, baseline R={baseline_r_text}, "
        f"fixed step={fixed_step}, baseline error={baseline_error:.3e}"
    )
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="best")

    for row in plot_rows:
        delta = float(row["delta_R"])
        ratio = float(row["normalized_inverse_depth"])
        ax.annotate(
            f"{delta:+.0e}",
            xy=(delta, ratio),
            xytext=(4, 5 if delta > 0 else -12),
            textcoords="offset points",
            fontsize=8,
        )

    fig.tight_layout()
    path = out_dir / "minimal_system_B_gray_direct_R_log_sensitivity_depth_v4.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--center-R", type=str, default=direct.DEFAULT_R_TEXT)
    parser.add_argument(
        "--baseline-R",
        type=str,
        help="R that fixes the evaluation step and normalization; defaults to center-R",
    )
    parser.add_argument("--steps", type=int, default=4096)
    parser.add_argument("--points-per-side", type=int, default=5)
    parser.add_argument("--target-depth-ratio", type=float, default=0.05)
    parser.add_argument("--start-delta", type=float, default=1.0e-12)
    parser.add_argument("--max-delta", type=float, default=1.0e-3)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("波の情報読出し/20260715/minimal_system_B_gray_bugcheck_result_v1/direct_R_log_sensitivity_v4"),
    )
    args = parser.parse_args()

    params = direct.Params(steps=args.steps)
    center_r = float(args.center_R)
    center_r_text = args.center_R
    baseline_r_text = args.baseline_R if args.baseline_R is not None else center_r_text
    baseline_r = float(baseline_r_text)
    baseline_checked = direct.aggregate_conditions(params, baseline_r, baseline_r_text)
    baseline_error, fixed_step, baseline_condition = global_prefix_min(baseline_checked)
    center_checked = direct.aggregate_conditions(params, center_r, center_r_text)
    center_error_at_baseline_step, center_condition_at_baseline_step = fixed_step_error(center_checked, fixed_step)
    center_global_error, center_global_step, center_global_condition = global_prefix_min(center_checked)

    left_edge = choose_edge_delta(
        params,
        center_r,
        baseline_error,
        fixed_step,
        -1.0,
        args.start_delta,
        args.max_delta,
        args.target_depth_ratio,
    )
    right_edge = choose_edge_delta(
        params,
        center_r,
        baseline_error,
        fixed_step,
        1.0,
        args.start_delta,
        args.max_delta,
        args.target_depth_ratio,
    )

    left_offsets = -build_offsets(left_edge, args.points_per_side)[::-1]
    right_offsets = build_offsets(right_edge, args.points_per_side)
    offsets = [*left_offsets.tolist(), 0.0, *right_offsets.tolist()]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: List[Dict[str, Any]] = []
    output_files: List[str] = []

    for delta in offsets:
        r_value = center_r + delta
        r_text = center_r_text if delta == 0.0 else format_r(r_value)
        checked = center_checked if delta == 0.0 else direct.aggregate_conditions(params, r_value, r_text)
        fixed_error, fixed_condition = fixed_step_error(checked, fixed_step)
        global_error, global_step, global_condition = global_prefix_min(checked)
        normalized_inverse_depth = baseline_error / fixed_error if fixed_error > 0.0 else float("inf")
        depth_log10 = -math.log10(max(fixed_error, 1.0e-300))
        case_name = "center" if delta == 0.0 else f"delta_{'p' if delta > 0 else 'm'}{format_delta(abs(delta))}"
        case_dir = args.out_dir / case_name
        output_files.extend(write_case_outputs(params, args.out_dir, case_dir, checked, r_text))
        summary_rows.append(
            {
                "case": case_name,
                "R_input_text": r_text,
                "R": r_value,
                "center_R": center_r,
                "baseline_R": baseline_r,
                "delta_R": delta,
                "abs_delta_R": abs(delta),
                "log10_abs_delta_R": "" if delta == 0.0 else math.log10(abs(delta)),
                "fixed_step": fixed_step,
                "baseline_condition": baseline_condition,
                "center_condition_at_baseline_step": center_condition_at_baseline_step,
                "fixed_step_condition": fixed_condition,
                "fixed_step_error": fixed_error,
                "fixed_step_depth_log10": depth_log10,
                "baseline_fixed_step_error": baseline_error,
                "center_fixed_step_error": center_error_at_baseline_step,
                "normalized_inverse_depth": normalized_inverse_depth,
                "center_global_min_prefix_error": center_global_error,
                "center_global_best_step": center_global_step,
                "center_global_best_condition": center_global_condition,
                "global_min_prefix_error": global_error,
                "global_best_step": global_step,
                "global_best_condition": global_condition,
                "tail_best_gray_error": checked["best_gray_error"],
                "tail_best_condition": checked["best_condition_id"],
            }
        )

    summary_csv = args.out_dir / "minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    depth_plot = plot_depth_summary(
        args.out_dir,
        summary_rows,
        center_r_text,
        baseline_r_text,
        fixed_step,
        baseline_error,
        args.target_depth_ratio,
    )
    output_files.append(str(summary_csv))
    output_files.append(str(depth_plot))

    payload = {
        "model": "minimal_system_B_gray_direct_R_log_sensitivity_v4",
        "center_R": center_r_text,
        "baseline_R": baseline_r_text,
        "steps": args.steps,
        "fixed_step": fixed_step,
        "baseline_condition": baseline_condition,
        "baseline_fixed_step_error": baseline_error,
        "center_condition_at_baseline_step": center_condition_at_baseline_step,
        "center_fixed_step_error": center_error_at_baseline_step,
        "center_global_min_prefix_error": center_global_error,
        "center_global_best_step": center_global_step,
        "center_global_best_condition": center_global_condition,
        "reuse_arguments": (
            f"--baseline-R {baseline_r_text} --steps {args.steps}"
        ),
        "left_edge_delta": left_edge,
        "right_edge_delta": right_edge,
        "target_depth_ratio": args.target_depth_ratio,
        "points_per_side": args.points_per_side,
        "summary_csv": str(summary_csv),
        "depth_plot": str(depth_plot),
        "output_files": output_files,
        "rows": summary_rows,
    }
    summary_json = args.out_dir / "minimal_system_B_gray_direct_R_log_sensitivity_summary_v4.json"
    summary_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k not in {"output_files", "rows"}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
