#!/usr/bin/env python3
"""Compute all locked f(t) regression, interval, and growth-end candidates."""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PACKAGE_ROOT / "config_candidates.json"
VERIFY_PATH = PACKAGE_ROOT / "logs" / "input_verification.json"
PROCESSED_DIR = PACKAGE_ROOT / "processed"
METRICS_PATH = PROCESSED_DIR / "f_growth_metrics.csv"
INTERVALS_CSV = PROCESSED_DIR / "growth_intervals_all_candidates.csv"
INTERVALS_MD = PROCESSED_DIR / "growth_intervals_all_candidates.md"
ENDS_CSV = PROCESSED_DIR / "growth_end_all_candidates.csv"
ENDS_MD = PROCESSED_DIR / "growth_end_all_candidates.md"
MANIFEST_PATH = PACKAGE_ROOT / "logs" / "growth_analysis_manifest.json"
TEXT_LOG_PATH = PACKAGE_ROOT / "logs" / "analyze_growth_candidates.log"


def resolve(relative: str) -> Path:
    return (PACKAGE_ROOT / relative).resolve()


def format_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "nan"
        if math.isinf(float(value)):
            return "inf" if value > 0 else "-inf"
        return f"{float(value):.17g}"
    return str(value)


def write_csv_and_markdown(
    csv_path: Path,
    markdown_path: Path,
    fieldnames: list[str],
    rows: list[dict],
    title: str,
) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in fieldnames})
    with markdown_path.open("w", encoding="utf-8") as handle:
        handle.write(f"# {title}\n\n")
        handle.write("| " + " | ".join(fieldnames) + " |\n")
        handle.write("|" + "|".join("---" for _ in fieldnames) + "|\n")
        for row in rows:
            values = [format_value(row.get(key)).replace("|", "\\|") for key in fieldnames]
            handle.write("| " + " | ".join(values) + " |\n")


def regression_metrics(
    steps: np.ndarray,
    log_f: np.ndarray,
    window: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    count = len(steps)
    slope = np.full(count, np.nan)
    r2 = np.full(count, np.nan)
    residual_std = np.full(count, np.nan)
    half = window // 2
    for center in range(half, count - half):
        x = steps[center - half : center + half + 1].astype(float)
        y = log_f[center - half : center + half + 1]
        if not np.all(np.isfinite(y)) or not np.all(np.diff(x) == 1):
            continue
        x_centered = x - np.mean(x)
        y_mean = float(np.mean(y))
        sxx = float(x_centered @ x_centered)
        if sxx == 0:
            continue
        current_slope = float((x_centered @ (y - y_mean)) / sxx)
        fitted = y_mean + current_slope * x_centered
        residual = y - fitted
        sse = float(residual @ residual)
        centered_y = y - y_mean
        sst = float(centered_y @ centered_y)
        slope[center] = current_slope
        r2[center] = 1.0 - sse / sst if sst > 0 else np.nan
        residual_std[center] = math.sqrt(sse / (window - 2))
    return slope, r2, residual_std


def maximal_true_runs(mask: np.ndarray, steps: np.ndarray) -> list[tuple[int, int]]:
    runs = []
    start = None
    for index, valid in enumerate(mask):
        continues = (
            start is not None
            and index > start
            and steps[index] - steps[index - 1] == 1
        )
        if valid and start is None:
            start = index
        elif valid and not continues:
            runs.append((start, index - 1))
            start = index
        elif not valid and start is not None:
            runs.append((start, index - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


def first_persistent_start(
    condition: np.ndarray,
    first_allowed_index: int,
    persistence: int,
) -> int | None:
    if persistence <= 0 or len(condition) < persistence:
        return None
    cumulative = np.concatenate(([0], np.cumsum(condition.astype(np.int64))))
    rolling_counts = cumulative[persistence:] - cumulative[:-persistence]
    starts = np.flatnonzero(rolling_counts == persistence)
    starts = starts[starts >= first_allowed_index]
    return int(starts[0]) if len(starts) else None


def main() -> int:
    started = time.perf_counter()
    planned = [
        METRICS_PATH,
        INTERVALS_CSV,
        INTERVALS_MD,
        ENDS_CSV,
        ENDS_MD,
        MANIFEST_PATH,
        TEXT_LOG_PATH,
    ]
    if any(path.exists() for path in planned):
        raise RuntimeError("成長候補解析出力の上書きを避けて停止")
    if not VERIFY_PATH.is_file():
        raise RuntimeError("verify_inputs.py の成功記録がない")
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    if verification.get("success") is not True:
        raise RuntimeError("入力検証が成功していない")

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    input_path = resolve(config["inputs"]["fcurve"]["path"])
    with input_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    steps = np.array([int(row["tau"]) for row in rows], dtype=np.int64)
    f_values = np.array([float(row["f"]) for row in rows], dtype=float)
    if not np.array_equal(steps, np.arange(steps[0], steps[-1] + 1)):
        raise RuntimeError("fcurveの1 step時間軸が連続していない")

    positive = f_values > 0
    log_f = np.full(len(f_values), np.nan)
    log10_f = np.full(len(f_values), np.nan)
    log_f[positive] = np.log(f_values[positive])
    log10_f[positive] = np.log10(f_values[positive])
    diff_1 = np.full(len(f_values), np.nan)
    central_diff = np.full(len(f_values), np.nan)
    adjacent = (
        np.isfinite(log_f[1:])
        & np.isfinite(log_f[:-1])
        & (np.diff(steps) == 1)
    )
    diff_1[1:][adjacent] = log_f[1:][adjacent] - log_f[:-1][adjacent]
    if len(f_values) >= 3:
        central_valid = (
            np.isfinite(log_f[2:])
            & np.isfinite(log_f[:-2])
            & (steps[2:] - steps[:-2] == 2)
        )
        central_diff[1:-1][central_valid] = (
            log_f[2:][central_valid] - log_f[:-2][central_valid]
        ) / (steps[2:][central_valid] - steps[:-2][central_valid])

    windows = [int(value) for value in config["regression_windows"]]
    metrics = {}
    for window in windows:
        if window % 2 != 1:
            raise RuntimeError(f"回帰窓が奇数でない: {window}")
        metrics[window] = regression_metrics(steps, log_f, window)

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    metric_fields = [
        "step",
        "f",
        "log_f",
        "log10_f",
        "log_f_diff_1",
        "log_f_central_diff",
    ]
    for window in windows:
        tag = f"w{window:03d}"
        metric_fields.extend(
            [f"slope_{tag}", f"r2_{tag}", f"residual_std_{tag}"]
        )
    with METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=metric_fields)
        writer.writeheader()
        for index, step in enumerate(steps):
            row = {
                "step": int(step),
                "f": f_values[index],
                "log_f": log_f[index],
                "log10_f": log10_f[index],
                "log_f_diff_1": diff_1[index],
                "log_f_central_diff": central_diff[index],
            }
            for window in windows:
                tag = f"w{window:03d}"
                slope, r2, residual_std = metrics[window]
                row[f"slope_{tag}"] = slope[index]
                row[f"r2_{tag}"] = r2[index]
                row[f"residual_std_{tag}"] = residual_std[index]
            writer.writerow({key: format_value(row[key]) for key in metric_fields})

    interval_fields = [
        "candidate_id",
        "window",
        "r2_threshold",
        "minimum_duration",
        "interval_start",
        "interval_end",
        "duration",
        "median_slope",
        "mean_slope",
        "minimum_r2",
        "maximum_r2",
        "log_f_start",
        "log_f_end",
    ]
    interval_rows = []
    interval_counter = 0
    counts_by_window = {str(window): 0 for window in windows}
    for window in windows:
        slope, r2, _ = metrics[window]
        for threshold in [float(value) for value in config["r2_thresholds"]]:
            mask = np.isfinite(slope) & np.isfinite(r2) & (slope > 0) & (r2 >= threshold)
            runs = maximal_true_runs(mask, steps)
            for minimum_duration in [
                int(value) for value in config["growth_minimum_durations"]
            ]:
                for start_index, end_index in runs:
                    duration = end_index - start_index + 1
                    if duration < minimum_duration:
                        continue
                    interval_counter += 1
                    interval_rows.append(
                        {
                            "candidate_id": f"GI{interval_counter:06d}",
                            "window": window,
                            "r2_threshold": threshold,
                            "minimum_duration": minimum_duration,
                            "interval_start": int(steps[start_index]),
                            "interval_end": int(steps[end_index]),
                            "duration": duration,
                            "median_slope": float(np.median(slope[start_index : end_index + 1])),
                            "mean_slope": float(np.mean(slope[start_index : end_index + 1])),
                            "minimum_r2": float(np.min(r2[start_index : end_index + 1])),
                            "maximum_r2": float(np.max(r2[start_index : end_index + 1])),
                            "log_f_start": float(log_f[start_index]),
                            "log_f_end": float(log_f[end_index]),
                            "_end_index": end_index,
                        }
                    )
                    counts_by_window[str(window)] += 1

    write_csv_and_markdown(
        INTERVALS_CSV,
        INTERVALS_MD,
        interval_fields,
        interval_rows,
        "growth_intervals_all_candidates",
    )

    end_fields = [
        "growth_end_candidate_id",
        "growth_interval_candidate_id",
        "window",
        "r2_threshold",
        "minimum_duration",
        "interval_start",
        "interval_end",
        "interval_duration",
        "interval_median_slope",
        "end_condition",
        "end_condition_definition",
        "slope_threshold",
        "end_persistence",
        "status",
        "growth_end_candidate",
        "confirmation_end",
        "observed_duration",
    ]
    end_rows = []
    end_counter = 0
    found_count = 0
    for interval in interval_rows:
        slope = metrics[int(interval["window"])][0]
        first_allowed = int(interval["_end_index"]) + 1
        median_slope = float(interval["median_slope"])
        for condition_name, condition_config in config["growth_end_conditions"].items():
            multiplier = float(condition_config["median_slope_multiplier"])
            threshold = 0.0 if condition_name == "A" else multiplier * median_slope
            condition = np.isfinite(slope) & (slope <= threshold)
            for persistence in [
                int(value) for value in config["growth_end_persistence"]
            ]:
                end_counter += 1
                candidate_index = first_persistent_start(
                    condition,
                    first_allowed,
                    persistence,
                )
                found = candidate_index is not None
                if found:
                    found_count += 1
                end_rows.append(
                    {
                        "growth_end_candidate_id": f"GE{end_counter:07d}",
                        "growth_interval_candidate_id": interval["candidate_id"],
                        "window": interval["window"],
                        "r2_threshold": interval["r2_threshold"],
                        "minimum_duration": interval["minimum_duration"],
                        "interval_start": interval["interval_start"],
                        "interval_end": interval["interval_end"],
                        "interval_duration": interval["duration"],
                        "interval_median_slope": median_slope,
                        "end_condition": condition_name,
                        "end_condition_definition": condition_config["description"],
                        "slope_threshold": threshold,
                        "end_persistence": persistence,
                        "status": "found" if found else "not_found",
                        "growth_end_candidate": (
                            int(steps[candidate_index]) if found else None
                        ),
                        "confirmation_end": (
                            int(steps[candidate_index + persistence - 1])
                            if found
                            else None
                        ),
                        "observed_duration": persistence if found else None,
                    }
                )

    write_csv_and_markdown(
        ENDS_CSV,
        ENDS_MD,
        end_fields,
        end_rows,
        "growth_end_all_candidates",
    )

    manifest = {
        "stage": "A1",
        "success": True,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "input": str(input_path),
        "input_rows": len(rows),
        "nonpositive_f_count": int(np.sum(~positive)),
        "regression_windows": windows,
        "r2_thresholds": config["r2_thresholds"],
        "growth_minimum_durations": config["growth_minimum_durations"],
        "growth_end_conditions": config["growth_end_conditions"],
        "growth_end_persistence": config["growth_end_persistence"],
        "growth_interval_candidate_count": len(interval_rows),
        "growth_interval_counts_by_window": counts_by_window,
        "growth_end_row_count": len(end_rows),
        "growth_end_found_count": found_count,
        "growth_end_not_found_count": len(end_rows) - found_count,
        "outputs": [
            str(METRICS_PATH),
            str(INTERVALS_CSV),
            str(INTERVALS_MD),
            str(ENDS_CSV),
            str(ENDS_MD),
        ],
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    log_lines = [
        f"growth_metrics_rows={len(rows)}",
        f"growth_interval_candidates={len(interval_rows)}",
        f"growth_end_rows={len(end_rows)}",
        f"growth_end_found={found_count}",
        f"duration_seconds={manifest['duration_seconds']:.6f}",
        "SUCCESS",
    ]
    TEXT_LOG_PATH.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print("\n".join(log_lines))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
