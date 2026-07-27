#!/usr/bin/env python3
"""Describe the first N=5 transition from the three verified Stage A0 CSVs."""

from __future__ import annotations

import bisect
import csv
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parent
LOG_DIR = PACKAGE_ROOT / "logs"
VERIFY_PATH = LOG_DIR / "input_verification.json"
MANIFEST_PATH = LOG_DIR / "analysis_manifest.json"
TEXT_LOG_PATH = LOG_DIR / "analyze_first_transition.log"
PROCESSED_DIR = PACKAGE_ROOT / "processed"

F_METRICS_PATH = PROCESSED_DIR / "transition_f_metrics_0_3000.csv"
OCC_ACTUAL_PATH = PROCESSED_DIR / "occupation_actual_records_0_3000.csv"
OCC_DISPLAY_PATH = (
    PROCESSED_DIR / "occupation_display_only_linear_interpolation_0_3000.csv"
)
Q_ACTUAL_PATH = PROCESSED_DIR / "q_actual_records_0_3000.csv"
DIRECT_PATH = PROCESSED_DIR / "transition_direct_observations.json"

FIRST_PASSAGE_CSV = PROCESSED_DIR / "f_first_passage_levels.csv"
FIRST_PASSAGE_MD = PROCESSED_DIR / "f_first_passage_levels.md"
RATES_CSV = PROCESSED_DIR / "f_decade_growth_rates.csv"
RATES_MD = PROCESSED_DIR / "f_decade_growth_rates.md"
OCC_NEAREST_CSV = (
    PROCESSED_DIR / "first_passage_nearest_occupation_records.csv"
)
OCC_NEAREST_MD = PROCESSED_DIR / "first_passage_nearest_occupation_records.md"
Q_NEAREST_CSV = PROCESSED_DIR / "first_passage_nearest_q_records.csv"
Q_NEAREST_MD = PROCESSED_DIR / "first_passage_nearest_q_records.md"
STATS_CSV = PROCESSED_DIR / "transition_window_descriptive_statistics.csv"
STATS_MD = PROCESSED_DIR / "transition_window_descriptive_statistics.md"

CROSSING = 1167
OBSERVATION_END = 3000
REGRESSION_WINDOWS = [11, 21, 41, 81, 161]
DISPLAY_WINDOWS = [
    ("basic_0_3000", 0, 3000),
    ("zoom_0_500", 0, 500),
    ("zoom_500_1000", 500, 1000),
    ("zoom_800_1400", 800, 1400),
    ("zoom_1000_1800", 1000, 1800),
    ("zoom_1400_2500", 1400, 2500),
]
OCCUPATION_COLUMNS = [
    "direction_1_occupation",
    "direction_2_occupation",
    "direction_3_occupation",
    "direction_4_occupation",
    "other_rotating_occupation",
    "kernel_occupation",
    "splitting_fraction",
]
Q_COLUMNS = ["q1", "q2", "q3", "q4", "rank_q"]


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


def load_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_rows(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in fieldnames})


def write_csv_and_markdown(
    csv_path: Path,
    md_path: Path,
    fieldnames: list[str],
    rows: list[dict],
    title: str,
) -> None:
    write_rows(csv_path, fieldnames, rows)
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write(f"# {title}\n\n")
        handle.write("| " + " | ".join(fieldnames) + " |\n")
        handle.write("|" + "|".join("---" for _ in fieldnames) + "|\n")
        for row in rows:
            values = [
                format_value(row.get(key)).replace("|", "\\|")
                for key in fieldnames
            ]
            handle.write("| " + " | ".join(values) + " |\n")


def regression_metrics(
    steps: np.ndarray,
    log_f: np.ndarray,
    window: int,
) -> tuple[np.ndarray, np.ndarray]:
    slope = np.full(len(steps), np.nan)
    r2 = np.full(len(steps), np.nan)
    half = window // 2
    for center in range(half, len(steps) - half):
        x = steps[center - half : center + half + 1].astype(float)
        y = log_f[center - half : center + half + 1]
        if not np.all(np.isfinite(y)) or not np.all(np.diff(x) == 1):
            continue
        xc = x - np.mean(x)
        ym = float(np.mean(y))
        sxx = float(xc @ xc)
        if sxx == 0:
            continue
        current_slope = float((xc @ (y - ym)) / sxx)
        fitted = ym + current_slope * xc
        residual = y - fitted
        sse = float(residual @ residual)
        centered = y - ym
        sst = float(centered @ centered)
        slope[center] = current_slope
        r2[center] = 1.0 - sse / sst if sst > 0 else np.nan
    return slope, r2


def bracket_indices(saved_steps: list[int], target: int) -> tuple[int | None, int | None]:
    before = bisect.bisect_right(saved_steps, target) - 1
    after = bisect.bisect_left(saved_steps, target)
    return (
        before if before >= 0 else None,
        after if after < len(saved_steps) else None,
    )


def safe_ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator != 0 else float("nan")


def descriptive_row(
    window_name: str,
    window_start: int,
    window_end: int,
    source: str,
    sampling: str,
    variable: str,
    steps: np.ndarray,
    values: np.ndarray,
) -> dict:
    mask = (
        (steps >= window_start)
        & (steps <= window_end)
        & np.isfinite(values)
    )
    selected_steps = steps[mask]
    selected_values = values[mask]
    if not len(selected_values):
        return {
            "window_name": window_name,
            "window_start": window_start,
            "window_end": window_end,
            "source": source,
            "sampling": sampling,
            "variable": variable,
            "record_count": 0,
        }
    differences = np.diff(selected_values)
    return {
        "window_name": window_name,
        "window_start": window_start,
        "window_end": window_end,
        "source": source,
        "sampling": sampling,
        "variable": variable,
        "record_count": len(selected_values),
        "minimum": float(np.min(selected_values)),
        "maximum": float(np.max(selected_values)),
        "mean": float(np.mean(selected_values)),
        "median": float(np.median(selected_values)),
        "first_record_step": int(selected_steps[0]),
        "first_value": float(selected_values[0]),
        "last_record_step": int(selected_steps[-1]),
        "last_value": float(selected_values[-1]),
        "net_change": float(selected_values[-1] - selected_values[0]),
        "positive_differences": int(np.sum(differences > 0)),
        "negative_differences": int(np.sum(differences < 0)),
        "zero_differences": int(np.sum(differences == 0)),
    }


def main() -> int:
    started = time.perf_counter()
    planned = [
        F_METRICS_PATH,
        OCC_ACTUAL_PATH,
        OCC_DISPLAY_PATH,
        Q_ACTUAL_PATH,
        DIRECT_PATH,
        FIRST_PASSAGE_CSV,
        FIRST_PASSAGE_MD,
        RATES_CSV,
        RATES_MD,
        OCC_NEAREST_CSV,
        OCC_NEAREST_MD,
        Q_NEAREST_CSV,
        Q_NEAREST_MD,
        STATS_CSV,
        STATS_MD,
        MANIFEST_PATH,
        TEXT_LOG_PATH,
    ]
    if any(path.exists() for path in planned):
        raise RuntimeError("Stage A1b解析出力の上書きを避けて停止")
    if not VERIFY_PATH.is_file():
        raise RuntimeError("verify_inputs.py の成功記録がない")
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    if verification.get("success") is not True:
        raise RuntimeError("入力検証が成功していない")

    input_checks = verification["input_checks"]
    _, f_all = load_csv(Path(input_checks["fcurve"]["path"]))
    _, q_all = load_csv(Path(input_checks["q_svd"]["path"]))
    _, occ_all = load_csv(Path(input_checks["paper7_long"]["path"]))
    f_rows = [row for row in f_all if int(row["tau"]) <= OBSERVATION_END]
    q_rows = [row for row in q_all if int(row["step"]) <= OBSERVATION_END]
    occ_rows = [row for row in occ_all if int(row["step"]) <= OBSERVATION_END]

    f_step = np.array([int(row["tau"]) for row in f_rows], dtype=np.int64)
    f_value = np.array([float(row["f"]) for row in f_rows], dtype=float)
    positive = f_value > 0
    log_f = np.full(len(f_value), np.nan)
    log10_f = np.full(len(f_value), np.nan)
    log_f[positive] = np.log(f_value[positive])
    log10_f[positive] = np.log10(f_value[positive])
    running_max_f = np.maximum.accumulate(f_value)
    running_max_log10_f = np.full(len(f_value), np.nan)
    running_positive = running_max_f > 0
    running_max_log10_f[running_positive] = np.log10(
        running_max_f[running_positive]
    )
    f_diff_1 = np.full(len(f_value), np.nan)
    f_diff_1[1:] = f_value[1:] - f_value[:-1]
    regression = {
        window: regression_metrics(f_step, log_f, window)
        for window in REGRESSION_WINDOWS
    }

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    f_metric_fields = [
        "step",
        "f",
        "log10_f",
        "running_max_f",
        "running_max_log10_f",
        "f_diff_1",
    ]
    for window in REGRESSION_WINDOWS:
        f_metric_fields.extend([f"slope_w{window:03d}", f"r2_w{window:03d}"])
    f_metric_rows = []
    for index, step in enumerate(f_step):
        row = {
            "step": int(step),
            "f": f_value[index],
            "log10_f": log10_f[index],
            "running_max_f": running_max_f[index],
            "running_max_log10_f": running_max_log10_f[index],
            "f_diff_1": f_diff_1[index],
        }
        for window in REGRESSION_WINDOWS:
            row[f"slope_w{window:03d}"] = regression[window][0][index]
            row[f"r2_w{window:03d}"] = regression[window][1][index]
        f_metric_rows.append(row)
    write_rows(F_METRICS_PATH, f_metric_fields, f_metric_rows)

    occ_fields = ["step", "time", "crossing_flag"] + OCCUPATION_COLUMNS
    occ_actual_rows = [
        {
            key: (
                int(row[key])
                if key in ("step", "time", "crossing_flag")
                else float(row[key])
            )
            for key in occ_fields
        }
        for row in occ_rows
    ]
    write_rows(OCC_ACTUAL_PATH, occ_fields, occ_actual_rows)

    occ_steps = np.array([int(row["step"]) for row in occ_rows], dtype=np.int64)
    display_grid = np.arange(0, OBSERVATION_END + 1, dtype=np.int64)
    before_indices = np.searchsorted(occ_steps, display_grid, side="right") - 1
    after_indices = np.searchsorted(occ_steps, display_grid, side="left")
    before_indices = np.clip(before_indices, 0, len(occ_steps) - 1)
    after_indices = np.clip(after_indices, 0, len(occ_steps) - 1)
    actual_step_set = set(int(value) for value in occ_steps)
    occ_display_fields = [
        "step",
        "display_only_linear_interpolation",
        "is_actual_saved_record",
        "source_before_step",
        "source_after_step",
    ] + [f"{column}_display_interp" for column in OCCUPATION_COLUMNS]
    occ_display_rows = []
    for grid_index, step in enumerate(display_grid):
        row = {
            "step": int(step),
            "display_only_linear_interpolation": 1,
            "is_actual_saved_record": int(int(step) in actual_step_set),
            "source_before_step": int(occ_steps[before_indices[grid_index]]),
            "source_after_step": int(occ_steps[after_indices[grid_index]]),
        }
        for column in OCCUPATION_COLUMNS:
            values = np.array([float(item[column]) for item in occ_rows])
            row[f"{column}_display_interp"] = float(
                np.interp(step, occ_steps, values)
            )
        occ_display_rows.append(row)
    write_rows(OCC_DISPLAY_PATH, occ_display_fields, occ_display_rows)

    q_steps = np.array([int(row["step"]) for row in q_rows], dtype=np.int64)
    q_actual_fields = [
        "step",
        "time",
        "relative_time",
        "q1",
        "q2",
        "q3",
        "q4",
        "rank_q",
        "q3_over_q1",
        "q4_over_q1",
        "min_q3_q4_over_q1",
    ]
    q_actual_rows = []
    for row in q_rows:
        q1 = float(row["q1"])
        q3 = float(row["q3"])
        q4 = float(row["q4"])
        q_actual_rows.append(
            {
                "step": int(row["step"]),
                "time": int(row["time"]),
                "relative_time": int(row["relative_time"]),
                "q1": q1,
                "q2": float(row["q2"]),
                "q3": q3,
                "q4": q4,
                "rank_q": int(row["rank_q"]),
                "q3_over_q1": safe_ratio(q3, q1),
                "q4_over_q1": safe_ratio(q4, q1),
                "min_q3_q4_over_q1": safe_ratio(min(q3, q4), q1),
            }
        )
    write_rows(Q_ACTUAL_PATH, q_actual_fields, q_actual_rows)

    positive_minimum = float(np.min(f_value[positive]))
    minimum_decade_exponent = int(math.floor(math.log10(positive_minimum)))
    level_map = {}

    def add_level(level: float, source: str, exponent: int | None = None) -> None:
        key = f"{level:.17g}"
        if key not in level_map:
            level_map[key] = {
                "level": level,
                "sources": set(),
                "decade_exponent": exponent,
            }
        level_map[key]["sources"].add(source)
        if exponent is not None:
            level_map[key]["decade_exponent"] = exponent

    for exponent in range(minimum_decade_exponent, 0):
        add_level(10.0**exponent, "decade", exponent)
    for level in [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.1]:
        add_level(level, "explicit", None)

    first_passage_fields = [
        "level_id",
        "level",
        "level_label",
        "log10_level",
        "level_source",
        "decade_exponent",
        "status",
        "first_passage_step",
        "f_at_first_passage",
        "previous_step",
        "f_at_previous_step",
        "relative_to_crossing",
    ]
    first_passage_rows = []
    for number, item in enumerate(
        sorted(level_map.values(), key=lambda value: value["level"]),
        1,
    ):
        level = float(item["level"])
        indices = np.flatnonzero(f_value >= level)
        found = bool(len(indices))
        index = int(indices[0]) if found else None
        exponent = item["decade_exponent"]
        label = f"1e{exponent}" if exponent is not None else format_value(level)
        passage_step = int(f_step[index]) if found else None
        first_passage_rows.append(
            {
                "level_id": f"FPL{number:03d}",
                "level": level,
                "level_label": label,
                "log10_level": math.log10(level),
                "level_source": "+".join(sorted(item["sources"])),
                "decade_exponent": exponent,
                "status": "found" if found else "not_found",
                "first_passage_step": passage_step,
                "f_at_first_passage": float(f_value[index]) if found else None,
                "previous_step": (
                    int(f_step[index - 1]) if found and index > 0 else None
                ),
                "f_at_previous_step": (
                    float(f_value[index - 1]) if found and index > 0 else None
                ),
                "relative_to_crossing": (
                    passage_step - CROSSING if found else None
                ),
            }
        )
    write_csv_and_markdown(
        FIRST_PASSAGE_CSV,
        FIRST_PASSAGE_MD,
        first_passage_fields,
        first_passage_rows,
        "f_first_passage_levels",
    )

    rate_fields = [
        "pair_id",
        "lower_level_id",
        "upper_level_id",
        "lower_level",
        "upper_level",
        "pair_type",
        "lower_first_passage_step",
        "upper_first_passage_step",
        "status",
        "step_difference",
        "threshold_log_amplitude_difference",
        "observed_log_f_difference",
        "mean_exponential_rate_per_step",
    ]
    rate_rows = []
    for number, (lower, upper) in enumerate(
        zip(first_passage_rows, first_passage_rows[1:]),
        1,
    ):
        both_found = lower["status"] == "found" and upper["status"] == "found"
        step_difference = (
            upper["first_passage_step"] - lower["first_passage_step"]
            if both_found
            else None
        )
        threshold_log_difference = math.log(
            float(upper["level"]) / float(lower["level"])
        )
        observed_log_difference = (
            math.log(float(upper["f_at_first_passage"]))
            - math.log(float(lower["f_at_first_passage"]))
            if both_found
            and float(lower["f_at_first_passage"]) > 0
            and float(upper["f_at_first_passage"]) > 0
            else None
        )
        if not both_found:
            status = "not_found"
            rate = None
        elif step_difference == 0:
            status = "zero_step_difference"
            rate = None
        else:
            status = "found"
            rate = threshold_log_difference / step_difference
        lower_exp = lower["decade_exponent"]
        upper_exp = upper["decade_exponent"]
        pair_type = (
            "decade_to_decade"
            if lower_exp is not None
            and upper_exp is not None
            and int(upper_exp) - int(lower_exp) == 1
            else "adjacent_coordinate_levels"
        )
        rate_rows.append(
            {
                "pair_id": f"FGR{number:03d}",
                "lower_level_id": lower["level_id"],
                "upper_level_id": upper["level_id"],
                "lower_level": lower["level"],
                "upper_level": upper["level"],
                "pair_type": pair_type,
                "lower_first_passage_step": lower["first_passage_step"],
                "upper_first_passage_step": upper["first_passage_step"],
                "status": status,
                "step_difference": step_difference,
                "threshold_log_amplitude_difference": threshold_log_difference,
                "observed_log_f_difference": observed_log_difference,
                "mean_exponential_rate_per_step": rate,
            }
        )
    write_csv_and_markdown(
        RATES_CSV,
        RATES_MD,
        rate_fields,
        rate_rows,
        "f_decade_growth_rates",
    )

    occ_nearest_fields = [
        "level_id",
        "level",
        "first_passage_step",
        "status",
        "before_or_at_step",
        "before_step_offset",
    ]
    for column in OCCUPATION_COLUMNS:
        occ_nearest_fields.append(f"before_{column}")
    occ_nearest_fields.extend(["after_or_at_step", "after_step_offset"])
    for column in OCCUPATION_COLUMNS:
        occ_nearest_fields.append(f"after_{column}")
    occ_saved_steps = [int(row["step"]) for row in occ_actual_rows]
    occ_nearest_rows = []
    for passage in first_passage_rows:
        target = passage["first_passage_step"]
        row = {
            "level_id": passage["level_id"],
            "level": passage["level"],
            "first_passage_step": target,
            "status": passage["status"],
        }
        if target is not None:
            before_index, after_index = bracket_indices(occ_saved_steps, int(target))
            if before_index is not None:
                record = occ_actual_rows[before_index]
                row["before_or_at_step"] = record["step"]
                row["before_step_offset"] = record["step"] - int(target)
                for column in OCCUPATION_COLUMNS:
                    row[f"before_{column}"] = record[column]
            if after_index is not None:
                record = occ_actual_rows[after_index]
                row["after_or_at_step"] = record["step"]
                row["after_step_offset"] = record["step"] - int(target)
                for column in OCCUPATION_COLUMNS:
                    row[f"after_{column}"] = record[column]
        occ_nearest_rows.append(row)
    write_csv_and_markdown(
        OCC_NEAREST_CSV,
        OCC_NEAREST_MD,
        occ_nearest_fields,
        occ_nearest_rows,
        "first_passage_nearest_occupation_records",
    )

    q_nearest_fields = [
        "level_id",
        "level",
        "first_passage_step",
        "status",
        "before_or_at_step",
        "before_step_offset",
    ]
    q_value_columns = [
        "q1",
        "q2",
        "q3",
        "q4",
        "rank_q",
        "q3_over_q1",
        "q4_over_q1",
        "min_q3_q4_over_q1",
    ]
    for column in q_value_columns:
        q_nearest_fields.append(f"before_{column}")
    q_nearest_fields.extend(["after_or_at_step", "after_step_offset"])
    for column in q_value_columns:
        q_nearest_fields.append(f"after_{column}")
    q_saved_steps = [int(row["step"]) for row in q_actual_rows]
    q_nearest_rows = []
    for passage in first_passage_rows:
        target = passage["first_passage_step"]
        row = {
            "level_id": passage["level_id"],
            "level": passage["level"],
            "first_passage_step": target,
            "status": passage["status"],
        }
        if target is not None:
            before_index, after_index = bracket_indices(q_saved_steps, int(target))
            if before_index is not None:
                record = q_actual_rows[before_index]
                row["before_or_at_step"] = record["step"]
                row["before_step_offset"] = record["step"] - int(target)
                for column in q_value_columns:
                    row[f"before_{column}"] = record[column]
            if after_index is not None:
                record = q_actual_rows[after_index]
                row["after_or_at_step"] = record["step"]
                row["after_step_offset"] = record["step"] - int(target)
                for column in q_value_columns:
                    row[f"after_{column}"] = record[column]
        q_nearest_rows.append(row)
    write_csv_and_markdown(
        Q_NEAREST_CSV,
        Q_NEAREST_MD,
        q_nearest_fields,
        q_nearest_rows,
        "first_passage_nearest_q_records",
    )

    stats_fields = [
        "window_name",
        "window_start",
        "window_end",
        "source",
        "sampling",
        "variable",
        "record_count",
        "minimum",
        "maximum",
        "mean",
        "median",
        "first_record_step",
        "first_value",
        "last_record_step",
        "last_value",
        "net_change",
        "positive_differences",
        "negative_differences",
        "zero_differences",
    ]
    stats_rows = []
    f_stat_values = {
        "f": f_value,
        "log10_f": log10_f,
        "running_max_f": running_max_f,
        "running_max_log10_f": running_max_log10_f,
        "f_diff_1": f_diff_1,
    }
    for window in REGRESSION_WINDOWS:
        f_stat_values[f"slope_w{window:03d}"] = regression[window][0]
        f_stat_values[f"r2_w{window:03d}"] = regression[window][1]
    occ_stat_values = {
        column: np.array([float(row[column]) for row in occ_actual_rows])
        for column in OCCUPATION_COLUMNS
    }
    q_stat_values = {
        column: np.array([float(row[column]) for row in q_actual_rows])
        for column in q_value_columns
    }
    for name, window_start, window_end in DISPLAY_WINDOWS:
        for variable, values in f_stat_values.items():
            stats_rows.append(
                descriptive_row(
                    name,
                    window_start,
                    window_end,
                    "fcurve",
                    "actual_every_step",
                    variable,
                    f_step,
                    values,
                )
            )
        for variable, values in occ_stat_values.items():
            stats_rows.append(
                descriptive_row(
                    name,
                    window_start,
                    window_end,
                    "paper7_long_timeseries",
                    "actual_saved_records_25_step",
                    variable,
                    occ_steps,
                    values,
                )
            )
        for variable, values in q_stat_values.items():
            stats_rows.append(
                descriptive_row(
                    name,
                    window_start,
                    window_end,
                    "q_svd",
                    "actual_saved_records_5_step_in_observation_range",
                    variable,
                    q_steps,
                    values,
                )
            )
    write_csv_and_markdown(
        STATS_CSV,
        STATS_MD,
        stats_fields,
        stats_rows,
        "transition_window_descriptive_statistics",
    )

    first_rank4 = next(
        (row for row in q_actual_rows if int(row["rank_q"]) == 4),
        None,
    )
    first_q4_positive = next(
        (row for row in q_actual_rows if float(row["q4"]) > 0),
        None,
    )
    crossing_f_index = int(np.where(f_step == CROSSING)[0][0])
    crossing_occ_indices = bracket_indices(occ_saved_steps, CROSSING)
    crossing_q_indices = bracket_indices(q_saved_steps, CROSSING)
    valid_decade_rates = [
        float(row["mean_exponential_rate_per_step"])
        for row in rate_rows
        if row["pair_type"] == "decade_to_decade" and row["status"] == "found"
    ]
    direct = {
        "stage": "A1b",
        "observation_range": [0, OBSERVATION_END],
        "existing_crossing": CROSSING,
        "first_rank_q_4_saved_record": first_rank4,
        "first_q4_positive_saved_record": first_q4_positive,
        "rank_q_4_saved_records_before_crossing": sum(
            int(row["rank_q"]) == 4 and int(row["step"]) < CROSSING
            for row in q_actual_rows
        ),
        "crossing_f_records": {
            "previous": {
                "step": int(f_step[crossing_f_index - 1]),
                "f": float(f_value[crossing_f_index - 1]),
            },
            "crossing": {
                "step": int(f_step[crossing_f_index]),
                "f": float(f_value[crossing_f_index]),
            },
        },
        "crossing_occupation_actual_bracket": {
            "before_or_at": occ_actual_rows[crossing_occ_indices[0]],
            "after_or_at": occ_actual_rows[crossing_occ_indices[1]],
        },
        "crossing_q_actual_bracket": {
            "before_or_at": q_actual_rows[crossing_q_indices[0]],
            "after_or_at": q_actual_rows[crossing_q_indices[1]],
        },
        "positive_minimum_f_in_observation_range": positive_minimum,
        "minimum_decade_exponent": minimum_decade_exponent,
        "first_passage_level_count": len(first_passage_rows),
        "first_passage_found_count": sum(
            row["status"] == "found" for row in first_passage_rows
        ),
        "valid_decade_rate_count": len(valid_decade_rates),
        "decade_rate_distribution": (
            {
                "minimum": float(np.min(valid_decade_rates)),
                "q25": float(np.quantile(valid_decade_rates, 0.25)),
                "median": float(np.median(valid_decade_rates)),
                "q75": float(np.quantile(valid_decade_rates, 0.75)),
                "maximum": float(np.max(valid_decade_rates)),
                "mean": float(np.mean(valid_decade_rates)),
                "std": float(np.std(valid_decade_rates)),
            }
            if valid_decade_rates
            else None
        ),
        "interpolation_policy": {
            "occupation_analysis": "actual_saved_records_only",
            "occupation_display": "separate_linear_interpolation_file",
            "q_analysis": "actual_saved_records_only",
            "q_display_interpolation": False,
        },
        "automatic_single_event_selection": False,
    }
    DIRECT_PATH.write_text(
        json.dumps(direct, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    manifest = {
        "stage": "A1b",
        "success": True,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "observation_range": [0, OBSERVATION_END],
        "display_windows": [
            {"name": name, "start": start, "end": end}
            for name, start, end in DISPLAY_WINDOWS
        ],
        "regression_windows": REGRESSION_WINDOWS,
        "input_rows_used": {
            "fcurve": len(f_rows),
            "q_svd": len(q_rows),
            "paper7_long": len(occ_rows),
        },
        "first_passage_rows": len(first_passage_rows),
        "growth_rate_rows": len(rate_rows),
        "nearest_occupation_rows": len(occ_nearest_rows),
        "nearest_q_rows": len(q_nearest_rows),
        "descriptive_statistic_rows": len(stats_rows),
        "full_time_candidate_cartesian_product_generated": False,
        "automatic_single_event_selection": False,
        "outputs": [str(path) for path in planned if path not in (MANIFEST_PATH, TEXT_LOG_PATH)],
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    log_lines = [
        f"f_rows_0_3000={len(f_rows)}",
        f"occupation_actual_rows_0_3000={len(occ_rows)}",
        f"q_actual_rows_0_3000={len(q_rows)}",
        f"first_passage_levels={len(first_passage_rows)}",
        f"growth_rate_pairs={len(rate_rows)}",
        f"descriptive_statistics_rows={len(stats_rows)}",
        "cartesian_candidate_table_generated=false",
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
