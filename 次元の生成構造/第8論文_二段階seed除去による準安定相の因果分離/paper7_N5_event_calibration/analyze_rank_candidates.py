#!/usr/bin/env python3
"""Compute all locked q-ratio, rank robustness, onset, and pair candidates."""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PACKAGE_ROOT / "config_candidates.json"
VERIFY_PATH = PACKAGE_ROOT / "logs" / "input_verification.json"
GROWTH_MANIFEST_PATH = PACKAGE_ROOT / "logs" / "growth_analysis_manifest.json"
PROCESSED_DIR = PACKAGE_ROOT / "processed"
GROWTH_INTERVALS_PATH = PROCESSED_DIR / "growth_intervals_all_candidates.csv"
GROWTH_ENDS_PATH = PROCESSED_DIR / "growth_end_all_candidates.csv"
RANK_METRICS_PATH = PROCESSED_DIR / "q_rank_candidate_metrics.csv"
RANK_ONSETS_CSV = PROCESSED_DIR / "rank4_onset_all_candidates.csv"
RANK_ONSETS_MD = PROCESSED_DIR / "rank4_onset_all_candidates.md"
PAIRS_CSV = PROCESSED_DIR / "growth_end_vs_rank4_onset_all_pairs.csv"
PAIRS_MD = PROCESSED_DIR / "growth_end_vs_rank4_onset_all_pairs.md"
SUMMARY_CSV = PROCESSED_DIR / "candidate_summary.csv"
SUMMARY_MD = PROCESSED_DIR / "candidate_summary.md"
MANIFEST_PATH = PACKAGE_ROOT / "logs" / "rank_analysis_manifest.json"
TEXT_LOG_PATH = PACKAGE_ROOT / "logs" / "analyze_rank_candidates.log"


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


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return float("nan")
    return numerator / denominator


def threshold_tag(threshold: float) -> str:
    return f"{threshold:.0e}"


def write_csv_and_markdown(
    csv_path: Path,
    markdown_path: Path,
    fieldnames: list[str],
    rows: list[dict],
    title: str,
) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as csv_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in fieldnames})
    with markdown_path.open("w", encoding="utf-8") as md_handle:
        md_handle.write(f"# {title}\n\n")
        md_handle.write("| " + " | ".join(fieldnames) + " |\n")
        md_handle.write("|" + "|".join("---" for _ in fieldnames) + "|\n")
        for row in rows:
            values = [format_value(row.get(key)).replace("|", "\\|") for key in fieldnames]
            md_handle.write("| " + " | ".join(values) + " |\n")


def first_true_record_run(mask: np.ndarray, length: int) -> tuple[int, int] | None:
    if length <= 0 or len(mask) < length:
        return None
    cumulative = np.concatenate(([0], np.cumsum(mask.astype(np.int64))))
    counts = cumulative[length:] - cumulative[:-length]
    starts = np.flatnonzero(counts == length)
    if not len(starts):
        return None
    start = int(starts[0])
    return start, start + length - 1


def numeric_summary(values: list[float]) -> dict:
    if not values:
        return {
            "count": 0,
            "minimum": None,
            "q25": None,
            "median": None,
            "q75": None,
            "maximum": None,
            "unique_count": 0,
        }
    array = np.asarray(values, dtype=float)
    return {
        "count": len(values),
        "minimum": float(np.min(array)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "q75": float(np.quantile(array, 0.75)),
        "maximum": float(np.max(array)),
        "unique_count": len(set(float(value) for value in values)),
    }


def main() -> int:
    started = time.perf_counter()
    planned = [
        RANK_METRICS_PATH,
        RANK_ONSETS_CSV,
        RANK_ONSETS_MD,
        PAIRS_CSV,
        PAIRS_MD,
        SUMMARY_CSV,
        SUMMARY_MD,
        MANIFEST_PATH,
        TEXT_LOG_PATH,
    ]
    if any(path.exists() for path in planned):
        raise RuntimeError("rank候補解析出力の上書きを避けて停止")
    if not VERIFY_PATH.is_file() or not GROWTH_MANIFEST_PATH.is_file():
        raise RuntimeError("入力検証または成長候補解析の成功記録がない")
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    growth_manifest = json.loads(GROWTH_MANIFEST_PATH.read_text(encoding="utf-8"))
    if verification.get("success") is not True or growth_manifest.get("success") is not True:
        raise RuntimeError("先行工程が成功していない")

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    q_path = resolve(config["inputs"]["q_svd"]["path"])
    with q_path.open(newline="", encoding="utf-8") as handle:
        q_rows = list(csv.DictReader(handle))
    steps = np.array([int(row["step"]) for row in q_rows], dtype=np.int64)
    times = np.array([int(row["time"]) for row in q_rows], dtype=np.int64)
    relative_times = np.array(
        [int(row["relative_time"]) for row in q_rows],
        dtype=np.int64,
    )
    q = np.array(
        [
            [
                float(row["q1"]),
                float(row["q2"]),
                float(row["q3"]),
                float(row["q4"]),
            ]
            for row in q_rows
        ],
        dtype=float,
    )
    rank_existing = np.array([int(row["rank_q"]) for row in q_rows], dtype=np.int64)
    if not np.array_equal(steps, times):
        raise RuntimeError("qのstepとtimeが一致しない")

    thresholds = [float(value) for value in config["rank_relative_thresholds"]]
    rank_candidates = {}
    for threshold in thresholds:
        rank_candidates[threshold] = np.sum(
            q > threshold * q[:, [0]],
            axis=1,
        ).astype(np.int64)
    existing_recomputed = rank_candidates[1e-8]
    mismatch_rows = np.flatnonzero(existing_recomputed != rank_existing)
    if len(mismatch_rows):
        raise RuntimeError(
            "既存rank_Q列が count(q_j > 1e-8 q1) と一致しない: "
            + ",".join(str(int(value)) for value in mismatch_rows[:20])
        )

    q3_over_q1 = np.array([safe_ratio(a, b) for a, b in zip(q[:, 2], q[:, 0])])
    q4_over_q1 = np.array([safe_ratio(a, b) for a, b in zip(q[:, 3], q[:, 0])])
    min_q34_over_q1 = np.array(
        [safe_ratio(min(a, b), c) for a, b, c in zip(q[:, 2], q[:, 3], q[:, 0])]
    )
    q3_minus_q4 = q[:, 2] - q[:, 3]
    q3_over_q4 = np.array([safe_ratio(a, b) for a, b in zip(q[:, 2], q[:, 3])])

    metric_fields = [
        "step",
        "time",
        "relative_time",
        "q1",
        "q2",
        "q3",
        "q4",
        "rank_q_existing_1e-8",
        "q3_over_q1",
        "q4_over_q1",
        "min_q3_q4_over_q1",
        "q3_minus_q4",
        "q3_over_q4",
    ]
    for threshold in thresholds:
        metric_fields.append(f"rank_candidate_{threshold_tag(threshold)}")
    with RANK_METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=metric_fields)
        writer.writeheader()
        for index in range(len(q_rows)):
            row = {
                "step": int(steps[index]),
                "time": int(times[index]),
                "relative_time": int(relative_times[index]),
                "q1": q[index, 0],
                "q2": q[index, 1],
                "q3": q[index, 2],
                "q4": q[index, 3],
                "rank_q_existing_1e-8": int(rank_existing[index]),
                "q3_over_q1": q3_over_q1[index],
                "q4_over_q1": q4_over_q1[index],
                "min_q3_q4_over_q1": min_q34_over_q1[index],
                "q3_minus_q4": q3_minus_q4[index],
                "q3_over_q4": q3_over_q4[index],
            }
            for threshold in thresholds:
                row[f"rank_candidate_{threshold_tag(threshold)}"] = int(
                    rank_candidates[threshold][index]
                )
            writer.writerow({key: format_value(row[key]) for key in metric_fields})

    onset_fields = [
        "rank4_onset_candidate_id",
        "relative_threshold",
        "persistence_records",
        "persistence_unit",
        "status",
        "rank4_onset_candidate",
        "confirmation_step",
        "observed_span_steps",
        "start_record_index",
        "confirmation_record_index",
    ]
    onset_rows = []
    onset_counter = 0
    for threshold in thresholds:
        rank4 = rank_candidates[threshold] == 4
        for persistence in [
            int(value) for value in config["rank_persistence_records"]
        ]:
            onset_counter += 1
            run = first_true_record_run(rank4, persistence)
            found = run is not None
            start_index, end_index = run if run is not None else (None, None)
            onset_rows.append(
                {
                    "rank4_onset_candidate_id": f"R4O{onset_counter:03d}",
                    "relative_threshold": threshold,
                    "persistence_records": persistence,
                    "persistence_unit": config["rank_persistence_unit"],
                    "status": "found" if found else "not_found",
                    "rank4_onset_candidate": (
                        int(steps[start_index]) if found else None
                    ),
                    "confirmation_step": int(steps[end_index]) if found else None,
                    "observed_span_steps": (
                        int(steps[end_index] - steps[start_index] + 1)
                        if found
                        else None
                    ),
                    "start_record_index": start_index,
                    "confirmation_record_index": end_index,
                }
            )
    write_csv_and_markdown(
        RANK_ONSETS_CSV,
        RANK_ONSETS_MD,
        onset_fields,
        onset_rows,
        "rank4_onset_all_candidates",
    )

    with GROWTH_INTERVALS_PATH.open(newline="", encoding="utf-8") as handle:
        growth_intervals = list(csv.DictReader(handle))
    with GROWTH_ENDS_PATH.open(newline="", encoding="utf-8") as handle:
        growth_ends = list(csv.DictReader(handle))
    found_ends = [row for row in growth_ends if row["status"] == "found"]
    found_onsets = [row for row in onset_rows if row["status"] == "found"]

    pair_fields = [
        "pair_id",
        "growth_end_candidate_id",
        "growth_end_candidate",
        "end_condition",
        "end_persistence",
        "rank4_onset_candidate_id",
        "relative_threshold",
        "rank_persistence_records",
        "rank4_onset_candidate",
        "time_difference_rank4_minus_growth_end",
    ]
    pair_counter = 0
    pair_differences = []
    pair_differences_by_condition = {
        name: [] for name in config["growth_end_conditions"]
    }
    with PAIRS_CSV.open("w", newline="", encoding="utf-8") as csv_handle, PAIRS_MD.open(
        "w",
        encoding="utf-8",
    ) as md_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=pair_fields)
        writer.writeheader()
        md_handle.write("# growth_end_vs_rank4_onset_all_pairs\n\n")
        md_handle.write("| " + " | ".join(pair_fields) + " |\n")
        md_handle.write("|" + "|".join("---" for _ in pair_fields) + "|\n")
        for end_row in found_ends:
            end_step = int(end_row["growth_end_candidate"])
            for onset_row in found_onsets:
                onset_step = int(onset_row["rank4_onset_candidate"])
                difference = onset_step - end_step
                pair_counter += 1
                pair = {
                    "pair_id": f"PAIR{pair_counter:08d}",
                    "growth_end_candidate_id": end_row[
                        "growth_end_candidate_id"
                    ],
                    "growth_end_candidate": end_step,
                    "end_condition": end_row["end_condition"],
                    "end_persistence": int(end_row["end_persistence"]),
                    "rank4_onset_candidate_id": onset_row[
                        "rank4_onset_candidate_id"
                    ],
                    "relative_threshold": onset_row["relative_threshold"],
                    "rank_persistence_records": onset_row[
                        "persistence_records"
                    ],
                    "rank4_onset_candidate": onset_step,
                    "time_difference_rank4_minus_growth_end": difference,
                }
                formatted = {
                    key: format_value(pair.get(key)) for key in pair_fields
                }
                writer.writerow(formatted)
                md_handle.write(
                    "| "
                    + " | ".join(
                        formatted[key].replace("|", "\\|") for key in pair_fields
                    )
                    + " |\n"
                )
                pair_differences.append(difference)
                pair_differences_by_condition[end_row["end_condition"]].append(
                    difference
                )

    summary_fields = [
        "category",
        "group",
        "parameter",
        "value",
        "count",
        "minimum",
        "q25",
        "median",
        "q75",
        "maximum",
        "unique_count",
        "notes",
    ]
    summary_rows = []

    def add_summary(
        category: str,
        group: str,
        parameter: str,
        value,
        values: list[float],
        notes: str = "",
    ) -> None:
        stats = numeric_summary(values)
        summary_rows.append(
            {
                "category": category,
                "group": group,
                "parameter": parameter,
                "value": value,
                **stats,
                "notes": notes,
            }
        )

    interval_starts = [int(row["interval_start"]) for row in growth_intervals]
    interval_ends = [int(row["interval_end"]) for row in growth_intervals]
    add_summary(
        "growth_interval",
        "all",
        "interval_start",
        "all_parameters",
        interval_starts,
    )
    add_summary(
        "growth_interval",
        "all",
        "interval_end",
        "all_parameters",
        interval_ends,
    )
    for window in config["regression_windows"]:
        subset = [
            int(row["interval_start"])
            for row in growth_intervals
            if int(row["window"]) == int(window)
        ]
        add_summary(
            "growth_interval",
            "by_window",
            "window",
            window,
            subset,
        )

    found_end_steps = [int(row["growth_end_candidate"]) for row in found_ends]
    add_summary(
        "growth_end",
        "all_found",
        "candidate_step",
        "all_parameters",
        found_end_steps,
        notes=f"not_found rows retained separately: {len(growth_ends) - len(found_ends)}",
    )
    for condition_name in config["growth_end_conditions"]:
        subset = [
            int(row["growth_end_candidate"])
            for row in found_ends
            if row["end_condition"] == condition_name
        ]
        add_summary(
            "growth_end",
            "by_condition",
            "end_condition",
            condition_name,
            subset,
        )

    onset_steps = [
        int(row["rank4_onset_candidate"]) for row in found_onsets
    ]
    add_summary(
        "rank4_onset",
        "all_found",
        "candidate_step",
        "all_parameters",
        onset_steps,
        notes="persistence counts consecutive saved q records; no interpolation",
    )
    for threshold in thresholds:
        subset = [
            int(row["rank4_onset_candidate"])
            for row in found_onsets
            if float(row["relative_threshold"]) == threshold
        ]
        add_summary(
            "rank4_onset",
            "by_threshold",
            "relative_threshold",
            threshold,
            subset,
        )

    add_summary(
        "time_difference",
        "all_pairs",
        "rank4_minus_growth_end",
        "all_parameters",
        pair_differences,
    )
    for condition_name, values in pair_differences_by_condition.items():
        add_summary(
            "time_difference",
            "by_end_condition",
            "end_condition",
            condition_name,
            values,
        )

    for category, values in [
        ("growth_interval_start", interval_starts),
        ("growth_end", found_end_steps),
        ("rank4_onset", onset_steps),
    ]:
        frequencies = Counter(values)
        for rank, (step, count) in enumerate(
            sorted(frequencies.items(), key=lambda item: (-item[1], item[0]))[:10],
            1,
        ):
            summary_rows.append(
                {
                    "category": category,
                    "group": "recurring_exact_step",
                    "parameter": "frequency_rank",
                    "value": rank,
                    "count": count,
                    "minimum": step,
                    "q25": step,
                    "median": step,
                    "q75": step,
                    "maximum": step,
                    "unique_count": 1,
                    "notes": "頻度順位であり採用順位ではない",
                }
            )

    write_csv_and_markdown(
        SUMMARY_CSV,
        SUMMARY_MD,
        summary_fields,
        summary_rows,
        "candidate_summary",
    )

    manifest = {
        "stage": "A1",
        "success": True,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "q_input": str(q_path),
        "q_input_rows": len(q_rows),
        "existing_rank_definition": "count(q_j > 1e-8*q1), j=1..4",
        "existing_rank_mismatch_count": int(len(mismatch_rows)),
        "rank_relative_thresholds": thresholds,
        "rank_persistence_records": config["rank_persistence_records"],
        "rank_persistence_unit": config["rank_persistence_unit"],
        "rank4_onset_row_count": len(onset_rows),
        "rank4_onset_found_count": len(found_onsets),
        "growth_end_found_count": len(found_ends),
        "all_pair_count": pair_counter,
        "candidate_summary_row_count": len(summary_rows),
        "outputs": [
            str(RANK_METRICS_PATH),
            str(RANK_ONSETS_CSV),
            str(RANK_ONSETS_MD),
            str(PAIRS_CSV),
            str(PAIRS_MD),
            str(SUMMARY_CSV),
            str(SUMMARY_MD),
        ],
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    log_lines = [
        f"q_rows={len(q_rows)}",
        "existing_rank_definition_mismatch=0",
        f"rank4_onset_candidates={len(onset_rows)}",
        f"rank4_onset_found={len(found_onsets)}",
        f"growth_end_vs_rank4_pairs={pair_counter}",
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
