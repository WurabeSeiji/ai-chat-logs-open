#!/usr/bin/env python3
"""Verify the three Stage A0 N=5 CSVs for read-only Stage A1b analysis."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
A0_ROOT = PACKAGE_ROOT.parent / "paper7_N5_reproduction"
REPORT_PATH = A0_ROOT / "reports" / "paper7_N5_reproduction_report.md"
COMPARISON_PATH = A0_ROOT / "comparison" / "csv_comparison.json"
LOG_DIR = PACKAGE_ROOT / "logs"
RESULT_PATH = LOG_DIR / "input_verification.json"
TEXT_LOG_PATH = LOG_DIR / "verify_inputs.log"
CROSSING = 1167
OBSERVATION_END = 3000

INPUTS = {
    "fcurve": {
        "path": A0_ROOT
        / "reproduced/metastable_series_result_v1/fcurve_N00005_delta1e-15_seed0.csv",
        "sha256": "9220c5f3c1f570c8a52ea24a3cdd95568354cea0943d9bee7d8ed20316d3a9d0",
        "comparison_key": "fcurve_N00005",
        "required_columns": ["tau", "f"],
    },
    "q_svd": {
        "path": A0_ROOT
        / "reproduced/exact_lowN_eigenspectrum_v2/raw/N00005_dimension_saturation_v2/q_svd_N00005.csv",
        "sha256": "7c16a364c6cc9145293c2625dfe4ebb1f9962655d212679188215e8fad5e5155",
        "comparison_key": "q_svd_N00005",
        "required_columns": [
            "step",
            "time",
            "relative_time",
            "q1",
            "q2",
            "q3",
            "q4",
            "rank_q",
        ],
    },
    "paper7_long": {
        "path": A0_ROOT
        / "reproduced/exact_lowN_eigenspectrum_v2/paper7_longtime/raw/N00005/paper7_long_timeseries.csv",
        "sha256": "efeaf9dab753c057ad0c6109b9e4a8919f8d8db1249da186658bfed9fda784e3",
        "comparison_key": "paper7_long_timeseries_N00005",
        "required_columns": [
            "step",
            "time",
            "crossing_flag",
            "splitting_fraction",
            "direction_1_occupation",
            "direction_2_occupation",
            "direction_3_occupation",
            "direction_4_occupation",
            "other_rotating_occupation",
            "kernel_occupation",
        ],
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def unique_positive_differences(values: list[int]) -> list[int]:
    return sorted(set(b - a for a, b in zip(values, values[1:])))


def main() -> int:
    started = time.perf_counter()
    if RESULT_PATH.exists() or TEXT_LOG_PATH.exists():
        raise RuntimeError("入力検証出力の上書きを避けて停止")
    if not REPORT_PATH.is_file() or not COMPARISON_PATH.is_file():
        raise RuntimeError("Stage A0報告書またはCSV比較記録が欠落")
    report_text = REPORT_PATH.read_text(encoding="utf-8")
    if "REPRODUCED_EXACTLY" not in report_text:
        raise RuntimeError("Stage A0が完全再現として承認された記録でない")
    comparison = json.loads(COMPARISON_PATH.read_text(encoding="utf-8"))

    loaded = {}
    checks = {}
    failures = []
    for name, item in INPUTS.items():
        path = item["path"].resolve()
        exists = path.is_file()
        actual_hash = sha256(path) if exists else None
        columns, rows = load_csv(path) if exists else ([], [])
        missing_columns = [
            column for column in item["required_columns"] if column not in columns
        ]
        comparison_row = comparison.get(item["comparison_key"], {})
        a0_match = (
            comparison_row.get("classification") == "exact"
            and comparison_row.get("reproduced_sha256") == item["sha256"]
        )
        report_match = f"`{item['comparison_key']}`" in report_text
        success = (
            exists
            and actual_hash == item["sha256"]
            and a0_match
            and report_match
            and not missing_columns
        )
        checks[name] = {
            "path": str(path),
            "expected_sha256": item["sha256"],
            "actual_sha256": actual_hash,
            "sha256_match": actual_hash == item["sha256"],
            "stage_a0_exact_comparison_match": a0_match,
            "stage_a0_report_row_present": report_match,
            "columns": columns,
            "missing_columns": missing_columns,
            "row_count": len(rows),
            "success": success,
        }
        loaded[name] = rows
        if not success:
            failures.append({"name": name, **checks[name]})

    time_axis = {"success": False}
    if not failures:
        f_rows = loaded["fcurve"]
        q_rows = loaded["q_svd"]
        long_rows = loaded["paper7_long"]
        f_steps = [int(row["tau"]) for row in f_rows]
        f_values = [float(row["f"]) for row in f_rows]
        q_steps = [int(row["step"]) for row in q_rows]
        q_times = [int(row["time"]) for row in q_rows]
        q_relative = [int(row["relative_time"]) for row in q_rows]
        long_steps = [int(row["step"]) for row in long_rows]
        long_times = [int(row["time"]) for row in long_rows]
        long_flags = [int(row["crossing_flag"]) for row in long_rows]

        f_contiguous = f_steps == list(range(f_steps[0], f_steps[-1] + 1))
        f_crossing = next(
            (step for step, value in zip(f_steps, f_values) if value > 0.05),
            None,
        )
        q_absolute = q_steps == q_times
        q_offset = sorted(
            set(step - relative for step, relative in zip(q_steps, q_relative))
        )
        long_absolute = long_steps == long_times
        flag_match = all(
            flag == int(step >= CROSSING)
            for step, flag in zip(long_steps, long_flags)
        )
        q_observation_steps = [step for step in q_steps if 0 <= step <= OBSERVATION_END]
        long_observation_steps = [
            step for step in long_steps if 0 <= step <= OBSERVATION_END
        ]
        q_intervals = unique_positive_differences(q_observation_steps)
        long_intervals = unique_positive_differences(long_observation_steps)
        success = all(
            [
                f_contiguous,
                f_crossing == CROSSING,
                q_absolute,
                q_offset == [CROSSING],
                long_absolute,
                flag_match,
                q_intervals == [5],
                long_intervals == [25],
            ]
        )
        time_axis = {
            "success": success,
            "axis": "absolute_step",
            "observation_range": [0, OBSERVATION_END],
            "f_full_range": [f_steps[0], f_steps[-1]],
            "f_interval_in_observation_range": [1],
            "f_contiguous": f_contiguous,
            "f_crossing_from_existing_rule": f_crossing,
            "q_full_range": [q_steps[0], q_steps[-1]],
            "q_step_equals_time": q_absolute,
            "q_step_minus_relative_time": q_offset,
            "q_intervals_in_observation_range": q_intervals,
            "paper7_full_range": [long_steps[0], long_steps[-1]],
            "paper7_step_equals_time": long_absolute,
            "paper7_intervals_in_observation_range": long_intervals,
            "paper7_crossing_flags_match": flag_match,
            "q_interpolation_for_analysis": False,
            "occupation_interpolation_for_analysis": False,
        }
        if not success:
            failures.append({"name": "time_axis", **time_axis})

    result = {
        "stage": "A1b",
        "locked_n": 5,
        "success": not failures,
        "status": (
            "TRANSITION_INPUT_VERIFIED" if not failures else "INPUT_MISMATCH"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "stage_a0_report": str(REPORT_PATH.resolve()),
        "stage_a0_csv_comparison": str(COMPARISON_PATH.resolve()),
        "input_checks": checks,
        "time_axis": time_axis,
        "failures": failures,
    }
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        f"[{'OK' if check['success'] else 'NG'}] {name} "
        f"sha256={check['actual_sha256']} rows={check['row_count']}"
        for name, check in checks.items()
    ]
    lines.append(
        f"[{'OK' if time_axis.get('success') else 'NG'}] "
        "absolute step and save intervals"
    )
    lines.append("SUCCESS" if result["success"] else "STOP: INPUT_MISMATCH")
    TEXT_LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if result["success"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
