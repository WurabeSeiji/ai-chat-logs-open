#!/usr/bin/env python3
"""Verify the three locked Stage A0 CSV inputs without running dynamics."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PACKAGE_ROOT / "config_candidates.json"
LOG_DIR = PACKAGE_ROOT / "logs"
RESULT_PATH = LOG_DIR / "input_verification.json"
TEXT_LOG_PATH = LOG_DIR / "verify_inputs.log"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(relative: str) -> Path:
    return (PACKAGE_ROOT / relative).resolve()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def main() -> int:
    started = time.perf_counter()
    if RESULT_PATH.exists() or TEXT_LOG_PATH.exists():
        raise RuntimeError("入力検証結果の上書きを避けて停止")

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    if config.get("locked_n") != 5:
        raise RuntimeError("configがN=5に固定されていない")

    report_path = resolve(config["stage_a0_report"])
    comparison_path = resolve(config["stage_a0_csv_comparison"])
    if not report_path.is_file() or not comparison_path.is_file():
        raise RuntimeError("Stage A0報告書またはCSV比較記録が欠落")
    report_text = report_path.read_text(encoding="utf-8")
    if "REPRODUCED_EXACTLY" not in report_text:
        raise RuntimeError("Stage A0報告書の総合判定が完全再現ではない")
    stage_a0_comparison = json.loads(comparison_path.read_text(encoding="utf-8"))

    checks = {}
    loaded = {}
    failures = []
    for name, item in config["inputs"].items():
        path = resolve(item["path"])
        key = item["stage_a0_comparison_key"]
        record = stage_a0_comparison.get(key)
        exists = path.is_file()
        actual_hash = sha256(path) if exists else None
        columns, rows = read_rows(path) if exists else ([], [])
        missing_columns = [column for column in item["required_columns"] if column not in columns]
        a0_record_ok = bool(
            record
            and record.get("classification") == "exact"
            and record.get("reproduced_sha256") == item["sha256"]
        )
        report_row_ok = f"`{key}`" in report_text and "`exact`" in report_text
        ok = (
            exists
            and actual_hash == item["sha256"]
            and a0_record_ok
            and report_row_ok
            and not missing_columns
        )
        check = {
            "name": name,
            "path": str(path),
            "exists": exists,
            "expected_sha256": item["sha256"],
            "actual_sha256": actual_hash,
            "sha256_match": actual_hash == item["sha256"],
            "stage_a0_comparison_record_match": a0_record_ok,
            "stage_a0_report_exact_row_present": report_row_ok,
            "columns": columns,
            "required_columns": item["required_columns"],
            "missing_columns": missing_columns,
            "row_count": len(rows),
            "success": ok,
        }
        checks[name] = check
        if not ok:
            failures.append(check)
        loaded[name] = rows

    time_axis = {"success": False}
    if not failures:
        crossing = int(config["existing_crossing"])
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
        q_step_time_match = q_steps == q_times
        q_crossing_offsets = sorted(
            set(step - relative for step, relative in zip(q_steps, q_relative))
        )
        q_strictly_increasing = all(b > a for a, b in zip(q_steps, q_steps[1:]))
        long_step_time_match = long_steps == long_times
        long_strictly_increasing = all(
            b > a for a, b in zip(long_steps, long_steps[1:])
        )
        long_flags_match = all(
            flag == int(step >= crossing)
            for step, flag in zip(long_steps, long_flags)
        )
        common_steps = sorted(set(f_steps) & set(long_steps))
        time_axis_success = all(
            [
                f_contiguous,
                f_crossing == crossing,
                q_step_time_match,
                q_crossing_offsets == [crossing],
                q_strictly_increasing,
                long_step_time_match,
                long_strictly_increasing,
                long_flags_match,
                bool(common_steps),
            ]
        )
        time_axis = {
            "success": time_axis_success,
            "shared_axis": "absolute_step",
            "f_step_range": [f_steps[0], f_steps[-1]],
            "f_contiguous_unit_steps": f_contiguous,
            "f_existing_crossing_from_f_gt_0_05": f_crossing,
            "q_step_range": [q_steps[0], q_steps[-1]],
            "q_step_equals_time": q_step_time_match,
            "q_step_minus_relative_time_unique": q_crossing_offsets,
            "q_steps_strictly_increasing": q_strictly_increasing,
            "paper7_step_range": [long_steps[0], long_steps[-1]],
            "paper7_step_equals_time": long_step_time_match,
            "paper7_steps_strictly_increasing": long_strictly_increasing,
            "paper7_crossing_flag_matches_saved_steps": long_flags_match,
            "common_f_paper7_step_count": len(common_steps),
            "rank_persistence_interpretation": config["rank_persistence_unit"],
            "unobserved_q_steps_interpolated": False,
        }
        if not time_axis_success:
            failures.append({"name": "time_axis", **time_axis})

    result = {
        "stage": "A1",
        "locked_n": 5,
        "success": not failures,
        "status": "CALIBRATION_INPUT_VERIFIED" if not failures else "INPUT_MISMATCH",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "stage_a0_report": str(report_path),
        "stage_a0_csv_comparison": str(comparison_path),
        "input_checks": checks,
        "time_axis": time_axis,
        "failures": failures,
    }
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = []
    for name, check in checks.items():
        lines.append(
            f"[{'OK' if check['success'] else 'NG'}] {name} "
            f"sha256={check['actual_sha256']} rows={check['row_count']}"
        )
    lines.append(
        f"[{'OK' if time_axis.get('success') else 'NG'}] time_axis absolute_step"
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
