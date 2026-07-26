#!/usr/bin/env python3
"""Stage A2aの同一条件2実行をbitwise/数値健全性について比較する。"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
PROCESSED = HERE / "processed"
LOGS = HERE / "logs"
CONFIG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
RUN1, RUN2 = CONFIG["run_ids"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader.fieldnames or []), list(reader)


def write_csv_md(stem: str, columns: list[str], rows: list[dict]) -> None:
    csv_path = PROCESSED / f"{stem}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    md = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in rows:
        md.append("| " + " | ".join(str(row.get(c, "")) for c in columns) + " |")
    (PROCESSED / f"{stem}.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def compare_csv(filename: str) -> tuple[list[dict], bool]:
    c1, r1 = read_csv(RAW / RUN1 / filename)
    c2, r2 = read_csv(RAW / RUN2 / filename)
    header_match = c1 == c2
    row_count_match = len(r1) == len(r2)
    rows: list[dict] = []
    exact_all = header_match and row_count_match
    for column in c1:
        strings1 = [r[column] for r in r1]
        strings2 = [r[column] for r in r2] if header_match and row_count_match else []
        bitwise = strings1 == strings2
        max_abs = ""
        if header_match and row_count_match:
            try:
                a = np.asarray([float(x) for x in strings1], dtype=np.float64)
                b = np.asarray([float(x) for x in strings2], dtype=np.float64)
                finite = np.isfinite(a) & np.isfinite(b)
                max_abs = format(float(np.max(np.abs(a[finite] - b[finite]))) if np.any(finite) else 0.0, ".17e")
            except ValueError:
                max_abs = ""
        rows.append({
            "artifact": filename,
            "column": column,
            "rows_exec1": len(r1),
            "rows_exec2": len(r2),
            "header_match": str(header_match),
            "bitwise_match": str(bitwise),
            "max_absolute_difference": max_abs,
            "sha256_exec1": sha256(RAW / RUN1 / filename),
            "sha256_exec2": sha256(RAW / RUN2 / filename),
        })
        exact_all = exact_all and bitwise
    return rows, exact_all


def main() -> None:
    manifest_path = LOGS / "execution_manifest.json"
    if not manifest_path.is_file() or json.loads(manifest_path.read_text(encoding="utf-8")).get("status") != "COMPLETED":
        raise SystemExit("EXECUTION_FAILED: run_seedless.pyの完了記録がない")
    for run in (RUN1, RUN2):
        summary = RAW / run / "run_summary.json"
        if not summary.is_file() or json.loads(summary.read_text(encoding="utf-8")).get("final_step") != 5000:
            raise SystemExit(f"EXECUTION_FAILED: {run}がstep 5000まで完了していない")

    PROCESSED.mkdir(parents=True, exist_ok=True)
    comparison_rows: list[dict] = []
    exact = True
    for filename in ("f_timeseries.csv", "q_timeseries.csv", "occupation_timeseries.csv",
                     "first_passage_measurements.csv"):
        rows, file_exact = compare_csv(filename)
        comparison_rows.extend(rows)
        exact = exact and file_exact

    for filename in ("dominant_plane_steps.npy", "dominant_plane_values.npy"):
        p1, p2 = RAW / RUN1 / filename, RAW / RUN2 / filename
        a, b = np.load(p1), np.load(p2)
        array_equal = np.array_equal(a, b, equal_nan=True)
        max_abs = float(np.max(np.abs(a - b))) if a.size and b.size and a.shape == b.shape else float("nan")
        hash1, hash2 = sha256(p1), sha256(p2)
        comparison_rows.append({
            "artifact": filename,
            "column": "(array)",
            "rows_exec1": a.shape[0] if a.ndim else 1,
            "rows_exec2": b.shape[0] if b.ndim else 1,
            "header_match": str(a.shape == b.shape),
            "bitwise_match": str(array_equal and hash1 == hash2),
            "max_absolute_difference": format(max_abs, ".17e"),
            "sha256_exec1": hash1,
            "sha256_exec2": hash2,
        })
        exact = exact and array_equal and hash1 == hash2

    compare_columns = [
        "artifact", "column", "rows_exec1", "rows_exec2", "header_match",
        "bitwise_match", "max_absolute_difference", "sha256_exec1", "sha256_exec2",
    ]
    write_csv_md("exec1_vs_exec2", compare_columns, comparison_rows)

    health_rows: list[dict] = []
    health_ok = True
    for run in (RUN1, RUN2):
        fcols, frows = read_csv(RAW / run / "f_timeseries.csv")
        qcols, qrows = read_csv(RAW / run / "q_timeseries.csv")
        ocols, orows = read_csv(RAW / run / "occupation_timeseries.csv")
        datasets = [("f_timeseries", fcols, frows), ("q_timeseries", qcols, qrows),
                    ("occupation_timeseries", ocols, orows)]
        for dataset, columns, rows in datasets:
            for column in columns:
                if column == "step":
                    continue
                values = np.asarray([float(r[column]) for r in rows], dtype=np.float64)
                nonfinite = int(np.sum(~np.isfinite(values)))
                health_rows.append({
                    "run_id": run,
                    "dataset": dataset,
                    "metric": column,
                    "row_count": len(values),
                    "nonfinite_count": nonfinite,
                    "minimum": format(float(np.nanmin(values)), ".17e"),
                    "maximum": format(float(np.nanmax(values)), ".17e"),
                })
                health_ok = health_ok and nonfinite == 0
        steps = np.asarray([int(r["step"]) for r in frows], dtype=int)
        step_ok = np.array_equal(steps, np.arange(5001))
        health_rows.append({
            "run_id": run,
            "dataset": "f_timeseries",
            "metric": "step_continuity_0_to_5000",
            "row_count": len(steps),
            "nonfinite_count": 0,
            "minimum": str(int(steps.min())),
            "maximum": str(int(steps.max())),
        })
        health_ok = health_ok and step_ok

    health_columns = ["run_id", "dataset", "metric", "row_count", "nonfinite_count", "minimum", "maximum"]
    write_csv_md("numerical_health", health_columns, health_rows)
    result = {
        "stage": "A2a",
        "status": "COMPARED",
        "execs_bitwise_identical": exact,
        "numerical_health_passed": health_ok,
        "run_ids": [RUN1, RUN2],
    }
    (PROCESSED / "exec_comparison_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"COMPARED: execs_bitwise_identical={exact}, numerical_health_passed={health_ok}")


if __name__ == "__main__":
    main()
