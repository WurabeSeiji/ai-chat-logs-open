#!/usr/bin/env python3
"""Compare isolated N=5 reproduction outputs with the locked Paper 7 baselines."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import struct
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
CONFIG_PATH = PACKAGE_ROOT / "config_locked.json"
EXPECTED_PATH = PACKAGE_ROOT / "expected_hashes.json"
RUN_MANIFEST_PATH = PACKAGE_ROOT / "logs" / "run_manifest.json"
COMPARISON_DIR = PACKAGE_ROOT / "comparison"
CSV_RESULT_PATH = COMPARISON_DIR / "csv_comparison.json"
JSON_RESULT_PATH = COMPARISON_DIR / "json_comparison.json"
FIGURE_RESULT_PATH = COMPARISON_DIR / "figure_comparison.json"
SUMMARY_PATH = COMPARISON_DIR / "comparison_summary.json"
LOG_PATH = PACKAGE_ROOT / "logs" / "compare_outputs.log"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bitwise_equal_float(a: float, b: float) -> bool:
    return struct.pack(">d", a) == struct.pack(">d", b)


def relative_error(expected: float, reproduced: float) -> float:
    absolute = abs(reproduced - expected)
    if expected == 0.0:
        return 0.0 if absolute == 0.0 else math.inf
    return absolute / abs(expected)


def display_number(value: float):
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "Infinity" if value > 0 else "-Infinity"
    return value


def resolve_expected(path: str) -> Path:
    return REPO_ROOT / path


def resolve_reproduced(path: str) -> Path:
    return PACKAGE_ROOT / path


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def can_parse_float(values: list[str]) -> bool:
    try:
        for value in values:
            float(value)
    except (TypeError, ValueError):
        return False
    return True


def compare_csv(
    name: str,
    expected_path: Path,
    reproduced_path: Path,
    exact_columns: set[str],
    abs_tol: float,
    rel_tol: float,
) -> dict:
    if not expected_path.is_file() or not reproduced_path.is_file():
        return {
            "name": name,
            "expected_path": str(expected_path),
            "reproduced_path": str(reproduced_path),
            "classification": "mismatch",
            "failure": "必須CSV欠落",
            "expected_exists": expected_path.is_file(),
            "reproduced_exists": reproduced_path.is_file(),
        }

    expected_header, expected_rows = read_csv(expected_path)
    reproduced_header, reproduced_rows = read_csv(reproduced_path)
    columns_match = expected_header == reproduced_header
    rows_match = len(expected_rows) == len(reproduced_rows)
    result = {
        "name": name,
        "expected_path": str(expected_path),
        "reproduced_path": str(reproduced_path),
        "expected_sha256": sha256(expected_path),
        "reproduced_sha256": sha256(reproduced_path),
        "column_names_match": columns_match,
        "expected_columns": expected_header,
        "reproduced_columns": reproduced_header,
        "expected_row_count": len(expected_rows),
        "reproduced_row_count": len(reproduced_rows),
        "row_count_match": rows_match,
        "column_results": {},
    }
    if not columns_match or not rows_match:
        result["classification"] = "mismatch"
        return result

    all_exact = True
    all_numerical = True
    for column in expected_header:
        expected_values = [row[column] for row in expected_rows]
        reproduced_values = [row[column] for row in reproduced_rows]
        if column in exact_columns:
            mismatch_rows = [
                index
                for index, (expected_value, reproduced_value) in enumerate(
                    zip(expected_values, reproduced_values)
                )
                if expected_value != reproduced_value
            ]
            ok = not mismatch_rows
            result["column_results"][column] = {
                "kind": "exact_string_or_integer",
                "exact_match": ok,
                "mismatch_count": len(mismatch_rows),
                "mismatch_rows": mismatch_rows,
            }
            all_exact = all_exact and ok
            all_numerical = all_numerical and ok
            continue

        if not can_parse_float(expected_values + reproduced_values):
            mismatch_rows = [
                index
                for index, (expected_value, reproduced_value) in enumerate(
                    zip(expected_values, reproduced_values)
                )
                if expected_value != reproduced_value
            ]
            ok = not mismatch_rows
            result["column_results"][column] = {
                "kind": "string",
                "exact_match": ok,
                "mismatch_count": len(mismatch_rows),
                "mismatch_rows": mismatch_rows,
            }
            all_exact = all_exact and ok
            all_numerical = all_numerical and ok
            continue

        expected_floats = [float(value) for value in expected_values]
        reproduced_floats = [float(value) for value in reproduced_values]
        expected_nan_rows = [index for index, value in enumerate(expected_floats) if math.isnan(value)]
        reproduced_nan_rows = [index for index, value in enumerate(reproduced_floats) if math.isnan(value)]
        nan_match = expected_nan_rows == reproduced_nan_rows
        bitwise_match = nan_match
        max_abs = 0.0
        max_rel = 0.0
        max_abs_row = None
        max_rel_row = None
        compared_count = 0
        for index, (expected_value, reproduced_value) in enumerate(
            zip(expected_floats, reproduced_floats)
        ):
            if math.isnan(expected_value) or math.isnan(reproduced_value):
                continue
            compared_count += 1
            if not bitwise_equal_float(expected_value, reproduced_value):
                bitwise_match = False
            if math.isinf(expected_value) or math.isinf(reproduced_value):
                absolute = 0.0 if expected_value == reproduced_value else math.inf
                relative = absolute
            else:
                absolute = abs(reproduced_value - expected_value)
                relative = relative_error(expected_value, reproduced_value)
            if absolute > max_abs:
                max_abs = absolute
                max_abs_row = index
            if relative > max_rel:
                max_rel = relative
                max_rel_row = index
        within_tolerance = nan_match and max_abs <= abs_tol and max_rel <= rel_tol
        result["column_results"][column] = {
            "kind": "float",
            "bitwise_match": bitwise_match,
            "nan_positions_match": nan_match,
            "expected_nan_rows": expected_nan_rows,
            "reproduced_nan_rows": reproduced_nan_rows,
            "compared_finite_or_infinite_count": compared_count,
            "max_absolute_error": display_number(max_abs),
            "max_absolute_error_row": max_abs_row,
            "max_relative_error": display_number(max_rel),
            "max_relative_error_row": max_rel_row,
            "within_tolerance": within_tolerance,
        }
        all_exact = all_exact and bitwise_match
        all_numerical = all_numerical and within_tolerance

    if all_exact:
        result["classification"] = "exact"
    elif all_numerical:
        result["classification"] = "numerical_match"
    else:
        result["classification"] = "mismatch"
    return result


def compare_json(
    name: str,
    expected_path: Path,
    reproduced_path: Path,
    abs_tol: float,
    rel_tol: float,
) -> dict:
    if not expected_path.is_file() or not reproduced_path.is_file():
        return {
            "name": name,
            "expected_path": str(expected_path),
            "reproduced_path": str(reproduced_path),
            "classification": "mismatch",
            "failure": "必須JSON欠落",
            "expected_exists": expected_path.is_file(),
            "reproduced_exists": reproduced_path.is_file(),
        }

    expected_data = json.loads(expected_path.read_text(encoding="utf-8"))
    reproduced_data = json.loads(reproduced_path.read_text(encoding="utf-8"))
    state = {
        "keys_match": True,
        "types_match": True,
        "exact_nonfloat_values_match": True,
        "all_float_bitwise": True,
        "all_float_numerical": True,
        "missing_keys": [],
        "extra_keys": [],
        "length_mismatches": [],
        "type_mismatches": [],
        "value_differences": [],
        "float_results": {},
    }

    def walk(expected_value, reproduced_value, path: str) -> None:
        if isinstance(expected_value, dict):
            if not isinstance(reproduced_value, dict):
                state["types_match"] = False
                state["type_mismatches"].append(
                    {"path": path, "expected": "dict", "reproduced": type(reproduced_value).__name__}
                )
                return
            expected_keys = set(expected_value)
            reproduced_keys = set(reproduced_value)
            missing = sorted(expected_keys - reproduced_keys)
            extra = sorted(reproduced_keys - expected_keys)
            if missing or extra:
                state["keys_match"] = False
                state["missing_keys"].extend(f"{path}.{key}" for key in missing)
                state["extra_keys"].extend(f"{path}.{key}" for key in extra)
            for key in sorted(expected_keys & reproduced_keys):
                walk(expected_value[key], reproduced_value[key], f"{path}.{key}")
            return

        if isinstance(expected_value, list):
            if not isinstance(reproduced_value, list):
                state["types_match"] = False
                state["type_mismatches"].append(
                    {"path": path, "expected": "list", "reproduced": type(reproduced_value).__name__}
                )
                return
            if len(expected_value) != len(reproduced_value):
                state["length_mismatches"].append(
                    {"path": path, "expected": len(expected_value), "reproduced": len(reproduced_value)}
                )
            for index, (expected_item, reproduced_item) in enumerate(
                zip(expected_value, reproduced_value)
            ):
                walk(expected_item, reproduced_item, f"{path}[{index}]")
            return

        if type(expected_value) is float:
            if type(reproduced_value) is not float:
                state["types_match"] = False
                state["type_mismatches"].append(
                    {"path": path, "expected": "float", "reproduced": type(reproduced_value).__name__}
                )
                return
            expected_nan = math.isnan(expected_value)
            reproduced_nan = math.isnan(reproduced_value)
            nan_match = expected_nan == reproduced_nan
            if expected_nan or reproduced_nan:
                bitwise = nan_match
                absolute = 0.0 if nan_match else math.inf
                relative = absolute
            elif math.isinf(expected_value) or math.isinf(reproduced_value):
                bitwise = bitwise_equal_float(expected_value, reproduced_value)
                absolute = 0.0 if expected_value == reproduced_value else math.inf
                relative = absolute
            else:
                bitwise = bitwise_equal_float(expected_value, reproduced_value)
                absolute = abs(reproduced_value - expected_value)
                relative = relative_error(expected_value, reproduced_value)
            within = nan_match and absolute <= abs_tol and relative <= rel_tol
            state["all_float_bitwise"] = state["all_float_bitwise"] and bitwise
            state["all_float_numerical"] = state["all_float_numerical"] and within
            state["float_results"][path] = {
                "expected": display_number(expected_value),
                "reproduced": display_number(reproduced_value),
                "bitwise_match": bitwise,
                "nan_position_match": nan_match,
                "absolute_error": display_number(absolute),
                "relative_error": display_number(relative),
                "within_tolerance": within,
            }
            if not within:
                state["value_differences"].append(
                    {
                        "path": path,
                        "expected": display_number(expected_value),
                        "reproduced": display_number(reproduced_value),
                    }
                )
            return

        if type(expected_value) is not type(reproduced_value):
            state["types_match"] = False
            state["type_mismatches"].append(
                {
                    "path": path,
                    "expected": type(expected_value).__name__,
                    "reproduced": type(reproduced_value).__name__,
                }
            )
            return
        if expected_value != reproduced_value:
            state["exact_nonfloat_values_match"] = False
            state["value_differences"].append(
                {"path": path, "expected": expected_value, "reproduced": reproduced_value}
            )

    walk(expected_data, reproduced_data, "$")
    structural = (
        state["keys_match"]
        and state["types_match"]
        and not state["length_mismatches"]
        and state["exact_nonfloat_values_match"]
    )
    if structural and state["all_float_bitwise"]:
        classification = "exact"
    elif structural and state["all_float_numerical"]:
        classification = "numerical_match"
    else:
        classification = "mismatch"
    return {
        "name": name,
        "expected_path": str(expected_path),
        "reproduced_path": str(reproduced_path),
        "expected_sha256": sha256(expected_path),
        "reproduced_sha256": sha256(reproduced_path),
        "classification": classification,
        **state,
    }


def png_dimensions(path: Path) -> list[int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError(f"PNG IHDRを読めない: {path}")
    width, height = struct.unpack(">II", header[16:24])
    return [width, height]


def compare_png(name: str, expected_path: Path, reproduced_path: Path) -> dict:
    expected_exists = expected_path.is_file()
    reproduced_exists = reproduced_path.is_file()
    result = {
        "name": name,
        "expected_path": str(expected_path),
        "reproduced_path": str(reproduced_path),
        "expected_exists": expected_exists,
        "reproduced_exists": reproduced_exists,
    }
    if not expected_exists or not reproduced_exists:
        result.update(
            {
                "dimensions_match": False,
                "sha256_match": False,
                "classification": "mismatch",
            }
        )
        return result
    expected_dimensions = png_dimensions(expected_path)
    reproduced_dimensions = png_dimensions(reproduced_path)
    expected_hash = sha256(expected_path)
    reproduced_hash = sha256(reproduced_path)
    dimensions_match = expected_dimensions == reproduced_dimensions
    sha_match = expected_hash == reproduced_hash
    result.update(
        {
            "expected_dimensions": expected_dimensions,
            "reproduced_dimensions": reproduced_dimensions,
            "dimensions_match": dimensions_match,
            "expected_sha256": expected_hash,
            "reproduced_sha256": reproduced_hash,
            "sha256_match": sha_match,
            "classification": (
                "exact"
                if dimensions_match and sha_match
                else "rendering_difference_only"
                if dimensions_match
                else "mismatch"
            ),
        }
    )
    return result


def main() -> int:
    started = time.perf_counter()
    if not RUN_MANIFEST_PATH.is_file():
        raise RuntimeError("run_reproduction.py のmanifestがない")
    run_manifest = json.loads(RUN_MANIFEST_PATH.read_text(encoding="utf-8"))
    if run_manifest.get("success") is not True or run_manifest.get("N_values_executed") != [5]:
        raise RuntimeError("run_reproduction.py がN=5だけで正常終了していない")

    for planned in (CSV_RESULT_PATH, JSON_RESULT_PATH, FIGURE_RESULT_PATH, SUMMARY_PATH):
        if planned.exists():
            raise RuntimeError(f"比較結果の上書きを避けて停止: {planned}")

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    abs_tol = float(config["absolute_tolerance"])
    rel_tol = float(config["relative_tolerance"])
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)

    csv_results = {}
    for name, item in expected["csv_baselines"].items():
        csv_results[name] = compare_csv(
            name,
            resolve_expected(item["expected"]),
            resolve_reproduced(item["reproduced"]),
            set(config["csv_exact_columns"][name]),
            abs_tol,
            rel_tol,
        )

    json_results = {}
    for name, item in expected["json_baselines"].items():
        json_results[name] = compare_json(
            name,
            resolve_expected(item["expected"]),
            resolve_reproduced(item["reproduced"]),
            abs_tol,
            rel_tol,
        )

    figure_results = {}
    for name, item in expected["png_baselines"].items():
        figure_results[name] = compare_png(
            name,
            resolve_expected(item["expected"]),
            resolve_reproduced(item["reproduced"]),
        )

    csv_classes = [item["classification"] for item in csv_results.values()]
    json_classes = [item["classification"] for item in json_results.values()]
    numeric_classes = csv_classes + json_classes
    if all(value == "exact" for value in numeric_classes):
        numerical_result = "exact"
    elif all(value in ("exact", "numerical_match") for value in numeric_classes):
        numerical_result = "numerical_match"
    else:
        numerical_result = "mismatch"
    figures_exist = all(item["reproduced_exists"] for item in figure_results.values())
    dimensions_match = all(item.get("dimensions_match", False) for item in figure_results.values())

    summary = {
        "stage": "A0",
        "locked_n": 5,
        "success": True,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "absolute_tolerance": abs_tol,
        "relative_tolerance": rel_tol,
        "csv_classifications": {
            name: item["classification"] for name, item in csv_results.items()
        },
        "json_classifications": {
            name: item["classification"] for name, item in json_results.items()
        },
        "figure_classifications": {
            name: item["classification"] for name, item in figure_results.items()
        },
        "numeric_comparison": numerical_result,
        "all_required_figures_exist": figures_exist,
        "all_figure_dimensions_match": dimensions_match,
        "image_sha256_alone_is_not_failure": True,
    }

    CSV_RESULT_PATH.write_text(
        json.dumps(csv_results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    JSON_RESULT_PATH.write_text(
        json.dumps(json_results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    FIGURE_RESULT_PATH.write_text(
        json.dumps(figure_results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    SUMMARY_PATH.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    lines = [
        f"CSV {name}: {item['classification']}" for name, item in csv_results.items()
    ]
    lines.extend(
        f"JSON {name}: {item['classification']}" for name, item in json_results.items()
    )
    lines.extend(
        f"PNG {name}: {item['classification']}" for name, item in figure_results.items()
    )
    lines.append(f"NUMERIC_COMPARISON: {numerical_result}")
    LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
