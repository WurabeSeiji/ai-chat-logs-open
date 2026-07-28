#!/usr/bin/env python3
"""Verify persisted Stage B rows and path arrays without accessing originals."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SEMANTIC_NOTICE = (
    "B_to_A_transfer is spectral cosine similarity of the A-channel state "
    "to the initial B spectrum; it is NOT a path-exchange norm."
)

AUDIT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = AUDIT_ROOT / "data" / "stage_B"
LOG_ROOT = AUDIT_ROOT / "logs"


def array_norm2(array: np.ndarray) -> np.ndarray:
    axes = tuple(range(1, array.ndim))
    return np.sum(np.abs(array) ** 2, axis=axes).real


def main() -> int:
    baseline_path = DATA_ROOT / "current_behavior_baseline.csv"
    state_path = DATA_ROOT / "state_parity_metrics.csv"
    controls_path = DATA_ROOT / "input_control_metrics.csv"
    diagnostics_path = DATA_ROOT / "current_behavior_diagnostics.json"
    with baseline_path.open(encoding="utf-8", newline="") as handle:
        baseline = list(csv.DictReader(handle))
    with state_path.open(encoding="utf-8", newline="") as handle:
        states = list(csv.DictReader(handle))
    with controls_path.open(encoding="utf-8", newline="") as handle:
        controls = list(csv.DictReader(handle))
    diagnostics = json.loads(diagnostics_path.read_text(encoding="utf-8"))

    if not baseline or not states or not controls:
        raise AssertionError("Stage B CSV output is empty")
    if any(row["B_to_A_transfer_semantics"] != SEMANTIC_NOTICE for row in baseline):
        raise AssertionError("baseline semantic notice mismatch")
    if any(row["B_to_A_transfer_semantics"] != SEMANTIC_NOTICE for row in states):
        raise AssertionError("state semantic notice mismatch")
    if any(row["B_to_A_transfer_semantics"] != SEMANTIC_NOTICE for row in controls):
        raise AssertionError("control semantic notice mismatch")
    if diagnostics["B_to_A_transfer_semantics"] != SEMANTIC_NOTICE:
        raise AssertionError("diagnostic semantic notice mismatch")
    if not all(diagnostics["assertions"].values()):
        raise AssertionError("reproduction diagnostic assertion is false")

    baseline_by_case: dict[str, list[dict]] = {}
    for row in baseline:
        baseline_by_case.setdefault(row["case_id"], []).append(row)

    maximums = {
        "persisted_path_norm_error": 0.0,
        "persisted_interference_sum_error": 0.0,
        "persisted_output_decomposition_error": 0.0,
    }
    array_shapes: dict[str, list[int]] = {}
    for case_id, rows in sorted(baseline_by_case.items()):
        rows.sort(key=lambda row: int(row["collision"]))
        path = DATA_ROOT / "path_arrays" / f"{case_id}_collision_paths.npz"
        with np.load(path, allow_pickle=False) as archive:
            if str(archive["B_to_A_transfer_semantics"]) != SEMANTIC_NOTICE:
                raise AssertionError(f"{case_id} NPZ semantic notice mismatch")
            collisions = archive["collisions"]
            expected_collisions = np.asarray(
                [int(row["collision"]) for row in rows], dtype=np.int64
            )
            if not np.array_equal(collisions, expected_collisions):
                raise AssertionError(f"{case_id} collision axis mismatch")
            array_shapes[case_id] = list(archive["path_a_to_a"].shape)
            path_a_to_a = archive["path_a_to_a"]
            path_b_to_a = archive["path_b_to_a"]
            path_b_to_b = archive["path_b_to_b"]
            path_a_to_b = archive["path_a_to_b"]
            interference_a = archive["interference_in_a_density"]
            interference_b = archive["interference_in_b_density"]
            stored_norms = {
                "path_a_to_a_norm_raw": array_norm2(path_a_to_a),
                "path_b_to_a_norm_raw": array_norm2(path_b_to_a),
                "path_b_to_b_norm_raw": array_norm2(path_b_to_b),
                "path_a_to_b_norm_raw": array_norm2(path_a_to_b),
            }
            for column, values in stored_norms.items():
                csv_values = np.asarray([float(row[column]) for row in rows])
                maximums["persisted_path_norm_error"] = max(
                    maximums["persisted_path_norm_error"],
                    float(np.max(np.abs(values - csv_values))),
                )
            for column, density in (
                ("interference_in_a_raw", interference_a),
                ("interference_in_b_raw", interference_b),
            ):
                density_sum = np.sum(density, axis=(1, 2))
                csv_values = np.asarray([float(row[column]) for row in rows])
                maximums["persisted_interference_sum_error"] = max(
                    maximums["persisted_interference_sum_error"],
                    float(np.max(np.abs(density_sum - csv_values))),
                )
            a_raw = path_a_to_a + path_b_to_a
            b_raw = path_b_to_b + path_a_to_b
            for column, values in (
                ("a_output_norm2_raw", array_norm2(a_raw)),
                ("b_output_norm2_raw", array_norm2(b_raw)),
            ):
                csv_values = np.asarray([float(row[column]) for row in rows])
                maximums["persisted_output_decomposition_error"] = max(
                    maximums["persisted_output_decomposition_error"],
                    float(np.max(np.abs(values - csv_values))),
                )

    tolerance = 3.0e-12
    assertions = {
        "all_expected_cases_have_path_arrays": set(baseline_by_case)
        == {"F1_x_F1", "FK_x_FK", "BK_x_BK", "FK_x_BK", "MIX_x_MIX"},
        "each_case_has_32_collision_arrays": all(
            shape == [32, 512, 16] for shape in array_shapes.values()
        ),
        "persisted_path_norms_match_csv": maximums["persisted_path_norm_error"]
        <= tolerance,
        "persisted_interference_matches_csv": maximums[
            "persisted_interference_sum_error"
        ]
        <= tolerance,
        "persisted_paths_reconstruct_raw_output_norms": maximums[
            "persisted_output_decomposition_error"
        ]
        <= tolerance,
        "semantic_notice_present_in_all_stage_B_data": True,
    }
    if not all(assertions.values()):
        failed = [name for name, passed in assertions.items() if not passed]
        raise AssertionError(f"Stage B persisted-output verification failed: {failed}")

    payload = {
        "schema": "stage_B_persisted_output_verification_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "B_to_A_transfer_semantics": SEMANTIC_NOTICE,
        "baseline_row_count": len(baseline),
        "state_row_count": len(states),
        "control_row_count": len(controls),
        "array_shapes": array_shapes,
        "maximum_errors": maximums,
        "assertions": assertions,
    }
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    output = LOG_ROOT / "stage_B_output_verification.json"
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"B_to_A_transfer semantics: {SEMANTIC_NOTICE}")
    print(output)
    print(json.dumps(assertions, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
