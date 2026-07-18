#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import statistics
from decimal import Decimal
from pathlib import Path

import run_minimal_system_B_gray_direct_check_v5 as v5


HERE = Path(__file__).resolve().parent
CONTROL_CSV = (
    HERE.parent
    / "working_layout"
    / "20260715"
    / "full_sweep_control"
    / "high_to_ext_full_delta1e-7_candidates_v5.csv"
)
SKIP_CSV = (
    HERE
    / "full_sweep_skip_normalize"
    / "high_to_ext_full_delta1e-7_candidates_skip_normalize_v1.csv"
)
OUTPUT_JSON = HERE / "normalize_skip_full_sweep_comparison_v1.json"

EXACT_ROOT_CASES = (
    (
        "R124_23",
        HERE / "R124_23_normalize_control_v1.json",
        HERE / "R124_23_skip_normalize_v1.json",
    ),
    (
        "R122_23",
        HERE / "R122_23_control_v1.json",
        HERE / "R122_23_skip_v1.json",
    ),
    (
        "R567_107",
        HERE / "R567_107_control_v1.json",
        HERE / "R567_107_skip_v1.json",
    ),
    (
        "R620_117",
        HERE / "R620_117_control_v1.json",
        HERE / "R620_117_skip_v1.json",
    ),
)

FLOAT_FIELDS = (
    "best_prefix_gray_error_no_phase",
    "best_prefix_gray_depth_no_phase",
    "best_S_mean",
    "best_S_amp",
    "best_S_drift",
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def bands(rows: list[dict[str, str]]) -> list[list[dict[str, str]]]:
    if not rows:
        return []
    result: list[list[dict[str, str]]] = [[rows[0]]]
    max_contiguous_gap = Decimal("0.0000001000000001")
    for row in rows[1:]:
        previous = result[-1][-1]
        gap = Decimal(row["R_input_text"]) - Decimal(previous["R_input_text"])
        if gap <= max_contiguous_gap:
            result[-1].append(row)
        else:
            result.append([row])
    return result


def peak(row_group: list[dict[str, str]]) -> dict[str, object]:
    best = max(
        row_group,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    return {
        "R_start": row_group[0]["R_input_text"],
        "R_end": row_group[-1]["R_input_text"],
        "point_count": len(row_group),
        "peak_R": best["R_input_text"],
        "peak_depth": float(best["best_prefix_gray_depth_no_phase"]),
        "peak_error": float(best["best_prefix_gray_error_no_phase"]),
        "peak_step": int(best["best_step"]),
    }


def norm_drift(
    reflection_rate: float,
    update_count: int,
    normalize_enabled: bool,
) -> dict[str, float]:
    a, b = v5.state_from_s_phi(0.01, 0.0, normalize_enabled)
    t, r, _t_power, _r_power = v5.scattering_coefficients(reflection_rate)
    q_values = [abs(a) ** 2 + abs(b) ** 2]
    for _ in range(update_count):
        a_next = r * a + t * b
        b_next = t * a + r * b
        if normalize_enabled:
            a_next, b_next = v5.normalize_pair(a_next, b_next)
        a, b = a_next, b_next
        q_values.append(abs(a) ** 2 + abs(b) ** 2)
    return {
        "initial_q": q_values[0],
        "final_q": q_values[-1],
        "max_abs_q_minus_1": max(abs(value - 1.0) for value in q_values),
    }


def main() -> None:
    control_rows = read_rows(CONTROL_CSV)
    skip_rows = read_rows(SKIP_CSV)
    control_by_r = {row["R_input_text"]: row for row in control_rows}
    skip_by_r = {row["R_input_text"]: row for row in skip_rows}
    control_keys = set(control_by_r)
    skip_keys = set(skip_by_r)
    shared_keys = sorted(control_keys & skip_keys, key=Decimal)

    field_differences: dict[str, dict[str, float]] = {}
    for field in FLOAT_FIELDS:
        differences = [
            abs(float(skip_by_r[key][field]) - float(control_by_r[key][field]))
            for key in shared_keys
        ]
        field_differences[field] = {
            "max_abs_difference": max(differences, default=0.0),
            "mean_abs_difference": statistics.fmean(differences) if differences else 0.0,
            "median_abs_difference": statistics.median(differences) if differences else 0.0,
        }

    control_bands = bands(control_rows)
    skip_bands = bands(skip_rows)
    band_comparison = []
    for index, (control_band, skip_band) in enumerate(
        zip(control_bands, skip_bands),
        start=1,
    ):
        control_peak = peak(control_band)
        skip_peak = peak(skip_band)
        band_comparison.append(
            {
                "band_index_by_R": index,
                "control": control_peak,
                "skip_normalize": skip_peak,
                "peak_R_same": control_peak["peak_R"] == skip_peak["peak_R"],
                "peak_depth_difference": (
                    float(skip_peak["peak_depth"])
                    - float(control_peak["peak_depth"])
                ),
            }
        )

    exact_root_comparison = []
    for label, control_path, skip_path in EXACT_ROOT_CASES:
        control_payload = json.loads(control_path.read_text(encoding="utf-8"))
        skip_payload = json.loads(skip_path.read_text(encoding="utf-8"))
        exact_root_comparison.append(
            {
                "root": label,
                "control_best_step": int(control_payload["best_step"]),
                "skip_best_step": int(skip_payload["best_step"]),
                "control_error": float(
                    control_payload["best_prefix_gray_error_no_phase"]
                ),
                "skip_error": float(
                    skip_payload["best_prefix_gray_error_no_phase"]
                ),
                "control_depth": float(
                    control_payload["best_prefix_gray_depth_no_phase"]
                ),
                "skip_depth": float(
                    skip_payload["best_prefix_gray_depth_no_phase"]
                ),
                "depth_difference": float(
                    skip_payload["best_prefix_gray_depth_no_phase"]
                )
                - float(control_payload["best_prefix_gray_depth_no_phase"]),
                "control_norm_drift": norm_drift(
                    float(control_payload["R"]),
                    int(control_payload["stopped_at_step"]),
                    True,
                ),
                "skip_norm_drift": norm_drift(
                    float(skip_payload["R"]),
                    int(skip_payload["stopped_at_step"]),
                    False,
                ),
            }
        )

    result = {
        "control_csv": str(CONTROL_CSV),
        "skip_csv": str(SKIP_CSV),
        "control_candidate_count": len(control_rows),
        "skip_candidate_count": len(skip_rows),
        "candidate_R_sets_identical": control_keys == skip_keys,
        "removed_candidate_count": len(control_keys - skip_keys),
        "added_candidate_count": len(skip_keys - control_keys),
        "shared_candidate_count": len(shared_keys),
        "best_step_changed_count": sum(
            control_by_r[key]["best_step"] != skip_by_r[key]["best_step"]
            for key in shared_keys
        ),
        "stop_step_changed_count": sum(
            control_by_r[key]["stopped_at_step"]
            != skip_by_r[key]["stopped_at_step"]
            for key in shared_keys
        ),
        "field_differences": field_differences,
        "control_band_count": len(control_bands),
        "skip_band_count": len(skip_bands),
        "band_comparison": band_comparison,
        "exact_root_comparison": exact_root_comparison,
    }
    OUTPUT_JSON.write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
