from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import Stage2Params
from run_abc_baseline_stationary_wave_stage3_cross_readout_preliminary_v1 import local_rows
from run_abc_baseline_stationary_wave_transition_protocol_preliminary_v1 import protocol_profiles


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_stage3_noise_robustness_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


GAUGE_COUNT = 97
ZERO_MEAN_LEVELS = [0.0, 0.05, 0.20, 0.50, 1.00]
COMMON_BIAS_LEVELS = [1.0e-4, 1.0e-3, 1.0e-2]
COMMON_BIAS_DETECTION_FLOOR = 1.0e-13


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def sign(value: float, floor: float = 1.0e-18) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


def max_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return max(abs(float(row[key])) for row in rows)


def gauge_wave(gauge_index: int, channel_index: int) -> float:
    theta = 2.0 * math.pi * gauge_index / GAUGE_COUNT
    return (
        math.sin(theta + 0.73 * channel_index)
        + 0.37 * math.sin(3.0 * theta + 0.41 * channel_index)
        - 0.19 * math.cos(5.0 * theta - 0.29 * channel_index)
    )


def zero_mean_patterns() -> Dict[str, List[float]]:
    patterns: Dict[str, List[float]] = {}
    for channel_index, channel in enumerate(["chi", "p", "g1"]):
        values = [gauge_wave(gauge, channel_index) for gauge in range(GAUGE_COUNT)]
        mean = float(np.mean(values))
        centered = [value - mean for value in values]
        max_value = max(abs(value) for value in centered)
        patterns[channel] = [value / max_value for value in centered]
    return patterns


def collect_true_rows(params: Stage2Params) -> Tuple[List[Dict[str, Any]], float]:
    rows: List[Dict[str, Any]] = []
    for protocol, (profile, metadata) in protocol_profiles(params).items():
        rows.extend(local_rows(protocol, profile, metadata, params))

    reference_scale = max(
        max_abs(rows, "chi_curvature_Ra"),
        max_abs(rows, "p_change_Ra"),
        max_abs(rows, "g1_slope_Ra"),
    )
    return rows, reference_scale


def noisy_value(
    true_value: float,
    mode: str,
    level: float,
    reference_scale: float,
    gauge_index: int,
    channel: str,
    patterns: Dict[str, List[float]],
) -> float:
    if mode == "zero_mean":
        return true_value + level * reference_scale * patterns[channel][gauge_index]
    if mode == "common_bias":
        channel_bias = {"chi": 1.0, "p": -0.50, "g1": 0.25}[channel]
        return true_value + level * reference_scale * channel_bias
    raise ValueError(f"unknown mode: {mode}")


def evaluate_case(
    true_rows: List[Dict[str, Any]],
    params: Stage2Params,
    mode: str,
    level: float,
    reference_scale: float,
    patterns: Dict[str, List[float]],
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    noisy_rows: List[Dict[str, Any]] = []
    eval_rows: List[Dict[str, Any]] = []
    for row_index, row in enumerate(true_rows):
        gauge_values: Dict[str, List[float]] = {"chi": [], "p": [], "g1": []}
        true_values = {
            "chi": float(row["chi_curvature_Ra"]),
            "p": float(row["p_change_Ra"]),
            "g1": float(row["g1_slope_Ra"]),
        }
        for gauge_index in range(GAUGE_COUNT):
            values = {
                channel: noisy_value(
                    true_values[channel],
                    mode,
                    level,
                    reference_scale,
                    gauge_index,
                    channel,
                    patterns,
                )
                for channel in ["chi", "p", "g1"]
            }
            for channel, value in values.items():
                gauge_values[channel].append(value)
            if abs(float(row["delta_C_memory_forward"])) > params.effect_floor:
                sign_match = (
                    sign(values["chi"]) == sign(values["p"]) == sign(values["g1"])
                    and sign(values["chi"]) == sign(true_values["chi"])
                )
            else:
                sign_match = True
            noisy_rows.append(
                {
                    "mode": mode,
                    "level": level,
                    "row_index": row_index,
                    "protocol": row["protocol"],
                    "kind": row["kind"],
                    "duration": row["duration"],
                    "step": row["step"],
                    "particle": row["particle"],
                    "gauge_index": gauge_index,
                    "chi_readout_noisy": values["chi"],
                    "p_readout_noisy": values["p"],
                    "g1_readout_noisy": values["g1"],
                    "true_chi_curvature_Ra": true_values["chi"],
                    "true_p_change_Ra": true_values["p"],
                    "true_g1_slope_Ra": true_values["g1"],
                    "single_gauge_sign_match": sign_match,
                }
            )

        chi_mean = float(np.mean(gauge_values["chi"]))
        p_mean = float(np.mean(gauge_values["p"]))
        g1_mean = float(np.mean(gauge_values["g1"]))
        eval_rows.append(
            {
                "mode": mode,
                "level": level,
                "row_index": row_index,
                "protocol": row["protocol"],
                "kind": row["kind"],
                "duration": row["duration"],
                "step": row["step"],
                "particle": row["particle"],
                "active": abs(float(row["delta_C_memory_forward"])) > params.effect_floor,
                "chi_mean": chi_mean,
                "p_mean": p_mean,
                "g1_mean": g1_mean,
                "chi_mean_error": chi_mean - true_values["chi"],
                "p_mean_error": p_mean - true_values["p"],
                "g1_mean_error": g1_mean - true_values["g1"],
                "chi_minus_p_mean_error": chi_mean - p_mean,
                "chi_minus_g1_mean_error": chi_mean - g1_mean,
                "p_minus_g1_mean_error": p_mean - g1_mean,
                "chi_std": float(np.std(gauge_values["chi"])),
                "p_std": float(np.std(gauge_values["p"])),
                "g1_std": float(np.std(gauge_values["g1"])),
            }
        )

    single_gauge_sign_failure_count = sum(
        1 for row in noisy_rows if not bool(row["single_gauge_sign_match"])
    )
    summary = {
        "mode": mode,
        "level": level,
        "gauge_count": GAUGE_COUNT,
        "row_count": len(eval_rows),
        "active_row_count": sum(1 for row in eval_rows if bool(row["active"])),
        "single_gauge_sign_failure_count": single_gauge_sign_failure_count,
        "max_chi_mean_error_abs": max_abs(eval_rows, "chi_mean_error"),
        "max_p_mean_error_abs": max_abs(eval_rows, "p_mean_error"),
        "max_g1_mean_error_abs": max_abs(eval_rows, "g1_mean_error"),
        "max_chi_minus_p_mean_error_abs": max_abs(eval_rows, "chi_minus_p_mean_error"),
        "max_chi_minus_g1_mean_error_abs": max_abs(eval_rows, "chi_minus_g1_mean_error"),
        "max_p_minus_g1_mean_error_abs": max_abs(eval_rows, "p_minus_g1_mean_error"),
        "max_channel_std": max(
            max_abs(eval_rows, "chi_std"),
            max_abs(eval_rows, "p_std"),
            max_abs(eval_rows, "g1_std"),
        ),
    }
    summary["zero_mean_valid"] = bool(
        mode == "zero_mean"
        and float(summary["max_chi_mean_error_abs"]) <= params.effect_floor
        and float(summary["max_p_mean_error_abs"]) <= params.effect_floor
        and float(summary["max_g1_mean_error_abs"]) <= params.effect_floor
        and float(summary["max_chi_minus_g1_mean_error_abs"]) <= params.effect_floor
    )
    summary["common_bias_detected"] = bool(
        mode == "common_bias"
        and max(
            float(summary["max_chi_minus_p_mean_error_abs"]),
            float(summary["max_chi_minus_g1_mean_error_abs"]),
            float(summary["max_p_minus_g1_mean_error_abs"]),
        )
        > COMMON_BIAS_DETECTION_FLOOR
    )
    return summary, eval_rows


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C Stage III 別読出しノイズ頑健性予備実験 v1",
        "",
        "## 目的",
        "",
        "Stage III 別読出し照合が、無雑音の同一式を見ているだけではないかを検査する。",
        "",
        "同一状態スナップショットから得る位置位相二階差分、`p_read` 差分、G1傾斜読出しに対して、ゼロ平均ゲージ揺らぎと全ゲージ共通バイアスを加える。",
        "",
        "ゼロ平均揺らぎは多ゲージ平均で相殺され、共通バイアスは相殺されず検出される必要がある。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## ノイズケース別サマリー",
            "",
            "| mode | level | gauge | single sign failures | max mean err | max cross err | max std | valid/detected |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        max_mean = max(
            float(row["max_chi_mean_error_abs"]),
            float(row["max_p_mean_error_abs"]),
            float(row["max_g1_mean_error_abs"]),
        )
        max_cross = max(
            float(row["max_chi_minus_p_mean_error_abs"]),
            float(row["max_chi_minus_g1_mean_error_abs"]),
            float(row["max_p_minus_g1_mean_error_abs"]),
        )
        valid_detected = bool(row["zero_mean_valid"]) or bool(row["common_bias_detected"])
        lines.append(
            f"| {row['mode']} | {row['level']:.1e} | {row['gauge_count']} | "
            f"{row['single_gauge_sign_failure_count']} | {max_mean:.16e} | "
            f"{max_cross:.16e} | {row['max_channel_std']:.16e} | `{valid_detected}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- ゼロ平均ゲージ揺らぎでは、個々のゲージ値が大きく揺れても、多ゲージ平均は元の Stage III 照合へ戻る。",
            "- 共通バイアスは多ゲージ平均では消えないため、読出し器由来の系統偏差として検出される。",
            "- これは標準重力の導出ではなく、準安定傾斜候補の Stage III 読出しが単一ゲージ依存ではないことを確認する予備検査である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_stage3_noise_robustness_preliminary_result_v1.json` |",
            "| cases CSV | `abc_baseline_stationary_wave_stage3_noise_robustness_cases_v1.csv` |",
            "| eval CSV | `abc_baseline_stationary_wave_stage3_noise_robustness_eval_rows_v1.csv` |",
            "| noisy gauge CSV | `abc_baseline_stationary_wave_stage3_noise_robustness_noisy_gauge_rows_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_stage3_noise_robustness_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (
        BASE_DIR
        / "ABCベースライン定常波C Stage III別読出しノイズ頑健性予備実験検証メモ_v1.md"
    ).write_text(report, encoding="utf-8")


def run() -> Dict[str, Any]:
    params = Stage2Params()
    true_rows, reference_scale = collect_true_rows(params)
    patterns = zero_mean_patterns()
    case_summaries: List[Dict[str, Any]] = []
    all_eval_rows: List[Dict[str, Any]] = []
    all_noisy_rows: List[Dict[str, Any]] = []

    for level in ZERO_MEAN_LEVELS:
        summary, eval_rows = evaluate_case(
            true_rows, params, "zero_mean", level, reference_scale, patterns
        )
        case_summaries.append(summary)
        all_eval_rows.extend(eval_rows)
        if level in [0.50, 1.00]:
            _, noisy_rows = evaluate_case(
                true_rows, params, "zero_mean", level, reference_scale, patterns
            )
            all_noisy_rows.extend(noisy_rows)

    for level in COMMON_BIAS_LEVELS:
        summary, eval_rows = evaluate_case(
            true_rows, params, "common_bias", level, reference_scale, patterns
        )
        case_summaries.append(summary)
        all_eval_rows.extend(eval_rows)
        _, noisy_rows = evaluate_case(
            true_rows, params, "common_bias", level, reference_scale, patterns
        )
        if level == COMMON_BIAS_LEVELS[-1]:
            all_noisy_rows.extend(noisy_rows)

    zero_mean_summaries = [row for row in case_summaries if row["mode"] == "zero_mean"]
    common_bias_summaries = [row for row in case_summaries if row["mode"] == "common_bias"]
    aggregate_verdict = {
        "case_count": len(case_summaries),
        "zero_mean_level_count": len(zero_mean_summaries),
        "common_bias_level_count": len(common_bias_summaries),
        "gauge_count": GAUGE_COUNT,
        "single_gauge_only_used": False,
        "zero_mean_multigauge_valid_all": bool_all(bool(row["zero_mean_valid"]) for row in zero_mean_summaries),
        "common_bias_detection_floor": COMMON_BIAS_DETECTION_FLOOR,
        "common_bias_detected_all": bool_all(bool(row["common_bias_detected"]) for row in common_bias_summaries),
        "zero_mean_single_gauge_failures_exist": any(
            int(row["single_gauge_sign_failure_count"]) > 0 for row in zero_mean_summaries
        ),
        "zero_mean_max_mean_error_abs": max(
            max(
                float(row["max_chi_mean_error_abs"]),
                float(row["max_p_mean_error_abs"]),
                float(row["max_g1_mean_error_abs"]),
            )
            for row in zero_mean_summaries
        ),
        "common_bias_min_detected_cross_error_abs": min(
            max(
                float(row["max_chi_minus_p_mean_error_abs"]),
                float(row["max_chi_minus_g1_mean_error_abs"]),
                float(row["max_p_minus_g1_mean_error_abs"]),
            )
            for row in common_bias_summaries
        ),
    }
    aggregate_verdict["stage3_noise_robustness_preliminary_valid"] = bool(
        aggregate_verdict["zero_mean_multigauge_valid_all"]
        and aggregate_verdict["common_bias_detected_all"]
        and aggregate_verdict["single_gauge_only_used"] is False
    )

    result = {
        "experiment": "abc_baseline_stationary_wave_stage3_noise_robustness_preliminary_v1",
        "params": asdict(params),
        "noise": {
            "gauge_count": GAUGE_COUNT,
            "zero_mean_levels": ZERO_MEAN_LEVELS,
            "common_bias_levels": COMMON_BIAS_LEVELS,
            "common_bias_detection_floor": COMMON_BIAS_DETECTION_FLOOR,
            "reference_scale": reference_scale,
        },
        "case_summaries": case_summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (
        OUT_DIR
        / "abc_baseline_stationary_wave_stage3_noise_robustness_preliminary_result_v1.json"
    ).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_stage3_noise_robustness_cases_v1.csv", case_summaries)
    write_csv(
        OUT_DIR / "abc_baseline_stationary_wave_stage3_noise_robustness_eval_rows_v1.csv",
        all_eval_rows,
    )
    write_csv(
        OUT_DIR / "abc_baseline_stationary_wave_stage3_noise_robustness_noisy_gauge_rows_v1.csv",
        all_noisy_rows,
    )
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
