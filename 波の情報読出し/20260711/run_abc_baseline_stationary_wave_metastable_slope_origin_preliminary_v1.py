from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_blind_metastable_window_preliminary_v1 import (
    blind_signal_rows,
    extract_blind_window,
    select_step_rows,
)
from run_abc_baseline_stationary_wave_g1_time_gradient_preliminary_v1 import (
    g1_rows_for_residual,
    g1_rows_for_single,
)
from run_abc_baseline_stationary_wave_metastable_window_sweep_preliminary_v1 import (
    WindowSweepCase,
    default_sweep_cases,
    params_for_case,
    timeline_summary,
)
from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import (
    Stage2Case,
    Stage2Params,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_metastable_slope_origin_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


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


def safe_corr(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) < 3:
        return None
    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)

    def nearly_constant(values: np.ndarray) -> bool:
        scale = max(1.0e-30, float(np.max(np.abs(values))))
        return bool(float(np.max(values) - np.min(values)) <= scale * 1.0e-9)

    if nearly_constant(x) or nearly_constant(y):
        return None
    return float(np.corrcoef(x, y)[0, 1])


def max_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return max(abs(float(row[key])) for row in rows)


def mean_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return float(np.mean([abs(float(row[key])) for row in rows]))


def profile_rows(label: str, profile: List[float], params: Stage2Params) -> List[Dict[str, Any]]:
    total = params.R_A + params.R_B
    reduced = params.R_A * params.R_B / total
    cumulative_A = 0.0
    deltas_A: List[float] = []
    for c_memory in profile:
        response = params.epsilon_C_return * c_memory
        cumulative_A += response * params.R_B / total
        deltas_A.append(cumulative_A)

    rows: List[Dict[str, Any]] = []
    for step in range(1, len(profile) - 1):
        actual_Ra = params.R_A * (deltas_A[step + 1] - 2.0 * deltas_A[step] + deltas_A[step - 1])
        delta_c_forward = profile[step + 1] - profile[step]
        predicted_Ra = reduced * params.epsilon_C_return * delta_c_forward
        rows.append(
            {
                "label": label,
                "step": step,
                "C_memory": profile[step],
                "delta_C_memory_forward": delta_c_forward,
                "actual_Ra": actual_Ra,
                "predicted_Ra_from_delta_C": predicted_Ra,
                "level_only_predictor": profile[step],
                "slope_predictor": delta_c_forward,
                "balance_error": actual_Ra - predicted_Ra,
            }
        )
    return rows


def controlled_profiles(params: Stage2Params) -> Dict[str, List[float]]:
    n = params.step_count
    return {
        "constant_nonzero_C": [2.5e-6 for _ in range(n)],
        "linear_increase_C": [2.0e-7 + 5.0e-8 * step for step in range(n)],
        "linear_decrease_C": [2.8e-6 - 5.0e-8 * step for step in range(n)],
        "exponential_approach_C": [3.8e-6 * (1.0 - 0.92 ** (step + 1)) for step in range(n)],
        "overshoot_relax_C": [2.4e-6 + 1.2e-6 * ((-0.82) ** step) for step in range(n)],
    }


def summarize_profile(label: str, rows: List[Dict[str, Any]], params: Stage2Params) -> Dict[str, Any]:
    actual_values = [float(row["actual_Ra"]) for row in rows]
    slope_values = [float(row["slope_predictor"]) for row in rows]
    level_values = [float(row["level_only_predictor"]) for row in rows]
    positive_count = sum(1 for value in actual_values if value > params.effect_floor)
    negative_count = sum(1 for value in actual_values if value < -params.effect_floor)
    near_zero_count = sum(1 for value in actual_values if abs(value) <= params.effect_floor)
    return {
        "label": label,
        "row_count": len(rows),
        "max_actual_Ra_abs": max_abs(rows, "actual_Ra"),
        "max_delta_C_abs": max_abs(rows, "delta_C_memory_forward"),
        "max_C_memory_abs": max_abs(rows, "C_memory"),
        "max_balance_error_abs": max_abs(rows, "balance_error"),
        "actual_slope_corr": safe_corr(actual_values, slope_values),
        "actual_level_corr": safe_corr(actual_values, level_values),
        "positive_actual_count": positive_count,
        "negative_actual_count": negative_count,
        "near_zero_actual_count": near_zero_count,
    }


def c_only_aging_rows(params: Stage2Params) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    timeline, _ = timeline_summary(Stage2Case("C_mediated_only_persistent", use_c=True), params)
    g1_rows = [row for row in g1_rows_for_single("C_only_aging", timeline, params, "single_C_memory_time_gradient") if row["particle"] == "A"]
    timeline_by_step = {int(row["step"]): row for row in timeline if row["particle"] == "A"}
    starts = [1, 8, 16, 24, 32, 40]
    window_rows: List[Dict[str, Any]] = []
    for start in starts:
        end = min(start + 5, params.step_count - 2)
        selected = [row for row in g1_rows if start <= int(row["step"]) <= end]
        if not selected:
            continue
        c_values = [abs(float(timeline_by_step[int(row["step"])]["C_memory"])) for row in selected]
        q_values = [abs(float(timeline_by_step[int(row["step"])]["Q_raw"])) for row in selected]
        window_rows.append(
            {
                "window_start": start,
                "window_end": end,
                "mean_abs_actual_Ra": mean_abs(selected, "actual_Ra"),
                "mean_abs_C_memory": float(np.mean(c_values)),
                "mean_abs_delta_C_memory": mean_abs(selected, "delta_C_memory_forward"),
                "mean_abs_Q_raw": float(np.mean(q_values)),
                "max_g1_balance_error_abs": max_abs(selected, "g1_balance_error"),
            }
        )
    return g1_rows, window_rows


def summarize_c_only_aging(rows: List[Dict[str, Any]], window_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    first = window_rows[0]
    last = window_rows[-1]
    actual_abs = [abs(float(row["actual_Ra"])) for row in rows]
    c_abs = [abs(float(row["C_memory"])) for row in rows]
    dc_abs = [abs(float(row["delta_C_memory_forward"])) for row in rows]
    return {
        "window_count": len(window_rows),
        "C_level_growth_ratio_last_over_first": float(last["mean_abs_C_memory"] / first["mean_abs_C_memory"]),
        "actual_Ra_decay_ratio_last_over_first": float(last["mean_abs_actual_Ra"] / first["mean_abs_actual_Ra"]),
        "delta_C_decay_ratio_last_over_first": float(last["mean_abs_delta_C_memory"] / first["mean_abs_delta_C_memory"]),
        "actual_abs_vs_delta_C_abs_corr": safe_corr(actual_abs, dc_abs),
        "actual_abs_vs_C_level_abs_corr": safe_corr(actual_abs, c_abs),
        "max_window_g1_balance_error_abs": max(float(row["max_g1_balance_error_abs"]) for row in window_rows),
    }


def residual_blind_rows_for_case(case: WindowSweepCase, base_params: Stage2Params) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    params = params_for_case(base_params, case)
    r2_timeline, _ = timeline_summary(Stage2Case("R2_only", use_r2=True), params)
    r3_timeline, _ = timeline_summary(Stage2Case("R3_only", use_r3=True), params)
    combined_timeline, _ = timeline_summary(
        Stage2Case("C_plus_R2_plus_R3", use_c=True, use_r2=True, use_r3=True), params
    )
    window = extract_blind_window(blind_signal_rows(case.name, combined_timeline, params))
    residual_rows = [
        row
        for row in select_step_rows(
            g1_rows_for_residual(case.name, combined_timeline, r2_timeline, r3_timeline, params),
            window["blind_window_start"],
            window["blind_window_end"],
        )
        if row["particle"] == "A"
    ]
    actual_abs = [abs(float(row["actual_Ra"])) for row in residual_rows]
    dc_abs = [abs(float(row["delta_C_memory_forward"])) for row in residual_rows]
    c_abs = [abs(float(row["C_memory"])) for row in residual_rows]
    summary = {
        "case": case.name,
        "blind_window_length": int(window["blind_window_length"]),
        "actual_abs_vs_delta_C_abs_corr": safe_corr(actual_abs, dc_abs),
        "actual_abs_vs_C_level_abs_corr": safe_corr(actual_abs, c_abs),
        "mean_abs_actual_Ra": float(np.mean(actual_abs)) if actual_abs else 0.0,
        "mean_abs_delta_C_memory": float(np.mean(dc_abs)) if dc_abs else 0.0,
        "mean_abs_C_memory": float(np.mean(c_abs)) if c_abs else 0.0,
        "max_g1_balance_error_abs": max_abs(residual_rows, "g1_balance_error"),
    }
    for row in residual_rows:
        row["sweep_case"] = case.name
        row["blind_window_length"] = int(window["blind_window_length"])
    return summary, residual_rows


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 準安定傾斜起源予備実験 v1",
        "",
        "## 目的",
        "",
        "引力的に見える加速度候補が、定常状態の場そのものから出ているのか、閉鎖へ戻る途中の準安定傾斜から出ているのかを検査する。",
        "",
        "本実験では、`C_memory` の定常レベルと `ΔC_memory` の傾斜を分ける。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## 制御Cプロファイル",
            "",
            "| profile | max |R*a| | max |ΔC| | corr(actual, ΔC) | corr(actual, C) | max error |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["controlled_profile_summaries"]:
        lines.append(
            f"| {row['label']} | {row['max_actual_Ra_abs']:.16e} | {row['max_delta_C_abs']:.16e} | "
            f"{row['actual_slope_corr']} | {row['actual_level_corr']} | {row['max_balance_error_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## C単独 aging",
            "",
            "| window | mean |R*a| | mean |C| | mean |ΔC| | max error |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in result["c_only_aging_windows"]:
        lines.append(
            f"| {row['window_start']}-{row['window_end']} | {row['mean_abs_actual_Ra']:.16e} | "
            f"{row['mean_abs_C_memory']:.16e} | {row['mean_abs_delta_C_memory']:.16e} | "
            f"{row['max_g1_balance_error_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## residual blind window",
            "",
            "| case | blind window | corr(|R*a|, |ΔC|) | corr(|R*a|, |C|) | mean |R*a| | max error |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["residual_blind_summaries"]:
        lines.append(
            f"| {row['case']} | {row['blind_window_length']} | {row['actual_abs_vs_delta_C_abs_corr']} | "
            f"{row['actual_abs_vs_C_level_abs_corr']} | {row['mean_abs_actual_Ra']:.16e} | "
            f"{row['max_g1_balance_error_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 非ゼロの一定 `C_memory` が加速度を生まない場合、加速度候補は定常レベルそのものではなく、応答の変化から出ている。",
            "- `C_memory` が増えているにもかかわらず `|R*a|` が減る場合、定常場レベル起源では説明しにくい。",
            "- `|R*a|` が `|ΔC_memory|` と強く対応する場合、力候補は準安定な位相決済傾斜に由来する。",
            "- これは標準重力の導出ではなく、引力的読出し候補が安定場か準安定遷移かを分けるための予備試験である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_metastable_slope_origin_preliminary_result_v1.json` |",
            "| controlled CSV | `abc_baseline_stationary_wave_metastable_slope_origin_controlled_profiles_v1.csv` |",
            "| aging CSV | `abc_baseline_stationary_wave_metastable_slope_origin_c_only_aging_v1.csv` |",
            "| residual CSV | `abc_baseline_stationary_wave_metastable_slope_origin_residual_blind_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_slope_origin_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 準安定傾斜起源予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params()
    controlled_all_rows: List[Dict[str, Any]] = []
    controlled_summaries: List[Dict[str, Any]] = []
    for label, profile in controlled_profiles(params).items():
        rows = profile_rows(label, profile, params)
        controlled_all_rows.extend(rows)
        controlled_summaries.append(summarize_profile(label, rows, params))

    c_only_rows, aging_windows = c_only_aging_rows(params)
    aging_summary = summarize_c_only_aging(c_only_rows, aging_windows)

    residual_summaries: List[Dict[str, Any]] = []
    residual_rows: List[Dict[str, Any]] = []
    for case in default_sweep_cases(params):
        summary, rows = residual_blind_rows_for_case(case, params)
        residual_summaries.append(summary)
        residual_rows.extend(rows)

    controlled_by_label = {row["label"]: row for row in controlled_summaries}
    residual_corrs = [
        abs(float(row["actual_abs_vs_delta_C_abs_corr"]))
        for row in residual_summaries
        if row["actual_abs_vs_delta_C_abs_corr"] is not None
    ]
    aggregate_verdict = {
        "controlled_profile_count": len(controlled_summaries),
        "residual_case_count": len(residual_summaries),
        "single_gauge_only_used": False,
        "constant_nonzero_C_acceleration_near_zero": bool(
            float(controlled_by_label["constant_nonzero_C"]["max_actual_Ra_abs"]) <= params.effect_floor
        ),
        "linear_increase_C_positive_acceleration": bool(
            int(controlled_by_label["linear_increase_C"]["positive_actual_count"])
            == int(controlled_by_label["linear_increase_C"]["row_count"])
        ),
        "linear_decrease_C_negative_acceleration": bool(
            int(controlled_by_label["linear_decrease_C"]["negative_actual_count"])
            == int(controlled_by_label["linear_decrease_C"]["row_count"])
        ),
        "controlled_max_balance_error": max(float(row["max_balance_error_abs"]) for row in controlled_summaries),
        "C_only_level_growth_ratio_last_over_first": aging_summary["C_level_growth_ratio_last_over_first"],
        "C_only_actual_decay_ratio_last_over_first": aging_summary["actual_Ra_decay_ratio_last_over_first"],
        "C_only_delta_C_decay_ratio_last_over_first": aging_summary["delta_C_decay_ratio_last_over_first"],
        "C_only_actual_tracks_delta_C": bool(abs(float(aging_summary["actual_abs_vs_delta_C_abs_corr"])) >= 0.999),
        "C_only_actual_anti_tracks_C_level": bool(float(aging_summary["actual_abs_vs_C_level_abs_corr"]) <= -0.99),
        "min_residual_blind_actual_delta_C_abs_corr": min(residual_corrs),
        "residual_blind_actual_tracks_delta_C_all_cases": bool(min(residual_corrs) >= 0.998),
        "metastable_slope_origin_preliminary_valid": bool(
            float(controlled_by_label["constant_nonzero_C"]["max_actual_Ra_abs"]) <= params.effect_floor
            and int(controlled_by_label["linear_increase_C"]["positive_actual_count"])
            == int(controlled_by_label["linear_increase_C"]["row_count"])
            and int(controlled_by_label["linear_decrease_C"]["negative_actual_count"])
            == int(controlled_by_label["linear_decrease_C"]["row_count"])
            and aging_summary["C_level_growth_ratio_last_over_first"] > 2.0
            and aging_summary["actual_Ra_decay_ratio_last_over_first"] < 0.1
            and aging_summary["delta_C_decay_ratio_last_over_first"] < 0.1
            and abs(float(aging_summary["actual_abs_vs_delta_C_abs_corr"])) >= 0.999
            and float(aging_summary["actual_abs_vs_C_level_abs_corr"]) <= -0.99
            and min(residual_corrs) >= 0.998
        ),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_metastable_slope_origin_preliminary_v1",
        "params": asdict(params),
        "controlled_profile_summaries": controlled_summaries,
        "c_only_aging_summary": aging_summary,
        "c_only_aging_windows": aging_windows,
        "residual_blind_summaries": residual_summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_slope_origin_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_slope_origin_controlled_profiles_v1.csv", controlled_all_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_slope_origin_c_only_aging_v1.csv", aging_windows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_slope_origin_residual_blind_v1.csv", residual_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
