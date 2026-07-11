from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_g1_time_gradient_preliminary_v1 import (
    g1_rows_for_residual,
    g1_rows_for_single,
    timeline_for,
)
from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import (
    Stage2Case,
    Stage2Params,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_metastable_window_sweep_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class WindowSweepCase:
    name: str
    memory_decay_C: float
    memory_decay_R: float
    epsilon_to_C: float
    epsilon_C_return: float
    reaction_scale: float


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
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def max_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return max(abs(float(row[key])) for row in rows)


def step_validity(rows: List[Dict[str, Any]], params: Stage2Params) -> Dict[int, bool]:
    grouped: Dict[int, List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(int(row["step"]), []).append(row)
    return {
        step: bool_all(
            bool(row["sign_match"]) and abs(float(row["g1_balance_error"])) <= params.residual_tol
            for row in selected
        )
        for step, selected in grouped.items()
    }


def first_step_break(rows: List[Dict[str, Any]], params: Stage2Params) -> Optional[int]:
    validity = step_validity(rows, params)
    for step in sorted(validity):
        if not validity[step]:
            return step
    return None


def first_sign_break(rows: List[Dict[str, Any]]) -> Optional[int]:
    grouped: Dict[int, List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(int(row["step"]), []).append(row)
    for step in sorted(grouped):
        if not bool_all(bool(row["sign_match"]) for row in grouped[step]):
            return step
    return None


def summarize_window(label: str, mode: str, rows: List[Dict[str, Any]], params: Stage2Params) -> Dict[str, Any]:
    steps = sorted({int(row["step"]) for row in rows})
    full_window_length = len(steps)
    break_step = first_step_break(rows, params)
    sign_break_step = first_sign_break(rows)
    metastable_window_length = full_window_length if break_step is None else max(0, break_step - min(steps))
    if break_step is None:
        inside_rows = rows
        outside_rows: List[Dict[str, Any]] = []
    else:
        inside_rows = [row for row in rows if int(row["step"]) < break_step]
        outside_rows = [row for row in rows if int(row["step"]) >= break_step]

    sign_count = sum(1 for row in rows if bool(row["sign_match"]))
    return {
        "label": label,
        "mode": mode,
        "row_count": len(rows),
        "full_window_length": full_window_length,
        "metastable_window_length": metastable_window_length,
        "first_validity_break_step": break_step,
        "first_sign_break_step": sign_break_step,
        "sign_match_ratio": sign_count / len(rows) if rows else 0.0,
        "max_actual_Ra_abs": max_abs(rows, "actual_Ra"),
        "max_predicted_minus_partial_chi_E_abs": max_abs(rows, "predicted_minus_partial_chi_E"),
        "max_g1_balance_error_abs": max_abs(rows, "g1_balance_error"),
        "max_g1_balance_error_inside_window": max_abs(inside_rows, "g1_balance_error"),
        "max_g1_balance_error_outside_window": max_abs(outside_rows, "g1_balance_error"),
        "full_window_valid": break_step is None,
        "classified": bool(rows and full_window_length > 0),
    }


def params_for_case(base: Stage2Params, case: WindowSweepCase) -> Stage2Params:
    return replace(
        base,
        memory_decay_C=case.memory_decay_C,
        memory_decay_R=case.memory_decay_R,
        epsilon_to_C=case.epsilon_to_C,
        epsilon_C_return=case.epsilon_C_return,
        epsilon_r2_pressure=base.epsilon_r2_pressure * case.reaction_scale,
        epsilon_r3_direct=base.epsilon_r3_direct * case.reaction_scale,
    )


def default_sweep_cases(base: Stage2Params) -> List[WindowSweepCase]:
    return [
        WindowSweepCase("base", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 1.0),
        WindowSweepCase("C_decay_0_78", 0.78, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 1.0),
        WindowSweepCase("C_decay_0_86", 0.86, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 1.0),
        WindowSweepCase("C_decay_0_97", 0.97, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 1.0),
        WindowSweepCase("R_decay_0_78", base.memory_decay_C, 0.78, base.epsilon_to_C, base.epsilon_C_return, 1.0),
        WindowSweepCase("R_decay_0_96", base.memory_decay_C, 0.96, base.epsilon_to_C, base.epsilon_C_return, 1.0),
        WindowSweepCase("C_source_half", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C * 0.5, base.epsilon_C_return, 1.0),
        WindowSweepCase("C_source_double", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C * 2.0, base.epsilon_C_return, 1.0),
        WindowSweepCase("C_return_half", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return * 0.5, 1.0),
        WindowSweepCase("C_return_double", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return * 2.0, 1.0),
        WindowSweepCase("reaction_quarter", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 0.25),
        WindowSweepCase("reaction_half", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 0.5),
        WindowSweepCase("reaction_double", base.memory_decay_C, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 2.0),
        WindowSweepCase("slow_C_weak_reaction", 0.97, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 0.25),
        WindowSweepCase("slow_C_strong_reaction", 0.97, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 2.0),
        WindowSweepCase("fast_C_weak_reaction", 0.78, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 0.25),
        WindowSweepCase("fast_C_strong_reaction", 0.78, base.memory_decay_R, base.epsilon_to_C, base.epsilon_C_return, 2.0),
    ]


def timeline_summary(case: Stage2Case, params: Stage2Params) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    return timeline_for(case, params)


def run_case(case: WindowSweepCase, base_params: Stage2Params) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    params = params_for_case(base_params, case)
    c_timeline, c_summary = timeline_summary(Stage2Case("C_mediated_only_persistent", use_c=True), params)
    r2_timeline, r2_summary = timeline_summary(Stage2Case("R2_only", use_r2=True), params)
    r3_timeline, r3_summary = timeline_summary(Stage2Case("R3_only", use_r3=True), params)
    combined_timeline, combined_summary = timeline_summary(
        Stage2Case("C_plus_R2_plus_R3", use_c=True, use_r2=True, use_r3=True), params
    )

    c_rows = g1_rows_for_single(
        f"{case.name}:C_only",
        c_timeline,
        params,
        "single_C_memory_time_gradient",
    )
    residual_rows = g1_rows_for_residual(
        f"{case.name}:combined_minus_R2_R3_residual",
        combined_timeline,
        r2_timeline,
        r3_timeline,
        params,
    )

    c_window = summarize_window(case.name, "C_only", c_rows, params)
    residual_window = summarize_window(case.name, "combined_minus_R2_R3_residual", residual_rows, params)

    c_distance_change = float(c_summary["distance_change"])
    r2_distance_change = float(r2_summary["distance_change"])
    r3_distance_change = float(r3_summary["distance_change"])
    combined_distance_change = float(combined_summary["distance_change"])
    residual_after_reactions = combined_distance_change - r2_distance_change - r3_distance_change
    reaction_distance_sum = r2_distance_change + r3_distance_change

    row = {
        "case": case.name,
        "memory_decay_C": params.memory_decay_C,
        "memory_decay_R": params.memory_decay_R,
        "epsilon_to_C": params.epsilon_to_C,
        "epsilon_C_return": params.epsilon_C_return,
        "epsilon_r2_pressure": params.epsilon_r2_pressure,
        "epsilon_r3_direct": params.epsilon_r3_direct,
        "reaction_scale": case.reaction_scale,
        "C_distance_change": c_distance_change,
        "R2_distance_change": r2_distance_change,
        "R3_distance_change": r3_distance_change,
        "reaction_distance_sum": reaction_distance_sum,
        "combined_distance_change": combined_distance_change,
        "residual_after_R2_R3": residual_after_reactions,
        "reaction_dominance_ratio_abs": abs(reaction_distance_sum / c_distance_change)
        if abs(c_distance_change) > 0.0
        else float("inf"),
        "C_final_memory": float(c_summary["final_C_memory"]),
        "combined_final_C_memory": float(combined_summary["final_C_memory"]),
        "combined_max_Q_raw_abs": float(combined_summary["max_Q_raw_abs"]),
        "combined_max_Q_closed_abs": float(combined_summary["max_Q_closed_abs"]),
        "C_only_window_length": c_window["metastable_window_length"],
        "C_only_full_window_valid": c_window["full_window_valid"],
        "C_only_max_error": c_window["max_g1_balance_error_abs"],
        "residual_window_length": residual_window["metastable_window_length"],
        "residual_full_window_valid": residual_window["full_window_valid"],
        "residual_first_validity_break_step": residual_window["first_validity_break_step"],
        "residual_first_sign_break_step": residual_window["first_sign_break_step"],
        "residual_sign_match_ratio": residual_window["sign_match_ratio"],
        "residual_max_error": residual_window["max_g1_balance_error_abs"],
        "residual_max_error_inside_window": residual_window["max_g1_balance_error_inside_window"],
        "residual_max_error_outside_window": residual_window["max_g1_balance_error_outside_window"],
        "residual_classified": residual_window["classified"],
    }

    for data_row in c_rows:
        data_row["sweep_case"] = case.name
        data_row["window_mode"] = "C_only"
    for data_row in residual_rows:
        data_row["sweep_case"] = case.name
        data_row["window_mode"] = "combined_minus_R2_R3_residual"
    return row, c_rows + residual_rows


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 準安定G1窓スイープ予備実験 v1",
        "",
        "## 目的",
        "",
        "`G1`: `R*aχ ≈ -∂χE_read` 候補が、定常的な恒等式ではなく準安定遷移窓として現れるかを調べる。",
        "",
        "前段の G1 予備実験では、C単独では整合し、R2/R3反力対照では棄却され、合成状態から R2/R3 を差し引いた residual は early window で整合したが full window では破れた。",
        "",
        "本予備実験では、C記憶減衰、反力記憶減衰、C源強度、C戻り強度、反力候補強度を掃引し、G1 residual の準安定窓長を分類する。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## ケース別サマリー",
            "",
            "| case | C decay | R decay | C source | C return | reaction scale | residual window | break | reaction/C | residual after R2/R3 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["case_summaries"]:
        break_step = "" if row["residual_first_validity_break_step"] is None else row["residual_first_validity_break_step"]
        lines.append(
            f"| {row['case']} | {row['memory_decay_C']:.3f} | {row['memory_decay_R']:.3f} | "
            f"{row['epsilon_to_C']:.2e} | {row['epsilon_C_return']:.2e} | {row['reaction_scale']:.3f} | "
            f"{row['residual_window_length']} | {break_step} | {row['reaction_dominance_ratio_abs']:.6e} | "
            f"{row['residual_after_R2_R3']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- `C_only_full_window_valid` が保たれる一方で、`combined_minus_R2_R3_residual` が有限窓で破れる場合、G1 は単純な常時成立式ではなく、相互作用後の準安定遷移で読まれる候補である。",
            "- 反力候補が強い条件で窓が短くなる場合、観測窓の切り方により G1 候補が隠れる。",
            "- C記憶が長い条件で窓が伸びる場合、G1 は C残渣遅延と結びついた候補である。",
            "- 本実験の `valid` は、全ケースで G1 が成立したという意味ではない。準安定窓を多条件で分類できたという意味である。",
            "- これは標準重力の導出ではなく、加速度読出し前の窓長・遅延・反力分離の予備計量である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_metastable_window_sweep_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_metastable_window_sweep_cases_v1.csv` |",
            "| rows CSV | `abc_baseline_stationary_wave_metastable_window_sweep_rows_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_window_sweep_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 準安定G1窓スイープ予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    base_params = Stage2Params()
    cases = default_sweep_cases(base_params)
    summaries: List[Dict[str, Any]] = []
    all_rows: List[Dict[str, Any]] = []
    for case in cases:
        summary, rows = run_case(case, base_params)
        summaries.append(summary)
        all_rows.extend(rows)

    residual_windows = [int(row["residual_window_length"]) for row in summaries]
    c_windows = [int(row["C_only_window_length"]) for row in summaries]
    base_row = next(row for row in summaries if row["case"] == "base")
    finite_break_count = sum(1 for row in summaries if row["residual_first_validity_break_step"] is not None)
    c_only_valid_all = bool_all(bool(row["C_only_full_window_valid"]) for row in summaries)
    all_classified = bool_all(bool(row["residual_classified"]) for row in summaries)

    aggregate_verdict = {
        "case_count": len(summaries),
        "all_cases_classified": all_classified,
        "single_gauge_only_used": False,
        "C_only_full_window_valid_all_cases": c_only_valid_all,
        "residual_cases_with_finite_break": finite_break_count,
        "residual_cases_without_break": len(summaries) - finite_break_count,
        "base_residual_window_length": int(base_row["residual_window_length"]),
        "base_residual_first_validity_break_step": base_row["residual_first_validity_break_step"],
        "min_residual_window_length": min(residual_windows),
        "median_residual_window_length": float(median(residual_windows)),
        "max_residual_window_length": max(residual_windows),
        "min_C_only_window_length": min(c_windows),
        "max_C_only_window_length": max(c_windows),
        "max_C_only_error": max(float(row["C_only_max_error"]) for row in summaries),
        "max_residual_inside_window_error": max(float(row["residual_max_error_inside_window"]) for row in summaries),
        "max_reaction_dominance_ratio_abs": max(float(row["reaction_dominance_ratio_abs"]) for row in summaries),
        "memory_decay_C_window_correlation": safe_corr(
            [float(row["memory_decay_C"]) for row in summaries],
            [float(row["residual_window_length"]) for row in summaries],
        ),
        "memory_decay_R_window_correlation": safe_corr(
            [float(row["memory_decay_R"]) for row in summaries],
            [float(row["residual_window_length"]) for row in summaries],
        ),
        "reaction_scale_window_correlation": safe_corr(
            [float(row["reaction_scale"]) for row in summaries],
            [float(row["residual_window_length"]) for row in summaries],
        ),
        "metastable_window_sweep_preliminary_valid": bool(
            all_classified
            and c_only_valid_all
            and int(base_row["residual_window_length"]) >= 8
            and base_row["residual_first_validity_break_step"] is not None
            and finite_break_count >= 1
        ),
    }

    result = {
        "experiment": "abc_baseline_stationary_wave_metastable_window_sweep_preliminary_v1",
        "base_params": asdict(base_params),
        "case_summaries": summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_window_sweep_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_window_sweep_cases_v1.csv", summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_window_sweep_rows_v1.csv", all_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
