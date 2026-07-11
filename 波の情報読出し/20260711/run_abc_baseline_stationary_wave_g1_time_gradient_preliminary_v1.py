from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import (
    Stage2Case,
    Stage2Params,
    default_gauges,
    simulate_case,
    summarize_fanout_rows,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_g1_time_gradient_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def sign(value: float, floor: float = 1.0e-18) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


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


def timeline_for(case: Stage2Case, params: Stage2Params) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    rows, summary = simulate_case(case, params, default_gauges())
    return summarize_fanout_rows(rows, params), summary


def row_map(timeline: List[Dict[str, Any]]) -> Dict[Tuple[int, str], Dict[str, Any]]:
    return {(int(row["step"]), str(row["particle"])): row for row in timeline}


def memory_series(timeline: List[Dict[str, Any]], key: str = "C_memory") -> Dict[int, float]:
    series: Dict[int, float] = {}
    for row in timeline:
        if row["particle"] == "A":
            series[int(row["step"])] = float(row[key])
    return series


def predicted_Ra_from_c_memory(delta_c: float, particle: str, params: Stage2Params) -> float:
    reduced = params.R_A * params.R_B / (params.R_A + params.R_B)
    value = reduced * params.epsilon_C_return * delta_c
    if particle == "B":
        value = -value
    return float(value)


def g1_rows_for_single(
    label: str,
    timeline: List[Dict[str, Any]],
    params: Stage2Params,
    mode: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    by_key = row_map(timeline)
    c_series = memory_series(timeline)
    for step in range(1, params.step_count - 1):
        delta_c = c_series[step + 1] - c_series[step]
        for particle, r_value in [("A", params.R_A), ("B", params.R_B)]:
            state = by_key[(step, particle)]
            actual_Ra = r_value * float(state["a_chi_read"])
            predicted_Ra = predicted_Ra_from_c_memory(delta_c, particle, params)
            partial_chi_E_read = -predicted_Ra
            error = actual_Ra + partial_chi_E_read
            rows.append(
                {
                    "label": label,
                    "mode": mode,
                    "step": step,
                    "particle": particle,
                    "actual_Ra": actual_Ra,
                    "predicted_minus_partial_chi_E": predicted_Ra,
                    "partial_chi_E_read": partial_chi_E_read,
                    "g1_balance_error": error,
                    "C_memory": c_series[step],
                    "delta_C_memory_forward": delta_c,
                    "sign_match": sign(actual_Ra) == sign(predicted_Ra),
                }
            )
    return rows


def g1_rows_for_residual(
    label: str,
    combined: List[Dict[str, Any]],
    r2: List[Dict[str, Any]],
    r3: List[Dict[str, Any]],
    params: Stage2Params,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    combined_map = row_map(combined)
    r2_map = row_map(r2)
    r3_map = row_map(r3)
    c_series = memory_series(combined)
    for step in range(1, params.step_count - 1):
        delta_c = c_series[step + 1] - c_series[step]
        for particle, r_value in [("A", params.R_A), ("B", params.R_B)]:
            actual_Ra = r_value * (
                float(combined_map[(step, particle)]["a_chi_read"])
                - float(r2_map[(step, particle)]["a_chi_read"])
                - float(r3_map[(step, particle)]["a_chi_read"])
            )
            predicted_Ra = predicted_Ra_from_c_memory(delta_c, particle, params)
            partial_chi_E_read = -predicted_Ra
            error = actual_Ra + partial_chi_E_read
            rows.append(
                {
                    "label": label,
                    "mode": "combined_minus_R2_R3_residual",
                    "step": step,
                    "particle": particle,
                    "actual_Ra": actual_Ra,
                    "predicted_minus_partial_chi_E": predicted_Ra,
                    "partial_chi_E_read": partial_chi_E_read,
                    "g1_balance_error": error,
                    "C_memory": c_series[step],
                    "delta_C_memory_forward": delta_c,
                    "sign_match": sign(actual_Ra) == sign(predicted_Ra),
                }
            )
    return rows


def summarize_g1(label: str, mode: str, rows: List[Dict[str, Any]], params: Stage2Params, expected: str) -> Dict[str, Any]:
    max_actual = max(abs(float(row["actual_Ra"])) for row in rows)
    max_predicted = max(abs(float(row["predicted_minus_partial_chi_E"])) for row in rows)
    max_error = max(abs(float(row["g1_balance_error"])) for row in rows)
    sign_count = sum(1 for row in rows if bool(row["sign_match"]))
    sign_ratio = sign_count / len(rows) if rows else 0.0
    first_sign_break_steps = [int(row["step"]) for row in rows if not bool(row["sign_match"])]
    first_sign_break_step = min(first_sign_break_steps) if first_sign_break_steps else None

    if expected == "match":
        valid = bool(max_error <= params.residual_tol and sign_ratio >= 0.98)
    elif expected == "reject":
        valid = bool(max_actual >= params.effect_floor and max_predicted <= params.effect_floor and max_error >= params.effect_floor)
    else:
        valid = False

    return {
        "label": label,
        "mode": mode,
        "expected": expected,
        "row_count": len(rows),
        "max_actual_Ra_abs": max_actual,
        "max_predicted_minus_partial_chi_E_abs": max_predicted,
        "max_g1_balance_error_abs": max_error,
        "sign_match_ratio": sign_ratio,
        "first_sign_break_step": first_sign_break_step,
        "valid": valid,
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C G1 時間位相勾配予備実験 v1",
        "",
        "## 目的",
        "",
        "`G1`: `τ`方向 `E_read` 勾配による `χ`方向補償候補を、C残渣の時間変化から読む。",
        "",
        "ここでは `E_read` 側を位置加速度から定義しない。C残渣の時間変化 `ΔC_memory` から `∂χE_read` proxy を構成し、`R*aχ ≈ -∂χE_read` と整合するかを見る。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## ケース別判定",
            "",
            "| label | mode | expected | max R*a | max -∂χE | max error | sign ratio | valid |",
            "|---|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        break_step = "" if row["first_sign_break_step"] is None else f", break={row['first_sign_break_step']}"
        lines.append(
            f"| {row['label']} | {row['mode']}{break_step} | {row['expected']} | "
            f"{row['max_actual_Ra_abs']:.16e} | {row['max_predicted_minus_partial_chi_E_abs']:.16e} | "
            f"{row['max_g1_balance_error_abs']:.16e} | {row['sign_match_ratio']:.6f} | `{row['valid']}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- C媒介のみのケースで `R*aχ ≈ -∂χE_read` が成立する場合、C残渣の時間変化は位置位相加速度候補と整合する。",
            "- R2/R3 単独で成立しないことは陰性対照である。反力候補を G1 と誤認しないために必要である。",
            "- 合成状態から R2/R3 を差し引いた residual は、full window と early window に分けて判定する。",
            "- full window で符号が崩れ、early window で整合する場合、G1 は定常的な同一視ではなく、準安定遷移窓に限って残る候補として扱う。",
            "- これは標準重力の導出ではなく、`G1` 候補の読出し系列整合性の予備検査である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_g1_time_gradient_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_g1_time_gradient_cases_v1.csv` |",
            "| rows CSV | `abc_baseline_stationary_wave_g1_time_gradient_rows_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_g1_time_gradient_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C G1時間位相勾配予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params()
    cases = {
        "C_mediated_only_persistent": Stage2Case("C_mediated_only_persistent", use_c=True),
        "C_mediated_only_reset": Stage2Case("C_mediated_only_reset", use_c=True, c_persistent=False),
        "R2_only": Stage2Case("R2_only", use_r2=True),
        "R3_only": Stage2Case("R3_only", use_r3=True),
        "C_plus_R2_plus_R3": Stage2Case("C_plus_R2_plus_R3", use_c=True, use_r2=True, use_r3=True),
        "C_plus_R2_plus_R3_mirrored": Stage2Case(
            "C_plus_R2_plus_R3_mirrored", use_c=True, use_r2=True, use_r3=True, mirrored=True
        ),
        "R2_only_mirrored": Stage2Case("R2_only_mirrored", use_r2=True, mirrored=True),
        "R3_only_mirrored": Stage2Case("R3_only_mirrored", use_r3=True, mirrored=True),
    }

    timelines: Dict[str, List[Dict[str, Any]]] = {}
    for label, case in cases.items():
        timelines[label], _ = timeline_for(case, params)

    grouped_rows: List[Tuple[str, str, str, List[Dict[str, Any]]]] = []
    grouped_rows.append(
        (
            "C_mediated_only_persistent",
            "single_C_memory_time_gradient",
            "match",
            g1_rows_for_single("C_mediated_only_persistent", timelines["C_mediated_only_persistent"], params, "single_C_memory_time_gradient"),
        )
    )
    grouped_rows.append(
        (
            "C_mediated_only_reset",
            "single_C_memory_time_gradient",
            "match",
            g1_rows_for_single("C_mediated_only_reset", timelines["C_mediated_only_reset"], params, "single_C_memory_time_gradient"),
        )
    )
    grouped_rows.append(
        (
            "R2_only_control",
            "reaction_control_without_C_memory",
            "reject",
            g1_rows_for_single("R2_only_control", timelines["R2_only"], params, "reaction_control_without_C_memory"),
        )
    )
    grouped_rows.append(
        (
            "R3_only_control",
            "reaction_control_without_C_memory",
            "reject",
            g1_rows_for_single("R3_only_control", timelines["R3_only"], params, "reaction_control_without_C_memory"),
        )
    )
    grouped_rows.append(
        (
            "combined_minus_R2_R3_residual",
            "combined_minus_R2_R3_residual",
            "match",
            g1_rows_for_residual(
                "combined_minus_R2_R3_residual",
                timelines["C_plus_R2_plus_R3"],
                timelines["R2_only"],
                timelines["R3_only"],
                params,
            ),
        )
    )
    grouped_rows.append(
        (
            "combined_minus_R2_R3_residual_early_window",
            "combined_minus_R2_R3_residual_steps_1_to_16",
            "match",
            [
                row
                for row in g1_rows_for_residual(
                    "combined_minus_R2_R3_residual_early_window",
                    timelines["C_plus_R2_plus_R3"],
                    timelines["R2_only"],
                    timelines["R3_only"],
                    params,
                )
                if 1 <= int(row["step"]) <= 16
            ],
        )
    )
    grouped_rows.append(
        (
            "combined_minus_R2_R3_residual_mirrored",
            "combined_minus_R2_R3_residual",
            "match",
            g1_rows_for_residual(
                "combined_minus_R2_R3_residual_mirrored",
                timelines["C_plus_R2_plus_R3_mirrored"],
                timelines["R2_only_mirrored"],
                timelines["R3_only_mirrored"],
                params,
            ),
        )
    )
    grouped_rows.append(
        (
            "combined_minus_R2_R3_residual_mirrored_early_window",
            "combined_minus_R2_R3_residual_steps_1_to_16",
            "match",
            [
                row
                for row in g1_rows_for_residual(
                    "combined_minus_R2_R3_residual_mirrored_early_window",
                    timelines["C_plus_R2_plus_R3_mirrored"],
                    timelines["R2_only_mirrored"],
                    timelines["R3_only_mirrored"],
                    params,
                )
                if 1 <= int(row["step"]) <= 16
            ],
        )
    )

    all_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for label, mode, expected, rows in grouped_rows:
        all_rows.extend(rows)
        summaries.append(summarize_g1(label, mode, rows, params, expected))

    aggregate_verdict = {
        "case_count": len(summaries),
        "all_cases_valid": bool_all(row["valid"] for row in summaries),
        "single_gauge_only_used": False,
        "g1_time_gradient_preliminary_valid": bool_all(row["valid"] for row in summaries),
        "max_match_error": max(
            float(row["max_g1_balance_error_abs"]) for row in summaries if row["expected"] == "match"
        ),
        "min_match_sign_ratio": min(float(row["sign_match_ratio"]) for row in summaries if row["expected"] == "match"),
        "reaction_controls_rejected": bool_all(row["valid"] for row in summaries if row["expected"] == "reject"),
        "early_window_residual_valid": bool_all(
            row["valid"] for row in summaries if "early_window" in str(row["label"])
        ),
        "full_window_residual_valid": bool_all(
            row["valid"]
            for row in summaries
            if row["label"] in {"combined_minus_R2_R3_residual", "combined_minus_R2_R3_residual_mirrored"}
        ),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_g1_time_gradient_preliminary_v1",
        "params": asdict(params),
        "case_summaries": summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_g1_time_gradient_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_g1_time_gradient_cases_v1.csv", summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_g1_time_gradient_rows_v1.csv", all_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
