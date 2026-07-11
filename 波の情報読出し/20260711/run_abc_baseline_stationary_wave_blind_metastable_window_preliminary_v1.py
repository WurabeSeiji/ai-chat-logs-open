from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from run_abc_baseline_stationary_wave_g1_time_gradient_preliminary_v1 import (
    g1_rows_for_residual,
    g1_rows_for_single,
)
from run_abc_baseline_stationary_wave_metastable_window_sweep_preliminary_v1 import (
    WindowSweepCase,
    default_sweep_cases,
    params_for_case,
    summarize_window,
    timeline_summary,
)
from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import (
    Stage2Case,
    Stage2Params,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_blind_metastable_window_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

BLIND_SIGNAL_RATIO_FLOOR = 1.0e-5


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


def max_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return max(abs(float(row[key])) for row in rows)


def rows_for_particle(timeline: List[Dict[str, Any]], particle: str = "A") -> List[Dict[str, Any]]:
    return sorted([row for row in timeline if row["particle"] == particle], key=lambda row: int(row["step"]))


def blind_signal_rows(
    label: str,
    combined_timeline: List[Dict[str, Any]],
    params: Stage2Params,
) -> List[Dict[str, Any]]:
    rows = rows_for_particle(combined_timeline, "A")
    by_step = {int(row["step"]): row for row in rows}
    signal_rows: List[Dict[str, Any]] = []
    for step in range(1, params.step_count - 1):
        current = by_step[step]
        nxt = by_step[step + 1]
        delta_c = float(nxt["C_memory"]) - float(current["C_memory"])
        reaction_memory = abs(float(current["R2_memory"])) + abs(float(current["R3_memory"]))
        c_effect_increment = abs(delta_c) * params.epsilon_C_return
        blind_signal_ratio = c_effect_increment / (reaction_memory + 1.0e-30)
        signal_rows.append(
            {
                "case": label,
                "step": step,
                "C_memory": float(current["C_memory"]),
                "delta_C_memory_forward": delta_c,
                "R2_memory": float(current["R2_memory"]),
                "R3_memory": float(current["R3_memory"]),
                "reaction_memory_abs_sum": reaction_memory,
                "c_effect_increment_abs": c_effect_increment,
                "blind_signal_ratio": blind_signal_ratio,
                "Q_raw": float(current["Q_raw"]),
                "Q_closed": float(current["Q_closed"]),
            }
        )
    return signal_rows


def extract_blind_window(signal_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not signal_rows:
        return {
            "blind_window_start": None,
            "blind_window_end": None,
            "blind_window_length": 0,
            "blind_break_step": None,
            "blind_signal_ratio_floor": BLIND_SIGNAL_RATIO_FLOOR,
        }

    start = min(int(row["step"]) for row in signal_rows)
    max_step = max(int(row["step"]) for row in signal_rows)
    end = max_step
    break_step: Optional[int] = None
    for row in sorted(signal_rows, key=lambda item: int(item["step"])):
        if float(row["blind_signal_ratio"]) < BLIND_SIGNAL_RATIO_FLOOR:
            break_step = int(row["step"])
            end = break_step - 1
            break

    if end < start:
        end = start - 1
    return {
        "blind_window_start": start,
        "blind_window_end": end if end >= start else None,
        "blind_window_length": max(0, end - start + 1),
        "blind_break_step": break_step,
        "blind_signal_ratio_floor": BLIND_SIGNAL_RATIO_FLOOR,
    }


def select_step_rows(rows: List[Dict[str, Any]], start: Optional[int], end: Optional[int]) -> List[Dict[str, Any]]:
    if start is None or end is None or end < start:
        return []
    return [row for row in rows if start <= int(row["step"]) <= end]


def evaluate_match_window(label: str, mode: str, rows: List[Dict[str, Any]], params: Stage2Params) -> Dict[str, Any]:
    sign_count = sum(1 for row in rows if bool(row["sign_match"]))
    max_error = max_abs(rows, "g1_balance_error")
    return {
        "label": label,
        "mode": mode,
        "row_count": len(rows),
        "step_count": len({int(row["step"]) for row in rows}),
        "sign_match_ratio": sign_count / len(rows) if rows else 0.0,
        "max_actual_Ra_abs": max_abs(rows, "actual_Ra"),
        "max_predicted_minus_partial_chi_E_abs": max_abs(rows, "predicted_minus_partial_chi_E"),
        "max_g1_balance_error_abs": max_error,
        "valid_match": bool(rows and sign_count == len(rows) and max_error <= params.residual_tol),
    }


def evaluate_reject_window(label: str, mode: str, rows: List[Dict[str, Any]], params: Stage2Params) -> Dict[str, Any]:
    max_actual = max_abs(rows, "actual_Ra")
    max_predicted = max_abs(rows, "predicted_minus_partial_chi_E")
    max_error = max_abs(rows, "g1_balance_error")
    return {
        "label": label,
        "mode": mode,
        "row_count": len(rows),
        "step_count": len({int(row["step"]) for row in rows}),
        "max_actual_Ra_abs": max_actual,
        "max_predicted_minus_partial_chi_E_abs": max_predicted,
        "max_g1_balance_error_abs": max_error,
        "valid_reject": bool(
            rows
            and max_actual >= params.effect_floor
            and max_predicted <= params.effect_floor
            and max_error >= params.effect_floor
        ),
    }


def shifted_window(window: Dict[str, Any], max_step: int) -> Dict[str, Any]:
    start = window["blind_window_start"]
    end = window["blind_window_end"]
    length = int(window["blind_window_length"])
    if start is None or end is None or length <= 0 or end >= max_step:
        return {
            "shifted_window_start": None,
            "shifted_window_end": None,
            "shifted_window_length": 0,
            "shifted_window_applicable": False,
        }
    shifted_start = end + 1
    shifted_end = min(max_step, shifted_start + length - 1)
    return {
        "shifted_window_start": shifted_start,
        "shifted_window_end": shifted_end,
        "shifted_window_length": max(0, shifted_end - shifted_start + 1),
        "shifted_window_applicable": shifted_end >= shifted_start,
    }


def late_window(window: Dict[str, Any], max_step: int) -> Dict[str, Any]:
    length = int(window["blind_window_length"])
    if length <= 0 or length >= max_step:
        return {
            "late_window_start": None,
            "late_window_end": None,
            "late_window_length": 0,
            "late_window_applicable": False,
        }
    late_start = max(1, max_step - length + 1)
    late_end = max_step
    same_as_blind = late_start == window["blind_window_start"] and late_end == window["blind_window_end"]
    return {
        "late_window_start": late_start,
        "late_window_end": late_end,
        "late_window_length": max(0, late_end - late_start + 1),
        "late_window_applicable": not same_as_blind and late_end >= late_start,
    }


def control_status(applicable: bool, rejected_or_worse: bool) -> str:
    if not applicable:
        return "not_applicable"
    if rejected_or_worse:
        return "rejected_or_worse"
    return "not_rejected"


def run_case(case: WindowSweepCase, base_params: Stage2Params) -> Tuple[Dict[str, Any], List[Dict[str, Any]], List[Dict[str, Any]]]:
    params = params_for_case(base_params, case)
    c_timeline, _ = timeline_summary(Stage2Case("C_mediated_only_persistent", use_c=True), params)
    r2_timeline, _ = timeline_summary(Stage2Case("R2_only", use_r2=True), params)
    r3_timeline, _ = timeline_summary(Stage2Case("R3_only", use_r3=True), params)
    combined_timeline, _ = timeline_summary(
        Stage2Case("C_plus_R2_plus_R3", use_c=True, use_r2=True, use_r3=True), params
    )

    signal_rows = blind_signal_rows(case.name, combined_timeline, params)
    window = extract_blind_window(signal_rows)
    max_step = params.step_count - 2
    shifted = shifted_window(window, max_step)
    late = late_window(window, max_step)

    residual_rows = g1_rows_for_residual(
        f"{case.name}:combined_minus_R2_R3_residual",
        combined_timeline,
        r2_timeline,
        r3_timeline,
        params,
    )
    c_rows = g1_rows_for_single(f"{case.name}:C_only", c_timeline, params, "single_C_memory_time_gradient")
    r2_rows = g1_rows_for_single(f"{case.name}:R2_only_control", r2_timeline, params, "reaction_control_without_C_memory")
    r3_rows = g1_rows_for_single(f"{case.name}:R3_only_control", r3_timeline, params, "reaction_control_without_C_memory")

    true_window = summarize_window(case.name, "combined_minus_R2_R3_residual", residual_rows, params)
    blind_rows = select_step_rows(residual_rows, window["blind_window_start"], window["blind_window_end"])
    shifted_rows = select_step_rows(residual_rows, shifted["shifted_window_start"], shifted["shifted_window_end"])
    late_rows = select_step_rows(residual_rows, late["late_window_start"], late["late_window_end"])
    c_blind_rows = select_step_rows(c_rows, window["blind_window_start"], window["blind_window_end"])
    r2_blind_rows = select_step_rows(r2_rows, window["blind_window_start"], window["blind_window_end"])
    r3_blind_rows = select_step_rows(r3_rows, window["blind_window_start"], window["blind_window_end"])

    blind_eval = evaluate_match_window(case.name, "blind_residual_window", blind_rows, params)
    shifted_eval = evaluate_match_window(case.name, "shifted_residual_window", shifted_rows, params)
    late_eval = evaluate_match_window(case.name, "late_residual_window", late_rows, params)
    c_eval = evaluate_match_window(case.name, "C_only_blind_window", c_blind_rows, params)
    r2_reject = evaluate_reject_window(case.name, "R2_only_blind_window", r2_blind_rows, params)
    r3_reject = evaluate_reject_window(case.name, "R3_only_blind_window", r3_blind_rows, params)

    shifted_rejected = (
        not shifted["shifted_window_applicable"]
        or not shifted_rows
        or not bool(shifted_eval["valid_match"])
        or float(shifted_eval["max_g1_balance_error_abs"]) > float(blind_eval["max_g1_balance_error_abs"])
    )
    late_rejected = (
        not late["late_window_applicable"]
        or not late_rows
        or not bool(late_eval["valid_match"])
        or float(late_eval["max_g1_balance_error_abs"]) > float(blind_eval["max_g1_balance_error_abs"])
    )

    case_summary = {
        "case": case.name,
        "memory_decay_C": params.memory_decay_C,
        "memory_decay_R": params.memory_decay_R,
        "epsilon_to_C": params.epsilon_to_C,
        "epsilon_C_return": params.epsilon_C_return,
        "reaction_scale": case.reaction_scale,
        **window,
        **shifted,
        **late,
        "true_residual_window_length": true_window["metastable_window_length"],
        "true_first_validity_break_step": true_window["first_validity_break_step"],
        "blind_window_inside_true_window": bool(
            int(window["blind_window_length"]) > 0
            and int(window["blind_window_length"]) <= int(true_window["metastable_window_length"])
        ),
        "blind_residual_valid": bool(blind_eval["valid_match"]),
        "blind_residual_max_error": blind_eval["max_g1_balance_error_abs"],
        "blind_residual_sign_match_ratio": blind_eval["sign_match_ratio"],
        "C_only_blind_window_valid": bool(c_eval["valid_match"]),
        "R2_blind_window_rejected": bool(r2_reject["valid_reject"]),
        "R3_blind_window_rejected": bool(r3_reject["valid_reject"]),
        "shifted_window_rejected_or_worse": bool(shifted_rejected),
        "shifted_window_control_status": control_status(
            bool(shifted["shifted_window_applicable"]), bool(shifted_rejected)
        ),
        "shifted_window_valid_match": bool(shifted_eval["valid_match"]),
        "shifted_window_max_error": shifted_eval["max_g1_balance_error_abs"],
        "late_window_rejected_or_worse": bool(late_rejected),
        "late_window_control_status": control_status(bool(late["late_window_applicable"]), bool(late_rejected)),
        "late_window_valid_match": bool(late_eval["valid_match"]),
        "late_window_max_error": late_eval["max_g1_balance_error_abs"],
    }

    eval_rows = [
        {**blind_eval, "case": case.name},
        {**shifted_eval, "case": case.name},
        {**late_eval, "case": case.name},
        {**c_eval, "case": case.name},
        {**r2_reject, "case": case.name},
        {**r3_reject, "case": case.name},
    ]
    return case_summary, signal_rows, eval_rows


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 準安定窓ブラインド抽出予備実験 v1",
        "",
        "## 目的",
        "",
        "`G1`: `R*aχ ≈ -∂χE_read` 候補の評価窓を、G1誤差や加速度を見ずに先に決める。",
        "",
        "窓抽出に使うのは、同一スナップショットから読まれる `C_memory`, `ΔC_memory`, `R2_memory`, `R3_memory`, `Q_raw` だけである。",
        "",
        "本実験では、C効果増分が反力記憶に対して十分大きい期間を blind window とし、その後で初めて G1 residual を評価する。",
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
            "| case | blind window | true window | inside | blind valid | shifted control | late control | R2/R3 rejected | max error |",
            "|---|---:|---:|---|---|---|---|---|---:|",
        ]
    )
    for row in result["case_summaries"]:
        r_reject = bool(row["R2_blind_window_rejected"]) and bool(row["R3_blind_window_rejected"])
        lines.append(
            f"| {row['case']} | {row['blind_window_length']} | {row['true_residual_window_length']} | "
            f"`{row['blind_window_inside_true_window']}` | `{row['blind_residual_valid']}` | "
            f"`{row['shifted_window_control_status']}` | `{row['late_window_control_status']}` | "
            f"`{r_reject}` | {row['blind_residual_max_error']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- blind window は、G1 の当たり外れを見ずに、C残渣遅延と反力記憶の比から決めた窓である。",
            "- blind window 内で G1 residual が成立し、shifted/late window が悪化するなら、窓選択は後付けではなく、C残渣遅延から先に指定できる。",
            "- R2/R3 単独が同じ blind window で棄却されるなら、反力候補を G1 と誤認していない。",
            "- 本実験の主張は、G1 が全期間成立することではない。G1 が、観測前に指定可能な準安定窓で読めることである。",
            "- これは標準重力の導出ではなく、加速度読出しへ進む前の窓選択バイアス検査である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_blind_metastable_window_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_blind_metastable_window_cases_v1.csv` |",
            "| signal CSV | `abc_baseline_stationary_wave_blind_metastable_window_signal_rows_v1.csv` |",
            "| eval CSV | `abc_baseline_stationary_wave_blind_metastable_window_eval_rows_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_blind_metastable_window_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 準安定窓ブラインド抽出予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    base_params = Stage2Params()
    cases = default_sweep_cases(base_params)
    summaries: List[Dict[str, Any]] = []
    signal_rows: List[Dict[str, Any]] = []
    eval_rows: List[Dict[str, Any]] = []
    for case in cases:
        summary, signals, evals = run_case(case, base_params)
        summaries.append(summary)
        signal_rows.extend(signals)
        eval_rows.extend(evals)

    shifted_applicable = [row for row in summaries if bool(row["shifted_window_applicable"])]
    late_applicable = [row for row in summaries if bool(row["late_window_applicable"])]
    aggregate_verdict = {
        "case_count": len(summaries),
        "single_gauge_only_used": False,
        "blind_signal_ratio_floor": BLIND_SIGNAL_RATIO_FLOOR,
        "all_blind_windows_nonempty": bool_all(int(row["blind_window_length"]) > 0 for row in summaries),
        "all_blind_windows_inside_true_windows": bool_all(row["blind_window_inside_true_window"] for row in summaries),
        "all_blind_residual_windows_valid": bool_all(row["blind_residual_valid"] for row in summaries),
        "all_C_only_blind_windows_valid": bool_all(row["C_only_blind_window_valid"] for row in summaries),
        "all_R2_R3_controls_rejected": bool_all(
            bool(row["R2_blind_window_rejected"]) and bool(row["R3_blind_window_rejected"]) for row in summaries
        ),
        "shifted_window_applicable_count": len(shifted_applicable),
        "shifted_windows_rejected_or_worse_all_applicable": bool_all(
            row["shifted_window_rejected_or_worse"] for row in shifted_applicable
        ),
        "late_window_applicable_count": len(late_applicable),
        "late_windows_rejected_or_worse_all_applicable": bool_all(
            row["late_window_rejected_or_worse"] for row in late_applicable
        ),
        "min_blind_window_length": min(int(row["blind_window_length"]) for row in summaries),
        "max_blind_window_length": max(int(row["blind_window_length"]) for row in summaries),
        "max_blind_residual_error": max(float(row["blind_residual_max_error"]) for row in summaries),
        "blind_metastable_window_preliminary_valid": bool(
            bool_all(int(row["blind_window_length"]) > 0 for row in summaries)
            and bool_all(row["blind_window_inside_true_window"] for row in summaries)
            and bool_all(row["blind_residual_valid"] for row in summaries)
            and bool_all(row["C_only_blind_window_valid"] for row in summaries)
            and bool_all(
                bool(row["R2_blind_window_rejected"]) and bool(row["R3_blind_window_rejected"]) for row in summaries
            )
            and bool_all(row["shifted_window_rejected_or_worse"] for row in shifted_applicable)
            and bool_all(row["late_window_rejected_or_worse"] for row in late_applicable)
        ),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_blind_metastable_window_preliminary_v1",
        "base_params": asdict(base_params),
        "case_summaries": summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_blind_metastable_window_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_blind_metastable_window_cases_v1.csv", summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_blind_metastable_window_signal_rows_v1.csv", signal_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_blind_metastable_window_eval_rows_v1.csv", eval_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
