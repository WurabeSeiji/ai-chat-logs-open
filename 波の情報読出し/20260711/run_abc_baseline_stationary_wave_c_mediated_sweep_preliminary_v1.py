from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_abc_baseline_stationary_wave_c_mediated_response_preliminary_v1 import (
    Params,
    add_derivatives,
    default_gauges,
    phase_distance,
    simulate_phase3_persistent,
    summarize_gauge_rows,
    summarize_phase,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_c_mediated_sweep_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

RATIO_ABS_TOL = 1.0e-9


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


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def run_case(params: Params, sweep: str, case: str, value: float) -> Dict[str, Any]:
    gauges = default_gauges()
    gauge_rows, c_rows = simulate_phase3_persistent(params, gauges)
    timeline_rows = summarize_gauge_rows(gauge_rows)
    add_derivatives(timeline_rows, params)
    summary = summarize_phase(
        timeline_rows,
        params,
        "phase3_c_return_persistent",
        "C_return_with_persistent_reembedding",
        f"{sweep}={value}",
    )

    final_a = next(
        row
        for row in sorted(timeline_rows, key=lambda item: (int(item["step"]), str(item["particle"])))
        if row["particle"] == "A" and int(row["step"]) == params.step_count - 1
    )
    final_b = next(
        row
        for row in sorted(timeline_rows, key=lambda item: (int(item["step"]), str(item["particle"])))
        if row["particle"] == "B" and int(row["step"]) == params.step_count - 1
    )
    final_delta_A = float(final_a["delta_chi_mean"])
    final_delta_B = float(final_b["delta_chi_mean"])
    ratio_observed = abs(final_delta_A / final_delta_B) if abs(final_delta_B) > 0.0 else float("inf")
    ratio_expected = params.R_B / params.R_A
    ratio_error = abs(ratio_observed - ratio_expected)
    initial_distance = abs(phase_distance(params.chi_A0, params.chi_B0))
    expected_shape = abs(math.sin(initial_distance))

    return {
        "sweep": sweep,
        "case": case,
        "value": value,
        "R_A": params.R_A,
        "R_B": params.R_B,
        "R_ratio_B_over_A": ratio_expected,
        "epsilon_to_C": params.epsilon_to_C,
        "epsilon_C_return": params.epsilon_C_return,
        "memory_decay": params.memory_decay,
        "initial_distance_AB_abs": initial_distance,
        "expected_abs_sin_distance_shape": expected_shape,
        "final_delta_A": final_delta_A,
        "final_delta_B": final_delta_B,
        "delta_ratio_observed_abs": ratio_observed,
        "delta_ratio_expected_abs": ratio_expected,
        "delta_ratio_error": ratio_error,
        "distance_change": float(summary["distance_change"]),
        "distance_change_abs": abs(float(summary["distance_change"])),
        "R_weighted_delta_balance_final": float(summary["R_weighted_delta_balance_final"]),
        "R_weighted_acceleration_balance_max": float(summary["R_weighted_acceleration_balance_max"]),
        "max_delta_chi_std": float(summary["max_delta_chi_std"]),
        "max_C_memory_abs": max(abs(float(row["C_memory"])) for row in c_rows),
        "gauge_count": len(gauges),
    }


def monotone_non_decreasing(values: List[float], rel_tol: float = 1.0e-9) -> bool:
    for left, right in zip(values, values[1:]):
        if right + rel_tol * max(1.0, abs(left), abs(right)) < left:
            return False
    return True


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    values = [abs(float(row[key])) for row in rows]
    return float(max(values)) if values else 0.0


def summarize_sweep(rows: List[Dict[str, Any]], sweep: str) -> Dict[str, Any]:
    selected = [row for row in rows if row["sweep"] == sweep]
    selected = sorted(selected, key=lambda row: float(row["value"]))
    distance_abs = [float(row["distance_change_abs"]) for row in selected]
    negative_all = all(float(row["distance_change"]) < 0.0 for row in selected)
    balance_max = max_abs(selected, "R_weighted_delta_balance_final")
    accel_balance_max = max_abs(selected, "R_weighted_acceleration_balance_max")
    ratio_error_max = max_abs(selected, "delta_ratio_error")
    std_max = max_abs(selected, "max_delta_chi_std")

    if sweep in {"epsilon_C_return", "epsilon_to_C", "memory_decay"}:
        shape_valid = monotone_non_decreasing(distance_abs)
    elif sweep == "R_ratio_B_over_A":
        shape_valid = monotone_non_decreasing(distance_abs)
    elif sweep == "phase_distance":
        x = np.array([float(row["expected_abs_sin_distance_shape"]) for row in selected])
        y = np.array(distance_abs)
        corr = float(np.corrcoef(x, y)[0, 1]) if len(selected) > 2 and np.std(x) > 0.0 and np.std(y) > 0.0 else 1.0
        shape_valid = corr >= 0.999
    else:
        shape_valid = True

    return {
        "sweep": sweep,
        "case_count": len(selected),
        "distance_change_negative_all_cases": negative_all,
        "distance_change_abs_monotone_or_shape_valid": bool(shape_valid),
        "max_R_weighted_delta_balance": balance_max,
        "max_R_weighted_acceleration_balance": accel_balance_max,
        "max_delta_ratio_error": ratio_error_max,
        "delta_ratio_abs_tolerance": RATIO_ABS_TOL,
        "max_gauge_delta_std": std_max,
        "valid": bool(
            selected
            and negative_all
            and shape_valid
            and balance_max <= 1.0e-12
            and accel_balance_max <= 1.0e-12
            and ratio_error_max <= RATIO_ABS_TOL
            and std_max <= 1.0e-14
        ),
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C媒介応答スイープ予備実験 v1",
        "",
        "## 目的",
        "",
        "C媒介 persistent 応答が、結合強度、R比、位相距離、C残渣減衰率に対して粗く整合したスケーリングを持つかを確認する。",
        "",
        "これは本命の重力的読出しの主張ではない。Stage I に入る前に、C媒介経路が符号、R重み付き収支、ゲージ安定性を保つかを調べる予備実験である。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## スイープ別判定",
            "",
            "| sweep | cases | negative distance | shape valid | max R balance | max R a balance | max ratio error | valid |",
            "|---|---:|---|---|---:|---:|---:|---|",
        ]
    )
    for row in result["sweep_summaries"]:
        lines.append(
            f"| {row['sweep']} | {row['case_count']} | `{row['distance_change_negative_all_cases']}` | "
            f"`{row['distance_change_abs_monotone_or_shape_valid']}` | "
            f"{row['max_R_weighted_delta_balance']:.16e} | "
            f"{row['max_R_weighted_acceleration_balance']:.16e} | "
            f"{row['max_delta_ratio_error']:.16e} | `{row['valid']}` |"
        )

    lines.extend(
        [
            "",
            "## ケース一覧",
            "",
            "| sweep | case | value | R_B/R_A | distance change | A delta | B delta | R balance |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["case_rows"]:
        lines.append(
            f"| {row['sweep']} | {row['case']} | {row['value']:.16e} | "
            f"{row['R_ratio_B_over_A']:.16e} | {row['distance_change']:.16e} | "
            f"{row['final_delta_A']:.16e} | {row['final_delta_B']:.16e} | "
            f"{row['R_weighted_delta_balance_final']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 結合強度と C 残渣減衰率に対して、距離位相縮小量が単調に増えるかを見る。",
            "- R比スイープでは、`|δχ_A/δχ_B| ≈ R_B/R_A` が保たれるかを見る。",
            f"- 比率検査は、今回の効果量が `1e-8` から `1e-7` 程度の微小差分であるため、絶対許容 `RATIO_ABS_TOL={RATIO_ABS_TOL}` として扱う。",
            "- R比スイープの形状判定では、比率誤差そのものではなく、R比に対する距離位相縮小量の単調性を見る。",
            "- 位相距離スイープでは、符号付き C 残渣を戻す今回の制御モデルで期待される `|sin(Δχ)|` 型の粗い形状と整合するかを見る。",
            "- すべてのケースで距離位相が縮み、R重み付き変位収支と加速度収支が小さい場合だけ、C媒介経路の Stage I 進行候補として扱う。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_c_mediated_sweep_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_c_mediated_sweep_cases_v1.csv` |",
            "| sweep CSV | `abc_baseline_stationary_wave_c_mediated_sweep_summary_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_c_mediated_sweep_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C媒介応答スイープ予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    base = Params()
    rows: List[Dict[str, Any]] = []

    for value in [2.5e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3]:
        params = replace(base, epsilon_C_return=value)
        rows.append(run_case(params, "epsilon_C_return", f"epsilon_C_return_{value:g}", value))

    for value in [2.5e-7, 5.0e-7, 1.0e-6, 2.0e-6, 4.0e-6]:
        params = replace(base, epsilon_to_C=value)
        rows.append(run_case(params, "epsilon_to_C", f"epsilon_to_C_{value:g}", value))

    for ratio in [0.25, 0.5, 1.0, 2.0, 4.0, 16.0]:
        params = replace(base, R_A=1.0, R_B=ratio)
        rows.append(run_case(params, "R_ratio_B_over_A", f"R_ratio_{ratio:g}", ratio))

    for distance in [0.2, 0.4, 0.8, 1.2, 1.6]:
        params = replace(base, chi_A0=-distance / 2.0, chi_B0=distance / 2.0)
        rows.append(run_case(params, "phase_distance", f"phase_distance_{distance:g}", distance))

    for value in [0.0, 0.5, 0.8, 0.92, 0.98]:
        params = replace(base, memory_decay=value)
        rows.append(run_case(params, "memory_decay", f"memory_decay_{value:g}", value))

    sweep_names = ["epsilon_C_return", "epsilon_to_C", "R_ratio_B_over_A", "phase_distance", "memory_decay"]
    sweep_summaries = [summarize_sweep(rows, name) for name in sweep_names]

    aggregate_verdict = {
        "sweep_count": len(sweep_summaries),
        "case_count": len(rows),
        "all_sweeps_valid": bool_all(row["valid"] for row in sweep_summaries),
        "single_gauge_only_used": False,
        "c_mediated_sweep_preliminary_valid": bool_all(row["valid"] for row in sweep_summaries),
        "max_R_weighted_delta_balance": max_abs(rows, "R_weighted_delta_balance_final"),
        "max_R_weighted_acceleration_balance": max_abs(rows, "R_weighted_acceleration_balance_max"),
        "max_delta_ratio_error": max_abs(rows, "delta_ratio_error"),
        "max_gauge_delta_std": max_abs(rows, "max_delta_chi_std"),
    }

    result = {
        "experiment": "abc_baseline_stationary_wave_c_mediated_sweep_preliminary_v1",
        "base_params": asdict(base),
        "case_rows": rows,
        "sweep_summaries": sweep_summaries,
        "aggregate_verdict": aggregate_verdict,
    }

    (OUT_DIR / "abc_baseline_stationary_wave_c_mediated_sweep_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_c_mediated_sweep_cases_v1.csv", rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_c_mediated_sweep_summary_v1.csv", sweep_summaries)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
