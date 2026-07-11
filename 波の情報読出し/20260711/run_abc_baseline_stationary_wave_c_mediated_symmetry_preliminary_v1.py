from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Dict, Iterable, List

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
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_c_mediated_symmetry_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

BALANCE_TOL = 1.0e-12
GAUGE_STD_TOL = 1.0e-14
EFFECT_FLOOR = 1.0e-12


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


def expected_attraction_signs(params: Params, inverted: bool = False) -> Dict[str, float]:
    d = phase_distance(params.chi_A0, params.chi_B0)
    direction_to_b = 1.0 if d > 0.0 else -1.0
    if abs(d) <= 0.0:
        direction_to_b = 0.0
    factor = -1.0 if inverted else 1.0
    return {
        "expected_delta_A_sign": factor * direction_to_b,
        "expected_delta_B_sign": -factor * direction_to_b,
        "expected_distance_change_sign": -factor,
    }


def sign(value: float, floor: float = EFFECT_FLOOR) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


def run_case(case: str, params: Params, inverted: bool = False) -> Dict[str, Any]:
    gauges = default_gauges()
    gauge_rows, _ = simulate_phase3_persistent(params, gauges)
    timeline_rows = summarize_gauge_rows(gauge_rows)
    add_derivatives(timeline_rows, params)
    summary = summarize_phase(
        timeline_rows,
        params,
        "phase3_c_return_persistent",
        "C_return_with_persistent_reembedding",
        case,
    )

    ordered = sorted(timeline_rows, key=lambda row: (int(row["step"]), str(row["particle"])))
    final_a = next(row for row in ordered if row["particle"] == "A" and int(row["step"]) == params.step_count - 1)
    final_b = next(row for row in ordered if row["particle"] == "B" and int(row["step"]) == params.step_count - 1)
    expected = expected_attraction_signs(params, inverted=inverted)
    delta_A = float(final_a["delta_chi_mean"])
    delta_B = float(final_b["delta_chi_mean"])
    distance_change = float(summary["distance_change"])
    observed = {
        "observed_delta_A_sign": sign(delta_A),
        "observed_delta_B_sign": sign(delta_B),
        "observed_distance_change_sign": sign(distance_change),
    }
    sign_valid = bool(
        observed["observed_delta_A_sign"] == expected["expected_delta_A_sign"]
        and observed["observed_delta_B_sign"] == expected["expected_delta_B_sign"]
        and observed["observed_distance_change_sign"] == expected["expected_distance_change_sign"]
    )
    valid = bool(
        sign_valid
        and abs(float(summary["R_weighted_delta_balance_final"])) <= BALANCE_TOL
        and float(summary["R_weighted_acceleration_balance_max"]) <= BALANCE_TOL
        and float(summary["max_delta_chi_std"]) <= GAUGE_STD_TOL
    )

    return {
        "case": case,
        "chi_A0": params.chi_A0,
        "chi_B0": params.chi_B0,
        "R_A": params.R_A,
        "R_B": params.R_B,
        "epsilon_to_C": params.epsilon_to_C,
        "epsilon_C_return": params.epsilon_C_return,
        "initial_signed_distance_B_minus_A": phase_distance(params.chi_A0, params.chi_B0),
        "final_delta_A": delta_A,
        "final_delta_B": delta_B,
        "distance_change": distance_change,
        "R_weighted_delta_balance_final": float(summary["R_weighted_delta_balance_final"]),
        "R_weighted_acceleration_balance_max": float(summary["R_weighted_acceleration_balance_max"]),
        "max_delta_chi_std": float(summary["max_delta_chi_std"]),
        **expected,
        **observed,
        "sign_valid": sign_valid,
        "valid": valid,
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C媒介応答の符号対称性予備実験 v1",
        "",
        "## 目的",
        "",
        "C媒介 persistent 応答が、A/B というラベルに固定された押し引きではなく、位置位相差の符号に従うかを確認する。",
        "",
        "これは、引力的に見える距離位相縮小が、単なる実装上のラベル依存反力ではないかを切り分けるための予備テストである。",
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
            "| case | d_AB sign | A delta | B delta | distance change | expected signs A/B/d | observed signs A/B/d | sign valid | valid |",
            "|---|---:|---:|---:|---:|---|---|---|---|",
        ]
    )
    for row in result["case_rows"]:
        lines.append(
            f"| {row['case']} | {row['initial_signed_distance_B_minus_A']:.16e} | "
            f"{row['final_delta_A']:.16e} | {row['final_delta_B']:.16e} | {row['distance_change']:.16e} | "
            f"{row['expected_delta_A_sign']:.0f}/{row['expected_delta_B_sign']:.0f}/{row['expected_distance_change_sign']:.0f} | "
            f"{row['observed_delta_A_sign']:.0f}/{row['observed_delta_B_sign']:.0f}/{row['observed_distance_change_sign']:.0f} | "
            f"`{row['sign_valid']}` | `{row['valid']}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- normal と mirrored_positions の両方で距離位相が縮む場合、応答は少なくともラベル固定の片側押しではない。",
            "- inverted_C_source と inverted_C_return では距離位相が広がることを要求する。これにより、C 残渣の符号枝を反転したときに応答符号も反転するかを見る。",
            "- R重み付き変位収支と加速度収支が小さいことを同時に要求する。",
            "- この予備テストが失敗する場合、C媒介応答は Stage I 候補ではなく、実装ラベル依存または読出し枝依存として保留する。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_c_mediated_symmetry_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_c_mediated_symmetry_cases_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_c_mediated_symmetry_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C媒介応答の符号対称性予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    base = Params()
    case_rows = [
        run_case("normal", base),
        run_case("mirrored_positions", replace(base, chi_A0=0.25, chi_B0=-0.25)),
        run_case("inverted_C_source", replace(base, epsilon_to_C=-base.epsilon_to_C), inverted=True),
        run_case("inverted_C_return", replace(base, epsilon_C_return=-base.epsilon_C_return), inverted=True),
    ]
    aggregate_verdict = {
        "case_count": len(case_rows),
        "all_cases_valid": bool_all(row["valid"] for row in case_rows),
        "single_gauge_only_used": False,
        "c_mediated_symmetry_preliminary_valid": bool_all(row["valid"] for row in case_rows),
        "max_R_weighted_delta_balance": max(abs(float(row["R_weighted_delta_balance_final"])) for row in case_rows),
        "max_R_weighted_acceleration_balance": max(
            abs(float(row["R_weighted_acceleration_balance_max"])) for row in case_rows
        ),
        "max_gauge_delta_std": max(abs(float(row["max_delta_chi_std"])) for row in case_rows),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_c_mediated_symmetry_preliminary_v1",
        "base_params": asdict(base),
        "case_rows": case_rows,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_c_mediated_symmetry_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_c_mediated_symmetry_cases_v1.csv", case_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
