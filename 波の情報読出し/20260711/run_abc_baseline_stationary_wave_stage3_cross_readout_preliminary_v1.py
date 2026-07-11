from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import Stage2Params
from run_abc_baseline_stationary_wave_transition_protocol_preliminary_v1 import protocol_profiles


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_stage3_cross_readout_preliminary_result_v1"
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


def reduced_R(params: Stage2Params) -> float:
    return float(params.R_A * params.R_B / (params.R_A + params.R_B))


def build_state(profile: List[float], params: Stage2Params) -> Dict[str, List[float]]:
    total = params.R_A + params.R_B
    cumulative_A = 0.0
    cumulative_B = 0.0
    delta_A: List[float] = []
    delta_B: List[float] = []
    increment_A: List[float] = []
    increment_B: List[float] = []
    for c_memory in profile:
        response = params.epsilon_C_return * c_memory
        inc_A = response * params.R_B / total
        inc_B = -response * params.R_A / total
        cumulative_A += inc_A
        cumulative_B += inc_B
        increment_A.append(inc_A)
        increment_B.append(inc_B)
        delta_A.append(cumulative_A)
        delta_B.append(cumulative_B)
    return {
        "delta_A": delta_A,
        "delta_B": delta_B,
        "increment_A": increment_A,
        "increment_B": increment_B,
    }


def local_rows(protocol: str, profile: List[float], metadata: Dict[str, Any], params: Stage2Params) -> List[Dict[str, Any]]:
    state = build_state(profile, params)
    reduced = reduced_R(params)
    rows: List[Dict[str, Any]] = []
    for step in range(1, len(profile) - 1):
        delta_c = profile[step + 1] - profile[step]
        for particle, r_value, sign_value in [("A", params.R_A, 1.0), ("B", params.R_B, -1.0)]:
            delta_values = state[f"delta_{particle}"]
            increment_values = state[f"increment_{particle}"]
            chi_curvature_Ra = r_value * (delta_values[step + 1] - 2.0 * delta_values[step] + delta_values[step - 1])
            p_change_Ra = r_value * (increment_values[step + 1] - increment_values[step])
            g1_slope_Ra = sign_value * reduced * params.epsilon_C_return * delta_c
            rows.append(
                {
                    "protocol": protocol,
                    "kind": metadata["kind"],
                    "duration": metadata["duration"],
                    "step": step,
                    "particle": particle,
                    "C_memory": profile[step],
                    "delta_C_memory_forward": delta_c,
                    "chi_curvature_Ra": chi_curvature_Ra,
                    "p_change_Ra": p_change_Ra,
                    "g1_slope_Ra": g1_slope_Ra,
                    "chi_minus_p_error": chi_curvature_Ra - p_change_Ra,
                    "chi_minus_g1_error": chi_curvature_Ra - g1_slope_Ra,
                    "p_minus_g1_error": p_change_Ra - g1_slope_Ra,
                }
            )
    return rows


def integrated_rows(protocol: str, profile: List[float], metadata: Dict[str, Any], params: Stage2Params) -> Dict[str, Any]:
    state = build_state(profile, params)
    total = params.R_A + params.R_B
    expected_A = sum(params.epsilon_C_return * c_memory * params.R_B / total for c_memory in profile)
    expected_B = sum(-params.epsilon_C_return * c_memory * params.R_A / total for c_memory in profile)
    final_A = state["delta_A"][-1]
    final_B = state["delta_B"][-1]
    return {
        "protocol": protocol,
        "kind": metadata["kind"],
        "duration": metadata["duration"],
        "final_delta_A": final_A,
        "final_delta_B": final_B,
        "integrated_from_C_A": expected_A,
        "integrated_from_C_B": expected_B,
        "integrated_error_A": final_A - expected_A,
        "integrated_error_B": final_B - expected_B,
        "R_weighted_integrated_balance": params.R_A * final_A + params.R_B * final_B,
    }


def max_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return max(abs(float(row[key])) for row in rows)


def sign(value: float, floor: float = 1.0e-18) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


def summarize_protocol(protocol: str, metadata: Dict[str, Any], rows: List[Dict[str, Any]], integrated: Dict[str, Any], params: Stage2Params) -> Dict[str, Any]:
    a_rows = [row for row in rows if row["particle"] == "A"]
    b_rows = [row for row in rows if row["particle"] == "B"]
    active_rows = [row for row in rows if abs(float(row["delta_C_memory_forward"])) > params.effect_floor]
    active_a = [row for row in active_rows if row["particle"] == "A"]
    active_b = [row for row in active_rows if row["particle"] == "B"]
    return {
        "protocol": protocol,
        "kind": metadata["kind"],
        "duration": metadata["duration"],
        "row_count": len(rows),
        "active_row_count": len(active_rows),
        "max_chi_curvature_Ra_abs": max_abs(rows, "chi_curvature_Ra"),
        "max_p_change_Ra_abs": max_abs(rows, "p_change_Ra"),
        "max_g1_slope_Ra_abs": max_abs(rows, "g1_slope_Ra"),
        "max_chi_minus_p_error_abs": max_abs(rows, "chi_minus_p_error"),
        "max_chi_minus_g1_error_abs": max_abs(rows, "chi_minus_g1_error"),
        "max_p_minus_g1_error_abs": max_abs(rows, "p_minus_g1_error"),
        "A_active_sign_match_chi_p_g1": bool_all(
            sign(float(row["chi_curvature_Ra"])) == sign(float(row["p_change_Ra"])) == sign(float(row["g1_slope_Ra"]))
            for row in active_a
        ),
        "B_active_sign_match_chi_p_g1": bool_all(
            sign(float(row["chi_curvature_Ra"])) == sign(float(row["p_change_Ra"])) == sign(float(row["g1_slope_Ra"]))
            for row in active_b
        ),
        "A_B_opposite_when_active": bool_all(
            sign(float(a_row["chi_curvature_Ra"])) == -sign(float(b_row["chi_curvature_Ra"]))
            for a_row, b_row in zip(active_a, active_b)
        ),
        "integrated_error_A_abs": abs(float(integrated["integrated_error_A"])),
        "integrated_error_B_abs": abs(float(integrated["integrated_error_B"])),
        "R_weighted_integrated_balance_abs": abs(float(integrated["R_weighted_integrated_balance"])),
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C Stage III 別読出し照合予備実験 v1",
        "",
        "## 目的",
        "",
        "準安定傾斜から読まれる加速度候補が、単一の観測方法に依存していないかを検査する。",
        "",
        "同じ遷移プロトコルに対して、位置位相二階差分、`p_read` 相当の速度差分、C傾斜による G1 予測、Cからの積分再構成を照合する。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## プロトコル別サマリー",
            "",
            "| protocol | kind | active | max chi | max p | max G1 | chi-p err | chi-G1 err | signs | integrated balance |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|---:|",
        ]
    )
    for row in result["protocol_summaries"]:
        signs = bool(row["A_active_sign_match_chi_p_g1"]) and bool(row["B_active_sign_match_chi_p_g1"]) and bool(row["A_B_opposite_when_active"])
        lines.append(
            f"| {row['protocol']} | {row['kind']} | {row['active_row_count']} | "
            f"{row['max_chi_curvature_Ra_abs']:.16e} | {row['max_p_change_Ra_abs']:.16e} | "
            f"{row['max_g1_slope_Ra_abs']:.16e} | {row['max_chi_minus_p_error_abs']:.16e} | "
            f"{row['max_chi_minus_g1_error_abs']:.16e} | `{signs}` | "
            f"{row['R_weighted_integrated_balance_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 位置位相二階差分と `p_read` 差分が一致する場合、加速度候補は位置読出しだけの人工物ではない。",
            "- それらが C傾斜による G1 予測とも一致する場合、時間位相側の補償読出しとも整合する。",
            "- Cからの積分再構成が最終変位と一致し、R重み付き積分バランスが保たれる場合、局所差分だけの偶然ではない。",
            "- 本実験は標準重力の導出ではなく、準安定傾斜候補の観測方法依存性を下げるための Stage III 予備照合である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_stage3_cross_readout_preliminary_result_v1.json` |",
            "| protocol CSV | `abc_baseline_stationary_wave_stage3_cross_readout_cases_v1.csv` |",
            "| rows CSV | `abc_baseline_stationary_wave_stage3_cross_readout_rows_v1.csv` |",
            "| integrated CSV | `abc_baseline_stationary_wave_stage3_cross_readout_integrated_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_stage3_cross_readout_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C Stage III別読出し照合予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params()
    all_rows: List[Dict[str, Any]] = []
    integrated_summaries: List[Dict[str, Any]] = []
    protocol_summaries: List[Dict[str, Any]] = []

    for protocol, (profile, metadata) in protocol_profiles(params).items():
        rows = local_rows(protocol, profile, metadata, params)
        integrated = integrated_rows(protocol, profile, metadata, params)
        all_rows.extend(rows)
        integrated_summaries.append(integrated)
        protocol_summaries.append(summarize_protocol(protocol, metadata, rows, integrated, params))

    aggregate_verdict = {
        "protocol_count": len(protocol_summaries),
        "single_gauge_only_used": False,
        "all_chi_p_readouts_match": bool_all(
            float(row["max_chi_minus_p_error_abs"]) <= params.effect_floor for row in protocol_summaries
        ),
        "all_chi_g1_readouts_match": bool_all(
            float(row["max_chi_minus_g1_error_abs"]) <= params.effect_floor for row in protocol_summaries
        ),
        "all_p_g1_readouts_match": bool_all(
            float(row["max_p_minus_g1_error_abs"]) <= params.effect_floor for row in protocol_summaries
        ),
        "all_active_signs_consistent": bool_all(
            bool(row["A_active_sign_match_chi_p_g1"])
            and bool(row["B_active_sign_match_chi_p_g1"])
            and bool(row["A_B_opposite_when_active"])
            for row in protocol_summaries
        ),
        "all_integrated_errors_near_zero": bool_all(
            float(row["integrated_error_A_abs"]) <= params.effect_floor
            and float(row["integrated_error_B_abs"]) <= params.effect_floor
            for row in protocol_summaries
        ),
        "all_R_weighted_integrated_balances_near_zero": bool_all(
            float(row["R_weighted_integrated_balance_abs"]) <= params.effect_floor for row in protocol_summaries
        ),
        "max_chi_minus_g1_error_abs": max(float(row["max_chi_minus_g1_error_abs"]) for row in protocol_summaries),
        "max_integrated_error_abs": max(
            max(float(row["integrated_error_A_abs"]), float(row["integrated_error_B_abs"]))
            for row in protocol_summaries
        ),
        "stage3_cross_readout_preliminary_valid": bool(
            bool_all(float(row["max_chi_minus_p_error_abs"]) <= params.effect_floor for row in protocol_summaries)
            and bool_all(float(row["max_chi_minus_g1_error_abs"]) <= params.effect_floor for row in protocol_summaries)
            and bool_all(
                bool(row["A_active_sign_match_chi_p_g1"])
                and bool(row["B_active_sign_match_chi_p_g1"])
                and bool(row["A_B_opposite_when_active"])
                for row in protocol_summaries
            )
            and bool_all(
                float(row["integrated_error_A_abs"]) <= params.effect_floor
                and float(row["integrated_error_B_abs"]) <= params.effect_floor
                for row in protocol_summaries
            )
        ),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_stage3_cross_readout_preliminary_v1",
        "params": asdict(params),
        "protocol_summaries": protocol_summaries,
        "integrated_summaries": integrated_summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_stage3_cross_readout_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_stage3_cross_readout_cases_v1.csv", protocol_summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_stage3_cross_readout_rows_v1.csv", all_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_stage3_cross_readout_integrated_v1.csv", integrated_summaries)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
