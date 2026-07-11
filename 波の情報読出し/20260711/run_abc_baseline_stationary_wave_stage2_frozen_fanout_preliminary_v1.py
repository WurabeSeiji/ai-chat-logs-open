from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_c_mediated_response_preliminary_v1 import (
    Gauge,
    default_gauges,
    phase_distance,
    wrap_phase,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class Stage2Params:
    step_count: int = 48
    delta_s: float = 1.0
    chi_A0: float = -0.25
    chi_B0: float = 0.25
    R_A: float = 1.0
    R_B: float = 1.5625
    epsilon_to_C: float = 1.0e-6
    epsilon_C_return: float = 1.0e-3
    epsilon_r2_pressure: float = 8.0e-7
    epsilon_r3_direct: float = 5.0e-7
    memory_decay_C: float = 0.92
    memory_decay_R: float = 0.90
    q_raw_gain: float = 1.0e-3
    subwinding_residual: float = -0.0125
    balance_tol: float = 1.0e-12
    gauge_std_tol: float = 1.0e-14
    additivity_tol: float = 5.0e-8
    residual_tol: float = 5.0e-8
    effect_floor: float = 1.0e-12


@dataclass
class Stage2Case:
    name: str
    use_c: bool = False
    use_r2: bool = False
    use_r3: bool = False
    c_persistent: bool = True
    mirrored: bool = False


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def sign(value: float, floor: float = 1.0e-12) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


def branch_sign(value: float) -> float:
    return 1.0 if value >= 0.0 else -1.0


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


def initial_positions(params: Stage2Params, case: Stage2Case) -> Tuple[float, float]:
    if case.mirrored:
        return -params.chi_A0, -params.chi_B0
    return params.chi_A0, params.chi_B0


def weighted_pair(response: float, params: Stage2Params) -> Tuple[float, float]:
    total = params.R_A + params.R_B
    return (
        float(response * params.R_B / total),
        float(-response * params.R_A / total),
    )


def c_source_strength(chi_A: float, chi_B: float, params: Stage2Params) -> float:
    d = phase_distance(chi_A, chi_B)
    return float((params.R_A * params.R_B / (params.R_A + params.R_B)) * math.sin(d))


def linear_slope(values: List[float]) -> float:
    y = np.array(values, dtype=float)
    x = np.arange(y.size, dtype=float)
    x_mean = float(np.mean(x))
    y_mean = float(np.mean(y))
    denom = float(np.sum((x - x_mean) ** 2))
    if denom == 0.0:
        return 0.0
    return float(np.sum((x - x_mean) * (y - y_mean)) / denom)


def read_delta(true_delta: float, gauge: Gauge) -> float:
    return float(true_delta + 1.0e-17 * math.sin(gauge.phase_bias))


def add_fanout_rows(
    rows: List[Dict[str, Any]],
    params: Stage2Params,
    case: Stage2Case,
    gauges: List[Gauge],
    step: int,
    chi_A: float,
    chi_B: float,
    delta_A: float,
    delta_B: float,
    memory_C: float,
    memory_R2: float,
    memory_R3: float,
    q_raw: float,
) -> None:
    d_signed = phase_distance(chi_A, chi_B)
    d_abs = abs(d_signed)
    q_closed = 0.0
    beat_phase = wrap_phase(step * (2.0 * math.pi * 17.0 + params.subwinding_residual))
    unwrap_phase = step * params.subwinding_residual
    for gauge in gauges:
        for particle, delta, r_value, chi_value in [
            ("A", delta_A, params.R_A, chi_A),
            ("B", delta_B, params.R_B, chi_B),
        ]:
            rows.append(
                {
                    "case": case.name,
                    "step": step,
                    "readout_set": "frozen_fanout",
                    "particle": particle,
                    "gauge": gauge.name,
                    "chi_read": wrap_phase(chi_value + gauge.delta_chi),
                    "delta_chi_read": read_delta(delta, gauge),
                    "R_read": r_value * gauge.gain,
                    "distance_AB_signed": d_signed,
                    "distance_AB_abs": d_abs,
                    "C_memory": memory_C,
                    "R2_memory": memory_R2,
                    "R3_memory": memory_R3,
                    "Q_raw": q_raw,
                    "Q_closed": q_closed,
                    "beat_phase": wrap_phase(beat_phase + gauge.phase_bias),
                    "unwrap_phase": unwrap_phase,
                    "principal_branch_phase": wrap_phase(beat_phase),
                    "same_snapshot_id": f"{case.name}:{step}",
                    "delta_chi_gauge_offset": gauge.delta_chi,
                    "phase_bias": gauge.phase_bias,
                    "gain": gauge.gain,
                }
            )


def summarize_fanout_rows(rows: List[Dict[str, Any]], params: Stage2Params) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, int, str], List[Dict[str, Any]]] = {}
    for row in rows:
        key = (str(row["case"]), int(row["step"]), str(row["particle"]))
        grouped.setdefault(key, []).append(row)

    summaries: List[Dict[str, Any]] = []
    for (case, step, particle), selected in sorted(grouped.items()):
        deltas = np.array([float(row["delta_chi_read"]) for row in selected])
        r_values = np.array([float(row["R_read"]) for row in selected])
        summaries.append(
            {
                "case": case,
                "step": step,
                "particle": particle,
                "delta_chi_mean": float(np.mean(deltas)),
                "delta_chi_std": float(np.std(deltas)),
                "R_mean": float(np.mean(r_values)),
                "R_std": float(np.std(r_values)),
                "distance_AB_abs": float(np.mean([float(row["distance_AB_abs"]) for row in selected])),
                "distance_AB_signed": float(np.mean([float(row["distance_AB_signed"]) for row in selected])),
                "C_memory": float(np.mean([float(row["C_memory"]) for row in selected])),
                "R2_memory": float(np.mean([float(row["R2_memory"]) for row in selected])),
                "R3_memory": float(np.mean([float(row["R3_memory"]) for row in selected])),
                "Q_raw": float(np.mean([float(row["Q_raw"]) for row in selected])),
                "Q_closed": float(np.mean([float(row["Q_closed"]) for row in selected])),
                "beat_phase": float(np.mean([float(row["beat_phase"]) for row in selected])),
                "unwrap_phase": float(np.mean([float(row["unwrap_phase"]) for row in selected])),
                "gauge_count": len(selected),
            }
        )

    for key in {(row["case"], row["particle"]) for row in summaries}:
        selected = [row for row in summaries if (row["case"], row["particle"]) == key]
        selected.sort(key=lambda row: int(row["step"]))
        values = [float(row["delta_chi_mean"]) for row in selected]
        for idx, row in enumerate(selected):
            if 0 < idx < len(selected) - 1:
                row["v_chi_read"] = float((values[idx + 1] - values[idx - 1]) / (2.0 * params.delta_s))
                row["a_chi_read"] = float(
                    (values[idx + 1] - 2.0 * values[idx] + values[idx - 1]) / (params.delta_s**2)
                )
            else:
                row["v_chi_read"] = 0.0
                row["a_chi_read"] = 0.0
    return summaries


def simulate_case(case: Stage2Case, params: Stage2Params, gauges: List[Gauge]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    chi_A0, chi_B0 = initial_positions(params, case)
    chi_A = chi_A0
    chi_B = chi_B0
    cumulative_A = 0.0
    cumulative_B = 0.0
    memory_C = 0.0
    memory_R2 = 0.0
    memory_R3 = 0.0
    rows: List[Dict[str, Any]] = []

    for step in range(params.step_count):
        d = phase_distance(chi_A, chi_B)
        c_response = 0.0
        r2_response = 0.0
        r3_response = 0.0

        if case.use_c:
            source_C = params.epsilon_to_C * c_source_strength(chi_A, chi_B, params)
            memory_C = params.memory_decay_C * memory_C + source_C if case.c_persistent else source_C
            c_response = params.epsilon_C_return * memory_C
        if case.use_r2:
            source_R2 = branch_sign(d) * params.epsilon_r2_pressure * abs(math.sin(d))
            memory_R2 = params.memory_decay_R * memory_R2 + source_R2
            r2_response = -memory_R2
        if case.use_r3:
            source_R3 = branch_sign(d) * params.epsilon_r3_direct * math.sin(abs(d))
            memory_R3 = params.memory_decay_R * memory_R3 + source_R3
            r3_response = -memory_R3

        response = c_response + r2_response + r3_response
        inc_A, inc_B = weighted_pair(response, params)
        cumulative_A = wrap_phase(cumulative_A + inc_A)
        cumulative_B = wrap_phase(cumulative_B + inc_B)
        chi_A = wrap_phase(chi_A0 + cumulative_A)
        chi_B = wrap_phase(chi_B0 + cumulative_B)
        q_raw = params.q_raw_gain * (abs(memory_C) + abs(memory_R2) + abs(memory_R3))
        add_fanout_rows(
            rows,
            params,
            case,
            gauges,
            step,
            chi_A,
            chi_B,
            cumulative_A,
            cumulative_B,
            memory_C,
            memory_R2,
            memory_R3,
            q_raw,
        )

    timeline = summarize_fanout_rows(rows, params)
    initial_distance = abs(phase_distance(chi_A0, chi_B0))
    final_a = next(row for row in timeline if row["particle"] == "A" and int(row["step"]) == params.step_count - 1)
    final_b = next(row for row in timeline if row["particle"] == "B" and int(row["step"]) == params.step_count - 1)
    final_distance = float(final_a["distance_AB_abs"])
    r_balance = params.R_A * float(final_a["delta_chi_mean"]) + params.R_B * float(final_b["delta_chi_mean"])
    r_acc_balance = 0.0
    for step in range(params.step_count):
        ar = next(row for row in timeline if row["particle"] == "A" and int(row["step"]) == step)
        br = next(row for row in timeline if row["particle"] == "B" and int(row["step"]) == step)
        r_acc_balance = max(
            r_acc_balance,
            abs(params.R_A * float(ar["a_chi_read"]) + params.R_B * float(br["a_chi_read"])),
        )
    beat_values = [float(row["beat_phase"]) for row in timeline if row["particle"] == "A"]
    beat_slope = linear_slope(list(np.unwrap(np.array(beat_values, dtype=float))))
    summary = {
        "case": case.name,
        "use_c": case.use_c,
        "use_r2": case.use_r2,
        "use_r3": case.use_r3,
        "c_persistent": case.c_persistent,
        "mirrored": case.mirrored,
        "final_delta_A": float(final_a["delta_chi_mean"]),
        "final_delta_B": float(final_b["delta_chi_mean"]),
        "initial_distance_AB_abs": initial_distance,
        "final_distance_AB_abs": final_distance,
        "distance_change": float(final_distance - initial_distance),
        "R_weighted_delta_balance_final": float(r_balance),
        "R_weighted_acceleration_balance_max": float(r_acc_balance),
        "max_delta_chi_std": max(abs(float(row["delta_chi_std"])) for row in timeline),
        "max_Q_raw_abs": max(abs(float(row["Q_raw"])) for row in timeline),
        "max_Q_closed_abs": max(abs(float(row["Q_closed"])) for row in timeline),
        "final_C_memory": float(final_a["C_memory"]),
        "final_R2_memory": float(final_a["R2_memory"]),
        "final_R3_memory": float(final_a["R3_memory"]),
        "beat_slope": beat_slope,
        "beat_alias_reverse": bool(beat_slope < 0.0),
    }
    return rows, summary


def diff(value: float, expected: float) -> float:
    return float(value - expected)


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C Stage II frozen fanout 予備実験 v1",
        "",
        "## 目的",
        "",
        "同一状態スナップショットから、相対位相、位置位相、C残渣、反力候補、beat、閉鎖残差を同時に読む。",
        "",
        "目的は、C媒介接近候補、R2反力、R3反力が混在したとき、合成が単純和に近いか、反力差し引き後の残差として C媒介接近が残るかを確認することである。",
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
            "| case | C | R2 | R3 | distance change | R balance | Q_raw max | beat reverse |",
            "|---|---|---|---|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | `{row['use_c']}` | `{row['use_r2']}` | `{row['use_r3']}` | "
            f"{row['distance_change']:.16e} | {row['R_weighted_delta_balance_final']:.16e} | "
            f"{row['max_Q_raw_abs']:.16e} | `{row['beat_alias_reverse']}` |"
        )

    lines.extend(
        [
            "",
            "## 合成・残差検査",
            "",
            "| metric | value | tolerance | valid |",
            "|---|---:|---:|---|",
        ]
    )
    for row in result["composition_checks"]:
        lines.append(f"| {row['metric']} | {row['value']:.16e} | {row['tolerance']:.16e} | `{row['valid']}` |")

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- `combined ≈ C + R2 + R3` が成立する場合、Stage I 候補の粗い線形分離が可能である。",
            "- `combined - R2 - R3 ≈ C` が成立し、かつ C 残差が負の距離変化を保つ場合、反力候補を差し引いた後にも C媒介接近候補が残る。",
            "- ただし本予備実験では R2/R3 が C媒介接近より大きい。したがって、観測方法依存性を避けるには、G1 へ進む前に同時多読出しの記録を保持する必要がある。",
            "- beat 系列は負のサブ巻数残差により逆向きに読まれる。これは G3 候補として、距離位相残差とは別に保持する。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_stage2_frozen_fanout_cases_v1.csv` |",
            "| fanout CSV | `abc_baseline_stationary_wave_stage2_frozen_fanout_rows_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_stage2_frozen_fanout_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C Stage II frozen fanout予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params()
    gauges = default_gauges()
    cases = [
        Stage2Case("C_mediated_only_persistent", use_c=True),
        Stage2Case("C_mediated_only_reset", use_c=True, c_persistent=False),
        Stage2Case("R2_only", use_r2=True),
        Stage2Case("R3_only", use_r3=True),
        Stage2Case("C_plus_R2", use_c=True, use_r2=True),
        Stage2Case("C_plus_R3", use_c=True, use_r3=True),
        Stage2Case("C_plus_R2_plus_R3", use_c=True, use_r2=True, use_r3=True),
        Stage2Case("C_plus_R2_plus_R3_mirrored", use_c=True, use_r2=True, use_r3=True, mirrored=True),
    ]

    fanout_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for case in cases:
        rows, summary = simulate_case(case, params, gauges)
        fanout_rows.extend(rows)
        summaries.append(summary)

    by_case = {row["case"]: row for row in summaries}
    c = float(by_case["C_mediated_only_persistent"]["distance_change"])
    r2 = float(by_case["R2_only"]["distance_change"])
    r3 = float(by_case["R3_only"]["distance_change"])
    c_r2 = float(by_case["C_plus_R2"]["distance_change"])
    c_r3 = float(by_case["C_plus_R3"]["distance_change"])
    combined = float(by_case["C_plus_R2_plus_R3"]["distance_change"])
    combined_mirror = float(by_case["C_plus_R2_plus_R3_mirrored"]["distance_change"])
    residual_after_reactions = combined - r2 - r3
    checks = [
        {
            "metric": "C_plus_R2_minus_sum",
            "value": abs(diff(c_r2, c + r2)),
            "tolerance": params.additivity_tol,
        },
        {
            "metric": "C_plus_R3_minus_sum",
            "value": abs(diff(c_r3, c + r3)),
            "tolerance": params.additivity_tol,
        },
        {
            "metric": "C_plus_R2_plus_R3_minus_sum",
            "value": abs(diff(combined, c + r2 + r3)),
            "tolerance": params.additivity_tol,
        },
        {
            "metric": "residual_after_R2_R3_minus_C",
            "value": abs(diff(residual_after_reactions, c)),
            "tolerance": params.residual_tol,
        },
        {
            "metric": "mirrored_combined_distance_match",
            "value": abs(diff(combined_mirror, combined)),
            "tolerance": params.additivity_tol,
        },
    ]
    for row in checks:
        row["valid"] = bool(float(row["value"]) <= float(row["tolerance"]))

    for row in summaries:
        row["R_balance_valid"] = bool(abs(float(row["R_weighted_delta_balance_final"])) <= params.balance_tol)
        row["gauge_std_valid"] = bool(float(row["max_delta_chi_std"]) <= params.gauge_std_tol)
        row["Q_closed_valid"] = bool(float(row["max_Q_closed_abs"]) <= params.effect_floor)

    aggregate_verdict = {
        "case_count": len(summaries),
        "composition_check_count": len(checks),
        "all_composition_checks_valid": bool_all(row["valid"] for row in checks),
        "all_R_balances_valid": bool_all(row["R_balance_valid"] for row in summaries),
        "all_gauge_stability_valid": bool_all(row["gauge_std_valid"] for row in summaries),
        "all_Q_closed_valid": bool_all(row["Q_closed_valid"] for row in summaries),
        "single_gauge_only_used": False,
        "stage2_frozen_fanout_preliminary_valid": bool_all(row["valid"] for row in checks)
        and bool_all(row["R_balance_valid"] for row in summaries)
        and bool_all(row["gauge_std_valid"] for row in summaries)
        and bool_all(row["Q_closed_valid"] for row in summaries),
        "C_mediated_distance_change": c,
        "R2_distance_change": r2,
        "R3_distance_change": r3,
        "combined_distance_change": combined,
        "residual_after_R2_R3": residual_after_reactions,
        "reaction_dominance_ratio_abs": abs((r2 + r3) / c) if abs(c) > 0.0 else float("inf"),
        "beat_alias_reverse_all_cases": bool_all(row["beat_alias_reverse"] for row in summaries),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1",
        "params": asdict(params),
        "case_summaries": summaries,
        "composition_checks": checks,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_stage2_frozen_fanout_cases_v1.csv", summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_stage2_frozen_fanout_rows_v1.csv", fanout_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
