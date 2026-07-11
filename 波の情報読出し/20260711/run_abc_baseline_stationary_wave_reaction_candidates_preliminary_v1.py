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
    Params,
    add_derivatives,
    default_gauges,
    phase_distance,
    summarize_gauge_rows,
    summarize_phase,
    wrap_phase,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_reaction_candidates_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class ReactionParams:
    step_count: int = 48
    delta_s: float = 1.0
    chi_A0: float = -0.25
    chi_B0: float = 0.25
    R_A: float = 1.0
    R_B: float = 1.5625
    epsilon_r1_readout: float = 2.0e-8
    epsilon_r2_pressure: float = 8.0e-7
    epsilon_r3_direct: float = 5.0e-7
    memory_decay: float = 0.90
    balance_tol: float = 1.0e-12
    gauge_std_tol: float = 1.0e-14
    effect_floor: float = 1.0e-12


@dataclass
class ReactionCase:
    candidate: str
    case: str
    reference_sign: float = 1.0
    mirror_positions: bool = False
    c_pressure_sign: float = 0.0
    direct_sign: float = 0.0
    persistent: bool = True


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


def sign(value: float, floor: float = 1.0e-12) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


def branch_sign(value: float) -> float:
    return 1.0 if value >= 0.0 else -1.0


def add_rows(
    rows: List[Dict[str, Any]],
    params: ReactionParams,
    gauges: List[Gauge],
    candidate: str,
    case: str,
    step: int,
    chi_A: float,
    chi_B: float,
    delta_A: float,
    delta_B: float,
    source_memory: float,
    reembedded: bool,
) -> None:
    distance = abs(phase_distance(chi_A, chi_B))
    for particle, delta, r_value, chi_value in [
        ("A", delta_A, params.R_A, chi_A),
        ("B", delta_B, params.R_B, chi_B),
    ]:
        for gauge in gauges:
            rows.append(
                {
                    "phase": candidate,
                    "case": case,
                    "step": step,
                    "particle": particle,
                    "gauge": gauge.name,
                    "chi_read": wrap_phase(chi_value + gauge.delta_chi),
                    "delta_chi_read": float(delta + 1.0e-17 * math.sin(gauge.phase_bias)),
                    "R_read": r_value * gauge.gain,
                    "C_memory": source_memory,
                    "distance_AB_abs": distance,
                    "reembedded": reembedded,
                    "delta_chi_gauge_offset": gauge.delta_chi,
                    "phase_bias": gauge.phase_bias,
                    "gain": gauge.gain,
                }
            )


def initial_positions(params: ReactionParams, case: ReactionCase) -> Tuple[float, float]:
    if case.mirror_positions:
        return -params.chi_A0, -params.chi_B0
    return params.chi_A0, params.chi_B0


def weighted_pair(response: float, params: ReactionParams) -> Tuple[float, float]:
    total = params.R_A + params.R_B
    return (
        float(response * params.R_B / total),
        float(-response * params.R_A / total),
    )


def simulate_r1(case: ReactionCase, params: ReactionParams, gauges: List[Gauge]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    chi_A, chi_B = initial_positions(params, case)
    # R1 is a readout-reference common mode: it should not change distance and
    # should not satisfy R-weighted external translation balance.
    common_bias = case.reference_sign * params.epsilon_r1_readout
    for step in range(params.step_count):
        add_rows(
            rows,
            params,
            gauges,
            case.candidate,
            case.case,
            step,
            chi_A,
            chi_B,
            common_bias,
            common_bias,
            common_bias,
            reembedded=False,
        )
    return rows


def simulate_r2(case: ReactionCase, params: ReactionParams, gauges: List[Gauge]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    chi_A, chi_B = initial_positions(params, case)
    cumulative_A = 0.0
    cumulative_B = 0.0
    memory = 0.0
    for step in range(params.step_count):
        d = phase_distance(chi_A, chi_B)
        pressure = case.c_pressure_sign * branch_sign(d) * params.epsilon_r2_pressure * abs(math.sin(d))
        memory = params.memory_decay * memory + pressure if case.persistent else pressure
        delta_A, delta_B = weighted_pair(-memory, params)
        cumulative_A = wrap_phase(cumulative_A + delta_A)
        cumulative_B = wrap_phase(cumulative_B + delta_B)
        chi_A = wrap_phase(initial_positions(params, case)[0] + cumulative_A)
        chi_B = wrap_phase(initial_positions(params, case)[1] + cumulative_B)
        add_rows(
            rows,
            params,
            gauges,
            case.candidate,
            case.case,
            step,
            chi_A,
            chi_B,
            cumulative_A,
            cumulative_B,
            memory,
            reembedded=True,
        )
    return rows


def simulate_r3(case: ReactionCase, params: ReactionParams, gauges: List[Gauge]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    chi_A, chi_B = initial_positions(params, case)
    cumulative_A = 0.0
    cumulative_B = 0.0
    memory = 0.0
    for step in range(params.step_count):
        d = phase_distance(chi_A, chi_B)
        source = case.direct_sign * branch_sign(d) * params.epsilon_r3_direct * math.sin(abs(d))
        memory = params.memory_decay * memory + source if case.persistent else source
        delta_A, delta_B = weighted_pair(-memory, params)
        cumulative_A = wrap_phase(cumulative_A + delta_A)
        cumulative_B = wrap_phase(cumulative_B + delta_B)
        chi_A = wrap_phase(initial_positions(params, case)[0] + cumulative_A)
        chi_B = wrap_phase(initial_positions(params, case)[1] + cumulative_B)
        add_rows(
            rows,
            params,
            gauges,
            case.candidate,
            case.case,
            step,
            chi_A,
            chi_B,
            cumulative_A,
            cumulative_B,
            memory,
            reembedded=True,
        )
    return rows


def run_case(case: ReactionCase, params: ReactionParams, gauges: List[Gauge]) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    if case.candidate == "R1":
        gauge_rows = simulate_r1(case, params, gauges)
    elif case.candidate == "R2":
        gauge_rows = simulate_r2(case, params, gauges)
    elif case.candidate == "R3":
        gauge_rows = simulate_r3(case, params, gauges)
    else:
        raise ValueError(f"unknown candidate: {case.candidate}")

    timeline_rows = summarize_gauge_rows(gauge_rows)
    add_derivatives(timeline_rows, Params(step_count=params.step_count, delta_s=params.delta_s))
    summary = summarize_phase(timeline_rows, Params(chi_A0=initial_positions(params, case)[0], chi_B0=initial_positions(params, case)[1], R_A=params.R_A, R_B=params.R_B), case.candidate, case.case, case.case)
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
    delta_A = float(final_a["delta_chi_mean"])
    delta_B = float(final_b["delta_chi_mean"])
    distance_change = float(summary["distance_change"])
    r_balance = float(summary["R_weighted_delta_balance_final"])

    if case.candidate == "R1":
        classified = bool(
            abs(delta_A) >= params.effect_floor
            and abs(delta_B) >= params.effect_floor
            and sign(delta_A, params.effect_floor) == sign(delta_B, params.effect_floor)
            and abs(distance_change) <= params.effect_floor
            and abs(r_balance) >= params.effect_floor
        )
        expected_class = "reference_common_mode"
    else:
        classified = bool(
            distance_change > params.effect_floor
            and abs(r_balance) <= params.balance_tol
            and float(summary["R_weighted_acceleration_balance_max"]) <= params.balance_tol
        )
        expected_class = "repulsive_distance_growth"

    valid = bool(classified and float(summary["max_delta_chi_std"]) <= params.gauge_std_tol)
    row = {
        "candidate": case.candidate,
        "case": case.case,
        "expected_class": expected_class,
        "reference_sign": case.reference_sign,
        "mirror_positions": case.mirror_positions,
        "c_pressure_sign": case.c_pressure_sign,
        "direct_sign": case.direct_sign,
        "final_delta_A": delta_A,
        "final_delta_B": delta_B,
        "distance_change": distance_change,
        "R_weighted_delta_balance_final": r_balance,
        "R_weighted_acceleration_balance_max": float(summary["R_weighted_acceleration_balance_max"]),
        "max_delta_chi_std": float(summary["max_delta_chi_std"]),
        "delta_A_sign": sign(delta_A, params.effect_floor),
        "delta_B_sign": sign(delta_B, params.effect_floor),
        "distance_change_sign": sign(distance_change, params.effect_floor),
        "classified_as_expected": classified,
        "valid": valid,
    }
    return row, gauge_rows


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波Cにおける反力候補粗計量予備実験 v1",
        "",
        "## 目的",
        "",
        "Stage I の前段として、引力様に見える距離位相縮小と混同し得る反力候補を、`R1`, `R2`, `R3` に分けて粗く計量する。",
        "",
        "本実験は、標準力や標準重力を導入しない。読出し器依存の共通モード、C定常波圧、A-B直接干渉が、それぞれどの符号と桁で現れるかを調べる。",
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
            "| candidate | case | class | A delta | B delta | distance change | R balance | R a balance | valid |",
            "|---|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_rows"]:
        lines.append(
            f"| {row['candidate']} | {row['case']} | {row['expected_class']} | "
            f"{row['final_delta_A']:.16e} | {row['final_delta_B']:.16e} | "
            f"{row['distance_change']:.16e} | {row['R_weighted_delta_balance_final']:.16e} | "
            f"{row['R_weighted_acceleration_balance_max']:.16e} | `{row['valid']}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- `R1` は C 参照に依存した共通モードとして分類する。A/B が同符号に読まれ、距離位相を変えず、R重み付き外部収支を満たさないなら、外部並進候補ではなく読出し器反力または読出し器バイアスとして保留する。",
            "- `R2` は C 定常波圧による離反候補として分類する。距離位相が増大し、R重み付き収支が小さいなら、引力様候補とは別の反力指紋として扱う。",
            "- `R3` は A-B 直接干渉による離反候補として分類する。C媒介を切っても距離位相が増大するなら、C残渣由来の接近候補とは別系列として扱う。",
            "- この予備実験で反力候補の桁が C媒介 persistent 応答と同程度なら、次は Stage II の同時多読出しへ進む必要がある。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_reaction_candidates_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_reaction_candidates_cases_v1.csv` |",
            "| gauge CSV | `abc_baseline_stationary_wave_reaction_candidates_gauge_rows_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_reaction_candidates_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波Cにおける反力候補粗計量予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = ReactionParams()
    gauges = default_gauges()
    cases = [
        ReactionCase("R1", "C_reference_plus", reference_sign=1.0),
        ReactionCase("R1", "C_reference_minus", reference_sign=-1.0),
        ReactionCase("R2", "C_pressure_repulsion", c_pressure_sign=1.0),
        ReactionCase("R2", "C_pressure_repulsion_mirrored", c_pressure_sign=1.0, mirror_positions=True),
        ReactionCase("R3", "AB_direct_repulsion", direct_sign=1.0),
        ReactionCase("R3", "AB_direct_repulsion_mirrored", direct_sign=1.0, mirror_positions=True),
    ]

    case_rows: List[Dict[str, Any]] = []
    gauge_rows: List[Dict[str, Any]] = []
    for case in cases:
        row, rows = run_case(case, params, gauges)
        case_rows.append(row)
        gauge_rows.extend(rows)

    by_candidate: Dict[str, List[Dict[str, Any]]] = {}
    for row in case_rows:
        by_candidate.setdefault(str(row["candidate"]), []).append(row)

    aggregate_verdict = {
        "candidate_count": len(by_candidate),
        "case_count": len(case_rows),
        "R1_valid": bool_all(row["valid"] for row in by_candidate.get("R1", [])),
        "R2_valid": bool_all(row["valid"] for row in by_candidate.get("R2", [])),
        "R3_valid": bool_all(row["valid"] for row in by_candidate.get("R3", [])),
        "all_candidates_valid": bool_all(row["valid"] for row in case_rows),
        "single_gauge_only_used": False,
        "reaction_candidates_preliminary_valid": bool_all(row["valid"] for row in case_rows),
        "max_R_weighted_delta_balance_for_R2_R3": max(
            abs(float(row["R_weighted_delta_balance_final"]))
            for row in case_rows
            if row["candidate"] in {"R2", "R3"}
        ),
        "max_R_weighted_acceleration_balance_for_R2_R3": max(
            abs(float(row["R_weighted_acceleration_balance_max"]))
            for row in case_rows
            if row["candidate"] in {"R2", "R3"}
        ),
        "max_gauge_delta_std": max(abs(float(row["max_delta_chi_std"])) for row in case_rows),
    }

    result = {
        "experiment": "abc_baseline_stationary_wave_reaction_candidates_preliminary_v1",
        "params": asdict(params),
        "case_rows": case_rows,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_reaction_candidates_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_reaction_candidates_cases_v1.csv", case_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_reaction_candidates_gauge_rows_v1.csv", gauge_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
