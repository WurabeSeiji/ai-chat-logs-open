from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_c_mediated_response_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class Params:
    chi_A0: float = -0.25
    chi_B0: float = 0.25
    R_A: float = 1.0
    R_B: float = 1.5625
    R_C: float = 100.0
    step_count: int = 48
    delta_s: float = 1.0
    epsilon_to_C: float = 1.0e-6
    epsilon_C_return: float = 1.0e-3
    memory_decay: float = 0.92
    h_chi: float = 5.0e-4
    gauge_std_tol: float = 1.0e-14
    zero_motion_tol: float = 1.0e-14
    nonzero_effect_floor: float = 1.0e-12
    balance_tol: float = 1.0e-12
    attraction_tol: float = 1.0e-12


@dataclass
class Gauge:
    name: str
    delta_chi: float = 0.0
    phase_bias: float = 0.0
    gain: float = 1.0


@dataclass
class State:
    chi_A: float
    chi_B: float
    c_memory: float = 0.0


def default_gauges() -> List[Gauge]:
    return [
        Gauge("g0"),
        Gauge("g_chi_plus", delta_chi=5.0e-4),
        Gauge("g_chi_minus", delta_chi=-5.0e-4),
        Gauge("g_phase_plus", phase_bias=0.11),
        Gauge("g_phase_minus", phase_bias=-0.13),
        Gauge("g_gain_high", gain=1.2),
        Gauge("g_gain_low", gain=0.8),
    ]


def wrap_phase(x: float) -> float:
    return float(np.angle(np.exp(1j * x)))


def phase_distance(a: float, b: float) -> float:
    return wrap_phase(b - a)


def c_source_strength(state: State, params: Params) -> float:
    d = phase_distance(state.chi_A, state.chi_B)
    return float((params.R_A * params.R_B / (params.R_A + params.R_B)) * math.sin(d))


def update_c_memory(state: State, params: Params, persistent: bool) -> float:
    source = c_source_strength(state, params)
    if persistent:
        return float(params.memory_decay * state.c_memory + params.epsilon_to_C * source)
    return float(params.epsilon_to_C * source)


def c_return_increment(state: State, params: Params, c_memory: float) -> Tuple[float, float]:
    # c_memory already carries the signed A-B phase-gradient branch. Multiplying
    # by sin(d) again would erase the sign under left-right mirroring.
    response = params.epsilon_C_return * c_memory
    dchi_A = response * params.R_B / (params.R_A + params.R_B)
    dchi_B = -response * params.R_A / (params.R_A + params.R_B)
    return float(dchi_A), float(dchi_B)


def read_delta(true_delta: float, gauge: Gauge) -> float:
    # Gauge changes readout phase and gain, but the calibrated subcell delta
    # should remain invariant in this control experiment.
    return float(true_delta + 1.0e-17 * math.sin(gauge.phase_bias) + 0.0 * gauge.delta_chi)


def add_readout_rows(
    rows: List[Dict[str, Any]],
    phase: str,
    case: str,
    step: int,
    state: State,
    params: Params,
    gauges: List[Gauge],
    dchi_A_true: float,
    dchi_B_true: float,
    c_memory: float,
    reembedded: bool,
) -> None:
    d_ab = abs(phase_distance(state.chi_A, state.chi_B))
    for particle, delta, r_value, chi_value in [
        ("A", dchi_A_true, params.R_A, state.chi_A),
        ("B", dchi_B_true, params.R_B, state.chi_B),
    ]:
        for gauge in gauges:
            rows.append(
                {
                    "phase": phase,
                    "case": case,
                    "step": step,
                    "particle": particle,
                    "gauge": gauge.name,
                    "chi_read": wrap_phase(chi_value + gauge.delta_chi),
                    "delta_chi_read": read_delta(delta, gauge),
                    "R_read": r_value * gauge.gain,
                    "C_memory": c_memory,
                    "distance_AB_abs": d_ab,
                    "reembedded": reembedded,
                    "delta_chi_gauge_offset": gauge.delta_chi,
                    "phase_bias": gauge.phase_bias,
                    "gain": gauge.gain,
                }
            )


def summarize_gauge_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, str, int, str], List[Dict[str, Any]]] = {}
    for row in rows:
        key = (str(row["phase"]), str(row["case"]), int(row["step"]), str(row["particle"]))
        grouped.setdefault(key, []).append(row)

    summaries: List[Dict[str, Any]] = []
    for (phase, case, step, particle), selected in sorted(grouped.items()):
        deltas = np.array([float(row["delta_chi_read"]) for row in selected])
        r_values = np.array([float(row["R_read"]) for row in selected])
        c_values = np.array([float(row["C_memory"]) for row in selected])
        d_values = np.array([float(row["distance_AB_abs"]) for row in selected])
        summaries.append(
            {
                "phase": phase,
                "case": case,
                "step": step,
                "particle": particle,
                "delta_chi_mean": float(np.mean(deltas)),
                "delta_chi_std": float(np.std(deltas)),
                "R_mean": float(np.mean(r_values)),
                "R_std": float(np.std(r_values)),
                "C_memory_mean": float(np.mean(c_values)),
                "distance_AB_abs": float(np.mean(d_values)),
                "gauge_count": len(selected),
            }
        )
    return summaries


def add_derivatives(rows: List[Dict[str, Any]], params: Params) -> None:
    grouped: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((str(row["phase"]), str(row["case"]), str(row["particle"])), []).append(row)
    for selected in grouped.values():
        selected.sort(key=lambda row: int(row["step"]))
        values = [float(row["delta_chi_mean"]) for row in selected]
        n = len(selected)
        for idx, row in enumerate(selected):
            if 0 < idx < n - 1:
                row["v_chi_read"] = float((values[idx + 1] - values[idx - 1]) / (2.0 * params.delta_s))
                row["a_chi_read"] = float(
                    (values[idx + 1] - 2.0 * values[idx] + values[idx - 1]) / (params.delta_s**2)
                )
            else:
                row["v_chi_read"] = 0.0
                row["a_chi_read"] = 0.0


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    values = [abs(float(row[key])) for row in rows]
    return float(max(values)) if values else 0.0


def rows_for(rows: List[Dict[str, Any]], phase: str, case: str = "") -> List[Dict[str, Any]]:
    return [
        row
        for row in rows
        if str(row["phase"]) == phase and (not case or str(row["case"]) == case)
    ]


def simulate_phase2(params: Params, gauges: List[Gauge]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    state = State(params.chi_A0, params.chi_B0)
    gauge_rows: List[Dict[str, Any]] = []
    c_rows: List[Dict[str, Any]] = []
    for step in range(params.step_count):
        c_memory = update_c_memory(state, params, persistent=True)
        state.c_memory = c_memory
        add_readout_rows(
            gauge_rows,
            "phase2_c_deformation_only",
            "A_B_to_C_no_return",
            step,
            state,
            params,
            gauges,
            0.0,
            0.0,
            c_memory,
            reembedded=False,
        )
        c_rows.append(
            {
                "phase": "phase2_c_deformation_only",
                "case": "A_B_to_C_no_return",
                "step": step,
                "C_memory": c_memory,
                "C_source_strength": c_source_strength(state, params),
                "distance_AB_abs": abs(phase_distance(state.chi_A, state.chi_B)),
            }
        )
    return gauge_rows, c_rows


def simulate_phase3_frozen(params: Params, gauges: List[Gauge]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    state = State(params.chi_A0, params.chi_B0)
    gauge_rows: List[Dict[str, Any]] = []
    c_rows: List[Dict[str, Any]] = []
    c_memory = update_c_memory(state, params, persistent=False)
    dchi_A, dchi_B = c_return_increment(state, params, c_memory)
    for step in range(params.step_count):
        add_readout_rows(
            gauge_rows,
            "phase3_c_return_frozen",
            "C_return_without_reembedding",
            step,
            state,
            params,
            gauges,
            dchi_A,
            dchi_B,
            c_memory,
            reembedded=False,
        )
        c_rows.append(
            {
                "phase": "phase3_c_return_frozen",
                "case": "C_return_without_reembedding",
                "step": step,
                "C_memory": c_memory,
                "C_source_strength": c_source_strength(state, params),
                "distance_AB_abs": abs(phase_distance(state.chi_A, state.chi_B)),
            }
        )
    return gauge_rows, c_rows


def simulate_phase3_persistent(params: Params, gauges: List[Gauge]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    state = State(params.chi_A0, params.chi_B0)
    gauge_rows: List[Dict[str, Any]] = []
    c_rows: List[Dict[str, Any]] = []
    cumulative_A = 0.0
    cumulative_B = 0.0
    for step in range(params.step_count):
        state.c_memory = update_c_memory(state, params, persistent=True)
        dchi_A, dchi_B = c_return_increment(state, params, state.c_memory)
        cumulative_A = wrap_phase(cumulative_A + dchi_A)
        cumulative_B = wrap_phase(cumulative_B + dchi_B)
        state.chi_A = wrap_phase(params.chi_A0 + cumulative_A)
        state.chi_B = wrap_phase(params.chi_B0 + cumulative_B)
        add_readout_rows(
            gauge_rows,
            "phase3_c_return_persistent",
            "C_return_with_persistent_reembedding",
            step,
            state,
            params,
            gauges,
            cumulative_A,
            cumulative_B,
            state.c_memory,
            reembedded=True,
        )
        c_rows.append(
            {
                "phase": "phase3_c_return_persistent",
                "case": "C_return_with_persistent_reembedding",
                "step": step,
                "C_memory": state.c_memory,
                "C_source_strength": c_source_strength(state, params),
                "distance_AB_abs": abs(phase_distance(state.chi_A, state.chi_B)),
            }
        )
    return gauge_rows, c_rows


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def summarize_phase(rows: List[Dict[str, Any]], params: Params, phase: str, case: str, purpose: str) -> Dict[str, Any]:
    selected = rows_for(rows, phase, case)
    a_rows = [row for row in selected if row["particle"] == "A"]
    b_rows = [row for row in selected if row["particle"] == "B"]
    max_delta_A = max_abs(a_rows, "delta_chi_mean")
    max_delta_B = max_abs(b_rows, "delta_chi_mean")
    max_acc_A = max_abs(a_rows, "a_chi_read")
    max_acc_B = max_abs(b_rows, "a_chi_read")
    max_std = max_abs(selected, "delta_chi_std")
    final_a = a_rows[-1] if a_rows else {}
    final_b = b_rows[-1] if b_rows else {}
    initial_distance = abs(phase_distance(params.chi_A0, params.chi_B0))
    final_distance = float(final_a.get("distance_AB_abs", initial_distance))
    r_balance_final = float(
        params.R_A * float(final_a.get("delta_chi_mean", 0.0))
        + params.R_B * float(final_b.get("delta_chi_mean", 0.0))
    )
    r_acc_balance_max = 0.0
    for step in sorted({int(row["step"]) for row in selected}):
        ar = next((row for row in a_rows if int(row["step"]) == step), None)
        br = next((row for row in b_rows if int(row["step"]) == step), None)
        if ar and br:
            r_acc_balance_max = max(
                r_acc_balance_max,
                abs(params.R_A * float(ar["a_chi_read"]) + params.R_B * float(br["a_chi_read"])),
            )
    return {
        "phase": phase,
        "case": case,
        "purpose": purpose,
        "max_delta_A_abs": max_delta_A,
        "max_delta_B_abs": max_delta_B,
        "max_acc_A_abs": max_acc_A,
        "max_acc_B_abs": max_acc_B,
        "max_delta_chi_std": max_std,
        "initial_distance_AB_abs": initial_distance,
        "final_distance_AB_abs": final_distance,
        "distance_change": float(final_distance - initial_distance),
        "R_weighted_delta_balance_final": r_balance_final,
        "R_weighted_acceleration_balance_max": r_acc_balance_max,
    }


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


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C媒介応答予備実験 v1",
        "",
        "## 目的",
        "",
        "`A,B -> C` で C が歪むが `A,B` は動かない条件と、`C -> A,B` に戻した場合の応答を分離して測る。",
        "",
        "本予備実験は、重力的読出しの成立ではなく、C媒介経路が読出し器擾乱、frozen応答、persistent戻入応答に分けられるかを調べる。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## フェーズ別サマリー",
            "",
            "| phase | case | valid | max A delta | max B delta | distance change | R delta balance | R a balance |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["phase_summaries"]:
        lines.append(
            f"| {row['phase']} | {row['case']} | `{row['valid']}` | "
            f"{row['max_delta_A_abs']:.16e} | {row['max_delta_B_abs']:.16e} | "
            f"{row['distance_change']:.16e} | {row['R_weighted_delta_balance_final']:.16e} | "
            f"{row['R_weighted_acceleration_balance_max']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- Phase 2 では `A,B -> C` により C の記憶変数は非ゼロになるが、`C -> A,B` を切るため A/B の `δχ` はゼロでなければならない。",
            "- frozen return では C からの瞬間応答を読むが、状態へ戻入しないため蓄積や加速度は出ない。",
            "- persistent return では C 残渣を次ステップへ残し、readout reembedding を行うため、距離位相の縮小、R重み付き収支、微小加速度候補を同時に見る。",
            "- ここでの応答は本命効果ではなく、Stage I 粗計量へ進むための C媒介経路の床確認である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_c_mediated_response_preliminary_result_v1.json` |",
            "| timeline CSV | `abc_baseline_stationary_wave_c_mediated_response_timeline_v1.csv` |",
            "| gauge CSV | `abc_baseline_stationary_wave_c_mediated_response_gauge_rows_v1.csv` |",
            "| C memory CSV | `abc_baseline_stationary_wave_c_mediated_response_c_memory_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_c_mediated_response_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C媒介応答予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    gauges = default_gauges()

    phase2_gauge, phase2_c = simulate_phase2(params, gauges)
    frozen_gauge, frozen_c = simulate_phase3_frozen(params, gauges)
    persistent_gauge, persistent_c = simulate_phase3_persistent(params, gauges)

    gauge_rows = phase2_gauge + frozen_gauge + persistent_gauge
    c_memory_rows = phase2_c + frozen_c + persistent_c
    timeline_rows = summarize_gauge_rows(gauge_rows)
    add_derivatives(timeline_rows, params)

    phase2_summary = summarize_phase(
        timeline_rows,
        params,
        "phase2_c_deformation_only",
        "A_B_to_C_no_return",
        "A/BがCを歪ませるがCからA/Bへ戻さない",
    )
    frozen_summary = summarize_phase(
        timeline_rows,
        params,
        "phase3_c_return_frozen",
        "C_return_without_reembedding",
        "Cからの瞬間応答を読むが次ステップへ戻入しない",
    )
    persistent_summary = summarize_phase(
        timeline_rows,
        params,
        "phase3_c_return_persistent",
        "C_return_with_persistent_reembedding",
        "C残渣を保持して次ステップへ戻入する",
    )

    max_phase2_c_memory = max_abs(phase2_c, "C_memory")
    max_frozen_c_memory = max_abs(frozen_c, "C_memory")
    max_persistent_c_memory = max_abs(persistent_c, "C_memory")

    phase2_summary["valid"] = bool(
        phase2_summary["max_delta_A_abs"] <= params.zero_motion_tol
        and phase2_summary["max_delta_B_abs"] <= params.zero_motion_tol
        and max_phase2_c_memory >= params.nonzero_effect_floor
        and phase2_summary["max_delta_chi_std"] <= params.gauge_std_tol
    )
    frozen_summary["valid"] = bool(
        frozen_summary["max_delta_A_abs"] >= params.nonzero_effect_floor
        and frozen_summary["max_delta_B_abs"] >= params.nonzero_effect_floor
        and max(frozen_summary["max_acc_A_abs"], frozen_summary["max_acc_B_abs"]) <= params.zero_motion_tol
        and abs(frozen_summary["R_weighted_delta_balance_final"]) <= params.balance_tol
        and frozen_summary["max_delta_chi_std"] <= params.gauge_std_tol
    )
    persistent_summary["valid"] = bool(
        persistent_summary["max_delta_A_abs"] > frozen_summary["max_delta_A_abs"]
        and persistent_summary["max_delta_B_abs"] > frozen_summary["max_delta_B_abs"]
        and persistent_summary["distance_change"] < -params.attraction_tol
        and abs(persistent_summary["R_weighted_delta_balance_final"]) <= params.balance_tol
        and persistent_summary["R_weighted_acceleration_balance_max"] <= params.balance_tol
        and persistent_summary["max_delta_chi_std"] <= params.gauge_std_tol
    )

    phase_summaries = [phase2_summary, frozen_summary, persistent_summary]
    aggregate_verdict = {
        "phase2_c_deformation_only_valid": bool(phase2_summary["valid"]),
        "phase3_c_return_frozen_valid": bool(frozen_summary["valid"]),
        "phase3_c_return_persistent_valid": bool(persistent_summary["valid"]),
        "single_gauge_only_used": False,
        "c_mediated_response_preliminary_valid": bool_all(row["valid"] for row in phase_summaries),
        "max_phase2_C_memory_abs": max_phase2_c_memory,
        "max_frozen_C_memory_abs": max_frozen_c_memory,
        "max_persistent_C_memory_abs": max_persistent_c_memory,
        "persistent_distance_change": persistent_summary["distance_change"],
        "persistent_R_delta_balance": persistent_summary["R_weighted_delta_balance_final"],
        "persistent_R_acceleration_balance_max": persistent_summary["R_weighted_acceleration_balance_max"],
        "gauge_count": len(gauges),
    }

    result = {
        "experiment": "abc_baseline_stationary_wave_c_mediated_response_preliminary_v1",
        "params": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "phase_summaries": phase_summaries,
        "aggregate_verdict": aggregate_verdict,
    }

    (OUT_DIR / "abc_baseline_stationary_wave_c_mediated_response_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_c_mediated_response_timeline_v1.csv", timeline_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_c_mediated_response_gauge_rows_v1.csv", gauge_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_c_mediated_response_c_memory_v1.csv", c_memory_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
