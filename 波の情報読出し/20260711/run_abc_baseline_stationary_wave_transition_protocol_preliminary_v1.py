from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import Stage2Params


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_transition_protocol_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


FINAL_C_MEMORY = 3.8e-6
ONSET_STEP = 4


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


def ramp_profile(step_count: int, onset: int, final_value: float, duration: int) -> List[float]:
    values: List[float] = []
    for step in range(step_count):
        if step < onset:
            values.append(0.0)
        elif step <= onset + duration:
            values.append(final_value * (step - onset) / duration)
        else:
            values.append(final_value)
    return values


def down_quench_profile(step_count: int, onset: int, initial_value: float) -> List[float]:
    return [initial_value if step <= onset else 0.0 for step in range(step_count)]


def overshoot_relax_profile(
    step_count: int,
    onset: int,
    final_value: float,
    rise_duration: int = 4,
    relax_duration: int = 12,
    overshoot_factor: float = 1.4,
) -> List[float]:
    peak = final_value * overshoot_factor
    peak_step = onset + rise_duration
    settle_step = peak_step + relax_duration
    values: List[float] = []
    for step in range(step_count):
        if step < onset:
            values.append(0.0)
        elif step <= peak_step:
            values.append(peak * (step - onset) / rise_duration)
        elif step <= settle_step:
            frac = (step - peak_step) / relax_duration
            values.append(peak + (final_value - peak) * frac)
        else:
            values.append(final_value)
    return values


def profile_rows(label: str, profile: List[float], params: Stage2Params) -> List[Dict[str, Any]]:
    total = params.R_A + params.R_B
    reduced = params.R_A * params.R_B / total
    cumulative_A = 0.0
    cumulative_B = 0.0
    delta_A_values: List[float] = []
    delta_B_values: List[float] = []
    for c_memory in profile:
        response = params.epsilon_C_return * c_memory
        cumulative_A += response * params.R_B / total
        cumulative_B -= response * params.R_A / total
        delta_A_values.append(cumulative_A)
        delta_B_values.append(cumulative_B)

    rows: List[Dict[str, Any]] = []
    for step in range(1, len(profile) - 1):
        delta_c_forward = profile[step + 1] - profile[step]
        predicted_ra_a = reduced * params.epsilon_C_return * delta_c_forward
        predicted_ra_b = -predicted_ra_a
        actual_ra_a = params.R_A * (delta_A_values[step + 1] - 2.0 * delta_A_values[step] + delta_A_values[step - 1])
        actual_ra_b = params.R_B * (delta_B_values[step + 1] - 2.0 * delta_B_values[step] + delta_B_values[step - 1])
        rows.append(
            {
                "protocol": label,
                "step": step,
                "particle": "A",
                "C_memory": profile[step],
                "delta_C_memory_forward": delta_c_forward,
                "actual_Ra": actual_ra_a,
                "predicted_Ra_from_delta_C": predicted_ra_a,
                "balance_error": actual_ra_a - predicted_ra_a,
            }
        )
        rows.append(
            {
                "protocol": label,
                "step": step,
                "particle": "B",
                "C_memory": profile[step],
                "delta_C_memory_forward": delta_c_forward,
                "actual_Ra": actual_ra_b,
                "predicted_Ra_from_delta_C": predicted_ra_b,
                "balance_error": actual_ra_b - predicted_ra_b,
            }
        )
    return rows


def protocol_profiles(params: Stage2Params) -> Dict[str, Tuple[List[float], Dict[str, Any]]]:
    profiles: Dict[str, Tuple[List[float], Dict[str, Any]]] = {}
    for duration in [1, 2, 4, 8, 16, 32]:
        profiles[f"ramp_duration_{duration}"] = (
            ramp_profile(params.step_count, ONSET_STEP, FINAL_C_MEMORY, duration),
            {"kind": "ramp", "duration": duration},
        )
    profiles["restart_stable_final_C"] = (
        [FINAL_C_MEMORY for _ in range(params.step_count)],
        {"kind": "stable_restart", "duration": 0},
    )
    profiles["down_quench_from_final_C"] = (
        down_quench_profile(params.step_count, ONSET_STEP, FINAL_C_MEMORY),
        {"kind": "down_quench", "duration": 1},
    )
    profiles["overshoot_relax_to_final_C"] = (
        overshoot_relax_profile(params.step_count, ONSET_STEP, FINAL_C_MEMORY),
        {"kind": "overshoot_relax", "duration": 16},
    )
    return profiles


def max_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return max(abs(float(row[key])) for row in rows)


def mean_abs(rows: List[Dict[str, Any]], key: str) -> float:
    if not rows:
        return 0.0
    return float(np.mean([abs(float(row[key])) for row in rows]))


def active_rows(rows: List[Dict[str, Any]], floor: float = 1.0e-18) -> List[Dict[str, Any]]:
    return [row for row in rows if abs(float(row["delta_C_memory_forward"])) > floor]


def post_settle_rows(rows: List[Dict[str, Any]], start_step: int) -> List[Dict[str, Any]]:
    return [row for row in rows if int(row["step"]) >= start_step]


def summarize_protocol(label: str, metadata: Dict[str, Any], profile: List[float], rows: List[Dict[str, Any]], params: Stage2Params) -> Dict[str, Any]:
    active = active_rows(rows)
    settle_start = ONSET_STEP + int(metadata["duration"]) + 3 if metadata["duration"] > 0 else 1
    post = post_settle_rows(rows, settle_start)
    a_rows = [row for row in rows if row["particle"] == "A"]
    positive_steps = sorted({int(row["step"]) for row in a_rows if float(row["actual_Ra"]) > params.effect_floor})
    negative_steps = sorted({int(row["step"]) for row in a_rows if float(row["actual_Ra"]) < -params.effect_floor})
    max_delta_c = max_abs(rows, "delta_C_memory_forward")
    max_actual = max_abs(rows, "actual_Ra")
    slope_gain = max_actual / max_delta_c if max_delta_c > 0.0 else None
    return {
        "protocol": label,
        "kind": metadata["kind"],
        "duration": metadata["duration"],
        "initial_C_memory": profile[0],
        "final_C_memory": profile[-1],
        "max_C_memory": max(abs(value) for value in profile),
        "max_delta_C_memory_abs": max_delta_c,
        "max_actual_Ra_abs": max_actual,
        "mean_active_actual_Ra_abs": mean_abs(active, "actual_Ra"),
        "post_settle_start_step": settle_start,
        "post_settle_mean_actual_Ra_abs": mean_abs(post, "actual_Ra"),
        "post_settle_max_actual_Ra_abs": max_abs(post, "actual_Ra"),
        "max_balance_error_abs": max_abs(rows, "balance_error"),
        "slope_gain": slope_gain,
        "positive_step_count_A": len(positive_steps),
        "negative_step_count_A": len(negative_steps),
        "positive_steps_A": ",".join(str(step) for step in positive_steps),
        "negative_steps_A": ",".join(str(step) for step in negative_steps),
    }


def strictly_decreasing(values: List[float]) -> bool:
    return all(values[idx] > values[idx + 1] for idx in range(len(values) - 1))


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 遷移プロトコル予備実験 v1",
        "",
        "## 目的",
        "",
        "同じ最終 `C_memory` に到達する場合でも、quench, adiabatic ramp, stable restart, overshoot により加速度候補がどう変わるかを調べる。",
        "",
        "準安定傾斜起源なら、加速度候補は最終 `C_memory` の大きさではなく、`ΔC_memory` の大きさと符号に従う。",
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
            "| protocol | kind | duration | max |ΔC| | max |R*a| | post max |R*a| | sign A | max error |",
            "|---|---|---:|---:|---:|---:|---|---:|",
        ]
    )
    for row in result["protocol_summaries"]:
        sign_text = f"+{row['positive_step_count_A']}/-{row['negative_step_count_A']}"
        lines.append(
            f"| {row['protocol']} | {row['kind']} | {row['duration']} | "
            f"{row['max_delta_C_memory_abs']:.16e} | {row['max_actual_Ra_abs']:.16e} | "
            f"{row['post_settle_max_actual_Ra_abs']:.16e} | {sign_text} | {row['max_balance_error_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 同じ最終 `C_memory` でも ramp が遅いほど `|R*a|` が小さくなるなら、加速度候補は定常レベルではなく遷移速度に支配される。",
            "- stable restart で非ゼロ `C_memory` があるにもかかわらず加速度候補が消えるなら、安定状態そのものは力候補を生まない。",
            "- down quench で符号が反転し、overshoot relax で正負の両符号が出るなら、候補はポテンシャルの絶対量ではなく、準安定傾斜の向きに従う。",
            "- これは標準重力の導出ではなく、引力的読出し候補が定常場ではなく遷移プロトコルに依存するかを調べる予備試験である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_transition_protocol_preliminary_result_v1.json` |",
            "| protocol CSV | `abc_baseline_stationary_wave_transition_protocol_cases_v1.csv` |",
            "| rows CSV | `abc_baseline_stationary_wave_transition_protocol_rows_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_transition_protocol_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 遷移プロトコル予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params()
    profiles = protocol_profiles(params)
    all_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for label, (profile, metadata) in profiles.items():
        rows = profile_rows(label, profile, params)
        all_rows.extend(rows)
        summaries.append(summarize_protocol(label, metadata, profile, rows, params))

    by_protocol = {row["protocol"]: row for row in summaries}
    ramp_rows = [row for row in summaries if row["kind"] == "ramp"]
    ramp_rows.sort(key=lambda row: int(row["duration"]))
    ramp_max_values = [float(row["max_actual_Ra_abs"]) for row in ramp_rows]
    quench = by_protocol["ramp_duration_1"]
    slowest = by_protocol["ramp_duration_32"]
    restart = by_protocol["restart_stable_final_C"]
    down = by_protocol["down_quench_from_final_C"]
    overshoot = by_protocol["overshoot_relax_to_final_C"]

    aggregate_verdict = {
        "protocol_count": len(summaries),
        "single_gauge_only_used": False,
        "final_C_memory_target": FINAL_C_MEMORY,
        "all_protocols_same_final_C_except_down_quench": bool_all(
            abs(float(row["final_C_memory"]) - FINAL_C_MEMORY) <= 1.0e-18
            for row in summaries
            if row["kind"] != "down_quench"
        ),
        "ramp_max_Ra_strictly_decreases_with_duration": strictly_decreasing(ramp_max_values),
        "quench_to_slowest_max_Ra_ratio": float(quench["max_actual_Ra_abs"] / slowest["max_actual_Ra_abs"]),
        "restart_stable_acceleration_near_zero": bool(float(restart["max_actual_Ra_abs"]) <= params.effect_floor),
        "down_quench_negative_only_A": bool(
            int(down["positive_step_count_A"]) == 0 and int(down["negative_step_count_A"]) > 0
        ),
        "overshoot_has_positive_and_negative_A": bool(
            int(overshoot["positive_step_count_A"]) > 0 and int(overshoot["negative_step_count_A"]) > 0
        ),
        "post_settle_acceleration_near_zero_all_protocols": bool_all(
            float(row["post_settle_max_actual_Ra_abs"]) <= params.effect_floor for row in summaries
        ),
        "max_balance_error_abs": max(float(row["max_balance_error_abs"]) for row in summaries),
        "transition_protocol_preliminary_valid": bool(
            strictly_decreasing(ramp_max_values)
            and float(quench["max_actual_Ra_abs"] / slowest["max_actual_Ra_abs"]) >= 31.0
            and float(restart["max_actual_Ra_abs"]) <= params.effect_floor
            and int(down["positive_step_count_A"]) == 0
            and int(down["negative_step_count_A"]) > 0
            and int(overshoot["positive_step_count_A"]) > 0
            and int(overshoot["negative_step_count_A"]) > 0
            and bool_all(float(row["post_settle_max_actual_Ra_abs"]) <= params.effect_floor for row in summaries)
        ),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_transition_protocol_preliminary_v1",
        "params": asdict(params),
        "protocol_summaries": summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_transition_protocol_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_transition_protocol_cases_v1.csv", summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_transition_protocol_rows_v1.csv", all_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
