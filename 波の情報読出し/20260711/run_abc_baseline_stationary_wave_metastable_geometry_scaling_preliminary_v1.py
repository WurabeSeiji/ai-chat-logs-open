from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np

from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import Stage2Params


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_metastable_geometry_scaling_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

RAMP_DURATION = 8


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


def reduced_R(params: Stage2Params) -> float:
    return float(params.R_A * params.R_B / (params.R_A + params.R_B))


def steady_c_memory_from_distance(distance: float, params: Stage2Params) -> float:
    source = params.epsilon_to_C * reduced_R(params) * math.sin(distance)
    return float(source / (1.0 - params.memory_decay_C))


def transition_profile(final_value: float, step_count: int, duration: int = RAMP_DURATION) -> List[float]:
    values: List[float] = []
    for step in range(step_count):
        if step <= duration:
            values.append(final_value * step / duration)
        else:
            values.append(final_value)
    return values


def profile_rows(label: str, distance: float, params: Stage2Params) -> List[Dict[str, Any]]:
    final_c = steady_c_memory_from_distance(distance, params)
    profile = transition_profile(final_c, params.step_count)
    total = params.R_A + params.R_B
    reduced = reduced_R(params)
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
    for step in range(1, params.step_count - 1):
        delta_c = profile[step + 1] - profile[step]
        predicted_a = reduced * params.epsilon_C_return * delta_c
        actual_a = params.R_A * (delta_A_values[step + 1] - 2.0 * delta_A_values[step] + delta_A_values[step - 1])
        actual_b = params.R_B * (delta_B_values[step + 1] - 2.0 * delta_B_values[step] + delta_B_values[step - 1])
        rows.append(
            {
                "case": label,
                "step": step,
                "particle": "A",
                "distance": distance,
                "sin_distance": math.sin(distance),
                "R_A": params.R_A,
                "R_B": params.R_B,
                "reduced_R": reduced,
                "C_memory": profile[step],
                "delta_C_memory_forward": delta_c,
                "actual_Ra": actual_a,
                "predicted_Ra": predicted_a,
                "balance_error": actual_a - predicted_a,
            }
        )
        rows.append(
            {
                "case": label,
                "step": step,
                "particle": "B",
                "distance": distance,
                "sin_distance": math.sin(distance),
                "R_A": params.R_A,
                "R_B": params.R_B,
                "reduced_R": reduced,
                "C_memory": profile[step],
                "delta_C_memory_forward": delta_c,
                "actual_Ra": actual_b,
                "predicted_Ra": -predicted_a,
                "balance_error": actual_b + predicted_a,
            }
        )
    return rows


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


def summarize_case(label: str, family: str, distance: float, params: Stage2Params, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    a_rows = [row for row in rows if row["particle"] == "A"]
    b_rows = [row for row in rows if row["particle"] == "B"]
    max_a = max_abs(a_rows, "actual_Ra")
    max_b = max_abs(b_rows, "actual_Ra")
    signed_a = max((float(row["actual_Ra"]) for row in a_rows), key=abs)
    signed_b = max((float(row["actual_Ra"]) for row in b_rows), key=abs)
    reduced = reduced_R(params)
    expected_scale = (reduced**2) * abs(math.sin(distance))
    signed_expected_scale = (reduced**2) * math.sin(distance)
    normalized_scale = max_a / expected_scale if expected_scale > 0.0 else None
    max_balance = max_abs(rows, "balance_error")
    return {
        "case": label,
        "family": family,
        "distance": distance,
        "sin_distance": math.sin(distance),
        "R_A": params.R_A,
        "R_B": params.R_B,
        "R_ratio_B_over_A": params.R_B / params.R_A,
        "reduced_R": reduced,
        "max_A_Ra_abs": max_a,
        "max_B_Ra_abs": max_b,
        "max_A_Ra_signed": signed_a,
        "max_B_Ra_signed": signed_b,
        "expected_scale_reduced2_abs_sin": expected_scale,
        "expected_scale_reduced2_signed_sin": signed_expected_scale,
        "normalized_scale": normalized_scale,
        "A_sign_matches_sin_distance": bool(sign(signed_a) == sign(math.sin(distance))),
        "B_opposes_A": bool(sign(signed_b) == -sign(signed_a)),
        "max_balance_error_abs": max_balance,
    }


def distance_cases(base: Stage2Params) -> List[tuple[str, float, Stage2Params]]:
    distances = [-1.0, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0]
    return [(f"distance_{distance:+.2f}", distance, base) for distance in distances]


def r_ratio_cases(base: Stage2Params) -> List[tuple[str, float, Stage2Params]]:
    ratios = [0.25, 0.5, 1.0, 1.5625, 2.0, 4.0, 8.0]
    distance = 0.5
    return [
        (
            f"R_ratio_{ratio:.4g}",
            distance,
            replace(base, R_A=1.0, R_B=ratio),
        )
        for ratio in ratios
    ]


def slope(values: List[float], targets: List[float]) -> float:
    x = np.array(values, dtype=float)
    y = np.array(targets, dtype=float)
    denom = float(np.sum(x**2))
    if denom == 0.0:
        return 0.0
    return float(np.sum(x * y) / denom)


def max_relative_error(values: List[float], targets: List[float], gain: float) -> float:
    errors: List[float] = []
    for value, target in zip(values, targets):
        denom = max(abs(target * gain), 1.0e-30)
        errors.append(abs(value - gain * target) / denom)
    return max(errors) if errors else 0.0


def summarize_family(family: str, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    selected = [row for row in rows if row["family"] == family]
    actual = [float(row["max_A_Ra_abs"]) for row in selected]
    expected = [float(row["expected_scale_reduced2_abs_sin"]) for row in selected]
    gain = slope(expected, actual)
    normalized = [float(row["normalized_scale"]) for row in selected if row["normalized_scale"] is not None]
    return {
        "family": family,
        "case_count": len(selected),
        "actual_expected_corr": safe_corr(actual, expected),
        "scale_gain": gain,
        "max_relative_error_after_gain": max_relative_error(actual, expected, gain),
        "normalized_scale_mean": float(np.mean(normalized)) if normalized else None,
        "normalized_scale_std": float(np.std(normalized)) if normalized else None,
        "all_A_sign_matches_sin_distance": bool_all(row["A_sign_matches_sin_distance"] for row in selected),
        "all_B_opposes_A": bool_all(row["B_opposes_A"] for row in selected),
        "max_balance_error_abs": max(float(row["max_balance_error_abs"]) for row in selected),
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 準安定傾斜の幾何スケーリング予備実験 v1",
        "",
        "## 目的",
        "",
        "準安定傾斜から生じる加速度候補が、距離位相 `sin(Δχ)` と R重みに従うかを調べる。",
        "",
        "ここでは定常場の大きさではなく、同じ ramp duration の `ΔC_memory` が、距離位相と R重みによりどう変わるかを読む。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## family summary",
            "",
            "| family | cases | corr(actual, expected) | rel error | norm std | signs | balance |",
            "|---|---:|---:|---:|---:|---|---:|",
        ]
    )
    for row in result["family_summaries"]:
        signs = bool(row["all_A_sign_matches_sin_distance"]) and bool(row["all_B_opposes_A"])
        lines.append(
            f"| {row['family']} | {row['case_count']} | {row['actual_expected_corr']} | "
            f"{row['max_relative_error_after_gain']:.16e} | {row['normalized_scale_std']:.16e} | "
            f"`{signs}` | {row['max_balance_error_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## case summary",
            "",
            "| case | family | distance | R_B/R_A | max |R*a| | expected | normalized | sign ok |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        sign_ok = bool(row["A_sign_matches_sin_distance"]) and bool(row["B_opposes_A"])
        lines.append(
            f"| {row['case']} | {row['family']} | {row['distance']:.6f} | {row['R_ratio_B_over_A']:.6f} | "
            f"{row['max_A_Ra_abs']:.16e} | {row['expected_scale_reduced2_abs_sin']:.16e} | "
            f"{row['normalized_scale']:.16e} | `{sign_ok}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 距離位相スイープで `|R*a|` が `|sin(Δχ)|` に比例するなら、準安定傾斜には幾何的な距離位相依存がある。",
            "- R比スイープで `|R*a|` が `R_red^2 |sin(Δχ)|` に比例するなら、Cへの書き込みとCからの戻りの両方に R重みが入っている。",
            "- `A` と `B` の符号が反対で、`A` の符号が `sin(Δχ)` と対応するなら、ラベル固定の押し込みではなく相対位相の向きに従っている。",
            "- これは標準重力の距離法則ではない。準安定傾斜候補に距離位相・R重みのスケーリングがあるかを見る予備試験である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_metastable_geometry_scaling_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_metastable_geometry_scaling_cases_v1.csv` |",
            "| rows CSV | `abc_baseline_stationary_wave_metastable_geometry_scaling_rows_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_geometry_scaling_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 準安定傾斜の幾何スケーリング予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    base = Stage2Params()
    all_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []

    for family, cases in [
        ("distance_phase", distance_cases(base)),
        ("R_weight", r_ratio_cases(base)),
    ]:
        for label, distance, params in cases:
            rows = profile_rows(label, distance, params)
            all_rows.extend(rows)
            summaries.append(summarize_case(label, family, distance, params, rows))

    family_summaries = [summarize_family("distance_phase", summaries), summarize_family("R_weight", summaries)]
    aggregate_verdict = {
        "case_count": len(summaries),
        "family_count": len(family_summaries),
        "single_gauge_only_used": False,
        "ramp_duration": RAMP_DURATION,
        "all_family_correlations_near_one": bool_all(
            float(row["actual_expected_corr"]) >= 0.999999999999 for row in family_summaries
        ),
        "all_family_relative_errors_small": bool_all(
            float(row["max_relative_error_after_gain"]) <= 1.0e-12 for row in family_summaries
        ),
        "all_family_signs_valid": bool_all(
            bool(row["all_A_sign_matches_sin_distance"]) and bool(row["all_B_opposes_A"]) for row in family_summaries
        ),
        "max_balance_error_abs": max(float(row["max_balance_error_abs"]) for row in family_summaries),
        "metastable_geometry_scaling_preliminary_valid": bool(
            bool_all(float(row["actual_expected_corr"]) >= 0.999999999999 for row in family_summaries)
            and bool_all(float(row["max_relative_error_after_gain"]) <= 1.0e-12 for row in family_summaries)
            and bool_all(
                bool(row["all_A_sign_matches_sin_distance"]) and bool(row["all_B_opposes_A"])
                for row in family_summaries
            )
        ),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_metastable_geometry_scaling_preliminary_v1",
        "params": asdict(base),
        "case_summaries": summaries,
        "family_summaries": family_summaries,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_geometry_scaling_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_geometry_scaling_cases_v1.csv", summaries)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_geometry_scaling_rows_v1.csv", all_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
