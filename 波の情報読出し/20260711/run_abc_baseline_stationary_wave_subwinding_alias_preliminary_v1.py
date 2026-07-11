from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_subwinding_alias_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class AliasCase:
    name: str
    winding_count: int
    residual_phase_per_step: float
    step_count: int = 192


@dataclass
class Gauge:
    name: str
    shutter_phase_offset: float = 0.0
    readout_bias: float = 0.0


def default_gauges() -> List[Gauge]:
    return [
        Gauge("g0"),
        Gauge("g_offset_plus", shutter_phase_offset=0.17),
        Gauge("g_offset_minus", shutter_phase_offset=-0.19),
        Gauge("g_bias_plus", readout_bias=1.0e-17),
        Gauge("g_bias_minus", readout_bias=-1.0e-17),
    ]


def wrap_phase(x: float) -> float:
    return float(np.angle(np.exp(1j * x)))


def sign(value: float, floor: float = 1.0e-14) -> float:
    if value > floor:
        return 1.0
    if value < -floor:
        return -1.0
    return 0.0


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


def linear_slope(values: np.ndarray) -> float:
    x = np.arange(values.size, dtype=float)
    x_mean = float(np.mean(x))
    y_mean = float(np.mean(values))
    denom = float(np.sum((x - x_mean) ** 2))
    if denom == 0.0:
        return 0.0
    return float(np.sum((x - x_mean) * (values - y_mean)) / denom)


def run_case(case: AliasCase) -> Dict[str, Any]:
    gauges = default_gauges()
    gauge_rows: List[Dict[str, Any]] = []
    slopes: List[float] = []
    true_step_phase = 2.0 * math.pi * case.winding_count + case.residual_phase_per_step
    true_velocity_sign = sign(true_step_phase)
    residual_sign = sign(case.residual_phase_per_step)

    for gauge in gauges:
        principal_values = []
        for step in range(case.step_count):
            true_phase = step * true_step_phase
            observed = wrap_phase(true_phase + gauge.shutter_phase_offset) + gauge.readout_bias
            principal_values.append(wrap_phase(observed))
            gauge_rows.append(
                {
                    "case": case.name,
                    "gauge": gauge.name,
                    "step": step,
                    "winding_count": case.winding_count,
                    "residual_phase_per_step": case.residual_phase_per_step,
                    "true_step_phase": true_step_phase,
                    "principal_phase": principal_values[-1],
                    "shutter_phase_offset": gauge.shutter_phase_offset,
                    "readout_bias": gauge.readout_bias,
                }
            )
        unwrapped = np.unwrap(np.array(principal_values, dtype=float))
        slopes.append(linear_slope(unwrapped))

    slope_mean = float(np.mean(slopes))
    slope_std = float(np.std(slopes))
    beat_sign = sign(slope_mean)
    apparent_reverse = bool(true_velocity_sign > 0.0 and beat_sign < 0.0)
    expected_reverse = bool(case.winding_count > 0 and residual_sign < 0.0)
    valid = bool(
        beat_sign == residual_sign
        and apparent_reverse == expected_reverse
        and abs(slope_mean - case.residual_phase_per_step) <= 1.0e-12
        and slope_std <= 1.0e-14
    )

    return {
        **asdict(case),
        "true_step_phase": true_step_phase,
        "true_velocity_sign": true_velocity_sign,
        "expected_beat_sign": residual_sign,
        "observed_beat_slope_mean": slope_mean,
        "observed_beat_slope_std": slope_std,
        "observed_beat_sign": beat_sign,
        "expected_apparent_reverse": expected_reverse,
        "observed_apparent_reverse": apparent_reverse,
        "gauge_count": len(gauges),
        "valid": valid,
        "_gauge_rows": gauge_rows,
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波Cにおけるサブ巻数ビート alias 予備実験 v1",
        "",
        "## 目的",
        "",
        "整数巻数に近い進行位相を、主値位相と unwrap 位相で読むと、巻数以下の残差が巨視的な逆向きビートとして見えるかを確認する。",
        "",
        "これは G3 候補、すなわちサブ巻数ビート/alias による見かけの接近・反転を、C媒介の引力様応答から分離するための予備実験である。",
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
            "| case | winding | residual | true step phase | beat slope | expected reverse | observed reverse | valid |",
            "|---|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for row in result["case_rows"]:
        lines.append(
            f"| {row['name']} | {row['winding_count']} | {row['residual_phase_per_step']:.16e} | "
            f"{row['true_step_phase']:.16e} | {row['observed_beat_slope_mean']:.16e} | "
            f"`{row['expected_apparent_reverse']}` | `{row['observed_apparent_reverse']}` | `{row['valid']}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 真の進行位相は全ケースで正方向である。",
            "- しかし整数巻数よりわずかに小さい場合、主値位相で見た beat は負方向へ進む。",
            "- したがって、後続の加速度実験で距離位相縮小が見えた場合、unwrap 系列と beat 系列を分離しないと、引力様ドリフトと alias 反転を混同する。",
            "- 本予備実験は重力的効果を主張しない。G3 を独立した観測方法依存候補として保持するための基準実験である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_subwinding_alias_preliminary_result_v1.json` |",
            "| case CSV | `abc_baseline_stationary_wave_subwinding_alias_cases_v1.csv` |",
            "| gauge CSV | `abc_baseline_stationary_wave_subwinding_alias_gauge_rows_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_subwinding_alias_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波Cにおけるサブ巻数ビートalias予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    cases = [
        AliasCase("above_integer_winding", winding_count=5, residual_phase_per_step=0.03),
        AliasCase("below_integer_winding", winding_count=5, residual_phase_per_step=-0.03),
        AliasCase("above_high_winding", winding_count=41, residual_phase_per_step=0.0125),
        AliasCase("below_high_winding", winding_count=41, residual_phase_per_step=-0.0125),
    ]
    raw_rows = [run_case(case) for case in cases]
    gauge_rows: List[Dict[str, Any]] = []
    case_rows: List[Dict[str, Any]] = []
    for row in raw_rows:
        gauge_rows.extend(row.pop("_gauge_rows"))
        case_rows.append(row)

    aggregate_verdict = {
        "case_count": len(case_rows),
        "all_cases_valid": bool_all(row["valid"] for row in case_rows),
        "single_gauge_only_used": False,
        "subwinding_alias_preliminary_valid": bool_all(row["valid"] for row in case_rows),
        "max_beat_slope_error": max(
            abs(float(row["observed_beat_slope_mean"]) - float(row["residual_phase_per_step"]))
            for row in case_rows
        ),
        "max_gauge_slope_std": max(abs(float(row["observed_beat_slope_std"])) for row in case_rows),
    }
    result = {
        "experiment": "abc_baseline_stationary_wave_subwinding_alias_preliminary_v1",
        "case_rows": case_rows,
        "aggregate_verdict": aggregate_verdict,
    }
    (OUT_DIR / "abc_baseline_stationary_wave_subwinding_alias_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_subwinding_alias_cases_v1.csv", case_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_subwinding_alias_gauge_rows_v1.csv", gauge_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
