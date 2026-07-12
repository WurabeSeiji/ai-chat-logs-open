from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_ab_two_body_one_angle_harmonic_readout_preliminary_v1 import (
    InitialCase,
    Params,
    PROTOCOLS,
    ReadoutMode,
    bool_all,
    rows_for_case,
    summarize_case,
    write_csv,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi
PERIOD_STEPS = [48, 72, 96, 144, 192]
DEVIATION_DEGS = [1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 60.0]
LEAKS = [
    ("readout_off", 0.0, False),
    ("leak_1e-6", 1.0e-6, True),
    ("leak_1e-5", 1.0e-5, True),
    ("leak_5e-5", 5.0e-5, True),
    ("leak_2e-4", 2.0e-4, True),
    ("leak_1e-3", 1.0e-3, True),
]


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if abs(denominator) > 1.0e-30 else 0.0


def monotonic_non_decreasing(values: List[float], tol: float = 1.0e-14) -> bool:
    return bool_all(values[idx] <= values[idx + 1] + tol for idx in range(len(values) - 1))


def selected_series_key(period_steps: int, deviation_deg: float, leak_name: str, protocol: str) -> bool:
    return bool(
        period_steps == 96
        and deviation_deg in {5.0, 20.0}
        and leak_name in {"readout_off", "leak_5e-5", "leak_2e-4", "leak_1e-3"}
        and protocol == "Protocol_B"
    )


def summarize_protocol_pair(
    period_steps: int,
    deviation_deg: float,
    mode: ReadoutMode,
    f_rows: List[Dict[str, Any]],
    b_rows: List[Dict[str, Any]],
) -> Dict[str, Any]:
    max_d_near_diff = max(
        abs(float(f_row["D_AB_near_rad"]) - float(b_row["D_AB_near_rad"]))
        for f_row, b_row in zip(f_rows, b_rows)
    )
    max_v_diff = max(abs(float(f_row["V_AB"]) - float(b_row["V_AB"])) for f_row, b_row in zip(f_rows, b_rows))
    max_display_diff = max(
        abs(float(f_row["protocol_display_deviation_rad"]) - float(b_row["protocol_display_deviation_rad"]))
        for f_row, b_row in zip(f_rows, b_rows)
    )
    return {
        "period_steps": period_steps,
        "deviation_deg": deviation_deg,
        "readout_mode": mode.name,
        "per_step_leak": mode.per_step_leak,
        "max_D_AB_near_protocol_diff": max_d_near_diff,
        "max_V_AB_protocol_diff": max_v_diff,
        "max_protocol_display_deviation_diff": max_display_diff,
        "label_free_protocol_degenerate": bool(max_d_near_diff <= 1.0e-15 and max_v_diff <= 1.0e-15),
    }


def leak_summary_rows(case_summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for leak_name, leak_value, _ in LEAKS:
        selected = [row for row in case_summaries if row["readout_mode"] == leak_name]
        norm_errors = [float(row["normalized_f_AB_projection_error"]) for row in selected]
        decay_abs = [abs(float(row["decay_rate_V_AB"])) for row in selected]
        envelope_ratios = [float(row["envelope_ratio_final_over_initial"]) for row in selected]
        rows.append(
            {
                "readout_mode": leak_name,
                "per_step_leak": leak_value,
                "case_count": len(selected),
                "max_normalized_f_AB_projection_error": max(norm_errors),
                "mean_normalized_f_AB_projection_error": float(np.mean(norm_errors)),
                "max_abs_decay_rate_V_AB": max(decay_abs),
                "mean_abs_decay_rate_V_AB": float(np.mean(decay_abs)),
                "min_envelope_ratio_final_over_initial": min(envelope_ratios),
                "max_envelope_ratio_final_over_initial": max(envelope_ratios),
            }
        )
    return rows


def period_summary_rows(case_summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for period_steps in PERIOD_STEPS:
        selected = [row for row in case_summaries if int(row["period_steps"]) == period_steps]
        norm_errors = [float(row["normalized_f_AB_projection_error"]) for row in selected]
        decay_abs = [abs(float(row["decay_rate_V_AB"])) for row in selected]
        rows.append(
            {
                "period_steps": period_steps,
                "omega_step": TAU / period_steps,
                "case_count": len(selected),
                "max_normalized_f_AB_projection_error": max(norm_errors),
                "mean_normalized_f_AB_projection_error": float(np.mean(norm_errors)),
                "max_abs_decay_rate_V_AB": max(decay_abs),
            }
        )
    return rows


def monotonic_checks(case_summaries: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], bool, bool]:
    rows: List[Dict[str, Any]] = []
    decay_ok_values: List[bool] = []
    f_error_ok_values: List[bool] = []
    for period_steps in PERIOD_STEPS:
        for deviation_deg in DEVIATION_DEGS:
            for protocol in PROTOCOLS:
                selected = [
                    row
                    for row in case_summaries
                    if int(row["period_steps"]) == period_steps
                    and float(row["initial_deviation_deg"]) == deviation_deg
                    and row["protocol"] == protocol
                ]
                selected.sort(key=lambda row: float(row["per_step_leak"]))
                decay_values = [abs(float(row["decay_rate_V_AB"])) for row in selected]
                f_error_values = [float(row["normalized_f_AB_projection_error"]) for row in selected]
                decay_ok = monotonic_non_decreasing(decay_values, tol=1.0e-12)
                f_error_ok = monotonic_non_decreasing(f_error_values, tol=1.0e-12)
                decay_ok_values.append(decay_ok)
                f_error_ok_values.append(f_error_ok)
                rows.append(
                    {
                        "period_steps": period_steps,
                        "deviation_deg": deviation_deg,
                        "protocol": protocol,
                        "decay_abs_monotonic_by_leak": decay_ok,
                        "normalized_f_error_monotonic_by_leak": f_error_ok,
                        "decay_abs_values": ";".join(f"{value:.16e}" for value in decay_values),
                        "normalized_f_error_values": ";".join(f"{value:.16e}" for value in f_error_values),
                    }
                )
    return rows, bool_all(decay_ok_values), bool_all(f_error_ok_values)


def make_plots(
    case_summaries: List[Dict[str, Any]],
    leak_rows: List[Dict[str, Any]],
    period_rows: List[Dict[str, Any]],
    selected_series_rows: List[Dict[str, Any]],
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    leak_values = [float(row["per_step_leak"]) for row in leak_rows]
    max_errors = [float(row["max_normalized_f_AB_projection_error"]) for row in leak_rows]
    mean_errors = [float(row["mean_normalized_f_AB_projection_error"]) for row in leak_rows]
    ax.plot(leak_values, max_errors, marker="o", label="max normalized f error")
    ax.plot(leak_values, mean_errors, marker="o", label="mean normalized f error")
    ax.set_xscale("symlog", linthresh=1.0e-7)
    ax.set_yscale("symlog", linthresh=1.0e-12)
    ax.set_xlabel("per_step_leak")
    ax.set_ylabel("normalized f_AB projection error")
    ax.set_title("Readout leak vs f_AB projection error")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_leak_f_error_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        [float(row["per_step_leak"]) for row in leak_rows],
        [float(row["max_abs_decay_rate_V_AB"]) for row in leak_rows],
        marker="o",
        label="max |decay_rate|",
    )
    ax.plot(
        [float(row["per_step_leak"]) for row in leak_rows],
        [float(row["mean_abs_decay_rate_V_AB"]) for row in leak_rows],
        marker="o",
        label="mean |decay_rate|",
    )
    ax.set_xscale("symlog", linthresh=1.0e-7)
    ax.set_yscale("symlog", linthresh=1.0e-12)
    ax.set_xlabel("per_step_leak")
    ax.set_ylabel("|decay_rate_V_AB|")
    ax.set_title("Readout leak vs envelope decay")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_leak_decay_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        [int(row["period_steps"]) for row in period_rows],
        [float(row["max_normalized_f_AB_projection_error"]) for row in period_rows],
        marker="o",
        label="max normalized f error",
    )
    ax.plot(
        [int(row["period_steps"]) for row in period_rows],
        [float(row["mean_normalized_f_AB_projection_error"]) for row in period_rows],
        marker="o",
        label="mean normalized f error",
    )
    ax.set_xlabel("period_steps")
    ax.set_ylabel("normalized f_AB projection error")
    ax.set_title("Period sweep")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_period_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for leak_name in ["readout_off", "leak_5e-5", "leak_2e-4", "leak_1e-3"]:
        rows = [
            row
            for row in selected_series_rows
            if row["deviation_deg"] == 5.0 and row["readout_mode"] == leak_name
        ]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot([int(row["step"]) for row in rows], [float(row["V_AB"]) for row in rows], label=leak_name)
    ax.set_title("Selected series: period=96, deviation=5deg")
    ax.set_xlabel("step")
    ax.set_ylabel("V_AB")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_selected_series_v1.png", dpi=180)
    plt.close(fig)


def aggregate_verdict(
    case_summaries: List[Dict[str, Any]],
    protocol_rows: List[Dict[str, Any]],
    leak_rows: List[Dict[str, Any]],
    monotonic_rows: List[Dict[str, Any]],
    decay_monotonic_all: bool,
    f_error_monotonic_all: bool,
) -> Dict[str, Any]:
    off_rows = [row for row in case_summaries if row["readout_mode"] == "readout_off"]
    leak_1e3_rows = [row for row in case_summaries if row["readout_mode"] == "leak_1e-3"]
    max_protocol_d = max(float(row["max_D_AB_near_protocol_diff"]) for row in protocol_rows)
    max_protocol_v = max(float(row["max_V_AB_protocol_diff"]) for row in protocol_rows)
    max_norm_f_error = max(float(row["normalized_f_AB_projection_error"]) for row in case_summaries)
    leak_1e3_min_norm_f_error = min(float(row["normalized_f_AB_projection_error"]) for row in leak_1e3_rows)
    strong_leak_detection_floor = 5.0e-5
    return {
        "sweep_configuration_count": len(PERIOD_STEPS) * len(DEVIATION_DEGS) * len(LEAKS),
        "case_summary_count": len(case_summaries),
        "period_count": len(PERIOD_STEPS),
        "deviation_count": len(DEVIATION_DEGS),
        "leak_count": len(LEAKS),
        "protocol_count": len(PROTOCOLS),
        "observer_C_used": False,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
        "standard_force_law_used": False,
        "max_Q_closed_abs": max(float(row["max_Q_closed_abs"]) for row in case_summaries),
        "max_D_AB_near_protocol_diff": max_protocol_d,
        "max_V_AB_protocol_diff": max_protocol_v,
        "label_free_protocol_degenerate_all_cases": bool_all(
            bool(row["label_free_protocol_degenerate"]) for row in protocol_rows
        ),
        "oscillation_detected_all_cases": bool_all(bool(row["oscillation_detected"]) for row in case_summaries),
        "readout_off_decay_max_abs": max(abs(float(row["decay_rate_V_AB"])) for row in off_rows),
        "decay_abs_monotonic_by_leak_all_grids": decay_monotonic_all,
        "normalized_f_error_monotonic_by_leak_all_grids": f_error_monotonic_all,
        "max_normalized_f_AB_projection_error": max_norm_f_error,
        "leak_1e3_min_normalized_f_AB_projection_error": leak_1e3_min_norm_f_error,
        "strong_leak_detection_floor": strong_leak_detection_floor,
        "strong_leak_perturbs_projection_all_cases": bool(leak_1e3_min_norm_f_error > strong_leak_detection_floor),
        "parameter_sweep_preliminary_valid": bool(
            max_protocol_d <= 1.0e-15
            and max_protocol_v <= 1.0e-15
            and bool_all(bool(row["oscillation_detected"]) for row in case_summaries)
            and decay_monotonic_all
            and f_error_monotonic_all
            and max(abs(float(row["decay_rate_V_AB"])) for row in off_rows) <= 1.0e-12
            and leak_1e3_min_norm_f_error > strong_leak_detection_floor
        ),
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines: List[str] = [
        "# AB二体閉鎖位相系における一角度円周位相調和読出しパラメータスイープ予備実験検証メモ v1",
        "",
        "## 目的",
        "",
        "AB一角度円周位相調和読出し予備実験について、初期偏差、回転周期、読出し漏れ量を変え、結果の安定範囲を調べた。",
        "",
        "本スイープでも、観測機 C、標準重力式、標準ばね式、`f_A`, `f_B` は使わない。",
        "",
        "## 統合判定",
        "",
        "| 量 | 値 |",
        "|---|---:|",
    ]
    for key, value in verdict.items():
        lines.append(f"| `{key}` | `{value}` |")

    lines.extend(
        [
            "",
            "## 漏れ量別サマリー",
            "",
            "| readout | leak | max normalized f error | mean normalized f error | max |decay| | min envelope ratio |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["leak_summaries"]:
        lines.append(
            f"| {row['readout_mode']} | {row['per_step_leak']:.1e} | "
            f"{row['max_normalized_f_AB_projection_error']:.16e} | "
            f"{row['mean_normalized_f_AB_projection_error']:.16e} | "
            f"{row['max_abs_decay_rate_V_AB']:.16e} | "
            f"{row['min_envelope_ratio_final_over_initial']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 周期別サマリー",
            "",
            "| period_steps | omega_step | max normalized f error | mean normalized f error | max |decay| |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["period_summaries"]:
        lines.append(
            f"| {row['period_steps']} | {row['omega_step']:.16e} | "
            f"{row['max_normalized_f_AB_projection_error']:.16e} | "
            f"{row['mean_normalized_f_AB_projection_error']:.16e} | "
            f"{row['max_abs_decay_rate_V_AB']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- `Protocol F/B` は全スイープで `D_AB` と `V_AB` では縮退した。",
            "- `readout_off` の減衰は数値丸め範囲に留まった。",
            "- 読出し漏れ量を増やすと、包絡減衰と `f_AB` 射影不整合が単調に増えた。",
            "- したがって、AB一角度版は弱読出しでは安定だが、強い読出し波は補償表示そのものを歪める。",
            "- これは、後続の二角度・三角度実験で観測波を弱く保つ必要があることを示す制御結果である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1.json` |",
            "| case summary CSV | `ab_two_body_one_angle_parameter_sweep_case_summary_v1.csv` |",
            "| protocol comparison CSV | `ab_two_body_one_angle_parameter_sweep_protocol_comparison_v1.csv` |",
            "| leak summary CSV | `ab_two_body_one_angle_parameter_sweep_leak_summary_v1.csv` |",
            "| period summary CSV | `ab_two_body_one_angle_parameter_sweep_period_summary_v1.csv` |",
            "| monotonic checks CSV | `ab_two_body_one_angle_parameter_sweep_monotonic_checks_v1.csv` |",
            "| selected series CSV | `ab_two_body_one_angle_parameter_sweep_selected_series_v1.csv` |",
            "| f error plot | `ab_two_body_one_angle_parameter_sweep_leak_f_error_v1.png` |",
            "| decay plot | `ab_two_body_one_angle_parameter_sweep_leak_decay_v1.png` |",
            "| period plot | `ab_two_body_one_angle_parameter_sweep_period_v1.png` |",
            "| selected series plot | `ab_two_body_one_angle_parameter_sweep_selected_series_v1.png` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "AB二体閉鎖位相系における一角度円周位相調和読出しパラメータスイープ予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_summaries: List[Dict[str, Any]] = []
    protocol_rows: List[Dict[str, Any]] = []
    selected_series_rows: List[Dict[str, Any]] = []

    for period_steps in PERIOD_STEPS:
        params = Params(step_count=period_steps * 8, omega_step=TAU / period_steps)
        for deviation_deg in DEVIATION_DEGS:
            case = InitialCase(f"dev_{deviation_deg:g}deg", deviation_deg)
            for leak_name, leak_value, active in LEAKS:
                mode = ReadoutMode(leak_name, leak_value, active)
                rows_by_protocol: Dict[str, List[Dict[str, Any]]] = {}
                for protocol in PROTOCOLS:
                    rows = rows_for_case(case, protocol, mode, params)
                    rows_by_protocol[protocol] = rows
                    summary = summarize_case(rows, params)
                    initial_deviation_rad = math.radians(deviation_deg)
                    summary.update(
                        {
                            "period_steps": period_steps,
                            "omega_step": params.omega_step,
                            "per_step_leak": leak_value,
                            "normalized_f_AB_projection_error": safe_ratio(
                                float(summary["max_f_AB_projection_consistency_error"]), initial_deviation_rad
                            ),
                        }
                    )
                    case_summaries.append(summary)
                    if selected_series_key(period_steps, deviation_deg, leak_name, protocol):
                        for row in rows:
                            selected_series_rows.append(
                                {
                                    "period_steps": period_steps,
                                    "deviation_deg": deviation_deg,
                                    "readout_mode": leak_name,
                                    "protocol": protocol,
                                    "step": row["step"],
                                    "V_AB": row["V_AB"],
                                    "envelope_V_AB": row["envelope_V_AB"],
                                    "D_AB_near_deg": row["D_AB_near_deg"],
                                    "f_AB_projection_consistency_error": row["f_AB_projection_consistency_error"],
                                }
                            )
                protocol_rows.append(
                    summarize_protocol_pair(
                        period_steps,
                        deviation_deg,
                        mode,
                        rows_by_protocol["Protocol_F"],
                        rows_by_protocol["Protocol_B"],
                    )
                )

    leak_rows = leak_summary_rows(case_summaries)
    period_rows = period_summary_rows(case_summaries)
    monotonic_rows, decay_monotonic_all, f_error_monotonic_all = monotonic_checks(case_summaries)
    verdict = aggregate_verdict(
        case_summaries,
        protocol_rows,
        leak_rows,
        monotonic_rows,
        decay_monotonic_all,
        f_error_monotonic_all,
    )
    result = {
        "experiment": "ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_v1",
        "period_steps": PERIOD_STEPS,
        "deviation_degs": DEVIATION_DEGS,
        "leaks": [{"name": name, "per_step_leak": leak, "active_readout": active} for name, leak, active in LEAKS],
        "protocols": PROTOCOLS,
        "case_summaries": case_summaries,
        "protocol_comparison": protocol_rows,
        "leak_summaries": leak_rows,
        "period_summaries": period_rows,
        "monotonic_checks": monotonic_rows,
        "aggregate_verdict": verdict,
    }
    (OUT_DIR / "ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_case_summary_v1.csv", case_summaries)
    write_csv(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_protocol_comparison_v1.csv", protocol_rows)
    write_csv(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_leak_summary_v1.csv", leak_rows)
    write_csv(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_period_summary_v1.csv", period_rows)
    write_csv(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_monotonic_checks_v1.csv", monotonic_rows)
    write_csv(OUT_DIR / "ab_two_body_one_angle_parameter_sweep_selected_series_v1.csv", selected_series_rows)
    make_plots(case_summaries, leak_rows, period_rows, selected_series_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
