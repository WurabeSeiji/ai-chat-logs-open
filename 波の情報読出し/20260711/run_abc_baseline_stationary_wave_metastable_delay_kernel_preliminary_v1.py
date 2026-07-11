from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import Stage2Params
from run_abc_baseline_stationary_wave_stage3_cross_readout_preliminary_v1 import reduced_R
from run_abc_baseline_stationary_wave_transition_protocol_preliminary_v1 import (
    FINAL_C_MEMORY,
    ONSET_STEP,
    ramp_profile,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


STEP_COUNT = 96
MAX_LAG = 8
FIT_FLOOR = 1.0e-30


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


def zeros() -> List[float]:
    return [0.0 for _ in range(STEP_COUNT)]


def pulse_profile(step: int, amplitude: float) -> List[float]:
    values = zeros()
    values[step] = amplitude
    return values


def hold_release_profile(onset: int, release: int, amplitude: float) -> List[float]:
    values = zeros()
    for step in range(STEP_COUNT):
        if onset <= step < release:
            values[step] = amplitude
    return values


def double_pulse_profile(first: int, second: int, amplitude: float) -> List[float]:
    values = zeros()
    values[first] = amplitude
    values[second] = -0.65 * amplitude
    return values


def delayed_step_pair_profile(onset: int, middle: int, second: int, amplitude: float) -> List[float]:
    values = zeros()
    for step in range(STEP_COUNT):
        if onset <= step < middle:
            values[step] = amplitude
        elif second <= step:
            values[step] = 0.55 * amplitude
    return values


def triangle_profile(onset: int, peak_step: int, end_step: int, amplitude: float) -> List[float]:
    values = zeros()
    for step in range(STEP_COUNT):
        if onset <= step <= peak_step:
            values[step] = amplitude * (step - onset) / (peak_step - onset)
        elif peak_step < step <= end_step:
            values[step] = amplitude * (1.0 - (step - peak_step) / (end_step - peak_step))
    return values


def sine_packet_profile(onset: int, end_step: int, amplitude: float, period: float) -> List[float]:
    values = zeros()
    for step in range(onset, min(end_step, STEP_COUNT)):
        phase = 2.0 * math.pi * (step - onset) / period
        envelope = math.sin(math.pi * (step - onset) / max(1.0, end_step - onset)) ** 2
        values[step] = amplitude * envelope * math.sin(phase)
    return values


def staircase_profile(onset: int, width: int, amplitude: float) -> List[float]:
    levels = [0.0, 0.20, 0.55, 0.35, 0.85, 0.10, 0.0]
    values = zeros()
    for step in range(STEP_COUNT):
        idx = (step - onset) // width
        if 0 <= idx < len(levels):
            values[step] = amplitude * levels[idx]
        elif idx >= len(levels):
            values[step] = 0.0
    return values


def protocol_profiles() -> Dict[str, Tuple[List[float], Dict[str, Any]]]:
    amp = FINAL_C_MEMORY
    return {
        "single_impulse": (pulse_profile(ONSET_STEP + 2, amp), {"kind": "impulse"}),
        "double_signed_pulse": (double_pulse_profile(ONSET_STEP + 2, ONSET_STEP + 18, amp), {"kind": "pulse"}),
        "hold_release": (hold_release_profile(ONSET_STEP + 2, ONSET_STEP + 28, amp), {"kind": "hold_release"}),
        "delayed_step_pair": (
            delayed_step_pair_profile(ONSET_STEP + 1, ONSET_STEP + 18, ONSET_STEP + 38, amp),
            {"kind": "delayed_step_pair"},
        ),
        "ramp_duration_8": (ramp_profile(STEP_COUNT, ONSET_STEP, amp, 8), {"kind": "ramp"}),
        "ramp_duration_24": (ramp_profile(STEP_COUNT, ONSET_STEP, amp, 24), {"kind": "ramp"}),
        "triangle": (triangle_profile(ONSET_STEP + 1, ONSET_STEP + 20, ONSET_STEP + 42, amp), {"kind": "triangle"}),
        "sine_packet": (sine_packet_profile(ONSET_STEP, ONSET_STEP + 56, amp, 15.0), {"kind": "sine_packet"}),
        "staircase": (staircase_profile(ONSET_STEP, 7, amp), {"kind": "staircase"}),
    }


def delta_at(profile: List[float], step: int, lag: int = 0) -> float:
    idx = step - lag
    if idx < 0 or idx >= len(profile) - 1:
        return 0.0
    return float(profile[idx + 1] - profile[idx])


def feature_values(profile: List[float], step: int, params: Stage2Params) -> Dict[str, float]:
    features: Dict[str, float] = {}
    for lag in range(MAX_LAG + 1):
        features[f"delta_C_lag_{lag}"] = delta_at(profile, step, lag)
    features["C_level"] = float(profile[step])
    features["C_integral"] = float(np.sum(profile[: step + 1]))
    decay = params.memory_decay_C
    features["delayed_delta_tail"] = float(
        sum((decay ** (lag - 1)) * delta_at(profile, step, lag) for lag in range(1, MAX_LAG + 1))
    )
    features["delta_tail_including_current"] = float(
        sum((decay**lag) * delta_at(profile, step, lag) for lag in range(0, MAX_LAG + 1))
    )
    features["Q_raw_abs_proxy"] = float(params.q_raw_gain * abs(delta_at(profile, step, 0)))
    features["Q_recovery_abs_tail_proxy"] = float(
        params.q_raw_gain
        * sum((decay ** (lag - 1)) * abs(delta_at(profile, step, lag)) for lag in range(1, MAX_LAG + 1))
    )
    features["Q_recovery_signed_tail_proxy"] = float(
        params.q_raw_gain
        * sum((decay ** (lag - 1)) * delta_at(profile, step, lag) for lag in range(1, MAX_LAG + 1))
    )
    return features


def build_rows(params: Stage2Params) -> List[Dict[str, Any]]:
    gain = reduced_R(params) * params.epsilon_C_return
    rows: List[Dict[str, Any]] = []
    for protocol, (profile, metadata) in protocol_profiles().items():
        for step in range(1, len(profile) - 1):
            target_a = gain * delta_at(profile, step, 0)
            features = feature_values(profile, step, params)
            for particle, sign_value in [("A", 1.0), ("B", -1.0)]:
                actual_ra = sign_value * target_a
                row: Dict[str, Any] = {
                    "protocol": protocol,
                    "kind": metadata["kind"],
                    "step": step,
                    "particle": particle,
                    "C_memory": profile[step],
                    "actual_Ra": actual_ra,
                    "signed_target_Ra": sign_value * actual_ra,
                    "sign_normalizer": sign_value,
                }
                row.update(features)
                rows.append(row)
    return rows


def fit_single_feature(rows: List[Dict[str, Any]], feature: str) -> Dict[str, Any]:
    x = np.array([float(row[feature]) for row in rows], dtype=float)
    y = np.array([float(row["signed_target_Ra"]) for row in rows], dtype=float)
    denom = float(np.dot(x, x))
    scale = float(np.dot(x, y) / denom) if denom > FIT_FLOOR else 0.0
    pred = scale * x
    residual = y - pred
    sse = float(np.dot(residual, residual))
    centered = y - float(np.mean(y))
    sst = float(np.dot(centered, centered))
    r2 = 1.0 - sse / sst if sst > FIT_FLOOR else 1.0
    nrmse = math.sqrt(sse / max(FIT_FLOOR, float(np.dot(y, y))))
    if float(np.std(x)) > FIT_FLOOR and float(np.std(y)) > FIT_FLOOR:
        corr = float(np.corrcoef(x, y)[0, 1])
    else:
        corr = 0.0
    return {
        "feature": feature,
        "scale": scale,
        "r2": r2,
        "nrmse": nrmse,
        "corr": corr,
        "max_residual_abs": float(np.max(np.abs(residual))) if residual.size else 0.0,
    }


def fit_lag_kernel(rows: List[Dict[str, Any]], lags: int = MAX_LAG) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    feature_names = [f"delta_C_lag_{lag}" for lag in range(lags + 1)]
    x = np.array([[float(row[name]) for name in feature_names] for row in rows], dtype=float)
    y = np.array([float(row["signed_target_Ra"]) for row in rows], dtype=float)
    scales = np.max(np.abs(x), axis=0)
    safe_scales = np.where(scales > FIT_FLOOR, scales, 1.0)
    x_scaled = x / safe_scales
    coef_scaled, *_ = np.linalg.lstsq(x_scaled, y, rcond=None)
    coef = coef_scaled / safe_scales
    pred = np.sum(x_scaled * coef_scaled, axis=1)
    residual = y - pred
    sse = float(np.dot(residual, residual))
    centered = y - float(np.mean(y))
    sst = float(np.dot(centered, centered))
    sum_abs = float(np.sum(np.abs(coef)))
    dominant_lag = int(np.argmax(np.abs(coef)))
    summary = {
        "lag_count": lags + 1,
        "dominant_lag": dominant_lag,
        "lag0_abs_fraction": float(abs(coef[0]) / sum_abs) if sum_abs > FIT_FLOOR else 0.0,
        "delayed_abs_fraction": float(np.sum(np.abs(coef[1:])) / sum_abs) if sum_abs > FIT_FLOOR else 0.0,
        "r2": 1.0 - sse / sst if sst > FIT_FLOOR else 1.0,
        "nrmse": math.sqrt(sse / max(FIT_FLOOR, float(np.dot(y, y)))),
        "max_residual_abs": float(np.max(np.abs(residual))) if residual.size else 0.0,
    }
    coef_rows = [
        {
            "lag": lag,
            "feature": feature_names[lag],
            "coefficient": float(coef[lag]),
            "abs_fraction": float(abs(coef[lag]) / sum_abs) if sum_abs > FIT_FLOOR else 0.0,
        }
        for lag in range(lags + 1)
    ]
    return summary, coef_rows


def summarize_by_protocol(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    summaries: List[Dict[str, Any]] = []
    for protocol in sorted({str(row["protocol"]) for row in rows}):
        selected = [row for row in rows if row["protocol"] == protocol and row["particle"] == "A"]
        target = np.array([float(row["signed_target_Ra"]) for row in selected], dtype=float)
        delta0 = np.array([float(row["delta_C_lag_0"]) for row in selected], dtype=float)
        active = [row for row in selected if abs(float(row["signed_target_Ra"])) > 1.0e-18]
        summaries.append(
            {
                "protocol": protocol,
                "kind": selected[0]["kind"] if selected else "",
                "active_step_count": len(active),
                "max_signed_target_Ra_abs": float(np.max(np.abs(target))) if target.size else 0.0,
                "max_delta_C_lag0_abs": float(np.max(np.abs(delta0))) if delta0.size else 0.0,
                "target_delta0_corr": float(np.corrcoef(target, delta0)[0, 1])
                if target.size > 1 and float(np.std(target)) > FIT_FLOOR and float(np.std(delta0)) > FIT_FLOOR
                else 0.0,
            }
        )
    return summaries


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 準安定傾斜遅延カーネル分解予備実験 v1",
        "",
        "## 目的",
        "",
        "準安定傾斜から読まれる `R*a` 候補が、現在の `ΔC_memory` による瞬時応答なのか、過去の `C` 変形残渣による遅延応答なのかを切り分ける。",
        "",
        "impulse, double pulse, hold-release, delayed step, ramp, triangle, sine packet, staircase を用い、`R*a` を複数の候補特徴量で回帰する。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## 単一特徴量フィット",
            "",
            "| feature | scale | R2 | nRMSE | corr | max residual |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["single_feature_fits"]:
        lines.append(
            f"| {row['feature']} | {row['scale']:.16e} | {row['r2']:.16e} | "
            f"{row['nrmse']:.16e} | {row['corr']:.16e} | {row['max_residual_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 遅延ラグカーネル",
            "",
            "| lag | coefficient | abs fraction |",
            "|---:|---:|---:|",
        ]
    )
    for row in result["lag_kernel_coefficients"]:
        lines.append(
            f"| {row['lag']} | {row['coefficient']:.16e} | {row['abs_fraction']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## プロファイル別サマリー",
            "",
            "| protocol | kind | active | max abs R*a | max abs ΔC | corr target-ΔC0 |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in result["protocol_summaries"]:
        lines.append(
            f"| {row['protocol']} | {row['kind']} | {row['active_step_count']} | "
            f"{row['max_signed_target_Ra_abs']:.16e} | {row['max_delta_C_lag0_abs']:.16e} | "
            f"{row['target_delta0_corr']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- `delta_C_lag_0` が支配的なら、現行の準安定傾斜候補は現在ステップの `ΔC_memory` による瞬時傾斜として読まれる。",
            "- `lag >= 1` に大きな係数が出る場合、過去の `C` 変形残渣が遅れて `R*a` に入る遅延カーネル候補になる。",
            "- `C_level` や `C_integral` が支配する場合、定常レベルまたは累積ポテンシャル型の候補になる。",
            "- 本実験は標準重力の導出ではなく、準安定傾斜候補の時間履歴構造を分解する予備検査である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_metastable_delay_kernel_preliminary_result_v1.json` |",
            "| profiles CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_profiles_v1.csv` |",
            "| rows CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_rows_v1.csv` |",
            "| feature fit CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_feature_fits_v1.csv` |",
            "| kernel CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_coefficients_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 準安定傾斜遅延カーネル分解予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params(step_count=STEP_COUNT)
    rows = build_rows(params)
    feature_names = [
        *[f"delta_C_lag_{lag}" for lag in range(MAX_LAG + 1)],
        "delayed_delta_tail",
        "delta_tail_including_current",
        "C_level",
        "C_integral",
        "Q_raw_abs_proxy",
        "Q_recovery_abs_tail_proxy",
        "Q_recovery_signed_tail_proxy",
    ]
    fits = [fit_single_feature(rows, feature) for feature in feature_names]
    fits.sort(key=lambda row: (float(row["r2"]), -float(row["nrmse"])), reverse=True)
    lag_summary, lag_coefficients = fit_lag_kernel(rows)
    protocol_summaries = summarize_by_protocol(rows)

    best_fit = fits[0]
    aggregate_verdict = {
        "protocol_count": len(protocol_summaries),
        "row_count": len(rows),
        "single_gauge_only_used": False,
        "best_single_feature": best_fit["feature"],
        "best_single_feature_r2": best_fit["r2"],
        "best_single_feature_nrmse": best_fit["nrmse"],
        "lag_kernel_dominant_lag": lag_summary["dominant_lag"],
        "lag_kernel_lag0_abs_fraction": lag_summary["lag0_abs_fraction"],
        "lag_kernel_delayed_abs_fraction": lag_summary["delayed_abs_fraction"],
        "lag_kernel_r2": lag_summary["r2"],
        "lag_kernel_max_residual_abs": lag_summary["max_residual_abs"],
        "C_level_r2": next(row["r2"] for row in fits if row["feature"] == "C_level"),
        "C_integral_r2": next(row["r2"] for row in fits if row["feature"] == "C_integral"),
        "Q_raw_abs_proxy_r2": next(row["r2"] for row in fits if row["feature"] == "Q_raw_abs_proxy"),
        "delayed_delta_tail_r2": next(row["r2"] for row in fits if row["feature"] == "delayed_delta_tail"),
    }
    aggregate_verdict["instant_delta_C_dominant"] = bool(
        aggregate_verdict["best_single_feature"] == "delta_C_lag_0"
        and float(aggregate_verdict["best_single_feature_r2"]) > 1.0 - 1.0e-14
        and int(aggregate_verdict["lag_kernel_dominant_lag"]) == 0
        and float(aggregate_verdict["lag_kernel_lag0_abs_fraction"]) > 1.0 - 1.0e-12
    )
    aggregate_verdict["delayed_kernel_detected"] = bool(
        int(aggregate_verdict["lag_kernel_dominant_lag"]) > 0
        or float(aggregate_verdict["lag_kernel_delayed_abs_fraction"]) > 1.0e-6
    )
    aggregate_verdict["level_or_integral_dominant"] = bool(
        aggregate_verdict["best_single_feature"] in ["C_level", "C_integral"]
    )
    aggregate_verdict["metastable_delay_kernel_preliminary_valid"] = bool(
        aggregate_verdict["instant_delta_C_dominant"]
        and not aggregate_verdict["delayed_kernel_detected"]
        and not aggregate_verdict["level_or_integral_dominant"]
    )

    profiles_rows: List[Dict[str, Any]] = []
    for protocol, (profile, metadata) in protocol_profiles().items():
        for step, value in enumerate(profile):
            profiles_rows.append({"protocol": protocol, "kind": metadata["kind"], "step": step, "C_memory": value})

    result = {
        "experiment": "abc_baseline_stationary_wave_metastable_delay_kernel_preliminary_v1",
        "params": asdict(params),
        "step_count": STEP_COUNT,
        "max_lag": MAX_LAG,
        "single_feature_fits": fits,
        "lag_kernel_summary": lag_summary,
        "lag_kernel_coefficients": lag_coefficients,
        "protocol_summaries": protocol_summaries,
        "aggregate_verdict": aggregate_verdict,
    }

    (
        OUT_DIR
        / "abc_baseline_stationary_wave_metastable_delay_kernel_preliminary_result_v1.json"
    ).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_profiles_v1.csv", profiles_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_rows_v1.csv", rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_feature_fits_v1.csv", fits)
    write_csv(
        OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_coefficients_v1.csv",
        lag_coefficients,
    )
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
