from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_baseline_stationary_wave_metastable_delay_kernel_preliminary_v1 import (
    MAX_LAG,
    delta_at,
    feature_values,
    fit_lag_kernel,
    fit_single_feature,
    protocol_profiles,
)
from run_abc_baseline_stationary_wave_stage2_frozen_fanout_preliminary_v1 import Stage2Params
from run_abc_baseline_stationary_wave_stage3_cross_readout_preliminary_v1 import reduced_R


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_injection_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

RECOVERY_TOL = 1.0e-12


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


def kernel_cases() -> Dict[str, Dict[str, Any]]:
    return {
        "instant_control": {
            "kind": "control",
            "weights": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "lag1_dominant": {
            "kind": "one_step_delay",
            "weights": [0.25, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "lag2_dominant": {
            "kind": "two_step_delay",
            "weights": [0.20, 0.0, 0.80, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "exponential_tail": {
            "kind": "tail",
            "weights": [0.25, 0.34, 0.22, 0.12, 0.07, 0.0, 0.0, 0.0, 0.0],
        },
        "pure_delayed_tail": {
            "kind": "tail_no_current",
            "weights": [0.0, 0.46, 0.28, 0.16, 0.10, 0.0, 0.0, 0.0, 0.0],
        },
        "alternating_delay": {
            "kind": "signed_tail",
            "weights": [0.30, -0.45, 0.25, -0.15, 0.08, 0.0, 0.0, 0.0, 0.0],
        },
    }


def padded_weights(weights: List[float]) -> List[float]:
    padded = list(weights[: MAX_LAG + 1])
    padded.extend([0.0 for _ in range(MAX_LAG + 1 - len(padded))])
    return padded


def expected_summary(weights: List[float]) -> Dict[str, Any]:
    padded = padded_weights(weights)
    abs_sum = float(sum(abs(value) for value in padded))
    dominant_lag = int(np.argmax(np.abs(np.array(padded, dtype=float))))
    return {
        "expected_dominant_lag": dominant_lag,
        "expected_lag0_abs_fraction": abs(padded[0]) / abs_sum if abs_sum else 0.0,
        "expected_delayed_abs_fraction": sum(abs(value) for value in padded[1:]) / abs_sum if abs_sum else 0.0,
    }


def build_rows(params: Stage2Params, kernel_name: str, kernel: Dict[str, Any]) -> List[Dict[str, Any]]:
    gain = reduced_R(params) * params.epsilon_C_return
    weights = padded_weights(kernel["weights"])
    rows: List[Dict[str, Any]] = []
    for protocol, (profile, metadata) in protocol_profiles().items():
        for step in range(1, len(profile) - 1):
            target_a = gain * sum(weights[lag] * delta_at(profile, step, lag) for lag in range(MAX_LAG + 1))
            features = feature_values(profile, step, params)
            for particle, sign_value in [("A", 1.0), ("B", -1.0)]:
                actual_ra = sign_value * target_a
                row: Dict[str, Any] = {
                    "kernel": kernel_name,
                    "kernel_kind": kernel["kind"],
                    "protocol": protocol,
                    "protocol_kind": metadata["kind"],
                    "step": step,
                    "particle": particle,
                    "C_memory": profile[step],
                    "actual_Ra": actual_ra,
                    "signed_target_Ra": sign_value * actual_ra,
                    "sign_normalizer": sign_value,
                }
                for lag, weight in enumerate(weights):
                    row[f"injected_weight_lag_{lag}"] = weight
                    row[f"expected_coefficient_lag_{lag}"] = gain * weight
                row.update(features)
                rows.append(row)
    return rows


def coefficient_rows(kernel_name: str, kernel: Dict[str, Any], params: Stage2Params, recovered: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    gain = reduced_R(params) * params.epsilon_C_return
    weights = padded_weights(kernel["weights"])
    expected = expected_summary(weights)
    rows: List[Dict[str, Any]] = []
    for row in recovered:
        lag = int(row["lag"])
        expected_coefficient = gain * weights[lag]
        coefficient = float(row["coefficient"])
        rows.append(
            {
                "kernel": kernel_name,
                "kernel_kind": kernel["kind"],
                "lag": lag,
                "injected_weight": weights[lag],
                "expected_coefficient": expected_coefficient,
                "recovered_coefficient": coefficient,
                "coefficient_error": coefficient - expected_coefficient,
                "recovered_abs_fraction": row["abs_fraction"],
                **expected,
            }
        )
    return rows


def summarize_kernel(
    kernel_name: str,
    kernel: Dict[str, Any],
    params: Stage2Params,
    rows: List[Dict[str, Any]],
    lag_summary: Dict[str, Any],
    coef_rows: List[Dict[str, Any]],
    fits: List[Dict[str, Any]],
) -> Dict[str, Any]:
    expected = expected_summary(kernel["weights"])
    max_coef_error = max(abs(float(row["coefficient_error"])) for row in coef_rows)
    expected_dominant = int(expected["expected_dominant_lag"])
    recovered_dominant = int(lag_summary["dominant_lag"])
    expected_delayed = float(expected["expected_delayed_abs_fraction"])
    recovered_delayed = float(lag_summary["delayed_abs_fraction"])
    best_fit = fits[0]
    return {
        "kernel": kernel_name,
        "kernel_kind": kernel["kind"],
        "row_count": len(rows),
        "expected_dominant_lag": expected_dominant,
        "recovered_dominant_lag": recovered_dominant,
        "expected_lag0_abs_fraction": expected["expected_lag0_abs_fraction"],
        "recovered_lag0_abs_fraction": lag_summary["lag0_abs_fraction"],
        "expected_delayed_abs_fraction": expected_delayed,
        "recovered_delayed_abs_fraction": recovered_delayed,
        "lag_kernel_r2": lag_summary["r2"],
        "lag_kernel_nrmse": lag_summary["nrmse"],
        "lag_kernel_max_residual_abs": lag_summary["max_residual_abs"],
        "max_coefficient_error_abs": max_coef_error,
        "best_single_feature": best_fit["feature"],
        "best_single_feature_r2": best_fit["r2"],
        "dominant_lag_recovered": recovered_dominant == expected_dominant,
        "delayed_fraction_recovered": abs(recovered_delayed - expected_delayed) <= RECOVERY_TOL,
        "coefficients_recovered": max_coef_error <= RECOVERY_TOL,
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C 準安定傾斜遅延カーネル注入予備実験 v1",
        "",
        "## 目的",
        "",
        "遅延カーネル分解器が本当に lag 1 以降の残渣を検出できるかを確認するため、既知の遅延重みを `R*a` 候補へ人工的に注入する。",
        "",
        "現行モデルでは lag 0 が支配的であった。そこで本実験では、lag 1, lag 2, 指数尾、純遅延尾、符号反転尾を陽性対照として与え、回帰で既知係数が回収されるかを調べる。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## カーネル別サマリー",
            "",
            "| kernel | kind | expected lag | recovered lag | expected delayed | recovered delayed | max coef err | R2 | valid |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["kernel_summaries"]:
        valid = bool(row["dominant_lag_recovered"]) and bool(row["delayed_fraction_recovered"]) and bool(row["coefficients_recovered"])
        lines.append(
            f"| {row['kernel']} | {row['kernel_kind']} | {row['expected_dominant_lag']} | "
            f"{row['recovered_dominant_lag']} | {row['expected_delayed_abs_fraction']:.16e} | "
            f"{row['recovered_delayed_abs_fraction']:.16e} | {row['max_coefficient_error_abs']:.16e} | "
            f"{row['lag_kernel_r2']:.16e} | `{valid}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- 既知の遅延尾を入れた場合に lag 1 以降が回収されるなら、前段の `lag 0` 支配は検出器の鈍さではなく現行モデル側の性質と読める。",
            "- 純遅延尾や符号反転尾も回収できるなら、単に正の指数尾だけを見ているのではなく、符号付きカーネルを識別できる。",
            "- 本実験は標準重力の導出ではなく、遅延カーネル分解法の陽性対照である。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_preliminary_result_v1.json` |",
            "| rows CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_rows_v1.csv` |",
            "| summaries CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_summaries_v1.csv` |",
            "| coefficients CSV | `abc_baseline_stationary_wave_metastable_delay_kernel_injection_coefficients_v1.csv` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_injection_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波C 準安定傾斜遅延カーネル注入予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Stage2Params()
    all_rows: List[Dict[str, Any]] = []
    all_coef_rows: List[Dict[str, Any]] = []
    kernel_summaries: List[Dict[str, Any]] = []
    feature_names = [f"delta_C_lag_{lag}" for lag in range(MAX_LAG + 1)]

    for kernel_name, kernel in kernel_cases().items():
        rows = build_rows(params, kernel_name, kernel)
        fits = [fit_single_feature(rows, feature) for feature in feature_names]
        fits.sort(key=lambda row: (float(row["r2"]), -float(row["nrmse"])), reverse=True)
        lag_summary, lag_coefficients = fit_lag_kernel(rows)
        coef_rows = coefficient_rows(kernel_name, kernel, params, lag_coefficients)
        all_rows.extend(rows)
        all_coef_rows.extend(coef_rows)
        kernel_summaries.append(
            summarize_kernel(kernel_name, kernel, params, rows, lag_summary, coef_rows, fits)
        )

    aggregate_verdict = {
        "kernel_count": len(kernel_summaries),
        "row_count": len(all_rows),
        "single_gauge_only_used": False,
        "all_dominant_lags_recovered": bool_all(bool(row["dominant_lag_recovered"]) for row in kernel_summaries),
        "all_delayed_fractions_recovered": bool_all(bool(row["delayed_fraction_recovered"]) for row in kernel_summaries),
        "all_coefficients_recovered": bool_all(bool(row["coefficients_recovered"]) for row in kernel_summaries),
        "max_coefficient_error_abs": max(float(row["max_coefficient_error_abs"]) for row in kernel_summaries),
        "min_lag_kernel_r2": min(float(row["lag_kernel_r2"]) for row in kernel_summaries),
        "delayed_positive_controls_detected": bool_all(
            int(row["recovered_dominant_lag"]) > 0
            or float(row["recovered_delayed_abs_fraction"]) > 0.50
            for row in kernel_summaries
            if row["kernel"] != "instant_control"
        ),
    }
    aggregate_verdict["metastable_delay_kernel_injection_preliminary_valid"] = bool(
        aggregate_verdict["all_dominant_lags_recovered"]
        and aggregate_verdict["all_delayed_fractions_recovered"]
        and aggregate_verdict["all_coefficients_recovered"]
        and aggregate_verdict["delayed_positive_controls_detected"]
    )

    result = {
        "experiment": "abc_baseline_stationary_wave_metastable_delay_kernel_injection_preliminary_v1",
        "params": asdict(params),
        "max_lag": MAX_LAG,
        "recovery_tol": RECOVERY_TOL,
        "kernel_summaries": kernel_summaries,
        "aggregate_verdict": aggregate_verdict,
    }

    (
        OUT_DIR
        / "abc_baseline_stationary_wave_metastable_delay_kernel_injection_preliminary_result_v1.json"
    ).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_injection_rows_v1.csv", all_rows)
    write_csv(
        OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_injection_summaries_v1.csv",
        kernel_summaries,
    )
    write_csv(
        OUT_DIR / "abc_baseline_stationary_wave_metastable_delay_kernel_injection_coefficients_v1.csv",
        all_coef_rows,
    )
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
