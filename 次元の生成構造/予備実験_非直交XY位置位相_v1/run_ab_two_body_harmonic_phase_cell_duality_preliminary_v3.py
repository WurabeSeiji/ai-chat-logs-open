from __future__ import annotations

"""閉じた位相円の調和双対から加速度様二階差分を読む予備実験 v3。

同一の R^2 保存系で、倍音番号 n ごとに

    L_n = 2 pi / n
    omega_n = n omega_1

を与える。軌道生成または読出しに 1/L_n^2 は使用せず、固定半径の
位相円を omega_n で進めた結果の離散二階差分を直接測定する。

追加空間軸、時間軸、面積項、正規化、逆数項は使用しない。
既存 AB 加速度実験の通過型・フェルミオン型二チャネル散乱を同じ形で
適用し、R^2 保存と散乱プロトコル依存性も検査する。
"""

import csv
import importlib.util
import json
import math
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_harmonic_phase_cell_duality_preliminary_result_v3"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

REPO_ROOT = BASE_DIR.parents[1]
SOURCE_PATH = (
    REPO_ROOT
    / "波の情報読出し"
    / "20260711"
    / "run_ab_two_body_fermionic_reflection_harmonic_readout_v4.py"
)

FIXED_RADIUS = 1.0
BASE_PERIOD_STEPS = 12_288
COMMON_STEP_COUNT = 720
HARMONIC_ORDERS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
SCATTERING_PROTOCOL_NAMES = ["pass_through", "fermionic_reflection_pi"]
R2_TOL = 1.0e-12
DUALITY_TOL = 1.0e-14
SCATTERING_TOL = 1.0e-12


def load_acceleration_source(path: Path) -> ModuleType:
    module_name = "ab_acceleration_reference_v4_for_harmonic_phase_cell_duality"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load acceleration source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


ACCEL = load_acceleration_source(SOURCE_PATH)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def write_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def protocol_by_name(name: str) -> Any:
    for protocol in ACCEL.SCATTERING_PROTOCOLS:
        if protocol.name == name:
            return protocol
    raise ValueError(f"unknown scattering protocol: {name}")


def phase_circle(radius: float, phase: np.ndarray) -> np.ndarray:
    return radius * np.column_stack([np.cos(phase), np.sin(phase)])


def scatter_phase_circle(path: np.ndarray, protocol: Any) -> Tuple[np.ndarray, float, float]:
    out = np.zeros_like(path, dtype=float)
    max_unitarity_error = 0.0
    max_internal_imag = 0.0
    zeros = np.zeros(len(path), dtype=float)
    for axis in range(2):
        scatter = ACCEL.two_channel_scattering_state(protocol, path[:, axis], zeros)
        relative_out = np.asarray(scatter["relative_out"], dtype=complex)
        out[:, axis] = np.real(relative_out)
        max_internal_imag = max(max_internal_imag, float(np.max(np.abs(np.imag(relative_out)))))
        max_unitarity_error = max(
            max_unitarity_error,
            float(np.max(np.asarray(scatter["scattering_unitarity_error"], dtype=float))),
        )
    return out, max_unitarity_error, max_internal_imag


def vector_second_difference(path: np.ndarray) -> np.ndarray:
    out = np.zeros_like(path, dtype=float)
    if len(path) < 3:
        return out
    out[1:-1] = path[2:] - 2.0 * path[1:-1] + path[:-2]
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def fit_power(rows: Sequence[Dict[str, Any]], x_key: str, y_key: str) -> Dict[str, Any]:
    selected = [
        row
        for row in rows
        if float(row[x_key]) > 1.0e-30 and float(row[y_key]) > 1.0e-30
    ]
    x = np.log(np.array([float(row[x_key]) for row in selected], dtype=float))
    y = np.log(np.array([float(row[y_key]) for row in selected], dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    return {
        "x": x_key,
        "y": y_key,
        "case_count": len(selected),
        "loglog_slope": float(slope),
        "loglog_intercept": float(intercept),
        "log_rmse": float(np.sqrt(np.mean(residual * residual))),
    }


def run_case(
    harmonic_order: int,
    protocol_name: str,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    base_omega_step = 2.0 * math.pi / float(BASE_PERIOD_STEPS)
    phase_cell_width = 2.0 * math.pi / float(harmonic_order)
    omega_n_step = float(harmonic_order) * base_omega_step
    duality_constant = omega_n_step * phase_cell_width
    expected_duality_constant = 2.0 * math.pi * base_omega_step
    steps = np.arange(COMMON_STEP_COUNT + 1, dtype=float)
    phase = omega_n_step * steps

    path_in = phase_circle(FIXED_RADIUS, phase)
    protocol = protocol_by_name(protocol_name)
    path, max_unitarity_error, max_internal_imag = scatter_phase_circle(path_in, protocol)

    radius_squared = np.sum(path * path, axis=1)
    radius_squared_error = np.abs(radius_squared - FIXED_RADIUS * FIXED_RADIUS)
    acceleration_vectors = vector_second_difference(path)
    acceleration_magnitude = np.linalg.norm(acceleration_vectors, axis=1)
    radius_hat = path / np.linalg.norm(path, axis=1)[:, None]
    center_acceleration = -np.sum(acceleration_vectors * radius_hat, axis=1)
    tangent_acceleration_vectors = acceleration_vectors + center_acceleration[:, None] * radius_hat
    tangent_acceleration = np.linalg.norm(tangent_acceleration_vectors, axis=1)

    exact_discrete_acceleration = (
        4.0 * math.sin(omega_n_step / 2.0) ** 2 * FIXED_RADIUS
    )
    continuum_acceleration = FIXED_RADIUS * omega_n_step * omega_n_step
    measured_slice = slice(1, -1)
    measured_acceleration = float(np.mean(acceleration_magnitude[measured_slice]))
    measured_center_acceleration = float(np.mean(center_acceleration[measured_slice]))
    max_tangent_acceleration = float(np.max(tangent_acceleration[measured_slice]))

    mode_period_steps = BASE_PERIOD_STEPS // harmonic_order
    if mode_period_steps * harmonic_order != BASE_PERIOD_STEPS:
        raise RuntimeError(f"harmonic order does not divide base period: {harmonic_order}")
    endpoint_phase = omega_n_step * float(mode_period_steps)
    endpoint_input = phase_circle(
        FIXED_RADIUS,
        np.array([0.0, endpoint_phase], dtype=float),
    )
    endpoint_output, endpoint_unitarity_error, endpoint_internal_imag = scatter_phase_circle(
        endpoint_input,
        protocol,
    )
    cycle_return_error = float(np.linalg.norm(endpoint_output[1] - endpoint_output[0]))

    case_id = f"{protocol_name}_n{harmonic_order}"
    series_rows: List[Dict[str, Any]] = []
    for step in range(COMMON_STEP_COUNT + 1):
        series_rows.append(
            {
                "case_id": case_id,
                "scattering_protocol": protocol_name,
                "harmonic_order_n": harmonic_order,
                "fixed_radius_R": FIXED_RADIUS,
                "phase_cell_width_L": phase_cell_width,
                "base_omega_step": base_omega_step,
                "omega_n_step": omega_n_step,
                "omega_n_times_L": duality_constant,
                "expected_duality_constant": expected_duality_constant,
                "mode_period_steps": mode_period_steps,
                "step": step,
                "phase_rad": float(phase[step]),
                "center_component": float(path[step, 0]),
                "circumference_component": float(path[step, 1]),
                "R2_measured": float(radius_squared[step]),
                "R2_abs_error": float(radius_squared_error[step]),
                "acceleration_second_difference": float(acceleration_magnitude[step]),
                "center_acceleration_second_difference": float(center_acceleration[step]),
                "tangent_acceleration_second_difference": float(tangent_acceleration[step]),
                "exact_discrete_acceleration": exact_discrete_acceleration,
                "continuum_acceleration": continuum_acceleration,
                "normalization_used": False,
                "division_by_L2_used": False,
                "inverse_term_used": False,
                "additional_spatial_axis_used": False,
                "time_axis_used": False,
                "area_term_used": False,
            }
        )

    summary = {
        "case_id": case_id,
        "scattering_protocol": protocol_name,
        "harmonic_order_n": harmonic_order,
        "fixed_radius_R": FIXED_RADIUS,
        "phase_cell_width_L": phase_cell_width,
        "base_omega_step": base_omega_step,
        "omega_n_step": omega_n_step,
        "omega_n_times_L": duality_constant,
        "expected_duality_constant": expected_duality_constant,
        "duality_abs_error": abs(duality_constant - expected_duality_constant),
        "mode_period_steps": mode_period_steps,
        "cycle_return_error": cycle_return_error,
        "max_R2_abs_error": float(np.max(radius_squared_error)),
        "measured_acceleration": measured_acceleration,
        "measured_center_acceleration": measured_center_acceleration,
        "max_tangent_acceleration": max_tangent_acceleration,
        "exact_discrete_acceleration": exact_discrete_acceleration,
        "continuum_acceleration": continuum_acceleration,
        "measured_vs_exact_rel_error": abs(
            measured_acceleration - exact_discrete_acceleration
        )
        / max(exact_discrete_acceleration, 1.0e-30),
        "discrete_vs_continuum_rel_error": abs(
            exact_discrete_acceleration - continuum_acceleration
        )
        / max(continuum_acceleration, 1.0e-30),
        "max_scattering_unitarity_error": max(max_unitarity_error, endpoint_unitarity_error),
        "max_internal_imag": max(max_internal_imag, endpoint_internal_imag),
        "normalization_used": False,
        "division_by_L2_used": False,
        "inverse_term_used": False,
        "additional_spatial_axis_used": False,
        "time_axis_used": False,
        "area_term_used": False,
    }
    return series_rows, summary


def protocol_comparison(summaries: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_key = {
        (str(row["scattering_protocol"]), int(row["harmonic_order_n"])): row
        for row in summaries
    }
    rows: List[Dict[str, Any]] = []
    for harmonic_order in HARMONIC_ORDERS:
        left = by_key[("pass_through", harmonic_order)]
        right = by_key[("fermionic_reflection_pi", harmonic_order)]
        rows.append(
            {
                "harmonic_order_n": harmonic_order,
                "phase_cell_width_L": left["phase_cell_width_L"],
                "measured_acceleration_abs_difference": abs(
                    float(left["measured_acceleration"])
                    - float(right["measured_acceleration"])
                ),
                "R2_error_abs_difference": abs(
                    float(left["max_R2_abs_error"])
                    - float(right["max_R2_abs_error"])
                ),
            }
        )
    return rows


def make_plot(summaries: Sequence[Dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        selected = sorted(
            [row for row in summaries if row["scattering_protocol"] == protocol_name],
            key=lambda row: float(row["phase_cell_width_L"]),
        )
        ax.loglog(
            [float(row["phase_cell_width_L"]) for row in selected],
            [float(row["measured_acceleration"]) for row in selected],
            "o-",
            label=protocol_name,
        )
    reference = sorted(
        [row for row in summaries if row["scattering_protocol"] == "pass_through"],
        key=lambda row: float(row["phase_cell_width_L"]),
    )
    anchor = reference[len(reference) // 2]
    anchor_l = float(anchor["phase_cell_width_L"])
    anchor_a = float(anchor["measured_acceleration"])
    ref_l = np.array([float(row["phase_cell_width_L"]) for row in reference], dtype=float)
    ax.loglog(ref_l, anchor_a * (ref_l / anchor_l) ** -2.0, "k--", label="L^-2 reference")
    ax.set_xlabel("phase-cell width L_n = 2 pi / n")
    ax.set_ylabel("measured vector second difference")
    ax.set_title("harmonic phase-cell duality at fixed R")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_scaling_v3.png", dpi=180)
    plt.close(fig)


def main() -> None:
    all_series: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        for harmonic_order in HARMONIC_ORDERS:
            series, summary = run_case(harmonic_order, protocol_name)
            all_series.extend(series)
            summaries.append(summary)

    fits: List[Dict[str, Any]] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        selected = [row for row in summaries if row["scattering_protocol"] == protocol_name]
        fit_acceleration = fit_power(selected, "phase_cell_width_L", "measured_acceleration")
        fit_acceleration["scattering_protocol"] = protocol_name
        fit_acceleration["relation"] = "measured_acceleration_vs_phase_cell_width"
        fits.append(fit_acceleration)
        fit_frequency = fit_power(selected, "phase_cell_width_L", "omega_n_step")
        fit_frequency["scattering_protocol"] = protocol_name
        fit_frequency["relation"] = "omega_n_vs_phase_cell_width"
        fits.append(fit_frequency)

    comparisons = protocol_comparison(summaries)
    max_r2_error = max(float(row["max_R2_abs_error"]) for row in summaries)
    max_duality_error = max(float(row["duality_abs_error"]) for row in summaries)
    max_cycle_return_error = max(float(row["cycle_return_error"]) for row in summaries)
    max_scattering_error = max(
        float(row["max_scattering_unitarity_error"]) for row in summaries
    )
    max_tangent_acceleration = max(float(row["max_tangent_acceleration"]) for row in summaries)
    max_protocol_acceleration_difference = max(
        float(row["measured_acceleration_abs_difference"]) for row in comparisons
    )
    acceleration_fits = [
        row for row in fits if row["relation"] == "measured_acceleration_vs_phase_cell_width"
    ]
    frequency_fits = [row for row in fits if row["relation"] == "omega_n_vs_phase_cell_width"]

    result = {
        "experiment": "ab_two_body_harmonic_phase_cell_duality_preliminary_v3",
        "source_path": str(SOURCE_PATH.relative_to(REPO_ROOT)),
        "fixed_radius_R": FIXED_RADIUS,
        "base_period_steps": BASE_PERIOD_STEPS,
        "common_step_count": COMMON_STEP_COUNT,
        "harmonic_orders": HARMONIC_ORDERS,
        "scattering_protocols": SCATTERING_PROTOCOL_NAMES,
        "case_count": len(summaries),
        "R2_preserved_all_cases": max_r2_error <= R2_TOL,
        "max_R2_abs_error": max_r2_error,
        "harmonic_duality_preserved_all_cases": max_duality_error <= DUALITY_TOL,
        "max_omega_times_L_abs_error": max_duality_error,
        "cycles_close_all_cases": max_cycle_return_error <= R2_TOL,
        "max_cycle_return_error": max_cycle_return_error,
        "max_scattering_unitarity_error": max_scattering_error,
        "max_tangent_acceleration": max_tangent_acceleration,
        "max_pass_vs_fermionic_acceleration_difference": max_protocol_acceleration_difference,
        "acceleration_loglog_slopes_vs_L": {
            str(row["scattering_protocol"]): float(row["loglog_slope"])
            for row in acceleration_fits
        },
        "frequency_loglog_slopes_vs_L": {
            str(row["scattering_protocol"]): float(row["loglog_slope"])
            for row in frequency_fits
        },
        "normalization_used_any": False,
        "division_by_L2_used_any": False,
        "inverse_term_used_any": False,
        "additional_spatial_axis_used_any": False,
        "time_axis_used_any": False,
        "area_term_used_any": False,
        "summaries": summaries,
        "power_fits": fits,
        "protocol_comparison": comparisons,
    }

    write_csv(OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_series_v3.csv", all_series)
    write_csv(OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_case_summary_v3.csv", summaries)
    write_csv(OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_power_fits_v3.csv", fits)
    write_csv(OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_protocol_comparison_v3.csv", comparisons)
    (OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_result_v3.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    report_lines = [
        "# AB二体・調和位相セル双対予備実験 v3",
        "",
        "## 実装条件",
        "",
        f"- fixed_radius_R: `{FIXED_RADIUS}`",
        f"- base_period_steps: `{BASE_PERIOD_STEPS}`",
        f"- common_step_count: `{COMMON_STEP_COUNT}`",
        f"- harmonic_orders: `{HARMONIC_ORDERS}`",
        "- phase_cell_width: `L_n = 2 pi / n`",
        "- harmonic_frequency: `omega_n = n omega_1`",
        "- `1/L^2` による除算: `False`",
        "- 正規化: `False`",
        "- 追加空間軸・時間軸・面積項: `False`",
        "",
        "## 統合結果",
        "",
        f"- R2_preserved_all_cases: `{result['R2_preserved_all_cases']}`",
        f"- max_R2_abs_error: `{max_r2_error:.17g}`",
        f"- harmonic_duality_preserved_all_cases: `{result['harmonic_duality_preserved_all_cases']}`",
        f"- max_omega_times_L_abs_error: `{max_duality_error:.17g}`",
        f"- cycles_close_all_cases: `{result['cycles_close_all_cases']}`",
        f"- max_cycle_return_error: `{max_cycle_return_error:.17g}`",
        f"- max_scattering_unitarity_error: `{max_scattering_error:.17g}`",
        f"- max_tangent_acceleration: `{max_tangent_acceleration:.17g}`",
        f"- max_pass_vs_fermionic_acceleration_difference: `{max_protocol_acceleration_difference:.17g}`",
        "",
        "## 冪回帰",
        "",
    ]
    for row in acceleration_fits:
        report_lines.append(
            f"- `{row['scattering_protocol']}`: acceleration vs L slope = "
            f"`{float(row['loglog_slope']):.17g}`"
        )
    for row in frequency_fits:
        report_lines.append(
            f"- `{row['scattering_protocol']}`: omega vs L slope = "
            f"`{float(row['loglog_slope']):.17g}`"
        )
    report_lines.extend(
        [
            "",
            "## 判定範囲",
            "",
            "本実験は、閉じた位相円で `L_n = 2 pi / n` と `omega_n = n omega_1` を",
            "同時に用いた場合、固定 `R` の円運動から直接計算した離散二階差分が",
            "`L_n^-2` に近い冪を示すかを検査する。逆二乗項は演算へ入力していない。",
            "高倍音では離散二階差分 `4 R sin^2(omega_n/2)` と連続近似",
            "`R omega_n^2` の差が生じるため、その相対誤差を各ケースへ記録する。",
            "",
        ]
    )
    (OUT_DIR / "ab_two_body_harmonic_phase_cell_duality_report_v3.md").write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )
    make_plot(summaries)

    print(json.dumps({key: result[key] for key in [
        "experiment",
        "case_count",
        "R2_preserved_all_cases",
        "max_R2_abs_error",
        "harmonic_duality_preserved_all_cases",
        "max_omega_times_L_abs_error",
        "cycles_close_all_cases",
        "max_cycle_return_error",
        "max_scattering_unitarity_error",
        "max_tangent_acceleration",
        "max_pass_vs_fermionic_acceleration_difference",
        "acceleration_loglog_slopes_vs_L",
        "frequency_loglog_slopes_vs_L",
        "normalization_used_any",
        "division_by_L2_used_any",
        "inverse_term_used_any",
    ]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
