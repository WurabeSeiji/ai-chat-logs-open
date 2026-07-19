from __future__ import annotations

"""経度 -2.5 度と +2.5 度の二大円による球面スイープ予備実験。

二つの大円を

    r_lambda(phi) = R (cos(phi), sin(phi) cos(lambda), sin(phi) sin(lambda))

で直接生成する。lambda=-2.5 deg と +2.5 deg の大円は phi=0, pi で
交差し、その間に球面リューネを張る。各大円は全ステップで
x^2+y^2+z^2=R^2 を満たす。

既存AB加速度実験の二チャネル散乱行列は各実空間成分へ同じ形で適用し、
離散二階差分から加速度様量を読む。逆数、逆面積、逆二乗項は使用しない。
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
OUT_DIR = BASE_DIR / "ab_two_body_meridian_lune_sweep_preliminary_result_v2"
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

RADIUS_LABEL_DEGS = [1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 60.0]
LONGITUDE_MINUS_DEG = -2.5
LONGITUDE_PLUS_DEG = 2.5
LONGITUDE_SEPARATION_DEG = LONGITUDE_PLUS_DEG - LONGITUDE_MINUS_DEG
SCATTERING_PROTOCOL_NAMES = ["pass_through", "fermionic_reflection_pi"]
INVERSE_POWER_TOL = 0.12


def load_acceleration_source(path: Path) -> ModuleType:
    module_name = "ab_acceleration_reference_v4_for_meridian_lune"
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


def great_circle(radius: float, longitude_rad: float, phase: np.ndarray) -> np.ndarray:
    return radius * np.column_stack(
        [
            np.cos(phase),
            np.sin(phase) * math.cos(longitude_rad),
            np.sin(phase) * math.sin(longitude_rad),
        ]
    )


def scatter_spatial_path(path: np.ndarray, protocol: Any) -> Tuple[np.ndarray, float, float]:
    out = np.zeros_like(path, dtype=float)
    max_unitarity_error = 0.0
    max_internal_imag = 0.0
    zeros = np.zeros(len(path), dtype=float)
    for axis in range(3):
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


def spherical_triangle_area(a: np.ndarray, b: np.ndarray, c: np.ndarray, radius: float) -> float:
    ua = a / radius
    ub = b / radius
    uc = c / radius
    numerator = abs(float(np.dot(ua, np.cross(ub, uc))))
    denominator = 1.0 + float(np.dot(ua, ub) + np.dot(ub, uc) + np.dot(uc, ua))
    return 2.0 * radius * radius * math.atan2(numerator, denominator)


def spherical_quad_area(
    minus_0: np.ndarray,
    minus_1: np.ndarray,
    plus_1: np.ndarray,
    plus_0: np.ndarray,
    radius: float,
) -> float:
    return spherical_triangle_area(minus_0, minus_1, plus_1, radius) + spherical_triangle_area(
        minus_0,
        plus_1,
        plus_0,
        radius,
    )


def one_cycle_lune_area(
    path_minus: np.ndarray,
    path_plus: np.ndarray,
    radius: float,
    period_steps: int,
) -> np.ndarray:
    half_steps = period_steps // 2
    area = np.zeros(period_steps + 1, dtype=float)
    for step in range(period_steps):
        strip = spherical_quad_area(
            path_minus[step],
            path_minus[step + 1],
            path_plus[step + 1],
            path_plus[step],
            radius,
        )
        direction = 1.0 if step < half_steps else -1.0
        area[step + 1] = area[step] + direction * strip
    return area


def repeated_cycle_values(one_cycle: np.ndarray, step_count: int, period_steps: int) -> np.ndarray:
    return np.array([one_cycle[step % period_steps] for step in range(step_count + 1)], dtype=float)


def fit_power(cases: Sequence[Dict[str, Any]], candidate: str) -> Dict[str, Any]:
    selected = [row for row in cases if float(row[candidate]) > 1.0e-30]
    selected.sort(key=lambda row: float(row["radius"]))
    if len(selected) < 4:
        return {
            "candidate": candidate,
            "fit_valid": False,
            "case_count": len(selected),
            "loglog_slope": 0.0,
            "inverse_power_alpha": 0.0,
            "log_rmse": 0.0,
        }
    x = np.log(np.array([float(row["radius"]) for row in selected], dtype=float))
    y = np.log(np.array([float(row[candidate]) for row in selected], dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    return {
        "candidate": candidate,
        "fit_valid": True,
        "case_count": len(selected),
        "loglog_slope": float(slope),
        "inverse_power_alpha": -float(slope),
        "loglog_intercept": float(intercept),
        "log_rmse": float(np.sqrt(np.mean(residual * residual))),
        "min_value": float(np.min(np.exp(y))),
        "max_value": float(np.max(np.exp(y))),
    }


def run_case(
    radius_label_deg: float,
    protocol_name: str,
    params: Any,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    radius = math.radians(radius_label_deg)
    phase = params.omega_step * np.arange(params.step_count + 1, dtype=float)
    longitude_minus_rad = math.radians(LONGITUDE_MINUS_DEG)
    longitude_plus_rad = math.radians(LONGITUDE_PLUS_DEG)
    longitude_separation_rad = longitude_plus_rad - longitude_minus_rad

    path_minus_in = great_circle(radius, longitude_minus_rad, phase)
    path_plus_in = great_circle(radius, longitude_plus_rad, phase)
    protocol = protocol_by_name(protocol_name)
    path_minus, minus_unitarity_error, minus_internal_imag = scatter_spatial_path(
        path_minus_in,
        protocol,
    )
    path_plus, plus_unitarity_error, plus_internal_imag = scatter_spatial_path(
        path_plus_in,
        protocol,
    )

    radius_squared = radius * radius
    norm_minus = np.sum(path_minus * path_minus, axis=1)
    norm_plus = np.sum(path_plus * path_plus, axis=1)
    norm_minus_abs_error = np.abs(norm_minus - radius_squared)
    norm_plus_abs_error = np.abs(norm_plus - radius_squared)
    norm_minus_rel_error = norm_minus_abs_error / max(radius_squared, 1.0e-30)
    norm_plus_rel_error = norm_plus_abs_error / max(radius_squared, 1.0e-30)

    dd_minus = vector_second_difference(path_minus)
    dd_plus = vector_second_difference(path_plus)
    f_minus = np.linalg.norm(dd_minus, axis=1)
    f_plus = np.linalg.norm(dd_plus, axis=1)
    f_native_mean = 0.5 * (f_minus + f_plus)
    omega_discrete_sq = 4.0 * math.sin(params.omega_step / 2.0) ** 2
    f_center_control = omega_discrete_sq * radius
    f_projection_error = np.abs(f_native_mean - f_center_control)

    period_minus = path_minus[: params.period_steps + 1]
    period_plus = path_plus[: params.period_steps + 1]
    area_one_cycle = one_cycle_lune_area(
        period_minus,
        period_plus,
        radius,
        params.period_steps,
    )
    area = repeated_cycle_values(area_one_cycle, params.step_count, params.period_steps)
    analytic_lune_area = 2.0 * longitude_separation_rad * radius_squared
    max_area = float(np.max(area_one_cycle))
    area_analytic_rel_error = abs(max_area - analytic_lune_area) / max(analytic_lune_area, 1.0e-30)
    area_cycle_return_error = abs(float(area_one_cycle[-1]))

    separation = np.linalg.norm(path_plus - path_minus, axis=1)
    expected_max_separation = 2.0 * radius * math.sin(longitude_separation_rad / 2.0)
    node_mask = np.isclose(np.mod(np.arange(params.step_count + 1), params.period_steps // 2), 0)
    node_intersection_error = float(np.max(separation[node_mask]))

    cycle_differences_minus = path_minus[params.period_steps :] - path_minus[: -params.period_steps]
    cycle_differences_plus = path_plus[params.period_steps :] - path_plus[: -params.period_steps]
    cycle_return_error = max(
        float(np.max(np.linalg.norm(cycle_differences_minus, axis=1))),
        float(np.max(np.linalg.norm(cycle_differences_plus, axis=1))),
    )

    k_area_f = np.abs(area) * f_native_mean
    case_id = f"{protocol_name}_R{radius_label_deg:g}deg".replace(".", "p")
    series_rows: List[Dict[str, Any]] = []
    for step in range(params.step_count + 1):
        series_rows.append(
            {
                "case_id": case_id,
                "scattering_protocol": protocol_name,
                "radius_label_deg": radius_label_deg,
                "radius": radius,
                "radius_squared_expected": radius_squared,
                "longitude_minus_deg": LONGITUDE_MINUS_DEG,
                "longitude_plus_deg": LONGITUDE_PLUS_DEG,
                "longitude_separation_deg": LONGITUDE_SEPARATION_DEG,
                "readout_mode": "readout_off",
                "per_step_leak": 0.0,
                "step": step,
                "phase_rad": float(phase[step]),
                "minus_x": float(path_minus[step, 0]),
                "minus_y": float(path_minus[step, 1]),
                "minus_z": float(path_minus[step, 2]),
                "plus_x": float(path_plus[step, 0]),
                "plus_y": float(path_plus[step, 1]),
                "plus_z": float(path_plus[step, 2]),
                "R2_minus": float(norm_minus[step]),
                "R2_plus": float(norm_plus[step]),
                "R2_minus_abs_error": float(norm_minus_abs_error[step]),
                "R2_plus_abs_error": float(norm_plus_abs_error[step]),
                "R2_minus_rel_error": float(norm_minus_rel_error[step]),
                "R2_plus_rel_error": float(norm_plus_rel_error[step]),
                "path_separation": float(separation[step]),
                "A_meridian_lune": float(area[step]),
                "f_minus_second_difference": float(f_minus[step]),
                "f_plus_second_difference": float(f_plus[step]),
                "f_native_mean": float(f_native_mean[step]),
                "f_center_control": f_center_control,
                "f_projection_error": float(f_projection_error[step]),
                "K_area_times_f": float(k_area_f[step]),
                "inverse_term_used": False,
                "readout_loss_used": False,
                "time_axis_area_used": False,
            }
        )

    measured_slice = slice(1, -1)
    summary = {
        "case_id": case_id,
        "scattering_protocol": protocol_name,
        "radius_label_deg": radius_label_deg,
        "radius": radius,
        "radius_squared_expected": radius_squared,
        "longitude_minus_deg": LONGITUDE_MINUS_DEG,
        "longitude_plus_deg": LONGITUDE_PLUS_DEG,
        "longitude_separation_deg": LONGITUDE_SEPARATION_DEG,
        "step_count": params.step_count,
        "period_steps": params.period_steps,
        "max_R2_minus_abs_error": float(np.max(norm_minus_abs_error)),
        "max_R2_plus_abs_error": float(np.max(norm_plus_abs_error)),
        "max_R2_minus_rel_error": float(np.max(norm_minus_rel_error)),
        "max_R2_plus_rel_error": float(np.max(norm_plus_rel_error)),
        "max_cycle_return_error": cycle_return_error,
        "max_node_intersection_error": node_intersection_error,
        "max_path_separation": float(np.max(separation)),
        "expected_max_path_separation": expected_max_separation,
        "max_path_separation_abs_error": abs(float(np.max(separation)) - expected_max_separation),
        "max_A_meridian_lune": max_area,
        "analytic_A_meridian_lune": analytic_lune_area,
        "A_meridian_lune_analytic_rel_error": area_analytic_rel_error,
        "A_meridian_lune_cycle_return_error": area_cycle_return_error,
        "max_f_native_mean": float(np.max(f_native_mean[measured_slice])),
        "rms_f_native_mean": float(np.sqrt(np.mean(f_native_mean[measured_slice] ** 2))),
        "max_K_area_times_f": float(np.max(k_area_f)),
        "max_f_projection_error": float(np.max(f_projection_error[measured_slice])),
        "max_scattering_unitarity_error": max(minus_unitarity_error, plus_unitarity_error),
        "max_scattering_internal_imag": max(minus_internal_imag, plus_internal_imag),
        "R2_preserved": bool(
            float(np.max(norm_minus_rel_error)) <= 1.0e-12
            and float(np.max(norm_plus_rel_error)) <= 1.0e-12
        ),
        "cycle_returns": bool(cycle_return_error <= 1.0e-12),
        "nodes_intersect": bool(node_intersection_error <= 1.0e-12),
        "lune_area_detected": bool(max_area > params.area_tol),
        "inverse_term_used": False,
        "readout_loss_used": False,
        "time_axis_area_used": False,
    }
    return series_rows, summary


def build_fits(case_summaries: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    candidates = [
        "max_f_native_mean",
        "rms_f_native_mean",
        "max_A_meridian_lune",
        "max_K_area_times_f",
    ]
    rows: List[Dict[str, Any]] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        selected = [row for row in case_summaries if row["scattering_protocol"] == protocol_name]
        for candidate in candidates:
            fit = fit_power(selected, candidate)
            fit.update(
                {
                    "scattering_protocol": protocol_name,
                    "inverse_square_match": bool(
                        fit["fit_valid"]
                        and candidate in {"max_f_native_mean", "rms_f_native_mean"}
                        and abs(float(fit["inverse_power_alpha"]) - 2.0) <= INVERSE_POWER_TOL
                    ),
                }
            )
            rows.append(fit)
    return rows


def build_validation(
    case_summaries: Sequence[Dict[str, Any]],
    fits: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    native_fits = [
        row
        for row in fits
        if row["candidate"] in {"max_f_native_mean", "rms_f_native_mean"}
        and bool(row["fit_valid"])
    ]
    primary_fit = next(
        row
        for row in native_fits
        if row["scattering_protocol"] == "fermionic_reflection_pi"
        and row["candidate"] == "max_f_native_mean"
    )
    return {
        "experiment": "ab_two_body_meridian_lune_sweep_preliminary_v2",
        "source_path": str(SOURCE_PATH.relative_to(REPO_ROOT)),
        "case_count": len(case_summaries),
        "fit_count": len(fits),
        "longitude_minus_deg": LONGITUDE_MINUS_DEG,
        "longitude_plus_deg": LONGITUDE_PLUS_DEG,
        "longitude_separation_deg": LONGITUDE_SEPARATION_DEG,
        "R2_preserved_all_cases": all(bool(row["R2_preserved"]) for row in case_summaries),
        "max_R2_abs_error": max(
            max(float(row["max_R2_minus_abs_error"]), float(row["max_R2_plus_abs_error"]))
            for row in case_summaries
        ),
        "max_R2_rel_error": max(
            max(float(row["max_R2_minus_rel_error"]), float(row["max_R2_plus_rel_error"]))
            for row in case_summaries
        ),
        "cycle_returns_all_cases": all(bool(row["cycle_returns"]) for row in case_summaries),
        "max_cycle_return_error": max(float(row["max_cycle_return_error"]) for row in case_summaries),
        "nodes_intersect_all_cases": all(bool(row["nodes_intersect"]) for row in case_summaries),
        "max_node_intersection_error": max(
            float(row["max_node_intersection_error"]) for row in case_summaries
        ),
        "lune_area_detected_all_cases": all(
            bool(row["lune_area_detected"]) for row in case_summaries
        ),
        "max_lune_area_analytic_rel_error": max(
            float(row["A_meridian_lune_analytic_rel_error"]) for row in case_summaries
        ),
        "max_lune_area_cycle_return_error": max(
            float(row["A_meridian_lune_cycle_return_error"]) for row in case_summaries
        ),
        "max_scattering_unitarity_error": max(
            float(row["max_scattering_unitarity_error"]) for row in case_summaries
        ),
        "inverse_term_used_any": any(bool(row["inverse_term_used"]) for row in case_summaries),
        "readout_loss_used_any": any(bool(row["readout_loss_used"]) for row in case_summaries),
        "time_axis_area_used_any": any(bool(row["time_axis_area_used"]) for row in case_summaries),
        "native_inverse_square_detected": any(
            bool(row["inverse_square_match"]) for row in native_fits
        ),
        "primary_native_acceleration_fit": primary_fit,
    }


def make_plots(
    series_rows: Sequence[Dict[str, Any]],
    case_summaries: Sequence[Dict[str, Any]],
) -> None:
    selected = [
        row
        for row in series_rows
        if row["scattering_protocol"] == "fermionic_reflection_pi"
        and float(row["radius_label_deg"]) == 10.0
        and int(row["step"]) <= 96
    ]
    selected.sort(key=lambda row: int(row["step"]))
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(
        [float(row["minus_x"]) for row in selected],
        [float(row["minus_y"]) for row in selected],
        [float(row["minus_z"]) for row in selected],
        label="longitude -2.5 deg",
    )
    ax.plot(
        [float(row["plus_x"]) for row in selected],
        [float(row["plus_y"]) for row in selected],
        [float(row["plus_z"]) for row in selected],
        label="longitude +2.5 deg",
    )
    ax.set_title("two meridian great circles")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_meridian_great_circles_v2.png", dpi=180)
    plt.close(fig)

    selected_cases = [
        row
        for row in case_summaries
        if row["scattering_protocol"] == "fermionic_reflection_pi"
    ]
    selected_cases.sort(key=lambda row: float(row["radius"]))
    radii = [float(row["radius"]) for row in selected_cases]
    fig, ax = plt.subplots(figsize=(8, 5))
    for key, label in [
        ("max_f_native_mean", "native |d2 r|"),
        ("max_A_meridian_lune", "meridian lune area"),
        ("max_K_area_times_f", "area x native |d2 r|"),
    ]:
        ax.plot(radii, [float(row[key]) for row in selected_cases], marker="o", label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("meridian-lune native scaling")
    ax.set_xlabel("R")
    ax.set_ylabel("measured candidate")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_meridian_lune_scaling_v2.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any], fits: Sequence[Dict[str, Any]]) -> None:
    fermionic_fits = [
        row for row in fits if row["scattering_protocol"] == "fermionic_reflection_pi"
    ]
    fit_lines = "\n".join(
        f'- `{row["candidate"]}`: slope={float(row["loglog_slope"]):.12g}, '
        f'alpha={float(row["inverse_power_alpha"]):.12g}, '
        f'log_rmse={float(row["log_rmse"]):.3e}'
        for row in fermionic_fits
    )
    report = f"""# AB二体・経度±2.5度大円リューネ予備実験 v2

## 1. R²保存の確認

各経度大円について、全ステップで

```text
x² + y² + z² = R²
```

を直接検査した。

- R2_preserved_all_cases: `{validation["R2_preserved_all_cases"]}`
- max_R2_abs_error: `{validation["max_R2_abs_error"]:.16e}`
- max_R2_rel_error: `{validation["max_R2_rel_error"]:.16e}`
- cycle_returns_all_cases: `{validation["cycle_returns_all_cases"]}`
- max_cycle_return_error: `{validation["max_cycle_return_error"]:.16e}`

## 2. 二大円と交点

- longitude_minus_deg: `{validation["longitude_minus_deg"]}`
- longitude_plus_deg: `{validation["longitude_plus_deg"]}`
- longitude_separation_deg: `{validation["longitude_separation_deg"]}`
- nodes_intersect_all_cases: `{validation["nodes_intersect_all_cases"]}`
- max_node_intersection_error: `{validation["max_node_intersection_error"]:.16e}`

二大円は位相0度と180度で交差する。

## 3. 空間的な球面スイープ

- lune_area_detected_all_cases: `{validation["lune_area_detected_all_cases"]}`
- max_lune_area_analytic_rel_error: `{validation["max_lune_area_analytic_rel_error"]:.16e}`
- max_lune_area_cycle_return_error: `{validation["max_lune_area_cycle_return_error"]:.16e}`
- time_axis_area_used_any: `{validation["time_axis_area_used_any"]}`

面積はXT・YT投影ではなく、二つの実三成分大円が挟む球面積として計算した。

## 4. 距離指数

{fit_lines}

- native_inverse_square_detected: `{validation["native_inverse_square_detected"]}`
- inverse_term_used_any: `{validation["inverse_term_used_any"]}`
- readout_loss_used_any: `{validation["readout_loss_used_any"]}`

主要なネイティブ加速度フィット：

```json
{json.dumps(validation["primary_native_acceleration_fit"], ensure_ascii=False, indent=2)}
```
"""
    (OUT_DIR / "ab_two_body_meridian_lune_sweep_preliminary_report_v2.md").write_text(
        report,
        encoding="utf-8",
    )


def main() -> None:
    params = ACCEL.Params()
    series_rows: List[Dict[str, Any]] = []
    case_summaries: List[Dict[str, Any]] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        for radius_label_deg in RADIUS_LABEL_DEGS:
            rows, summary = run_case(radius_label_deg, protocol_name, params)
            series_rows.extend(rows)
            case_summaries.append(summary)

    fits = build_fits(case_summaries)
    validation = build_validation(case_summaries, fits)
    result = {
        "params": {
            "step_count": params.step_count,
            "period_steps": params.period_steps,
            "omega_step": params.omega_step,
            "radius_label_degs": RADIUS_LABEL_DEGS,
            "longitude_minus_deg": LONGITUDE_MINUS_DEG,
            "longitude_plus_deg": LONGITUDE_PLUS_DEG,
            "longitude_separation_deg": LONGITUDE_SEPARATION_DEG,
            "scattering_protocols": SCATTERING_PROTOCOL_NAMES,
        },
        "validation": validation,
        "case_summaries": case_summaries,
        "power_fits": fits,
        "notes": [
            "The area is a spatial spherical lune between two meridian great circles.",
            "No XT or YT projected area is used.",
            "No reciprocal or inverse-square term is used.",
        ],
    }

    write_csv(OUT_DIR / "ab_two_body_meridian_lune_series_v2.csv", series_rows)
    write_csv(OUT_DIR / "ab_two_body_meridian_lune_case_summary_v2.csv", case_summaries)
    write_csv(OUT_DIR / "ab_two_body_meridian_lune_power_fits_v2.csv", fits)
    (OUT_DIR / "ab_two_body_meridian_lune_sweep_preliminary_result_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(series_rows, case_summaries)
    write_report(validation, fits)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
