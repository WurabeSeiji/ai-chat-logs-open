from __future__ import annotations

"""AB二体加速度実験へ非直交な第二位置位相関係を加える予備実験。

既存のAB散乱行列と離散二階差分は変更しない。対称な関係を delta=0 deg
とし、そこから外れた delta=+/-5 deg の第二位置位相 Y を加える。

    Z_X = Z / sqrt(2)
    Z_Y = exp(i delta) Z / sqrt(2)

この分配は |Z_X|^2 + |Z_Y|^2 = |Z|^2 を保存する。逆数、逆面積、
逆二乗項は運動にも候補量にも入力しない。XY有向面積、実ランク、
二階差分、保存誤差を出力し、初期位置位相偏差 L に対する冪指数を測る。
"""

import csv
import importlib.util
import json
import math
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, Iterable, List, Sequence

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_nonorthogonal_xy_phase_preliminary_result_v1"
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

INITIAL_DEVIATION_DEGS = [1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 60.0]
RELATION_OFFSET_DEGS = [-5.0, 0.0, 5.0]
SCATTERING_PROTOCOL_NAMES = ["pass_through", "fermionic_reflection_pi"]
POWER_ALPHA_TOL = 0.12


def load_acceleration_source(path: Path) -> ModuleType:
    module_name = "ab_acceleration_reference_v4"
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


def signed_second_difference(values: np.ndarray) -> np.ndarray:
    out = np.zeros_like(values, dtype=float)
    if len(values) < 3:
        return out
    out[1:-1] = values[2:] - 2.0 * values[1:-1] + values[:-2]
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def fit_power(cases: Sequence[Dict[str, Any]], candidate: str) -> Dict[str, Any]:
    selected = [row for row in cases if float(row[candidate]) > 1.0e-30]
    selected.sort(key=lambda row: float(row["initial_deviation_rad"]))
    if len(selected) < 4:
        return {
            "candidate": candidate,
            "fit_valid": False,
            "case_count": len(selected),
            "loglog_slope": 0.0,
            "inverse_power_alpha": 0.0,
            "log_rmse": 0.0,
        }
    x = np.log(np.array([float(row["initial_deviation_rad"]) for row in selected]))
    y = np.log(np.array([float(row[candidate]) for row in selected]))
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


def protocol_by_name(name: str) -> Any:
    for protocol in ACCEL.SCATTERING_PROTOCOLS:
        if protocol.name == name:
            return protocol
    raise ValueError(f"unknown scattering protocol: {name}")


def run_case(
    initial_deviation_deg: float,
    relation_offset_deg: float,
    protocol_name: str,
    params: Any,
    readout_off: Any,
) -> tuple[List[Dict[str, Any]], Dict[str, Any]]:
    initial_deviation_rad = math.radians(initial_deviation_deg)
    relation_offset_rad = math.radians(relation_offset_deg)
    protocol = protocol_by_name(protocol_name)

    chi_base, complement_base = ACCEL.base_rotation(
        initial_deviation_rad,
        readout_off,
        params,
    )
    z_base = chi_base + 1j * complement_base
    z_x_in = z_base / math.sqrt(2.0)
    z_y_in = np.exp(1j * relation_offset_rad) * z_base / math.sqrt(2.0)

    scatter_x = ACCEL.two_channel_scattering_state(
        protocol,
        np.real(z_x_in),
        np.imag(z_x_in),
    )
    scatter_y = ACCEL.two_channel_scattering_state(
        protocol,
        np.real(z_y_in),
        np.imag(z_y_in),
    )
    z_x = np.asarray(scatter_x["relative_out"], dtype=complex)
    z_y = np.asarray(scatter_y["relative_out"], dtype=complex)
    x_read = np.real(z_x)
    y_read = np.real(z_y)

    ddx = signed_second_difference(x_read)
    ddy = signed_second_difference(y_read)
    f_xy_second_difference = np.hypot(ddx, ddy)

    omega_discrete_sq = 4.0 * math.sin(params.omega_step / 2.0) ** 2
    f_x_center = omega_discrete_sq * np.abs(x_read)
    f_y_center = omega_discrete_sq * np.abs(y_read)
    f_xy_center = np.hypot(f_x_center, f_y_center)

    area = ACCEL.rolling_closed_area(x_read, y_read, params.period_steps)
    rank = ACCEL.rolling_rank(
        x_read,
        y_read,
        params.period_steps,
        True,
        params.rank_tol,
    )
    relation_norm = np.abs(z_x) ** 2 + np.abs(z_y) ** 2
    expected_norm = initial_deviation_rad**2
    norm_abs_error = np.abs(relation_norm - expected_norm)
    norm_rel_error = norm_abs_error / max(expected_norm, 1.0e-30)
    scattering_unitarity_error = np.maximum(
        np.asarray(scatter_x["scattering_unitarity_error"], dtype=float),
        np.asarray(scatter_y["scattering_unitarity_error"], dtype=float),
    )
    projection_error = np.abs(f_xy_second_difference - f_xy_center)
    k_area_f = np.abs(area) * f_xy_second_difference

    case_id = (
        f"{protocol_name}_L{initial_deviation_deg:g}_delta{relation_offset_deg:+g}"
        .replace("+", "p")
        .replace("-", "m")
        .replace(".", "p")
    )
    series_rows: List[Dict[str, Any]] = []
    for step in range(params.step_count + 1):
        series_rows.append(
            {
                "case_id": case_id,
                "scattering_protocol": protocol_name,
                "initial_deviation_deg": initial_deviation_deg,
                "initial_deviation_rad": initial_deviation_rad,
                "relation_offset_deg": relation_offset_deg,
                "relation_offset_rad": relation_offset_rad,
                "readout_mode": "readout_off",
                "per_step_leak": 0.0,
                "step": step,
                "x_read": float(x_read[step]),
                "y_read": float(y_read[step]),
                "x_internal_imag": float(np.imag(z_x[step])),
                "y_internal_imag": float(np.imag(z_y[step])),
                "ddx": float(ddx[step]),
                "ddy": float(ddy[step]),
                "f_xy_second_difference": float(f_xy_second_difference[step]),
                "f_xy_center_control": float(f_xy_center[step]),
                "f_projection_error": float(projection_error[step]),
                "A_XY": float(area[step]),
                "abs_A_XY": abs(float(area[step])),
                "rank_XY": int(rank[step]),
                "relation_norm": float(relation_norm[step]),
                "relation_norm_abs_error": float(norm_abs_error[step]),
                "relation_norm_rel_error": float(norm_rel_error[step]),
                "scattering_unitarity_error": float(scattering_unitarity_error[step]),
                "K_area_times_f": float(k_area_f[step]),
                "inverse_term_used": False,
                "readout_loss_used": False,
                "external_axis_used": False,
            }
        )

    measurement_slice = slice(params.period_steps, None)
    measured_area = area[measurement_slice]
    measured_rank = rank[measurement_slice]
    measured_f = f_xy_second_difference[measurement_slice]
    measured_k = k_area_f[measurement_slice]
    signed_area_median = float(np.median(measured_area))
    summary = {
        "case_id": case_id,
        "scattering_protocol": protocol_name,
        "initial_deviation_deg": initial_deviation_deg,
        "initial_deviation_rad": initial_deviation_rad,
        "relation_offset_deg": relation_offset_deg,
        "relation_offset_rad": relation_offset_rad,
        "readout_mode": "readout_off",
        "per_step_leak": 0.0,
        "step_count": params.step_count,
        "rank_XY_max": int(np.max(measured_rank)),
        "max_abs_A_XY": float(np.max(np.abs(measured_area))),
        "signed_A_XY_median": signed_area_median,
        "signed_A_XY_orientation": int(np.sign(signed_area_median)),
        "max_f_xy_second_difference": float(np.max(measured_f)),
        "rms_f_xy_second_difference": float(np.sqrt(np.mean(measured_f**2))),
        "max_K_area_times_f": float(np.max(measured_k)),
        "max_relation_norm_abs_error": float(np.max(norm_abs_error)),
        "max_relation_norm_rel_error": float(np.max(norm_rel_error)),
        "max_scattering_unitarity_error": float(np.max(scattering_unitarity_error)),
        "max_f_projection_error": float(np.max(projection_error[1:-1])),
        "area_detected": bool(
            float(np.max(np.abs(measured_area))) > params.area_tol
            and int(np.max(measured_rank)) >= 2
        ),
        "inverse_term_used": False,
        "readout_loss_used": False,
        "external_axis_used": False,
    }
    return series_rows, summary


def build_fits(case_summaries: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    candidates = [
        "max_f_xy_second_difference",
        "rms_f_xy_second_difference",
        "max_abs_A_XY",
        "max_K_area_times_f",
    ]
    rows: List[Dict[str, Any]] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        for offset_deg in RELATION_OFFSET_DEGS:
            selected = [
                row
                for row in case_summaries
                if row["scattering_protocol"] == protocol_name
                and float(row["relation_offset_deg"]) == offset_deg
            ]
            for candidate in candidates:
                fit = fit_power(selected, candidate)
                fit.update(
                    {
                        "scattering_protocol": protocol_name,
                        "relation_offset_deg": offset_deg,
                        "inverse_square_match": bool(
                            fit["fit_valid"]
                            and candidate
                            in {
                                "max_f_xy_second_difference",
                                "rms_f_xy_second_difference",
                            }
                            and abs(float(fit["inverse_power_alpha"]) - 2.0)
                            <= POWER_ALPHA_TOL
                        ),
                    }
                )
                rows.append(fit)
    return rows


def antisymmetry_error(case_summaries: Sequence[Dict[str, Any]]) -> float:
    errors: List[float] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        for deviation_deg in INITIAL_DEVIATION_DEGS:
            plus = next(
                row
                for row in case_summaries
                if row["scattering_protocol"] == protocol_name
                and float(row["initial_deviation_deg"]) == deviation_deg
                and float(row["relation_offset_deg"]) == 5.0
            )
            minus = next(
                row
                for row in case_summaries
                if row["scattering_protocol"] == protocol_name
                and float(row["initial_deviation_deg"]) == deviation_deg
                and float(row["relation_offset_deg"]) == -5.0
            )
            a_plus = float(plus["signed_A_XY_median"])
            a_minus = float(minus["signed_A_XY_median"])
            scale = max(abs(a_plus), abs(a_minus), 1.0e-30)
            errors.append(abs(a_plus + a_minus) / scale)
    return max(errors) if errors else 0.0


def build_validation(
    case_summaries: Sequence[Dict[str, Any]],
    fits: Sequence[Dict[str, Any]],
    params: Any,
) -> Dict[str, Any]:
    zero = [row for row in case_summaries if float(row["relation_offset_deg"]) == 0.0]
    nonzero = [row for row in case_summaries if float(row["relation_offset_deg"]) != 0.0]
    native_fits = [
        row
        for row in fits
        if row["candidate"]
        in {"max_f_xy_second_difference", "rms_f_xy_second_difference"}
        and bool(row["fit_valid"])
    ]
    closest = min(
        native_fits,
        key=lambda row: abs(float(row["inverse_power_alpha"]) - 2.0),
    )
    return {
        "experiment": "ab_two_body_nonorthogonal_xy_phase_preliminary_v1",
        "source_path": str(SOURCE_PATH.relative_to(REPO_ROOT)),
        "case_count": len(case_summaries),
        "fit_count": len(fits),
        "readout_loss_used_any": any(bool(row["readout_loss_used"]) for row in case_summaries),
        "inverse_term_used_any": any(bool(row["inverse_term_used"]) for row in case_summaries),
        "external_axis_used_any": any(bool(row["external_axis_used"]) for row in case_summaries),
        "zero_offset_rank1_all": all(int(row["rank_XY_max"]) == 1 for row in zero),
        "zero_offset_area_zero_all": all(
            float(row["max_abs_A_XY"]) <= params.area_tol for row in zero
        ),
        "nonzero_offset_rank2_all": all(int(row["rank_XY_max"]) == 2 for row in nonzero),
        "nonzero_offset_area_detected_all": all(bool(row["area_detected"]) for row in nonzero),
        "signed_area_plus_minus_antisymmetry_max_rel_error": antisymmetry_error(
            case_summaries
        ),
        "max_relation_norm_rel_error": max(
            float(row["max_relation_norm_rel_error"]) for row in case_summaries
        ),
        "max_scattering_unitarity_error": max(
            float(row["max_scattering_unitarity_error"]) for row in case_summaries
        ),
        "native_inverse_square_detected": any(
            bool(row["inverse_square_match"]) for row in native_fits
        ),
        "closest_native_inverse_square_fit": closest,
    }


def make_plots(
    series_rows: Sequence[Dict[str, Any]],
    case_summaries: Sequence[Dict[str, Any]],
    fits: Sequence[Dict[str, Any]],
) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    for offset_deg in RELATION_OFFSET_DEGS:
        rows = [
            row
            for row in series_rows
            if row["scattering_protocol"] == "fermionic_reflection_pi"
            and float(row["initial_deviation_deg"]) == 10.0
            and float(row["relation_offset_deg"]) == offset_deg
        ]
        rows.sort(key=lambda row: int(row["step"]))
        ax.plot(
            [float(row["x_read"]) for row in rows],
            [float(row["y_read"]) for row in rows],
            label=f"delta={offset_deg:+g} deg",
        )
    ax.set_title("nonorthogonal XY position-phase paths")
    ax.set_xlabel("X position-phase readout")
    ax.set_ylabel("Y position-phase readout")
    ax.axis("equal")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_nonorthogonal_xy_paths_v1.png", dpi=180)
    plt.close(fig)

    selected = [
        row
        for row in case_summaries
        if row["scattering_protocol"] == "fermionic_reflection_pi"
        and float(row["relation_offset_deg"]) == 5.0
    ]
    selected.sort(key=lambda row: float(row["initial_deviation_rad"]))
    x = [float(row["initial_deviation_rad"]) for row in selected]
    fig, ax = plt.subplots(figsize=(8, 5))
    for key, label in [
        ("max_f_xy_second_difference", "native |d2 XY|"),
        ("max_abs_A_XY", "|A_XY|"),
        ("max_K_area_times_f", "|A_XY| |d2 XY|"),
    ]:
        ax.plot(x, [float(row[key]) for row in selected], marker="o", label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("native scaling at delta=+5 deg")
    ax.set_xlabel("initial position-phase deviation L [rad]")
    ax.set_ylabel("measured candidate")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_nonorthogonal_xy_scaling_v1.png", dpi=180)
    plt.close(fig)

    native = [
        row
        for row in fits
        if row["candidate"] == "max_f_xy_second_difference" and bool(row["fit_valid"])
    ]
    labels = [
        f'{row["scattering_protocol"]}\ndelta={float(row["relation_offset_deg"]):+g}'
        for row in native
    ]
    alphas = [float(row["inverse_power_alpha"]) for row in native]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(range(len(native)), alphas)
    ax.axhline(2.0, color="black", linestyle="--", linewidth=1.0, label="inverse square")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xticks(range(len(native)), labels, rotation=30, ha="right")
    ax.set_ylabel("inverse-power alpha")
    ax.set_title("native acceleration power scan")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_nonorthogonal_xy_native_alpha_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any], fits: Sequence[Dict[str, Any]]) -> None:
    closest = validation["closest_native_inverse_square_fit"]
    selected_fits = [
        row
        for row in fits
        if row["scattering_protocol"] == "fermionic_reflection_pi"
        and float(row["relation_offset_deg"]) == 5.0
        and row["candidate"]
        in {
            "max_f_xy_second_difference",
            "max_abs_A_XY",
            "max_K_area_times_f",
        }
    ]
    fit_lines = "\n".join(
        f'- `{row["candidate"]}`: slope={float(row["loglog_slope"]):.12g}, '
        f'alpha={float(row["inverse_power_alpha"]):.12g}, '
        f'log_rmse={float(row["log_rmse"]):.3e}'
        for row in selected_fits
    )
    report = f"""# AB二体・非直交XY位置位相予備実験 v1

## 実装

既存の `run_ab_two_body_fermionic_reflection_harmonic_readout_v4.py` から、
AB基礎回転、二チャネル散乱行列、面積・実ランク計算をそのまま使用した。

対称関係を `delta=0 deg` とし、第二位置位相を
`delta=-5 deg, +5 deg` だけ外した。二位置位相には総振幅を等分し、
読出し損失、外部軸、逆数、逆面積、逆二乗項は使用していない。

## 判別結果

- case_count: `{validation["case_count"]}`
- zero_offset_rank1_all: `{validation["zero_offset_rank1_all"]}`
- zero_offset_area_zero_all: `{validation["zero_offset_area_zero_all"]}`
- nonzero_offset_rank2_all: `{validation["nonzero_offset_rank2_all"]}`
- nonzero_offset_area_detected_all: `{validation["nonzero_offset_area_detected_all"]}`
- signed_area_plus_minus_antisymmetry_max_rel_error: `{validation["signed_area_plus_minus_antisymmetry_max_rel_error"]:.16e}`
- max_relation_norm_rel_error: `{validation["max_relation_norm_rel_error"]:.16e}`
- max_scattering_unitarity_error: `{validation["max_scattering_unitarity_error"]:.16e}`
- native_inverse_square_detected: `{validation["native_inverse_square_detected"]}`

## fermionic reflection, delta=+5 deg の冪指数

{fit_lines}

## 読み

`delta=0 deg` は一次元へ退化し、XY有向面積を生成しない。
`delta=+/-5 deg` は実ランク2と非零の有向面積を生成し、その符号は反転する。

この最小追加だけで既存の二階差分が逆二乗へ変化したかどうかは、
`native_inverse_square_detected` で判定する。最も指数2に近いネイティブ結果は次である。

```json
{json.dumps(closest, ensure_ascii=False, indent=2)}
```
"""
    (OUT_DIR / "ab_two_body_nonorthogonal_xy_phase_preliminary_report_v1.md").write_text(
        report,
        encoding="utf-8",
    )


def main() -> None:
    params = ACCEL.Params()
    readout_off = next(mode for mode in ACCEL.READOUT_MODES if mode.name == "readout_off")
    series_rows: List[Dict[str, Any]] = []
    case_summaries: List[Dict[str, Any]] = []
    for protocol_name in SCATTERING_PROTOCOL_NAMES:
        for relation_offset_deg in RELATION_OFFSET_DEGS:
            for initial_deviation_deg in INITIAL_DEVIATION_DEGS:
                rows, summary = run_case(
                    initial_deviation_deg,
                    relation_offset_deg,
                    protocol_name,
                    params,
                    readout_off,
                )
                series_rows.extend(rows)
                case_summaries.append(summary)

    fits = build_fits(case_summaries)
    validation = build_validation(case_summaries, fits, params)
    result = {
        "params": {
            "step_count": params.step_count,
            "period_steps": params.period_steps,
            "omega_step": params.omega_step,
            "rank_tol": params.rank_tol,
            "area_tol": params.area_tol,
            "initial_deviation_degs": INITIAL_DEVIATION_DEGS,
            "relation_offset_degs": RELATION_OFFSET_DEGS,
            "scattering_protocols": SCATTERING_PROTOCOL_NAMES,
        },
        "validation": validation,
        "case_summaries": case_summaries,
        "power_fits": fits,
        "notes": [
            "delta is the relational phase offset from the symmetric X=Y state, not a pre-existing metric angle.",
            "X and Y receive the same unchanged AB scattering map.",
            "No reciprocal or inverse-square term is used.",
        ],
    }

    write_csv(
        OUT_DIR / "ab_two_body_nonorthogonal_xy_phase_series_v1.csv",
        series_rows,
    )
    write_csv(
        OUT_DIR / "ab_two_body_nonorthogonal_xy_phase_case_summary_v1.csv",
        case_summaries,
    )
    write_csv(
        OUT_DIR / "ab_two_body_nonorthogonal_xy_phase_power_fits_v1.csv",
        fits,
    )
    (OUT_DIR / "ab_two_body_nonorthogonal_xy_phase_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(series_rows, case_summaries, fits)
    write_report(validation, fits)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
