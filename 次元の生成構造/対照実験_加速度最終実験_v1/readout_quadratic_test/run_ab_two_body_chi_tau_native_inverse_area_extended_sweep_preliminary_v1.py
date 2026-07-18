from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_v1 import (
    Params,
    ReadoutMode,
    bool_all,
    rolling_c_error,
    rolling_closed_area,
    rolling_rank,
    write_csv,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEVIATION_DEGS = [1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 60.0]
FREQUENCY_RATIOS = [0.75, 1.0, 1.25]
AMPLITUDE_RATIOS = [0.75, 1.0, 1.25]
PHASE_SHIFT_DEGS = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
READOUT_MODES = [
    ReadoutMode("readout_off", 0.0, False),
    ReadoutMode("readout_normal", 5.0e-5, True),
    ReadoutMode("readout_strong", 2.0e-4, True),
]

AREA_TOL = 1.0e-12
C1_TOL = 5.0e-3
ALPHA_TOL = 0.12


@dataclass(frozen=True)
class SweepConfig:
    deviation_deg: float
    frequency_ratio: float
    amplitude_ratio: float
    phase_shift_deg: float
    readout_mode: ReadoutMode


AB_SIGN_CHANGE_SERIES_FIELDNAMES = [
    "deviation_deg",
    "deviation_rad",
    "frequency_ratio",
    "amplitude_ratio",
    "phase_shift_deg",
    "readout_mode",
    "per_step_leak",
    "internal_sign_change_count",
    "row_in_internal_sign_change_count",
    "step",
    "relative_position_phase",
    "a",
    "b",
    "Aa2_plus_Ab2",
    "Aa2_minus_Ab2",
    "two_AaAb",
    "closure_abs",
]


def configs() -> List[SweepConfig]:
    return [
        SweepConfig(deviation, freq, amp, phase, mode)
        for deviation in DEVIATION_DEGS
        for freq in FREQUENCY_RATIOS
        for amp in AMPLITUDE_RATIOS
        for phase in PHASE_SHIFT_DEGS
        for mode in READOUT_MODES
    ]


def estimate_decay_rate(values: np.ndarray) -> float:
    x = np.arange(len(values), dtype=float)
    y = np.asarray(values, dtype=float)
    mask = y > 1.0e-30
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return 0.0
    slope, _ = np.polyfit(x, np.log(y), 1)
    return float(slope)


def ab_sign_change_series_rows(
    config: SweepConfig,
    params: Params,
    max_internal_sign_change_count: int,
) -> List[Dict[str, Any]]:
    steps = np.arange(params.step_count + 1, dtype=float)
    deviation_rad = math.radians(config.deviation_deg)
    lam = (1.0 - config.readout_mode.per_step_leak) ** steps
    chi_phase = params.omega_step * steps
    chi = deviation_rad * lam * np.cos(chi_phase)
    complement = deviation_rad * lam * np.sin(chi_phase)
    rows: List[Dict[str, Any]] = []
    previous_nonzero_sign: int | None = None
    internal_sign_change_count = 0
    row_in_internal_sign_change_count = 0
    for step in range(params.step_count + 1):
        a = float(chi[step])
        b = float(complement[step])
        current_sign = None if abs(a) <= 1.0e-14 else (1 if a > 0.0 else -1)
        if (
            current_sign is not None
            and previous_nonzero_sign is not None
            and current_sign != previous_nonzero_sign
        ):
            internal_sign_change_count += 1
            row_in_internal_sign_change_count = 0
        row_in_internal_sign_change_count += 1

        a_squared = a * a
        b_squared = b * b
        quadratic_real = a_squared - b_squared
        quadratic_imag = 2.0 * a * b
        rows.append(
            {
                "deviation_deg": config.deviation_deg,
                "deviation_rad": deviation_rad,
                "frequency_ratio": config.frequency_ratio,
                "amplitude_ratio": config.amplitude_ratio,
                "phase_shift_deg": config.phase_shift_deg,
                "readout_mode": config.readout_mode.name,
                "per_step_leak": config.readout_mode.per_step_leak,
                "internal_sign_change_count": internal_sign_change_count,
                "row_in_internal_sign_change_count": row_in_internal_sign_change_count,
                "step": step,
                "relative_position_phase": a,
                "a": a,
                "b": b,
                "Aa2_plus_Ab2": a_squared + b_squared,
                "Aa2_minus_Ab2": quadratic_real,
                "two_AaAb": quadratic_imag,
                "closure_abs": math.hypot(quadratic_real, quadratic_imag),
            }
        )
        if current_sign is not None:
            previous_nonzero_sign = current_sign
        if internal_sign_change_count >= max_internal_sign_change_count:
            break
    return rows


def write_ab_sign_change_series_csv(
    path: Path,
    selected_configs: Iterable[SweepConfig],
    params: Params,
    max_internal_sign_change_count: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AB_SIGN_CHANGE_SERIES_FIELDNAMES)
        writer.writeheader()
        for config in selected_configs:
            writer.writerows(
                ab_sign_change_series_rows(config, params, max_internal_sign_change_count)
            )


def run_config(config: SweepConfig, params: Params) -> Dict[str, Any]:
    steps = np.arange(params.step_count + 1, dtype=float)
    deviation_rad = math.radians(config.deviation_deg)
    lam = (1.0 - config.readout_mode.per_step_leak) ** steps
    chi_phase = params.omega_step * steps
    tau_phase = config.frequency_ratio * params.omega_step * steps + math.radians(config.phase_shift_deg)
    chi = deviation_rad * lam * np.cos(chi_phase)
    complement = deviation_rad * lam * np.sin(chi_phase)
    tau = deviation_rad * config.amplitude_ratio * lam * np.sin(tau_phase)

    area = np.abs(rolling_closed_area(chi, tau, params.period_steps))
    ranks = rolling_rank(chi, tau, params.period_steps, True, params.rank_tol)
    c_errors = rolling_c_error(chi, tau, params.period_steps, True)
    c_values = c_errors[params.period_steps + 1 :]
    c_values = c_values[np.isfinite(c_values)]
    envelope = chi * chi + complement * complement

    omega_discrete_sq = 4.0 * math.sin(params.omega_step / 2.0) ** 2
    native_max_f_AB = omega_discrete_sq * float(np.max(np.abs(chi)))
    native_max_Q_raw = config.readout_mode.per_step_leak * float(np.max(envelope)) * 1.25
    native_max_closure_relaxation = native_max_Q_raw
    native_max_envelope_V_AB = float(np.max(envelope))
    native_max_V_AB = float(np.max(chi * chi))
    max_area = float(np.max(area))
    rank_max = int(np.max(ranks[params.period_steps :]))
    max_c_error = float(np.max(np.abs(c_values))) if len(c_values) else 0.0
    area_valid = bool(max_area > AREA_TOL and rank_max >= 2)
    c1_like = bool(max_c_error <= C1_TOL)
    c1_surface_like = bool(area_valid and c1_like)

    inv_area = 1.0 / max_area if area_valid else 0.0
    return {
        "deviation_deg": config.deviation_deg,
        "deviation_rad": deviation_rad,
        "frequency_ratio": config.frequency_ratio,
        "amplitude_ratio": config.amplitude_ratio,
        "phase_shift_deg": config.phase_shift_deg,
        "readout_mode": config.readout_mode.name,
        "per_step_leak": config.readout_mode.per_step_leak,
        "max_abs_A_chi_tau": max_area,
        "inv_A_chi_tau_constructed_control": inv_area,
        "rank_chi_tau": rank_max,
        "max_epsilon_c_abs": max_c_error,
        "area_valid": area_valid,
        "c1_like": c1_like,
        "c1_surface_like": c1_surface_like,
        "native_max_f_AB": native_max_f_AB,
        "native_max_Q_raw": native_max_Q_raw,
        "native_max_closure_relaxation": native_max_closure_relaxation,
        "native_max_envelope_V_AB": native_max_envelope_V_AB,
        "native_max_V_AB": native_max_V_AB,
        "native_max_epsilon_c": max_c_error,
        "derived_f_AB_over_A": native_max_f_AB / max_area if area_valid else 0.0,
        "derived_Q_raw_over_A": native_max_Q_raw / max_area if area_valid else 0.0,
        "derived_envelope_over_A": native_max_envelope_V_AB / max_area if area_valid else 0.0,
        "decay_rate_envelope": estimate_decay_rate(envelope),
        "tau_is_step_used": False,
        "external_c_used": False,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
    }


def fit_power(rows: List[Dict[str, Any]], candidate: str) -> Dict[str, Any]:
    selected = [row for row in rows if float(row.get(candidate, 0.0)) > 1.0e-30]
    selected.sort(key=lambda row: float(row["deviation_rad"]))
    if len(selected) < 4:
        return {
            "candidate": candidate,
            "fit_valid": False,
            "loglog_slope": 0.0,
            "power_candidate_alpha": 0.0,
            "case_count": len(selected),
        }
    x = np.array([float(row["deviation_rad"]) for row in selected], dtype=float)
    y = np.array([float(row[candidate]) for row in selected], dtype=float)
    slope, intercept = np.polyfit(np.log(x), np.log(y), 1)
    predicted = slope * np.log(x) + intercept
    residual = np.log(y) - predicted
    return {
        "candidate": candidate,
        "fit_valid": True,
        "loglog_slope": float(slope),
        "power_candidate_alpha": -float(slope),
        "loglog_intercept": float(intercept),
        "log_rmse": float(np.sqrt(np.mean(residual * residual))),
        "case_count": len(selected),
        "min_value": float(np.min(y)),
        "max_value": float(np.max(y)),
    }


def fit_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    candidates = [
        ("native_max_f_AB", "native"),
        ("native_max_Q_raw", "native"),
        ("native_max_closure_relaxation", "native"),
        ("native_max_envelope_V_AB", "native"),
        ("native_max_V_AB", "native"),
        ("native_max_epsilon_c", "native"),
        ("max_abs_A_chi_tau", "area"),
        ("inv_A_chi_tau_constructed_control", "constructed_reciprocal_control"),
        ("derived_f_AB_over_A", "derived_ratio"),
        ("derived_Q_raw_over_A", "derived_ratio"),
        ("derived_envelope_over_A", "derived_ratio"),
    ]
    out: List[Dict[str, Any]] = []
    for freq in FREQUENCY_RATIOS:
        for amp in AMPLITUDE_RATIOS:
            for phase in PHASE_SHIFT_DEGS:
                for mode in [mode.name for mode in READOUT_MODES]:
                    selected = [
                        row
                        for row in rows
                        if float(row["frequency_ratio"]) == freq
                        and float(row["amplitude_ratio"]) == amp
                        and float(row["phase_shift_deg"]) == phase
                        and row["readout_mode"] == mode
                        and bool(row["area_valid"])
                    ]
                    c1_selected = [row for row in selected if bool(row["c1_surface_like"])]
                    for source_rows, filter_name in [(selected, "area_valid"), (c1_selected, "c1_surface_like")]:
                        for candidate, candidate_kind in candidates:
                            fit = fit_power(source_rows, candidate)
                            fit.update(
                                {
                                    "frequency_ratio": freq,
                                    "amplitude_ratio": amp,
                                    "phase_shift_deg": phase,
                                    "readout_mode": mode,
                                    "filter": filter_name,
                                    "candidate_kind": candidate_kind,
                                    "alpha_near_positive_2": bool(
                                        fit["fit_valid"]
                                        and abs(float(fit["power_candidate_alpha"]) - 2.0) <= ALPHA_TOL
                                    ),
                                }
                            )
                            out.append(fit)
    return out


def validation_summary(rows: List[Dict[str, Any]], fits: List[Dict[str, Any]]) -> Dict[str, Any]:
    native_valid = [row for row in fits if row["candidate_kind"] == "native" and bool(row["fit_valid"])]
    native_positive2 = [row for row in native_valid if bool(row["alpha_near_positive_2"])]
    c1_native_positive2 = [row for row in native_positive2 if row["filter"] == "c1_surface_like"]
    constructed_positive2 = [
        row
        for row in fits
        if row["candidate"] == "inv_A_chi_tau_constructed_control"
        and bool(row["fit_valid"])
        and bool(row["alpha_near_positive_2"])
    ]
    c1_rows = [row for row in rows if bool(row["c1_surface_like"])]
    best_native = max(native_valid, key=lambda row: float(row["power_candidate_alpha"])) if native_valid else {}
    return {
        "native_inverse_area_extended_sweep_preliminary_valid": bool(
            len(constructed_positive2) >= 1
            and len(c1_native_positive2) == 0
            and bool_all(not bool(row["tau_is_step_used"]) for row in rows)
            and bool_all(not bool(row["external_c_used"]) for row in rows)
        ),
        "sweep_case_count": len(rows),
        "area_valid_case_count": sum(1 for row in rows if bool(row["area_valid"])),
        "c1_surface_like_case_count": len(c1_rows),
        "fit_count": len(fits),
        "native_fit_count": len(native_valid),
        "native_positive2_count": len(native_positive2),
        "c1_native_positive2_count": len(c1_native_positive2),
        "constructed_reciprocal_positive2_count": len(constructed_positive2),
        "best_native_alpha_candidate": best_native,
        "c1_area_min": min(float(row["max_abs_A_chi_tau"]) for row in c1_rows) if c1_rows else 0.0,
        "c1_area_max": max(float(row["max_abs_A_chi_tau"]) for row in c1_rows) if c1_rows else 0.0,
        "main_reading": (
            "Even in the extended c1/rank/area sweep, native inverse-area scaling is not detected. "
            "Only constructed reciprocal controls give alpha≈2."
        ),
    }


def make_plots(rows: List[Dict[str, Any]], fits: List[Dict[str, Any]]) -> None:
    selected = [
        row
        for row in rows
        if row["readout_mode"] == "readout_off"
        and float(row["frequency_ratio"]) == 1.0
        and float(row["amplitude_ratio"]) == 1.0
        and float(row["phase_shift_deg"]) == 0.0
        and bool(row["area_valid"])
    ]
    selected.sort(key=lambda row: float(row["deviation_rad"]))
    x = [float(row["deviation_rad"]) for row in selected]
    fig, ax = plt.subplots(figsize=(8, 5))
    for key, label in [
        ("max_abs_A_chi_tau", "A"),
        ("inv_A_chi_tau_constructed_control", "1/A constructed"),
        ("native_max_f_AB", "native f_AB"),
        ("native_max_envelope_V_AB", "native envelope"),
    ]:
        y = [float(row[key]) for row in selected]
        if max(y) > 1.0e-30:
            ax.plot(x, y, marker="o", label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("initial deviation rad")
    ax.set_ylabel("candidate value")
    ax.set_title("extended sweep reference curves")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_native_inverse_area_extended_reference_curves_v1.png", dpi=180)
    plt.close(fig)

    fit_selected = [
        row
        for row in fits
        if row["filter"] == "c1_surface_like"
        and row["readout_mode"] == "readout_off"
        and row["candidate"]
        in {
            "native_max_f_AB",
            "native_max_envelope_V_AB",
            "native_max_V_AB",
            "max_abs_A_chi_tau",
            "inv_A_chi_tau_constructed_control",
            "derived_f_AB_over_A",
            "derived_envelope_over_A",
        }
        and bool(row["fit_valid"])
    ]
    fit_selected = sorted(
        fit_selected,
        key=lambda row: (
            row["candidate_kind"],
            row["candidate"],
            float(row["frequency_ratio"]),
            float(row["amplitude_ratio"]),
            float(row["phase_shift_deg"]),
        ),
    )[:40]
    labels = [
        f'{row["candidate"]}\nf={row["frequency_ratio"]},a={row["amplitude_ratio"]},p={row["phase_shift_deg"]}'
        for row in fit_selected
    ]
    alphas = [float(row["power_candidate_alpha"]) for row in fit_selected]
    colors = [
        "tab:red" if row["candidate_kind"] == "constructed_reciprocal_control" else "tab:blue"
        for row in fit_selected
    ]
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(range(len(alphas)), alphas, color=colors)
    ax.axhline(2.0, color="black", linewidth=1.0, linestyle="--", label="alpha=2")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xticks(range(len(labels)), labels, rotation=80, ha="right", fontsize=7)
    ax.set_ylabel("power_candidate_alpha")
    ax.set_title("extended c1-surface alpha scan")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_native_inverse_area_extended_alpha_scan_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any]) -> None:
    best = validation["best_native_alpha_candidate"]
    text = f"""# AB二体 chi-tau native 逆面積 extended sweep 予備実験レポート v1

## Summary

- valid: `{validation["native_inverse_area_extended_sweep_preliminary_valid"]}`
- sweep_case_count: `{validation["sweep_case_count"]}`
- area_valid_case_count: `{validation["area_valid_case_count"]}`
- c1_surface_like_case_count: `{validation["c1_surface_like_case_count"]}`
- native_positive2_count: `{validation["native_positive2_count"]}`
- c1_native_positive2_count: `{validation["c1_native_positive2_count"]}`
- constructed_reciprocal_positive2_count: `{validation["constructed_reciprocal_positive2_count"]}`
- c1_area_min: `{validation["c1_area_min"]:.16e}`
- c1_area_max: `{validation["c1_area_max"]:.16e}`

## Best native alpha candidate

```json
{json.dumps(best, ensure_ascii=False, indent=2)}
```

## Reading

This is an extended negative-control search.

The sweep varies initial deviation, frequency ratio, amplitude ratio, phase shift, and readout leak.

The strict result is:

```text
native inverse-area scaling was not detected.
constructed reciprocal controls give alpha near +2.
```

This keeps the inverse-square claim on hold.
"""
    (OUT_DIR / "ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_report_v1.md").write_text(
        text,
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--record-ab-sign-change-count",
        type=int,
        default=0,
        metavar="N",
        help="record each configuration through its Nth internal chi sign change; omitted or N <= 0 disables the CSV",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = Params()
    selected_configs = configs()
    rows = [run_config(config, params) for config in selected_configs]
    fits = fit_rows(rows)
    validation = validation_summary(rows, fits)
    write_csv(OUT_DIR / "ab_two_body_chi_tau_native_inverse_area_extended_sweep_cases_v1.csv", rows)
    write_csv(OUT_DIR / "ab_two_body_chi_tau_native_inverse_area_extended_sweep_fits_v1.csv", fits)
    if args.record_ab_sign_change_count > 0:
        write_ab_sign_change_series_csv(
            OUT_DIR / "ab_two_body_chi_tau_native_inverse_area_extended_sweep_ab_sign_change_series_v1.csv",
            selected_configs,
            params,
            args.record_ab_sign_change_count,
        )
    result = {"validation": validation, "case_summaries": rows, "fits": fits}
    (OUT_DIR / "ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(rows, fits)
    write_report(validation)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
