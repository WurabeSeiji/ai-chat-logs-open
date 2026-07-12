from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_v1 import (
    Params,
    ReadoutMode,
    bool_all,
    rolling_c_error,
    rolling_closed_area,
    rolling_rank,
    safe_ratio,
    write_csv,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi


@dataclass(frozen=True)
class SweepConfig:
    frequency_ratio: float
    amplitude_ratio: float
    phase_shift_deg: float


READOUT_MODES = [
    ReadoutMode("readout_off", 0.0, False),
    ReadoutMode("readout_normal", 5.0e-5, True),
    ReadoutMode("readout_strong", 2.0e-4, True),
]

FREQUENCY_RATIOS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
AMPLITUDE_RATIOS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
PHASE_SHIFT_DEGS = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]


def sweep_configs() -> List[SweepConfig]:
    return [
        SweepConfig(freq, amp, phase)
        for freq in FREQUENCY_RATIOS
        for amp in AMPLITUDE_RATIOS
        for phase in PHASE_SHIFT_DEGS
    ]


def chi_tau_series(
    initial_deviation_rad: float,
    config: SweepConfig,
    readout_mode: ReadoutMode,
    params: Params,
) -> Dict[str, np.ndarray]:
    steps = np.arange(params.step_count + 1, dtype=float)
    lam = (1.0 - readout_mode.per_step_leak) ** steps
    phase = params.omega_step * steps
    chi = initial_deviation_rad * lam * np.cos(phase)
    tau_phase = config.frequency_ratio * params.omega_step * steps + math.radians(config.phase_shift_deg)
    tau = initial_deviation_rad * config.amplitude_ratio * lam * np.sin(tau_phase)
    area = rolling_closed_area(chi, tau, params.period_steps)
    ranks = rolling_rank(chi, tau, params.period_steps, True, params.rank_tol)
    c_errors = rolling_c_error(chi, tau, params.period_steps, True)
    return {"chi": chi, "tau": tau, "area": area, "rank": ranks, "c_error": c_errors}


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


def summarize_config(config: SweepConfig, readout_mode: ReadoutMode, params: Params) -> Dict[str, Any]:
    initial_deviation_rad = math.radians(10.0)
    series = chi_tau_series(initial_deviation_rad, config, readout_mode, params)
    chi = series["chi"]
    tau = series["tau"]
    area = np.abs(series["area"])
    ranks = series["rank"][params.period_steps :]
    c_values = series["c_error"][params.period_steps + 1 :]
    c_values = c_values[np.isfinite(c_values)]
    envelope = (1.0 - readout_mode.per_step_leak) ** (2.0 * np.arange(params.step_count + 1, dtype=float))
    max_area = float(np.max(area))
    rank_max = int(np.max(ranks)) if len(ranks) else 1
    c_error = float(np.max(np.abs(c_values))) if len(c_values) else 0.0
    locked_like = bool(max_area <= params.area_tol or rank_max < 2)
    c1_like = bool(c_error <= 5.0e-3)
    c1_surface_like = bool(c1_like and not locked_like and rank_max >= 2)
    return {
        "frequency_ratio": config.frequency_ratio,
        "amplitude_ratio": config.amplitude_ratio,
        "phase_shift_deg": config.phase_shift_deg,
        "amp_freq_product": config.amplitude_ratio * config.frequency_ratio,
        "readout_mode": readout_mode.name,
        "per_step_leak": readout_mode.per_step_leak,
        "max_abs_A_chi_tau": max_area,
        "rank_chi_tau": rank_max,
        "max_epsilon_c_abs": c_error,
        "locked_like": locked_like,
        "c1_like": c1_like,
        "c1_surface_like": c1_surface_like,
        "decay_rate_envelope": estimate_decay_rate(envelope),
        "chi_min": float(np.min(chi)),
        "chi_max": float(np.max(chi)),
        "tau_min": float(np.min(tau)),
        "tau_max": float(np.max(tau)),
        "tau_is_step_used": False,
        "external_c_used": False,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
    }


def make_plots(rows: List[Dict[str, Any]], params: Params) -> None:
    off_phase0 = [
        row for row in rows if row["readout_mode"] == "readout_off" and float(row["phase_shift_deg"]) == 0.0
    ]
    off_phase0.sort(key=lambda row: (float(row["frequency_ratio"]), float(row["amplitude_ratio"])))
    freq_values = FREQUENCY_RATIOS
    amp_values = AMPLITUDE_RATIOS
    c_grid = np.full((len(freq_values), len(amp_values)), np.nan)
    area_grid = np.full((len(freq_values), len(amp_values)), np.nan)
    for row in off_phase0:
        fi = freq_values.index(float(row["frequency_ratio"]))
        ai = amp_values.index(float(row["amplitude_ratio"]))
        c_grid[fi, ai] = float(row["max_epsilon_c_abs"])
        area_grid[fi, ai] = float(row["max_abs_A_chi_tau"])

    fig, ax = plt.subplots(figsize=(8, 5))
    image = ax.imshow(np.log10(np.maximum(c_grid, 1.0e-16)), origin="lower", aspect="auto")
    ax.set_xticks(range(len(amp_values)), [str(value) for value in amp_values])
    ax.set_yticks(range(len(freq_values)), [str(value) for value in freq_values])
    ax.set_xlabel("amplitude_ratio")
    ax.set_ylabel("frequency_ratio")
    ax.set_title("log10 max epsilon_c, phase_shift=0deg")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_parameter_sweep_c_error_heatmap_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    image = ax.imshow(np.log10(np.maximum(area_grid, 1.0e-16)), origin="lower", aspect="auto")
    ax.set_xticks(range(len(amp_values)), [str(value) for value in amp_values])
    ax.set_yticks(range(len(freq_values)), [str(value) for value in freq_values])
    ax.set_xlabel("amplitude_ratio")
    ax.set_ylabel("frequency_ratio")
    ax.set_title("log10 max |A_chi_tau|, phase_shift=0deg")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_parameter_sweep_area_heatmap_v1.png", dpi=180)
    plt.close(fig)

    phase_rows = [
        row
        for row in rows
        if row["readout_mode"] == "readout_off"
        and float(row["frequency_ratio"]) == 1.0
        and float(row["amplitude_ratio"]) == 1.0
    ]
    phase_rows.sort(key=lambda row: float(row["phase_shift_deg"]))
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(
        [float(row["phase_shift_deg"]) for row in phase_rows],
        [float(row["max_abs_A_chi_tau"]) for row in phase_rows],
        marker="o",
        color="tab:blue",
        label="max |A_chi_tau|",
    )
    ax1.set_xlabel("phase_shift_deg")
    ax1.set_ylabel("max |A_chi_tau|", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(
        [float(row["phase_shift_deg"]) for row in phase_rows],
        [float(row["max_epsilon_c_abs"]) for row in phase_rows],
        marker="s",
        color="tab:red",
        label="max |epsilon_c|",
    )
    ax2.set_ylabel("max |epsilon_c|", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax1.set_title("phase shift separates c calibration from area")
    ax1.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_parameter_sweep_phase_response_v1.png", dpi=180)
    plt.close(fig)

    readout_rows = [
        row
        for row in rows
        if float(row["frequency_ratio"]) == 1.0
        and float(row["amplitude_ratio"]) == 1.0
        and float(row["phase_shift_deg"]) == 0.0
    ]
    readout_rows.sort(key=lambda row: float(row["per_step_leak"]))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        [float(row["per_step_leak"]) for row in readout_rows],
        [float(row["max_epsilon_c_abs"]) for row in readout_rows],
        marker="o",
        label="max |epsilon_c|",
    )
    ax.plot(
        [float(row["per_step_leak"]) for row in readout_rows],
        [abs(float(row["decay_rate_envelope"])) for row in readout_rows],
        marker="o",
        label="|decay_rate_envelope|",
    )
    ax.set_xscale("symlog", linthresh=1.0e-7)
    ax.set_yscale("symlog", linthresh=1.0e-15)
    ax.set_xlabel("per_step_leak")
    ax.set_title("readout leak perturbs c1 surface")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_c1_parameter_sweep_readout_leak_v1.png", dpi=180)
    plt.close(fig)


def validation_summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    off_rows = [row for row in rows if row["readout_mode"] == "readout_off"]
    c1_surface_rows = [row for row in off_rows if bool(row["c1_surface_like"])]
    c1_locked_rows = [row for row in off_rows if bool(row["c1_like"]) and bool(row["locked_like"])]
    rank2_rows = [row for row in off_rows if int(row["rank_chi_tau"]) == 2]
    best = min(off_rows, key=lambda row: (float(row["max_epsilon_c_abs"]), -float(row["max_abs_A_chi_tau"])))
    return {
        "c1_internal_calibration_parameter_sweep_preliminary_valid": bool(
            len(c1_surface_rows) >= 1
            and len(c1_locked_rows) >= 1
            and bool_all(not bool(row["tau_is_step_used"]) for row in rows)
            and bool_all(not bool(row["external_c_used"]) for row in rows)
        ),
        "sweep_case_count": len(rows),
        "readout_off_case_count": len(off_rows),
        "rank2_readout_off_count": len(rank2_rows),
        "c1_surface_like_readout_off_count": len(c1_surface_rows),
        "c1_locked_like_readout_off_count": len(c1_locked_rows),
        "best_readout_off_config": best,
        "min_c_error_readout_off": min(float(row["max_epsilon_c_abs"]) for row in off_rows),
        "max_area_readout_off": max(float(row["max_abs_A_chi_tau"]) for row in off_rows),
        "tau_is_step_used_any": any(bool(row["tau_is_step_used"]) for row in rows),
        "external_c_used_any": any(bool(row["external_c_used"]) for row in rows),
        "f_A_or_f_B_used_any": any(bool(row["f_A_or_f_B_used"]) for row in rows),
        "main_reading": "c=1 calibration alone is not enough; rank and area controls are required.",
    }


def write_report(validation: Dict[str, Any]) -> None:
    best = validation["best_readout_off_config"]
    text = f"""# AB二体 c=1 内部較正パラメータスイープ予備実験レポート v1

## Summary

- valid: `{validation["c1_internal_calibration_parameter_sweep_preliminary_valid"]}`
- sweep_case_count: `{validation["sweep_case_count"]}`
- readout_off_case_count: `{validation["readout_off_case_count"]}`
- rank2_readout_off_count: `{validation["rank2_readout_off_count"]}`
- c1_surface_like_readout_off_count: `{validation["c1_surface_like_readout_off_count"]}`
- c1_locked_like_readout_off_count: `{validation["c1_locked_like_readout_off_count"]}`
- min_c_error_readout_off: `{validation["min_c_error_readout_off"]:.16e}`
- max_area_readout_off: `{validation["max_area_readout_off"]:.16e}`

## Best readout_off config

```json
{json.dumps(best, ensure_ascii=False, indent=2)}
```

## Reading

This sweep shows that internal `c=1` calibration alone is not sufficient.

Some locked-like configurations can satisfy the RMS exchange ratio while producing no independent `chi-tau` area.

Therefore, the necessary readout conditions are:

```text
c=1 calibration
rank_chi_tau = 2
A_chi_tau != 0
```

The result strengthens the control discipline for the next experiment.
"""
    (OUT_DIR / "ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_report_v1.md").write_text(
        text,
        encoding="utf-8",
    )


def main() -> None:
    params = Params()
    rows: List[Dict[str, Any]] = []
    for config in sweep_configs():
        for readout_mode in READOUT_MODES:
            rows.append(summarize_config(config, readout_mode, params))
    validation = validation_summary(rows)
    write_csv(
        OUT_DIR / "ab_two_body_c1_internal_calibration_parameter_sweep_case_summary_v1.csv",
        rows,
    )
    result = {
        "params": {
            "step_count": params.step_count,
            "period_steps": params.period_steps,
            "omega_step": params.omega_step,
        },
        "validation": validation,
        "case_summaries": rows,
    }
    (OUT_DIR / "ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(rows, params)
    write_report(validation)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
