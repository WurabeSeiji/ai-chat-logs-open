#!/usr/bin/env python3
"""Blind control sweep around the previously found high-energy candidate band.

The grid phase, System B v5 dynamics, depth definition, stopping rule, and plot
style are inherited from the saved experiments.  The only observational guide
drawn on the sweep plot is the black dashed line at N_obs,high = 128.946.
The finite-order root is evaluated only after the grid sweep and is kept out of
the control figure.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import sys
import tempfile
import time
from decimal import Decimal
from pathlib import Path
from typing import Any

PLOT_CACHE = Path(tempfile.gettempdir()) / "wave_readout_matplotlib_cache"
PLOT_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(PLOT_CACHE / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
V5_SOURCE = (
    HERE.parent / "20260715" / "run_minimal_system_B_gray_direct_check_v5.py"
)
OUTPUT_DIR = HERE / "high_energy_control_delta1e-7_v1"

# This range preserves the phase of the original full-range delta_R=1e-7 grid.
# Its first point is the nearest old-grid point below R_obs,high, and its last
# point is the saved upper boundary of candidate band rank 2.
MIN_R = "0.68782280255614798"
MAX_R = "0.68853950255614798"
DELTA_R = "0.0000001"

N_OBS_HIGH = 128.946
OLD_PEAK_R = 0.68836390255614798
OLD_PEAK_N = 129.394062925467
R122_23 = math.cos(23.0 * math.pi / 122.0) ** 2

STEPS = 1024
PHI_MODE = "zero"
MIN_STEPS = 256
EARLY_STOP_PATIENCE = 20
Y_LIMITS = (7.8, 15.0)


def load_v5_module() -> Any:
    spec = importlib.util.spec_from_file_location("system_b_direct_v5", V5_SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load v5 source: {V5_SOURCE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def n_of_r(reflection_rate: float) -> float:
    return 4.0 * math.pi / (1.0 - reflection_rate) ** 2


def r_of_n(capacity: float) -> float:
    return 1.0 - math.sqrt(4.0 * math.pi / capacity)


def decimal_range(start_text: str, stop_text: str, step_text: str) -> list[str]:
    start = Decimal(start_text)
    stop = Decimal(stop_text)
    step = Decimal(step_text)
    values: list[str] = []
    current = start
    while current <= stop:
        values.append(format(current, "f"))
        current += step
    return values


def evaluate(v5: Any, r_text: str, kind: str, index: int | str) -> dict[str, object]:
    payload = v5.probe_r(
        float(r_text),
        r_text,
        STEPS,
        PHI_MODE,
        MIN_STEPS,
        EARLY_STOP_PATIENCE,
    )
    reflection_rate = float(payload["R"])
    return {
        "sample_kind": kind,
        "grid_index": index,
        "R": f"{reflection_rate:.18f}",
        "R_input_text": payload["R_input_text"],
        "N_of_R": f"{n_of_r(reflection_rate):.15f}",
        "offset_from_R_obs_high": (
            f"{reflection_rate - r_of_n(N_OBS_HIGH):+.18e}"
        ),
        "offset_from_old_peak": f"{reflection_rate - OLD_PEAK_R:+.18e}",
        "offset_from_R_122_23": f"{reflection_rate - R122_23:+.18e}",
        "best_step": int(payload["best_step"]),
        "best_prefix_gray_error_no_phase": (
            f"{float(payload['best_prefix_gray_error_no_phase']):.18e}"
        ),
        "best_prefix_gray_depth_no_phase": (
            f"{float(payload['best_prefix_gray_depth_no_phase']):.15f}"
        ),
        "best_condition_id": payload["best_condition_id"],
        "best_S_mean": f"{float(payload['best_S_mean']):+.18e}",
        "best_S_amp": f"{float(payload['best_S_amp']):.18e}",
        "best_S_drift": f"{float(payload['best_S_drift']):.18e}",
        "stopped_at_step": int(payload["stopped_at_step"]),
        "stop_reason": payload["stop_reason"],
        "is_v5_candidate": int(v5.stable_fixed_point_candidate(payload)),
        "T": f"{float(payload['T']):.18e}",
        "reflection_power": f"{float(payload['reflection_power']):.18e}",
        "phi_mode": payload["phi_mode"],
        "steps": int(payload["steps"]),
        "min_steps": int(payload["min_steps"]),
        "early_stop_patience": int(payload["early_stop_patience"]),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def local_maxima(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    maxima = [
        rows[index]
        for index in range(1, len(rows) - 1)
        if depths[index] > depths[index - 1]
        and depths[index] >= depths[index + 1]
    ]
    return sorted(
        maxima,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
        reverse=True,
    )


def make_plot(
    rows: list[dict[str, object]],
    obs_row: dict[str, object],
    output_png: Path,
    output_svg: Path,
) -> None:
    r_obs = r_of_n(N_OBS_HIGH)
    r_values = np.array([float(row["R"]) for row in rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    x_values = (r_values - r_obs) * 1.0e9
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )

    fig, axis = plt.subplots(figsize=(10.5, 6.4))
    axis.plot(
        x_values,
        depths,
        color="#2b6cb0",
        linewidth=1.6,
        marker="o",
        markersize=2.2,
        label=r"grid $\Delta R=10^{-7}$",
    )
    axis.axvline(
        0.0,
        color="#111111",
        linestyle="--",
        linewidth=1.6,
    )
    axis.scatter(
        [0.0],
        [float(obs_row["best_prefix_gray_depth_no_phase"])],
        color="#111111",
        edgecolor="white",
        linewidth=0.7,
        s=62,
        zorder=6,
        label=r"$R_{\rm obs,high}$ ($N=128.946$)",
    )
    axis.scatter(
        [(float(best["R"]) - r_obs) * 1.0e9],
        [float(best["best_prefix_gray_depth_no_phase"])],
        marker="x",
        color="#ff7f0e",
        linewidths=2.2,
        s=80,
        zorder=7,
        label="best grid point",
    )
    span = x_values[-1] - x_values[0]
    axis.set_xlim(x_values[0] - 0.015 * span, x_values[-1] + 0.015 * span)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlabel(r"$(R-R_{\rm obs,high})\times10^9$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(r"System B high-energy control sweep: $\Delta R=10^{-7}$")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_svg, bbox_inches="tight")
    plt.close(fig)


def make_markdown(
    path: Path,
    rows: list[dict[str, object]],
    obs_row: dict[str, object],
    root_row: dict[str, object],
    maxima: list[dict[str, object]],
    summary: dict[str, object],
) -> None:
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    lines = [
        "# 高エネルギー側ピーク対照掃引（delta R = 1e-7）",
        "",
        "## 条件",
        "",
        "System B v5 の計算核、深度定義、停止条件、旧全域掃引の格子位相を保持した。",
        "図中の理論根は伏せ、観測値 N_obs,high = 128.946 に対応する黒破線だけを置いた。",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| min R | {MIN_R} |",
        f"| max R | {MAX_R} |",
        f"| delta R | {DELTA_R} |",
        f"| grid points | {len(rows)} |",
        f"| R_obs,high | {summary['R_obs_high']} |",
        f"| old peak R | {OLD_PEAK_R:.18f} |",
        "",
        "## 盲検格子掃引の最深点",
        "",
        "| R | N(R) | depth | error | best step |",
        "|---:|---:|---:|---:|---:|",
        f"| {best['R']} | {best['N_of_R']} | "
        f"{best['best_prefix_gray_depth_no_phase']} | "
        f"{best['best_prefix_gray_error_no_phase']} | {best['best_step']} |",
        "",
        "## 掃引後の有限位数根照合",
        "",
        "| 点 | R | N(R) | depth | error | best step |",
        "|---|---:|---:|---:|---:|---:|",
        f"| observation | {obs_row['R']} | {obs_row['N_of_R']} | "
        f"{obs_row['best_prefix_gray_depth_no_phase']} | "
        f"{obs_row['best_prefix_gray_error_no_phase']} | {obs_row['best_step']} |",
        f"| exact R_122,23 | {root_row['R']} | {root_row['N_of_R']} | "
        f"{root_row['best_prefix_gray_depth_no_phase']} | "
        f"{root_row['best_prefix_gray_error_no_phase']} | {root_row['best_step']} |",
        "",
        f"格子最深点と R_122,23 の差は {best['offset_from_R_122_23']} である。",
        "",
        "## 局所極大（上位10点）",
        "",
        "| rank | R | N(R) | depth | best step |",
        "|---:|---:|---:|---:|---:|",
    ]
    for rank, row in enumerate(maxima[:10], 1):
        lines.append(
            f"| {rank} | {row['R']} | {row['N_of_R']} | "
            f"{row['best_prefix_gray_depth_no_phase']} | {row['best_step']} |"
        )
    lines.extend(
        [
            "",
            "## 図",
            "",
            "![high energy control delta R 1e-7](high_energy_control_delta1e-7_v1.png)",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    v5 = load_v5_module()
    grid = decimal_range(MIN_R, MAX_R, DELTA_R)
    rows: list[dict[str, object]] = []
    started = time.perf_counter()

    for index, r_text in enumerate(grid):
        rows.append(evaluate(v5, r_text, "grid", index))

    # These two direct evaluations are controls and do not alter the grid search.
    r_obs = r_of_n(N_OBS_HIGH)
    obs_row = evaluate(v5, f"{r_obs:.18f}", "observation", "")
    root_row = evaluate(v5, f"{R122_23:.18f}", "posthoc_exact_root", "")
    elapsed = time.perf_counter() - started

    maxima = local_maxima(rows)
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    summary: dict[str, object] = {
        "experiment_class": "blind_grid_control_then_posthoc_root_check",
        "min_R": MIN_R,
        "max_R": MAX_R,
        "delta_R": DELTA_R,
        "grid_point_count": len(rows),
        "R_obs_high": f"{r_obs:.18f}",
        "N_obs_high": N_OBS_HIGH,
        "old_peak_R": f"{OLD_PEAK_R:.18f}",
        "old_peak_N": OLD_PEAK_N,
        "grid_best": best,
        "local_maximum_count": len(maxima),
        "top_local_maxima": maxima[:10],
        "posthoc_exact_R_122_23": root_row,
        "elapsed_sec": elapsed,
        "plot_format": {
            "x_quantity": "(R-R_obs_high)*1e9",
            "y_quantity": "best_prefix_gray_depth_no_phase",
            "y_limits": Y_LIMITS,
            "reference_lines": ["R_obs_high only"],
            "figure_size": [10.5, 6.4],
            "dpi": 220,
        },
        "numerics": {
            "steps": STEPS,
            "phi_mode": PHI_MODE,
            "min_steps": MIN_STEPS,
            "early_stop_patience": EARLY_STOP_PATIENCE,
        },
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "high_energy_control_delta1e-7_all_v1.csv", rows)
    write_csv(
        OUTPUT_DIR / "high_energy_control_delta1e-7_direct_checks_v1.csv",
        [obs_row, root_row],
    )
    write_csv(
        OUTPUT_DIR / "high_energy_control_delta1e-7_local_maxima_v1.csv",
        maxima,
    )
    make_plot(
        rows,
        obs_row,
        OUTPUT_DIR / "high_energy_control_delta1e-7_v1.png",
        OUTPUT_DIR / "high_energy_control_delta1e-7_v1.svg",
    )
    (OUTPUT_DIR / "high_energy_control_delta1e-7_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "high_energy_control_delta1e-7_result_v1.md",
        rows,
        obs_row,
        root_row,
        maxima,
        summary,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
