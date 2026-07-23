#!/usr/bin/env python3
"""Multiprecision root-centered sweeps for the two observation-linked roots."""

from __future__ import annotations

import csv
import json
import math
import os
import tempfile
import time
from pathlib import Path
from typing import Any

PLOT_CACHE = Path(tempfile.gettempdir()) / "wave_readout_matplotlib_cache"
PLOT_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(PLOT_CACHE / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "two_physical_roots_multiprecision_v1"

ROOTS = [
    {
        "stem": "R124_23",
        "latex": r"R_{124,23}",
        "n": 124,
        "m": 23,
        "physical_label": r"low-energy readout",
    },
    {
        "stem": "R620_117",
        "latex": r"R_{620,117}",
        "n": 620,
        "m": 117,
        "physical_label": r"high-energy readout",
    },
]

CONFIGS = [
    {"name": "50digit_delta1e-14", "dps": 50, "delta": "1e-14", "exponent": 14},
    {"name": "80digit_delta1e-16", "dps": 80, "delta": "1e-16", "exponent": 16},
]

HALF_STEPS = 1242
POINT_COUNT = 2 * HALF_STEPS + 1
X_LIMITS = (-1250.0, 1250.0)
Y_LIMITS = (7.8, 85.0)
AMP_TARGET_TEXT = "0.02"


def n_of_r(reflection_rate: mp.mpf) -> mp.mpf:
    return 4 * mp.pi / (1 - reflection_rate) ** 2


def resonant_prefix_metrics(
    reflection_rate: mp.mpf,
    n: int,
) -> dict[str, mp.mpf | int]:
    """Evaluate the v5 gray metric at the verified resonant prefix 2n."""
    phase = mp.pi + 2 * mp.asin(mp.sqrt(reflection_rate))
    phase_cos = mp.cos(phase)
    amplitude = mp.mpf(AMP_TARGET_TEXT)
    sample_count = 2 * n

    cosine_previous = mp.mpf(1)
    cosine_current = phase_cos
    total = mp.mpf(0)
    first_half_total = mp.mpf(0)
    s_min = mp.inf
    s_max = -mp.inf

    for sample_index in range(sample_count):
        if sample_index == 0:
            cosine_value = cosine_previous
        elif sample_index == 1:
            cosine_value = cosine_current
        else:
            cosine_previous, cosine_current = (
                cosine_current,
                2 * phase_cos * cosine_current - cosine_previous,
            )
            cosine_value = cosine_current

        s_value = amplitude * cosine_value
        total += s_value
        s_min = min(s_min, s_value)
        s_max = max(s_max, s_value)
        if sample_index == n - 1:
            first_half_total = total

    s_mean = total / sample_count
    s_amp = (s_max - s_min) / 2
    first_mean = first_half_total / n
    second_mean = (total - first_half_total) / n
    s_drift = abs(second_mean - first_mean)
    gray_error = abs(s_mean) + abs(s_amp - amplitude) + s_drift
    gray_depth = -mp.log10(gray_error)
    return {
        "best_step": sample_count - 1,
        "gray_error": gray_error,
        "gray_depth": gray_depth,
        "S_mean": s_mean,
        "S_amp": s_amp,
        "S_drift": s_drift,
    }


def value_text(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def run_sweep(root: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    dps = int(config["dps"])
    n = int(root["n"])
    m = int(root["m"])
    with mp.workdps(dps):
        exact_root = mp.cos(mp.mpf(m) * mp.pi / mp.mpf(n)) ** 2
        delta = mp.mpf(str(config["delta"]))
        rows: list[dict[str, object]] = []
        for grid_index, offset in enumerate(range(-HALF_STEPS, HALF_STEPS + 1)):
            reflection_rate = exact_root + mp.mpf(offset) * delta
            metrics = resonant_prefix_metrics(reflection_rate, n)
            rows.append(
                {
                    "grid_index": grid_index,
                    "grid_offset": offset,
                    "R": value_text(reflection_rate, dps + 8),
                    "offset_from_root": value_text(mp.mpf(offset) * delta, dps + 4),
                    "N_of_R": value_text(n_of_r(reflection_rate), dps),
                    "best_step": int(metrics["best_step"]),
                    "best_prefix_gray_error_no_phase": value_text(
                        mp.mpf(metrics["gray_error"]), dps + 8
                    ),
                    "best_prefix_gray_depth_no_phase": value_text(
                        mp.mpf(metrics["gray_depth"]), 24
                    ),
                    "best_S_mean": value_text(mp.mpf(metrics["S_mean"]), dps + 8),
                    "best_S_amp": value_text(mp.mpf(metrics["S_amp"]), dps + 8),
                    "best_S_drift": value_text(mp.mpf(metrics["S_drift"]), dps + 8),
                    "working_decimal_digits": dps,
                    "delta_R": str(config["delta"]),
                }
            )

        best = max(
            rows,
            key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
        )
        center = rows[HALF_STEPS]
        return {
            "root_label": root["stem"],
            "root_latex": root["latex"],
            "n": n,
            "m": m,
            "physical_label": root["physical_label"],
            "working_decimal_digits": dps,
            "delta_R": str(config["delta"]),
            "root_high_precision_at_working_dps": value_text(exact_root, dps),
            "point_count": len(rows),
            "half_steps": HALF_STEPS,
            "grid_best": best,
            "center": center,
            "center_is_best": int(best["grid_index"]) == HALF_STEPS,
            "expected_best_step": 2 * n - 1,
            "rows": rows,
        }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_one(axis: plt.Axes, result: dict[str, Any], exponent: int) -> None:
    rows = result["rows"]
    x_values = np.array([float(row["grid_offset"]) for row in rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    center = result["center"]
    best = result["grid_best"]
    root_latex = str(result["root_latex"])
    dps = int(result["working_decimal_digits"])

    axis.plot(
        x_values,
        depths,
        color="#2b6cb0",
        linewidth=1.6,
        marker="o",
        markersize=2.2,
        label=rf"{dps}-digit grid $\Delta R=10^{{-{exponent}}}$",
    )
    axis.axvline(0.0, color="#d627a5", linestyle="--", linewidth=1.6)
    axis.scatter(
        [0.0],
        [float(center["best_prefix_gray_depth_no_phase"])],
        color="#d627a5",
        edgecolor="white",
        linewidth=0.7,
        s=72,
        zorder=6,
        label=rf"exact ${root_latex}$",
    )
    axis.scatter(
        [float(best["grid_offset"])],
        [float(best["best_prefix_gray_depth_no_phase"])],
        marker="x",
        color="#ff7f0e",
        linewidths=2.2,
        s=85,
        zorder=7,
        label="best grid point",
    )
    axis.annotate(
        f"center depth = {float(center['best_prefix_gray_depth_no_phase']):.6f}",
        xy=(0.0, float(center["best_prefix_gray_depth_no_phase"])),
        xytext=(100.0, min(81.0, float(center["best_prefix_gray_depth_no_phase"]) + 4.0)),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=10,
    )
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlabel(rf"$(R-{root_latex})\times10^{{{exponent}}}$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(
        rf"System B multiprecision: ${root_latex}$, "
        rf"{dps} digits, $\Delta R=10^{{-{exponent}}}$"
    )
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")


def save_individual(result: dict[str, Any], config: dict[str, Any]) -> None:
    exponent = int(config["exponent"])
    fig, axis = plt.subplots(figsize=(10.5, 6.4))
    plot_one(axis, result, exponent)
    fig.tight_layout()
    stem = OUTPUT_DIR / f"{result['root_label']}_{config['name']}_v1"
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def make_markdown(path: Path, compact_results: list[dict[str, Any]]) -> None:
    lines = [
        "# 二つの観測対応根：多倍長精度掃引",
        "",
        "対象は低エネルギー側 R124,23 と高エネルギー側 R620,117 である。",
        "50桁・delta R=1e-14、および80桁・delta R=1e-16で、各厳密根を中央格子点に置いた。",
        "",
        "| root | precision | delta R | center depth | center error | best step | center is best |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for result in compact_results:
        center = result["center"]
        lines.append(
            f"| {result['root_label']} | {result['working_decimal_digits']} | "
            f"{result['delta_R']} | {center['best_prefix_gray_depth_no_phase']} | "
            f"{center['best_prefix_gray_error_no_phase']} | "
            f"{center['best_step']} | {result['center_is_best']} |"
        )
    for config in CONFIGS:
        lines.extend(
            [
                "",
                f"## {config['dps']}桁・delta R={config['delta']}",
                "",
                f"![{config['name']} comparison]({config['name']}_comparison_v1.png)",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def compact_result(result: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in result.items() if key != "rows"}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_compact: list[dict[str, Any]] = []
    started = time.perf_counter()

    for config in CONFIGS:
        config_results: list[dict[str, Any]] = []
        for root in ROOTS:
            result = run_sweep(root, config)
            config_results.append(result)
            compact = compact_result(result)
            all_compact.append(compact)
            write_csv(
                OUTPUT_DIR / f"{root['stem']}_{config['name']}_all_v1.csv",
                result["rows"],
            )
            save_individual(result, config)
            print(
                json.dumps(
                    {
                        "root": root["stem"],
                        "precision": config["dps"],
                        "delta_R": config["delta"],
                        "center_depth": result["center"]["best_prefix_gray_depth_no_phase"],
                    },
                    ensure_ascii=False,
                ),
                flush=True,
            )

        fig, axes = plt.subplots(1, 2, figsize=(17.0, 6.3), sharex=True, sharey=True)
        for axis, result in zip(axes, config_results):
            plot_one(axis, result, int(config["exponent"]))
        fig.suptitle(
            f"System B multiprecision root-centered sweeps: {config['dps']} digits",
            fontsize=17,
            fontweight="bold",
        )
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
        fig.savefig(
            OUTPUT_DIR / f"{config['name']}_comparison_v1.png",
            dpi=220,
            bbox_inches="tight",
        )
        fig.savefig(
            OUTPUT_DIR / f"{config['name']}_comparison_v1.svg",
            bbox_inches="tight",
        )
        plt.close(fig)

    summary = {
        "experiment": "multiprecision precision-floor test",
        "roots": [root["stem"] for root in ROOTS],
        "configs": CONFIGS,
        "half_steps": HALF_STEPS,
        "point_count_per_sweep": POINT_COUNT,
        "x_limits_in_grid_units": X_LIMITS,
        "y_limits": Y_LIMITS,
        "method": (
            "exact eigenphase recurrence S_j=0.02*cos(j*(pi+2*asin(sqrt(R)))) "
            "evaluated at the endpoint-validated resonant prefix 2n"
        ),
        "results": all_compact,
        "elapsed_sec": time.perf_counter() - started,
    }
    (OUTPUT_DIR / "two_physical_roots_multiprecision_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "two_physical_roots_multiprecision_result_v1.md",
        all_compact,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
