#!/usr/bin/env python3
"""Compare root-centered System B sweeps at delta_R=1e-10 and 1e-12."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import sys
import tempfile
import time
from decimal import Decimal, localcontext
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
OUTPUT_DIR = HERE / "R122_23_root_centered_resolution_comparison_v1"

# 80-digit evaluation of cos^2(23*pi/122).  The v5 calculation itself remains
# IEEE-754 double precision, as in all preceding System B experiments.
ROOT_TEXT = (
    "0.6883639468175925497192711956024063181992978311211591523570398745029067266126481208"
)
DELTAS = {
    "1e-10": "0.0000000001",
    "1e-12": "0.000000000001",
}
HALF_STEPS = 1242
POINT_COUNT = 2 * HALF_STEPS + 1

STEPS = 1024
PHI_MODE = "zero"
MIN_STEPS = 256
EARLY_STOP_PATIENCE = 20

Y_LIMITS = (7.8, 15.0)
X_LIMITS = (-1250.0, 1250.0)


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


def centered_grid(delta_text: str) -> list[str]:
    with localcontext() as context:
        context.prec = 100
        root = Decimal(ROOT_TEXT)
        delta = Decimal(delta_text)
        return [
            format(root + Decimal(offset) * delta, "f")
            for offset in range(-HALF_STEPS, HALF_STEPS + 1)
        ]


def evaluate(
    v5: Any,
    r_text: str,
    kind: str,
    index: int | str,
) -> dict[str, object]:
    payload = v5.probe_r(
        float(r_text),
        r_text,
        STEPS,
        PHI_MODE,
        MIN_STEPS,
        EARLY_STOP_PATIENCE,
    )
    reflection_rate = float(payload["R"])
    root_float = float(ROOT_TEXT)
    return {
        "sample_kind": kind,
        "grid_index": index,
        "R": f"{reflection_rate:.18f}",
        "R_input_text": payload["R_input_text"],
        "N_of_R": f"{n_of_r(reflection_rate):.15f}",
        "offset_from_R_122_23_float": (
            f"{reflection_rate - root_float:+.18e}"
        ),
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


def plot_one(
    axis: plt.Axes,
    rows: list[dict[str, object]],
    root_row: dict[str, object],
    exponent: int,
) -> None:
    root = float(ROOT_TEXT)
    scale = 10.0**exponent
    r_values = np.array([float(row["R"]) for row in rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    x_values = (r_values - root) * scale
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )

    axis.plot(
        x_values,
        depths,
        color="#2b6cb0",
        linewidth=1.6,
        marker="o",
        markersize=2.2,
        label=rf"grid $\Delta R=10^{{-{exponent}}}$",
    )
    axis.axvline(
        0.0,
        color="#d627a5",
        linestyle="--",
        linewidth=1.6,
    )
    axis.scatter(
        [0.0],
        [float(root_row["best_prefix_gray_depth_no_phase"])],
        color="#d627a5",
        edgecolor="white",
        linewidth=0.7,
        s=72,
        zorder=6,
        label=r"exact $R_{122,23}$",
    )
    axis.scatter(
        [(float(best["R"]) - root) * scale],
        [float(best["best_prefix_gray_depth_no_phase"])],
        marker="x",
        color="#ff7f0e",
        linewidths=2.2,
        s=85,
        zorder=7,
        label="best grid point",
    )
    axis.annotate(
        "exact root = center grid point",
        xy=(0.0, float(root_row["best_prefix_gray_depth_no_phase"])),
        xytext=(120.0, 14.35),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=10,
    )
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlabel(rf"$(R-R_{{122,23}})\times10^{{{exponent}}}$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(rf"System B root-centered sweep: $\Delta R=10^{{-{exponent}}}$")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")


def mirrored_stats(rows: list[dict[str, object]]) -> dict[str, float]:
    differences = [
        abs(
            float(rows[HALF_STEPS - offset]["best_prefix_gray_depth_no_phase"])
            - float(rows[HALF_STEPS + offset]["best_prefix_gray_depth_no_phase"])
        )
        for offset in range(1, HALF_STEPS + 1)
    ]
    return {
        "max_mirrored_depth_difference": max(differences),
        "mean_mirrored_depth_difference": float(np.mean(differences)),
    }


def make_markdown(
    path: Path,
    results: dict[str, dict[str, object]],
) -> None:
    lines = [
        "# R122,23 厳密根中心：1e-10 / 1e-12 解像度比較",
        "",
        "両掃引は R122,23 を中央格子点とし、点数、左右のステップ数、縦軸範囲、描画条件を共通化した。",
        "",
        f"高精度根: `{ROOT_TEXT}`",
        "",
        "| delta R | points | min R | max R | best grid index | best R input | depth | error | best step |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label in ("1e-10", "1e-12"):
        item = results[label]
        best = item["grid_best"]
        lines.append(
            f"| {label} | {item['point_count']} | {item['min_R']} | {item['max_R']} | "
            f"{best['grid_index']} | {best['R_input_text']} | "
            f"{best['best_prefix_gray_depth_no_phase']} | "
            f"{best['best_prefix_gray_error_no_phase']} | {best['best_step']} |"
        )
    lines.extend(
        [
            "",
            "## delta R = 1e-10",
            "",
            "![R122 23 root centered 1e-10](R122_23_root_centered_delta1e-10_v1.png)",
            "",
            "## delta R = 1e-12",
            "",
            "![R122 23 root centered 1e-12](R122_23_root_centered_delta1e-12_v1.png)",
            "",
            "## 同一軸比較",
            "",
            "![R122 23 root centered resolution comparison](R122_23_root_centered_resolution_comparison_v1.png)",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    v5 = load_v5_module()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, object]] = {}
    rows_by_label: dict[str, list[dict[str, object]]] = {}
    roots_by_label: dict[str, dict[str, object]] = {}
    started = time.perf_counter()

    for label, delta_text in DELTAS.items():
        grid = centered_grid(delta_text)
        rows = [
            evaluate(v5, r_text, "grid", index)
            for index, r_text in enumerate(grid)
        ]
        root_row = evaluate(v5, ROOT_TEXT, "exact_reference", "")
        best = max(
            rows,
            key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
        )
        stats: dict[str, object] = {
            "delta_R": delta_text,
            "half_steps": HALF_STEPS,
            "point_count": len(rows),
            "min_R": grid[0],
            "max_R": grid[-1],
            "center_grid_index": HALF_STEPS,
            "center_input_text": rows[HALF_STEPS]["R_input_text"],
            "center_is_exact_root_text": (
                rows[HALF_STEPS]["R_input_text"] == ROOT_TEXT
            ),
            "grid_best": best,
            "exact_root": root_row,
            **mirrored_stats(rows),
        }
        results[label] = stats
        rows_by_label[label] = rows
        roots_by_label[label] = root_row
        write_csv(
            OUTPUT_DIR / f"R122_23_root_centered_delta{label}_all_v1.csv",
            rows,
        )
        write_csv(
            OUTPUT_DIR / f"R122_23_root_centered_delta{label}_exact_root_v1.csv",
            [root_row],
        )

        exponent = 10 if label == "1e-10" else 12
        fig, axis = plt.subplots(figsize=(10.5, 6.4))
        plot_one(axis, rows, root_row, exponent)
        fig.tight_layout()
        stem = OUTPUT_DIR / f"R122_23_root_centered_delta{label}_v1"
        fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(17.0, 6.3), sharex=True, sharey=True)
    plot_one(axes[0], rows_by_label["1e-10"], roots_by_label["1e-10"], 10)
    plot_one(axes[1], rows_by_label["1e-12"], roots_by_label["1e-12"], 12)
    fig.suptitle(
        r"System B: $R_{122,23}$ root-centered resolution comparison",
        fontsize=17,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(
        OUTPUT_DIR / "R122_23_root_centered_resolution_comparison_v1.png",
        dpi=220,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTPUT_DIR / "R122_23_root_centered_resolution_comparison_v1.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    summary = {
        "root_definition": "cos^2(23*pi/122)",
        "root_high_precision": ROOT_TEXT,
        "root_float": f"{float(ROOT_TEXT):.18f}",
        "same_plot_conditions": {
            "half_steps": HALF_STEPS,
            "point_count": POINT_COUNT,
            "x_limits_in_grid_units": X_LIMITS,
            "y_quantity": "best_prefix_gray_depth_no_phase",
            "y_limits": Y_LIMITS,
            "figure_size": [10.5, 6.4],
            "dpi": 220,
        },
        "numerics": {
            "steps": STEPS,
            "phi_mode": PHI_MODE,
            "min_steps": MIN_STEPS,
            "early_stop_patience": EARLY_STOP_PATIENCE,
        },
        "results": results,
        "elapsed_sec": time.perf_counter() - started,
    }
    (OUTPUT_DIR / "R122_23_root_centered_resolution_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "R122_23_root_centered_resolution_comparison_v1.md",
        results,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
