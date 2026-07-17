#!/usr/bin/env python3
"""Run delta_R=1e-8 and plot it in exactly the same form as delta_R=1e-10."""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import sys
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
import numpy as np


HERE = Path(__file__).resolve().parent
PRECISION_RUNNER = HERE / "run_three_point_precision_sweep_v1.py"
PRECISION_OUTPUT = HERE / "three_point_precision_sweep_v1"
OUTPUT_DIR = HERE / "three_point_matched_resolution_comparison_v1"

MIN_R = "0.697177779231003050"
MAX_R = "0.697178027556659305"
DELTA_1E8 = "0.00000001"
DELTA_1E10 = "0.0000000001"

X_LIMITS = (-100.0, 150.0)
Y_LIMITS = (7.8, 15.0)


def load_precision_runner() -> Any:
    spec = importlib.util.spec_from_file_location(
        "three_point_precision_runner", PRECISION_RUNNER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load: {PRECISION_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_1e8(
    runner: Any,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    v5 = runner.load_v5_module()
    refs = runner.references()
    r_grid = runner.decimal_range(MIN_R, MAX_R, DELTA_1E8)
    grid_rows: list[dict[str, object]] = []
    exact_rows: list[dict[str, object]] = []
    total_loops = 0
    started = time.perf_counter()

    for index, r_text in enumerate(r_grid):
        row, loops = runner.evaluate(v5, r_text, "grid", "", index, refs)
        grid_rows.append(row)
        total_loops += loops

    for ref in refs:
        row, loops = runner.evaluate(
            v5,
            f"{float(ref['R']):.18f}",
            "exact_reference",
            str(ref["label"]),
            "",
            refs,
        )
        exact_rows.append(row)
        total_loops += loops

    elapsed = time.perf_counter() - started
    best = max(
        grid_rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    stats = {
        "min_R": MIN_R,
        "max_R": MAX_R,
        "delta_R": DELTA_1E8,
        "grid_point_count": len(grid_rows),
        "exact_reference_count": len(exact_rows),
        "elapsed_sec": elapsed,
        "total_loop_count": total_loops,
        "grid_best": best,
    }
    return grid_rows, exact_rows, stats


def plot_one(
    axis: plt.Axes,
    grid_rows: list[dict[str, object]],
    exact_rows: list[dict[str, object]],
    refs: list[dict[str, float | str]],
    delta_label: str,
) -> None:
    r_obs = float(refs[0]["R"])
    r_values = np.array([float(row["R"]) for row in grid_rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in grid_rows]
    )
    colors = {"R_obs": "#111111", "R_g": "#2ca02c", "R_124_23": "#d627a5"}
    labels = {
        "R_obs": r"$R_{\rm obs}$",
        "R_g": r"$R_g$",
        "R_124_23": r"$R_{124,23}$",
    }

    axis.plot(
        (r_values - r_obs) * 1.0e9,
        depths,
        color="#2b6cb0",
        linewidth=1.6,
        marker="o",
        markersize=2.2,
        label=rf"grid $\Delta R={delta_label}$",
    )
    for ref, exact in zip(refs, exact_rows):
        offset = (float(ref["R"]) - r_obs) * 1.0e9
        axis.axvline(
            offset,
            color=colors[str(ref["label"])],
            linestyle="--",
            linewidth=1.6,
        )
        axis.scatter(
            [offset],
            [float(exact["best_prefix_gray_depth_no_phase"])],
            color=colors[str(ref["label"])],
            edgecolor="white",
            linewidth=0.7,
            s=62,
            zorder=6,
            label=labels[str(ref["label"])],
        )

    best = max(
        grid_rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
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
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlabel(r"$(R-R_{\rm obs})\times10^9$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(rf"System B sweep: $\Delta R={delta_label}$")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")


def save_individual(
    grid_rows: list[dict[str, object]],
    exact_rows: list[dict[str, object]],
    refs: list[dict[str, float | str]],
    delta_label: str,
    stem: str,
) -> None:
    fig, axis = plt.subplots(figsize=(10.5, 6.4))
    plot_one(axis, grid_rows, exact_rows, refs, delta_label)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def make_markdown(
    path: Path,
    rows_1e8: list[dict[str, object]],
    rows_1e10: list[dict[str, object]],
    exact_rows: list[dict[str, object]],
) -> None:
    best_1e8 = max(
        rows_1e8,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    best_1e10 = max(
        rows_1e10,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    lines = [
        "# 低エネルギー3点：同一形式による解像度比較",
        "",
        "両図は横軸、縦軸、軸範囲、基準線、凡例を完全に共通化している。",
        "",
        "| 刻み | 点数 | grid best R | rootとの差 | depth | error |",
        "|---|---:|---:|---:|---:|---:|",
        f"| 1e-8 | {len(rows_1e8)} | {best_1e8['R']} | "
        f"{best_1e8['offset_from_R_124_23']} | "
        f"{best_1e8['best_prefix_gray_depth_no_phase']} | "
        f"{best_1e8['best_prefix_gray_error_no_phase']} |",
        f"| 1e-10 | {len(rows_1e10)} | {best_1e10['R']} | "
        f"{best_1e10['offset_from_R_124_23']} | "
        f"{best_1e10['best_prefix_gray_depth_no_phase']} | "
        f"{best_1e10['best_prefix_gray_error_no_phase']} |",
        f"| exact root | 1 | {exact_rows[2]['R']} | 0 | "
        f"{exact_rows[2]['best_prefix_gray_depth_no_phase']} | "
        f"{exact_rows[2]['best_prefix_gray_error_no_phase']} |",
        "",
        "## 1e-8",
        "",
        "![delta R 1e-8](three_point_sweep_delta1e-8_matched_v1.png)",
        "",
        "## 1e-10",
        "",
        "![delta R 1e-10](three_point_sweep_delta1e-10_matched_v1.png)",
        "",
        "## 左右比較",
        "",
        "![matched resolution comparison](three_point_matched_resolution_comparison_v1.png)",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    runner = load_precision_runner()
    refs = runner.references()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows_1e8, exact_1e8, stats_1e8 = run_1e8(runner)
    rows_1e10 = read_csv(
        PRECISION_OUTPUT / "three_point_precision_sweep_all_v1.csv"
    )
    exact_1e10 = read_csv(
        PRECISION_OUTPUT / "three_point_precision_exact_references_v1.csv"
    )

    write_csv(OUTPUT_DIR / "three_point_sweep_delta1e-8_all_v1.csv", rows_1e8)
    write_csv(
        OUTPUT_DIR / "three_point_sweep_delta1e-8_exact_references_v1.csv",
        exact_1e8,
    )

    save_individual(
        rows_1e8,
        exact_1e8,
        refs,
        "10^{-8}",
        "three_point_sweep_delta1e-8_matched_v1",
    )
    save_individual(
        rows_1e10,
        exact_1e10,
        refs,
        "10^{-10}",
        "three_point_sweep_delta1e-10_matched_v1",
    )

    fig, axes = plt.subplots(1, 2, figsize=(17.0, 6.3), sharex=True, sharey=True)
    plot_one(axes[0], rows_1e8, exact_1e8, refs, "10^{-8}")
    plot_one(axes[1], rows_1e10, exact_1e10, refs, "10^{-10}")
    fig.suptitle(
        "System B matched-format resolution comparison",
        fontsize=17,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(
        OUTPUT_DIR / "three_point_matched_resolution_comparison_v1.png",
        dpi=220,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTPUT_DIR / "three_point_matched_resolution_comparison_v1.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    summary = {
        "format": {
            "x_quantity": "(R-R_obs)*1e9",
            "x_limits": X_LIMITS,
            "y_quantity": "best_prefix_gray_depth_no_phase",
            "y_limits": Y_LIMITS,
            "reference_lines": ["R_obs", "R_g", "R_124_23"],
        },
        "delta_1e8": stats_1e8,
        "delta_1e10": {
            "delta_R": DELTA_1E10,
            "grid_point_count": len(rows_1e10),
            "grid_best": max(
                rows_1e10,
                key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
            ),
        },
    }
    (OUTPUT_DIR / "three_point_matched_resolution_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "three_point_matched_resolution_comparison_v1.md",
        rows_1e8,
        rows_1e10,
        exact_1e10,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
