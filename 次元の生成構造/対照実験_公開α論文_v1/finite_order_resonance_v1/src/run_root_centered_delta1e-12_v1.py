#!/usr/bin/env python3
"""Run a delta_R=1e-12 sweep centered exactly on R_(124,23)."""

from __future__ import annotations

import csv
import importlib.util
import json
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
PRECISION_RUNNER = HERE / "run_three_point_precision_sweep_v1.py"
OUTPUT_DIR = HERE / "root_centered_delta1e-12_v1"

ROOT_TEXT = "0.697177927556659305"
DELTA_TEXT = "0.000000000001"
HALF_STEPS = 1242
POINT_COUNT = 2 * HALF_STEPS + 1

Y_LIMITS = (7.8, 15.0)
X_LIMITS = (-1250.0, 1250.0)


def load_precision_runner() -> Any:
    spec = importlib.util.spec_from_file_location(
        "root_centered_precision_runner", PRECISION_RUNNER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load: {PRECISION_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def centered_grid() -> list[str]:
    root = Decimal(ROOT_TEXT)
    delta = Decimal(DELTA_TEXT)
    return [
        format(root + Decimal(offset) * delta, "f")
        for offset in range(-HALF_STEPS, HALF_STEPS + 1)
    ]


def make_plot(
    rows: list[dict[str, object]],
    root_row: dict[str, object],
    output_png: Path,
    output_svg: Path,
) -> None:
    root = float(ROOT_TEXT)
    r_values = np.array([float(row["R"]) for row in rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    x_values = (r_values - root) * 1.0e12
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
        label=r"grid $\Delta R=10^{-12}$",
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
        label=r"exact $R_{124,23}$",
    )
    axis.scatter(
        [(float(best["R"]) - root) * 1.0e12],
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
    axis.set_xlabel(r"$(R-R_{124,23})\times10^{12}$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(r"System B root-centered sweep: $\Delta R=10^{-12}$")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_svg, bbox_inches="tight")
    plt.close(fig)


def make_markdown(
    path: Path,
    rows: list[dict[str, object]],
    root_row: dict[str, object],
    stats: dict[str, object],
) -> None:
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    lines = [
        "# 有限位数根中心 \(10^{-12}\) 拡大スイープ",
        "",
        "## 条件",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| center | {ROOT_TEXT} |",
        f"| delta_R | {DELTA_TEXT} |",
        f"| half steps | {HALF_STEPS} |",
        f"| total points | {POINT_COUNT} |",
        f"| min_R | {stats['min_R']} |",
        f"| max_R | {stats['max_R']} |",
        "",
        "中央格子点を解析的有限位数根 \(R_{124,23}\) に一致させた。",
        "",
        "## 結果",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| best grid index | {best['grid_index']} |",
        f"| best grid R | {best['R_input_text']} |",
        f"| best − exact root | {best['offset_from_R_124_23']} |",
        f"| depth | {best['best_prefix_gray_depth_no_phase']} |",
        f"| error | {best['best_prefix_gray_error_no_phase']} |",
        f"| exact-root depth | {root_row['best_prefix_gray_depth_no_phase']} |",
        f"| exact-root error | {root_row['best_prefix_gray_error_no_phase']} |",
        f"| maximum mirrored depth difference | {stats['max_mirrored_depth_difference']} |",
        "",
        "最深点は中央格子点そのものであり、有限位数根と一致する。",
        "",
        "## 図",
        "",
        "![root centered delta R 1e-12](root_centered_delta1e-12_v1.png)",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    runner = load_precision_runner()
    v5 = runner.load_v5_module()
    refs = runner.references()
    grid = centered_grid()
    rows: list[dict[str, object]] = []
    total_loops = 0
    started = time.perf_counter()

    for index, r_text in enumerate(grid):
        row, loops = runner.evaluate(v5, r_text, "grid", "", index, refs)
        rows.append(row)
        total_loops += loops

    root_row, root_loops = runner.evaluate(
        v5,
        ROOT_TEXT,
        "exact_reference",
        "R_124_23",
        "",
        refs,
    )
    total_loops += root_loops
    elapsed = time.perf_counter() - started

    center = rows[HALF_STEPS]
    mirrored_differences = [
        abs(
            float(rows[HALF_STEPS - offset]["best_prefix_gray_depth_no_phase"])
            - float(rows[HALF_STEPS + offset]["best_prefix_gray_depth_no_phase"])
        )
        for offset in range(1, HALF_STEPS + 1)
    ]
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    stats: dict[str, object] = {
        "center_R": ROOT_TEXT,
        "delta_R": DELTA_TEXT,
        "half_steps": HALF_STEPS,
        "point_count": len(rows),
        "min_R": grid[0],
        "max_R": grid[-1],
        "center_grid_index": HALF_STEPS,
        "center_input_text": center["R_input_text"],
        "center_is_exact_root": center["R_input_text"] == ROOT_TEXT,
        "grid_best": best,
        "exact_root": root_row,
        "max_mirrored_depth_difference": max(mirrored_differences),
        "mean_mirrored_depth_difference": float(np.mean(mirrored_differences)),
        "elapsed_sec": elapsed,
        "total_loop_count": total_loops,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "root_centered_delta1e-12_all_v1.csv", rows)
    write_csv(OUTPUT_DIR / "root_centered_delta1e-12_exact_root_v1.csv", [root_row])
    make_plot(
        rows,
        root_row,
        OUTPUT_DIR / "root_centered_delta1e-12_v1.png",
        OUTPUT_DIR / "root_centered_delta1e-12_v1.svg",
    )
    (OUTPUT_DIR / "root_centered_delta1e-12_summary_v1.json").write_text(
        json.dumps(stats, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "root_centered_delta1e-12_result_v1.md",
        rows,
        root_row,
        stats,
    )
    print(json.dumps(stats, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
