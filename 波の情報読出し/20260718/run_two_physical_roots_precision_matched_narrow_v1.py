#!/usr/bin/env python3
"""Precision-matched narrow sweeps for R_(124,23) and R_(620,117)."""

from __future__ import annotations

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
BASE_SCRIPT = HERE / "run_two_physical_roots_multiprecision_v1.py"
OUTPUT_DIR = HERE / "two_physical_roots_precision_matched_narrow_v1"

CONFIGS = [
    {
        "name": "50digit_delta1e-50",
        "dps": 50,
        "delta": "1e-50",
        "exponent": 50,
        "y_limits": (47.8, 54.0),
    },
    {
        "name": "80digit_delta1e-80",
        "dps": 80,
        "delta": "1e-80",
        "exponent": 80,
        "y_limits": (77.8, 84.0),
    },
]


def load_base() -> Any:
    spec = importlib.util.spec_from_file_location("multiprecision_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def plot_one(
    axis: plt.Axes,
    result: dict[str, Any],
    config: dict[str, Any],
) -> None:
    rows = result["rows"]
    x_values = np.array([float(row["grid_offset"]) for row in rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    center = result["center"]
    best = result["grid_best"]
    root_latex = str(result["root_latex"])
    exponent = int(config["exponent"])
    dps = int(config["dps"])

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
    center_depth = float(center["best_prefix_gray_depth_no_phase"])
    axis.annotate(
        f"center depth = {center_depth:.6f}",
        xy=(0.0, center_depth),
        xytext=(110.0, center_depth + 0.45),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=10,
    )
    axis.set_xlim(-1250.0, 1250.0)
    axis.set_ylim(*config["y_limits"])
    axis.set_xlabel(rf"$(R-{root_latex})\times10^{{{exponent}}}$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(
        rf"System B precision-matched zoom: ${root_latex}$, "
        rf"{dps} digits, $\Delta R=10^{{-{exponent}}}$"
    )
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")


def save_individual(
    result: dict[str, Any],
    config: dict[str, Any],
) -> None:
    fig, axis = plt.subplots(figsize=(10.5, 6.4))
    plot_one(axis, result, config)
    fig.tight_layout()
    stem = OUTPUT_DIR / f"{result['root_label']}_{config['name']}_narrow_v1"
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def compact(result: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in result.items() if key != "rows"}


def make_markdown(path: Path, results: list[dict[str, Any]]) -> None:
    lines = [
        "# 二つの観測対応根：精度適合局所掃引",
        "",
        "50桁では delta R=1e-50、80桁では delta R=1e-80 とし、両方とも2,485点に統一した。",
        "",
        "| root | precision | delta R | full width | points | edge depth left | center depth | edge depth right |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        rows = result["_edge_rows"]
        lines.append(
            f"| {result['root_label']} | {result['working_decimal_digits']} | "
            f"{result['delta_R']} | {result['_full_width']} | {result['point_count']} | "
            f"{rows[0]['best_prefix_gray_depth_no_phase']} | "
            f"{result['center']['best_prefix_gray_depth_no_phase']} | "
            f"{rows[1]['best_prefix_gray_depth_no_phase']} |"
        )
    for config in CONFIGS:
        lines.extend(
            [
                "",
                f"## {config['dps']}桁・delta R={config['delta']}",
                "",
                f"![{config['name']} narrow comparison]({config['name']}_narrow_comparison_v1.png)",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    base = load_base()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_results: list[dict[str, Any]] = []
    markdown_results: list[dict[str, Any]] = []
    started = time.perf_counter()

    for config in CONFIGS:
        config_results: list[dict[str, Any]] = []
        for root in base.ROOTS:
            result = base.run_sweep(root, config)
            config_results.append(result)
            base.write_csv(
                OUTPUT_DIR / f"{root['stem']}_{config['name']}_all_v1.csv",
                result["rows"],
            )
            save_individual(result, config)
            compact_result = compact(result)
            summary_results.append(compact_result)
            markdown_result = dict(compact_result)
            markdown_result["_edge_rows"] = [result["rows"][0], result["rows"][-1]]
            markdown_result["_full_width"] = (
                f"{2 * base.HALF_STEPS}e-{int(config['exponent'])}"
            )
            markdown_results.append(markdown_result)
            print(
                json.dumps(
                    {
                        "root": root["stem"],
                        "dps": config["dps"],
                        "delta_R": config["delta"],
                        "center_depth": result["center"]["best_prefix_gray_depth_no_phase"],
                    },
                    ensure_ascii=False,
                ),
                flush=True,
            )

        fig, axes = plt.subplots(1, 2, figsize=(17.0, 6.3), sharex=True, sharey=True)
        for axis, result in zip(axes, config_results):
            plot_one(axis, result, config)
        fig.suptitle(
            f"System B precision-matched root zooms: {config['dps']} digits",
            fontsize=17,
            fontweight="bold",
        )
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
        fig.savefig(
            OUTPUT_DIR / f"{config['name']}_narrow_comparison_v1.png",
            dpi=220,
            bbox_inches="tight",
        )
        fig.savefig(
            OUTPUT_DIR / f"{config['name']}_narrow_comparison_v1.svg",
            bbox_inches="tight",
        )
        plt.close(fig)

    summary = {
        "experiment": "precision-matched narrow multiprecision sweeps",
        "sampling_count_per_sweep": base.POINT_COUNT,
        "half_steps": base.HALF_STEPS,
        "configs": CONFIGS,
        "results": summary_results,
        "elapsed_sec": time.perf_counter() - started,
    }
    (OUTPUT_DIR / "precision_matched_narrow_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "precision_matched_narrow_result_v1.md",
        markdown_results,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
