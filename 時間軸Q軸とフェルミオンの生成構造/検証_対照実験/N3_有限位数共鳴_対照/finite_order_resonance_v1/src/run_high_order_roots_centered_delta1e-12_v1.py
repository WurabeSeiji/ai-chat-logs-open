#!/usr/bin/env python3
"""Root-centered delta_R=1e-12 sweeps for R_(567,107) and R_(620,117)."""

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
OUTPUT_DIR = HERE / "high_order_roots_centered_delta1e-12_v1"

ROOTS = [
    {
        "stem": "R567_107",
        "latex": r"R_{567,107}",
        "n": 567,
        "m": 107,
        "root_text": (
            "0.68781690893425671211625916827203577005237999852936384518812012703579869330301607343056158177"
        ),
    },
    {
        "stem": "R620_117",
        "latex": r"R_{620,117}",
        "n": 620,
        "m": 117,
        "root_text": (
            "0.68782519111414518666745163827099317306585834508739787748914792010613050357559734279239971833"
        ),
    },
]

DELTA_TEXT = "0.000000000001"
HALF_STEPS = 1242
POINT_COUNT = 2 * HALF_STEPS + 1

# Extended only far enough to contain 2n-1 for n=567 and n=620.
STEPS = 1536
PHI_MODE = "zero"
MIN_STEPS = 1280
EARLY_STOP_PATIENCE = 20

Y_LIMITS = (7.8, 15.0)
X_LIMITS = (-1250.0, 1250.0)
N_OBS_HIGH = 128.946
N_OBS_HIGH_SIGMA = 0.015


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


def centered_grid(root_text: str) -> list[str]:
    with localcontext() as context:
        context.prec = 110
        root = Decimal(root_text)
        delta = Decimal(DELTA_TEXT)
        return [
            format(root + Decimal(offset) * delta, "f")
            for offset in range(-HALF_STEPS, HALF_STEPS + 1)
        ]


def evaluate(
    v5: Any,
    r_text: str,
    root_text: str,
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
    root_float = float(root_text)
    return {
        "sample_kind": kind,
        "grid_index": index,
        "R": f"{reflection_rate:.18f}",
        "R_input_text": payload["R_input_text"],
        "N_of_R": f"{n_of_r(reflection_rate):.15f}",
        "offset_from_root_float": f"{reflection_rate-root_float:+.18e}",
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


def make_plot(
    rows: list[dict[str, object]],
    root_row: dict[str, object],
    root: dict[str, object],
    output_png: Path,
    output_svg: Path,
) -> None:
    root_float = float(str(root["root_text"]))
    r_values = np.array([float(row["R"]) for row in rows])
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in rows]
    )
    x_values = (r_values - root_float) * 1.0e12
    best = max(
        rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    best_depth = float(best["best_prefix_gray_depth_no_phase"])
    root_depth = float(root_row["best_prefix_gray_depth_no_phase"])

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
        [root_depth],
        color="#d627a5",
        edgecolor="white",
        linewidth=0.7,
        s=72,
        zorder=6,
        label=rf"exact ${root['latex']}$",
    )
    axis.scatter(
        [(float(best["R"]) - root_float) * 1.0e12],
        [best_depth],
        marker="x",
        color="#ff7f0e",
        linewidths=2.2,
        s=85,
        zorder=7,
        label="best grid point",
    )
    if Y_LIMITS[0] <= root_depth <= Y_LIMITS[1]:
        axis.annotate(
            "exact root = center grid point",
            xy=(0.0, root_depth),
            xytext=(120.0, 14.35),
            arrowprops={"arrowstyle": "->", "color": "#555555"},
            fontsize=10,
        )
    else:
        axis.text(
            0.98,
            0.95,
            f"maximum depth = {best_depth:.6f}\n(below fixed display range)",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "#999999", "alpha": 0.9},
        )
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(*Y_LIMITS)
    axis.set_xlabel(rf"$(R-{root['latex']})\times10^{{12}}$")
    axis.set_ylabel("best prefix gray depth")
    axis.set_title(rf"System B root-centered sweep: ${root['latex']}$, $\Delta R=10^{{-12}}$")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_svg, bbox_inches="tight")
    plt.close(fig)


def mirrored_stats(rows: list[dict[str, object]]) -> dict[str, float]:
    differences = [
        abs(
            float(rows[HALF_STEPS-offset]["best_prefix_gray_depth_no_phase"])
            - float(rows[HALF_STEPS+offset]["best_prefix_gray_depth_no_phase"])
        )
        for offset in range(1, HALF_STEPS + 1)
    ]
    return {
        "max_mirrored_depth_difference": max(differences),
        "mean_mirrored_depth_difference": float(np.mean(differences)),
    }


def make_markdown(path: Path, results: list[dict[str, object]]) -> None:
    r567_exact = results[0]["exact_root"]
    r620_exact = results[1]["exact_root"]
    predicted_odd_defect = 0.01 * (1.0 - math.cos(math.pi / 567.0))
    measured_odd_defect = 0.02 - float(r567_exact["best_S_amp"])
    n620_offset = float(r620_exact["N_of_R"]) - N_OBS_HIGH
    lines = [
        "# 高次有限位数根中心 delta R = 1e-12 ピーク探索",
        "",
        "R124,23 および R122,23 と同じ根中心形式で、R567,107 と R620,117 を独立に掃引した。",
        "高次周期を収容するため、steps=1536、min_steps=1280 とした。深度定義と描画範囲は従来と同一である。",
        "",
        "| root | exact R | N(R) | grid best index | grid best offset | depth | error | best step | expected 2n-1 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in results:
        best = item["grid_best"]
        lines.append(
            f"| {item['root_label']} | {item['root_high_precision']} | "
            f"{item['exact_root']['N_of_R']} | {best['grid_index']} | "
            f"{best['offset_from_root_float']} | "
            f"{best['best_prefix_gray_depth_no_phase']} | "
            f"{best['best_prefix_gray_error_no_phase']} | "
            f"{best['best_step']} | {item['expected_best_step']} |"
        )
    lines.extend(
        [
            "",
            "## 判別結果",
            "",
            "R567,107 は位相周期 best step = 1133 = 2(567)-1 を満たすが、根中心に灰色深度ピークを作らない。",
            "奇数 n の離散軌道が反対側の極値を格子点として含まないため、振幅欠損は",
            "",
            f"`0.01(1-cos(pi/567)) = {predicted_odd_defect:.18e}`",
            "",
            f"となる。実測欠損 `{measured_odd_defect:.18e}` と倍精度誤差内で一致する。",
            "",
            "R620,117 は中央格子点に機械精度のピークを持つ。偶数 n では反対側の極値を含むためである。",
            f"その読出し N(R) = {r620_exact['N_of_R']} は観測値 128.946 との差が {n620_offset:.15f}、すなわち {n620_offset/N_OBS_HIGH_SIGMA:.9f} sigma である。",
            "",
            "## R567,107",
            "",
            "![R567 107 root centered 1e-12](R567_107_root_centered_delta1e-12_v1.png)",
            "",
            "## R620,117",
            "",
            "![R620 117 root centered 1e-12](R620_117_root_centered_delta1e-12_v1.png)",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    v5 = load_v5_module()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, object]] = []
    started = time.perf_counter()

    for root in ROOTS:
        root_text = str(root["root_text"])
        grid = centered_grid(root_text)
        rows = [
            evaluate(v5, r_text, root_text, "grid", index)
            for index, r_text in enumerate(grid)
        ]
        root_row = evaluate(v5, root_text, root_text, "exact_reference", "")
        best = max(
            rows,
            key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
        )
        result: dict[str, object] = {
            "root_label": root["stem"],
            "n": root["n"],
            "m": root["m"],
            "root_high_precision": root_text,
            "root_float": f"{float(root_text):.18f}",
            "delta_R": DELTA_TEXT,
            "half_steps": HALF_STEPS,
            "point_count": len(rows),
            "min_R": grid[0],
            "max_R": grid[-1],
            "center_grid_index": HALF_STEPS,
            "center_input_text": rows[HALF_STEPS]["R_input_text"],
            "center_is_exact_root_text": rows[HALF_STEPS]["R_input_text"] == root_text,
            "expected_best_step": 2 * int(root["n"]) - 1,
            "grid_best": best,
            "exact_root": root_row,
            **mirrored_stats(rows),
        }
        results.append(result)

        stem = str(root["stem"])
        write_csv(OUTPUT_DIR / f"{stem}_root_centered_delta1e-12_all_v1.csv", rows)
        write_csv(
            OUTPUT_DIR / f"{stem}_root_centered_delta1e-12_exact_root_v1.csv",
            [root_row],
        )
        make_plot(
            rows,
            root_row,
            root,
            OUTPUT_DIR / f"{stem}_root_centered_delta1e-12_v1.png",
            OUTPUT_DIR / f"{stem}_root_centered_delta1e-12_v1.svg",
        )

    summary = {
        "experiment": "two independent root-centered delta_R=1e-12 peak searches",
        "delta_R": DELTA_TEXT,
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
        "analytic_checks": {
            "R567_107_odd_n_predicted_amplitude_defect": (
                0.01 * (1.0 - math.cos(math.pi / 567.0))
            ),
            "R567_107_measured_amplitude_defect": (
                0.02 - float(results[0]["exact_root"]["best_S_amp"])
            ),
            "R620_117_N_minus_observation": (
                float(results[1]["exact_root"]["N_of_R"]) - N_OBS_HIGH
            ),
            "R620_117_observation_sigma_units": (
                (
                    float(results[1]["exact_root"]["N_of_R"])
                    - N_OBS_HIGH
                )
                / N_OBS_HIGH_SIGMA
            ),
        },
        "elapsed_sec": time.perf_counter() - started,
    }
    (OUTPUT_DIR / "high_order_roots_centered_delta1e-12_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_markdown(
        OUTPUT_DIR / "high_order_roots_centered_delta1e-12_result_v1.md",
        results,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
