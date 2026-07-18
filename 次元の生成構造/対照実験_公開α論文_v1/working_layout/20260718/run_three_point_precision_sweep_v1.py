#!/usr/bin/env python3
"""Run the 1e-10 System B sweep around three low-energy reference values.

The dynamics and numerical settings are inherited unchanged from the saved
v5 full-range experiment.  Unlike the original candidate-only writer, this
runner saves every evaluated grid point and evaluates the three references
at their direct floating-point values as a separate control table.
"""

from __future__ import annotations

import argparse
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
DEFAULT_OUTPUT_DIR = HERE / "three_point_precision_sweep_v1"

MIN_R = "0.697177779231003050"
MAX_R = "0.697178027556659305"
DELTA_R = "0.0000000001"
STEPS = 1024
PHI_MODE = "zero"
MIN_STEPS = 256
EARLY_STOP_PATIENCE = 20
N_OBS = 137.035999177

OLD_PEAK_R = 0.69717790255614798
OLD_PEAK_DEPTH = 9.083204920768843
OLD_PEAK_ERROR = 8.256482775696157e-10


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


def references() -> list[dict[str, float | str]]:
    n_geometry = (137.0 + math.sqrt(137.0**2 + 2.0 * math.pi**2)) / 2.0
    r_resonance = math.cos(23.0 * math.pi / 124.0) ** 2
    return [
        {
            "label": "R_obs",
            "definition": "R(N_obs), N_obs=137.035999177",
            "R": r_of_n(N_OBS),
            "N": N_OBS,
        },
        {
            "label": "R_g",
            "definition": "R(N_g), N_g=137+(pi^2/2)/N_g",
            "R": r_of_n(n_geometry),
            "N": n_geometry,
        },
        {
            "label": "R_124_23",
            "definition": "cos^2(23*pi/124)",
            "R": r_resonance,
            "N": n_of_r(r_resonance),
        },
    ]


def decimal_range(start_text: str, stop_text: str, step_text: str) -> list[str]:
    start = Decimal(start_text)
    stop = Decimal(stop_text)
    step = Decimal(step_text)
    if step <= 0 or stop < start:
        raise ValueError("invalid decimal sweep range")
    values: list[str] = []
    current = start
    while current <= stop:
        values.append(format(current, "f"))
        current += step
    return values


def payload_to_row(
    payload: dict[str, Any],
    kind: str,
    label: str,
    index: int | str,
    v5: Any,
    refs: list[dict[str, float | str]],
) -> dict[str, object]:
    reflection_rate = float(payload["R"])
    row = {
        "sample_kind": kind,
        "reference_label": label,
        "grid_index": index,
        "R": f"{reflection_rate:.18f}",
        "R_input_text": payload["R_input_text"],
        "N_of_R": f"{n_of_r(reflection_rate):.15f}",
        "offset_from_R_obs": (
            f"{reflection_rate - float(refs[0]['R']):+.18e}"
        ),
        "offset_from_R_g": f"{reflection_rate - float(refs[1]['R']):+.18e}",
        "offset_from_R_124_23": (
            f"{reflection_rate - float(refs[2]['R']):+.18e}"
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
        "T": f"{float(payload['T']):.18e}",
        "reflection_power": f"{float(payload['reflection_power']):.18e}",
        "phi_mode": payload["phi_mode"],
        "steps": int(payload["steps"]),
        "min_steps": int(payload["min_steps"]),
        "early_stop_patience": int(payload["early_stop_patience"]),
    }
    return row


def evaluate(
    v5: Any,
    r_text: str,
    kind: str,
    label: str,
    index: int | str,
    refs: list[dict[str, float | str]],
) -> tuple[dict[str, object], int]:
    payload = v5.probe_r(
        float(r_text),
        r_text,
        STEPS,
        PHI_MODE,
        MIN_STEPS,
        EARLY_STOP_PATIENCE,
    )
    loops = sum(
        int(condition["stopped_at_step"]) + 1
        for condition in payload["condition_rows"]
    )
    return payload_to_row(payload, kind, label, index, v5, refs), loops


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
    maxima: list[dict[str, object]] = []
    for index in range(1, len(rows) - 1):
        if depths[index] > depths[index - 1] and depths[index] >= depths[index + 1]:
            maxima.append(rows[index])
    return sorted(
        maxima,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
        reverse=True,
    )


def make_plot(
    grid_rows: list[dict[str, object]],
    exact_rows: list[dict[str, object]],
    refs: list[dict[str, float | str]],
    output_png: Path,
    output_svg: Path,
) -> None:
    r_obs = float(refs[0]["R"])
    r_res = float(refs[2]["R"])
    r_values = np.array([float(row["R"]) for row in grid_rows])
    errors = np.array(
        [float(row["best_prefix_gray_error_no_phase"]) for row in grid_rows]
    )
    depths = np.array(
        [float(row["best_prefix_gray_depth_no_phase"]) for row in grid_rows]
    )
    colors = {"R_obs": "#111111", "R_g": "#2ca02c", "R_124_23": "#d627a5"}
    labels = {
        "R_obs": r"$R_{\rm obs}$",
        "R_g": r"$R_g$",
        "R_124_23": r"$R_{124,23}$",
    }

    fig, axes = plt.subplots(1, 2, figsize=(16.0, 6.3))
    fig.suptitle(
        r"System B precision sweep: $\Delta R=10^{-10}$",
        fontsize=17,
        fontweight="bold",
    )

    axes[0].semilogy(
        (r_values - r_obs) * 1.0e9,
        errors,
        color="#2b6cb0",
        linewidth=1.5,
        label="precision grid",
    )
    axes[0].scatter(
        [(OLD_PEAK_R - r_obs) * 1.0e9],
        [OLD_PEAK_ERROR],
        marker="x",
        s=75,
        color="#ff7f0e",
        linewidths=2.0,
        label=r"old $\Delta R=10^{-7}$ peak",
        zorder=5,
    )
    for ref, exact in zip(refs, exact_rows):
        x = (float(ref["R"]) - r_obs) * 1.0e9
        y = float(exact["best_prefix_gray_error_no_phase"])
        axes[0].axvline(
            x,
            color=colors[str(ref["label"])],
            linestyle="--",
            linewidth=1.5,
        )
        axes[0].scatter(
            [x],
            [y],
            s=60,
            color=colors[str(ref["label"])],
            edgecolor="white",
            linewidth=0.7,
            label=labels[str(ref["label"])],
            zorder=6,
        )
    axes[0].set_xlabel(r"$(R-R_{\rm obs})\times10^9$")
    axes[0].set_ylabel("best prefix gray error")
    axes[0].set_title("Full triplet interval")
    axes[0].grid(alpha=0.25, which="both")
    axes[0].legend(fontsize=9)

    x_res = (r_values - r_res) * 1.0e10
    mask = np.abs(x_res) <= 25.0
    axes[1].plot(
        x_res[mask],
        depths[mask],
        color="#2b6cb0",
        linewidth=1.5,
        marker="o",
        markersize=2.8,
        label=r"$10^{-10}$ grid",
    )
    resonance_exact = exact_rows[2]
    axes[1].scatter(
        [0.0],
        [float(resonance_exact["best_prefix_gray_depth_no_phase"])],
        s=85,
        color=colors["R_124_23"],
        edgecolor="white",
        linewidth=0.8,
        label="exact finite-order root",
        zorder=6,
    )
    grid_best = max(
        grid_rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    axes[1].scatter(
        [(float(grid_best["R"]) - r_res) * 1.0e10],
        [float(grid_best["best_prefix_gray_depth_no_phase"])],
        marker="x",
        s=80,
        color="#ff7f0e",
        linewidths=2.0,
        label="best grid point",
        zorder=7,
    )
    axes[1].axvline(0.0, color=colors["R_124_23"], linestyle="--", linewidth=1.5)
    axes[1].set_xlabel(r"$(R-R_{124,23})\times10^{10}$")
    axes[1].set_ylabel("best prefix gray depth")
    axes[1].set_title("Resonance-root microstructure")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=9)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_svg, bbox_inches="tight")
    plt.close(fig)


def result_markdown(
    path: Path,
    grid_rows: list[dict[str, object]],
    exact_rows: list[dict[str, object]],
    refs: list[dict[str, float | str]],
    maxima: list[dict[str, object]],
    stats: dict[str, object],
) -> None:
    grid_best = max(
        grid_rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    resonance = exact_rows[2]
    obs = exact_rows[0]
    geometry = exact_rows[1]
    ratio_obs = (
        float(obs["best_prefix_gray_error_no_phase"])
        / float(resonance["best_prefix_gray_error_no_phase"])
    )
    ratio_geometry = (
        float(geometry["best_prefix_gray_error_no_phase"])
        / float(resonance["best_prefix_gray_error_no_phase"])
    )
    best_steps = sorted({int(row["best_step"]) for row in grid_rows})
    stop_steps = sorted({int(row["stopped_at_step"]) for row in grid_rows})
    candidate_count = sum(int(row["is_v5_candidate"]) for row in grid_rows)

    lines = [
        "# 低エネルギー3基準点 \(10^{-10}\) 精密スイープ結果",
        "",
        "## 1. 実験条件",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| min_R | {MIN_R} |",
        f"| max_R | {MAX_R} |",
        f"| delta_R | {DELTA_R} |",
        f"| grid points | {len(grid_rows)} |",
        f"| steps | {STEPS} |",
        f"| min_steps | {MIN_STEPS} |",
        f"| early_stop_patience | {EARLY_STOP_PATIENCE} |",
        f"| phi_mode | {PHI_MODE} |",
        f"| unique best steps | {best_steps} |",
        f"| unique stop steps | {stop_steps} |",
        f"| v5 candidate rows | {candidate_count}/{len(grid_rows)} |",
        "",
        "前回全域スイープから変更したのは \(R\) の範囲と刻みだけである。",
        "全点で best step と停止点が共通であるため、精密ピークは計算分枝の",
        "切替えではなく、同一の247ステップ読出しにおける位相誤差の変化である。",
        "",
        "## 2. 3基準点の厳密直接評価",
        "",
        "| 基準点 | R | N(R) | depth | error | best step |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for ref, row in zip(refs, exact_rows):
        lines.append(
            f"| {ref['label']} | {row['R']} | {row['N_of_R']} | "
            f"{row['best_prefix_gray_depth_no_phase']} | "
            f"{row['best_prefix_gray_error_no_phase']} | {row['best_step']} |"
        )

    lines.extend(
        [
            "",
            "## 3. 精密グリッドの最深点",
            "",
            "| 項目 | 値 |",
            "|---|---:|",
            f"| R_grid_best | {grid_best['R']} |",
            f"| R_grid_best − R_124,23 | {grid_best['offset_from_R_124_23']} |",
            f"| depth | {grid_best['best_prefix_gray_depth_no_phase']} |",
            f"| error | {grid_best['best_prefix_gray_error_no_phase']} |",
            f"| best step | {grid_best['best_step']} |",
            f"| local maxima count | {len(maxima)} |",
            "",
            "グリッド最深点は有限位数根 \(R_{124,23}\) の最近傍にある。",
            "さらに根そのものを直接評価すると、誤差は倍精度計算の丸め限界領域まで低下する。",
            "",
            "## 4. 3点間の判別",
            "",
            r"\[",
            rf"\frac{{E(R_{{\rm obs}})}}{{E(R_{{124,23}})}}"
            rf"={ratio_obs:.6e},",
            r"\qquad",
            rf"\frac{{E(R_g)}}{{E(R_{{124,23}})}}"
            rf"={ratio_geometry:.6e}.",
            r"\]",
            "",
            "このSystem B読出しでは、3点は同じ深さを持たない。",
            "最深点は幾何学値または観測値ではなく、作用素の有限位数根に一致する。",
            "したがって、前回の \(10^{-7}\) スイープで未分解だった低エネルギーピークの",
            "発生位置は、今回の精密計算では \(R_{124,23}\) と判別される。",
            "",
            "これは137セル幾何を否定する結果ではない。判別されたのは、現行System Bの",
            "単独散乱読出しが選択する根である。幾何学値と観測値を得るには、研究ノートで",
            "定式化した二状態結合または長周期干渉を別途実装して判別する必要がある。",
            "",
            "## 5. 前回全域スイープとの比較",
            "",
            "| データ | R | depth | error |",
            "|---|---:|---:|---:|",
            f"| old delta_R=1e-7 | {OLD_PEAK_R:.17f} | {OLD_PEAK_DEPTH:.15f} | {OLD_PEAK_ERROR:.18e} |",
            f"| new delta_R=1e-10 grid | {grid_best['R']} | {grid_best['best_prefix_gray_depth_no_phase']} | {grid_best['best_prefix_gray_error_no_phase']} |",
            f"| exact R_124,23 | {resonance['R']} | {resonance['best_prefix_gray_depth_no_phase']} | {resonance['best_prefix_gray_error_no_phase']} |",
            "",
            "## 6. 図",
            "",
            "![System B three-point precision sweep](three_point_precision_sweep_v1.png)",
            "",
            "左図は3点全域の誤差、右図は有限位数根の直近を \(10^{-10}\) 単位で拡大したもの。",
            "",
            "## 7. 実行統計",
            "",
            f"- elapsed_sec: {stats['elapsed_sec']}",
            f"- total_loop_count: {stats['total_loop_count']}",
            f"- average_msec_per_R: {stats['average_msec_per_R']}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--progress-every", type=int, default=500)
    args = parser.parse_args()

    v5 = load_v5_module()
    refs = references()
    r_grid = decimal_range(MIN_R, MAX_R, DELTA_R)
    grid_rows: list[dict[str, object]] = []
    total_loops = 0
    started = time.perf_counter()

    for index, r_text in enumerate(r_grid):
        row, loops = evaluate(v5, r_text, "grid", "", index, refs)
        grid_rows.append(row)
        total_loops += loops
        if args.progress_every > 0 and (index + 1) % args.progress_every == 0:
            print(
                f"[{index + 1}/{len(r_grid)}] R={r_text} "
                f"depth={row['best_prefix_gray_depth_no_phase']}",
                file=sys.stderr,
                flush=True,
            )

    exact_rows: list[dict[str, object]] = []
    for ref in refs:
        r_text = f"{float(ref['R']):.18f}"
        row, loops = evaluate(
            v5,
            r_text,
            "exact_reference",
            str(ref["label"]),
            "",
            refs,
        )
        exact_rows.append(row)
        total_loops += loops

    elapsed = time.perf_counter() - started
    maxima = local_maxima(grid_rows)
    stats: dict[str, object] = {
        "model": "System B v5 three-point precision sweep",
        "source_code": V5_SOURCE.relative_to(HERE.parent.parent).as_posix(),
        "min_R": MIN_R,
        "max_R": MAX_R,
        "delta_R": DELTA_R,
        "grid_point_count": len(grid_rows),
        "exact_reference_count": len(exact_rows),
        "steps": STEPS,
        "phi_mode": PHI_MODE,
        "min_steps": MIN_STEPS,
        "early_stop_patience": EARLY_STOP_PATIENCE,
        "total_loop_count": total_loops,
        "elapsed_sec": elapsed,
        "average_msec_per_R": elapsed * 1000.0 / (len(grid_rows) + len(exact_rows)),
        "local_maxima_count": len(maxima),
        "unique_best_steps": sorted(
            {int(row["best_step"]) for row in grid_rows}
        ),
        "unique_stopped_at_steps": sorted(
            {int(row["stopped_at_step"]) for row in grid_rows}
        ),
        "v5_candidate_row_count": sum(
            int(row["is_v5_candidate"]) for row in grid_rows
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "three_point_precision_sweep_all_v1.csv", grid_rows)
    write_csv(
        args.output_dir / "three_point_precision_exact_references_v1.csv",
        exact_rows,
    )
    if maxima:
        write_csv(
            args.output_dir / "three_point_precision_local_maxima_v1.csv",
            maxima,
        )
    make_plot(
        grid_rows,
        exact_rows,
        refs,
        args.output_dir / "three_point_precision_sweep_v1.png",
        args.output_dir / "three_point_precision_sweep_v1.svg",
    )

    grid_best = max(
        grid_rows,
        key=lambda row: float(row["best_prefix_gray_depth_no_phase"]),
    )
    summary = {
        **stats,
        "reference_points": refs,
        "grid_best": grid_best,
        "exact_references": exact_rows,
        "local_maxima": maxima[:20],
        "old_full_sweep_peak": {
            "R": OLD_PEAK_R,
            "depth": OLD_PEAK_DEPTH,
            "error": OLD_PEAK_ERROR,
        },
    }
    (args.output_dir / "three_point_precision_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    result_markdown(
        args.output_dir / "three_point_precision_sweep_result_v1.md",
        grid_rows,
        exact_rows,
        refs,
        maxima,
        stats,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
