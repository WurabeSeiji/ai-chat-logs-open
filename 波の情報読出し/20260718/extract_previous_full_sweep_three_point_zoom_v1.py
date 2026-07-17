#!/usr/bin/env python3
"""Prepare old full-sweep control data around three low-energy R values.

No dynamics are recomputed here.  The saved 2026-07-15 candidate CSV is
extracted, plotted, and compared with R_obs, R_g, and R_(124,23).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from decimal import Decimal
from pathlib import Path
from statistics import median

PLOT_CACHE = Path(tempfile.gettempdir()) / "wave_readout_matplotlib_cache"
PLOT_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(PLOT_CACHE / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = (
    HERE.parent
    / "20260715"
    / "minimal_system_B_gray_bugcheck_result_v1"
    / "direct_depth_probe_v5_sweep_control"
    / "high_to_ext_full_delta1e-7_candidates_v5.csv"
)
DEFAULT_OUTPUT_DIR = HERE / "three_point_precision_sweep_preparation_v1"

N_OBS = 137.035999177
LOCAL_HALF_WINDOW = 3.0e-6
MICRO_HALF_WINDOW = 1.5e-7
PRECISION_DELTA_R = 1.0e-10
PRECISION_DELTA_TEXT = "0.0000000001"
PRECISION_MARGIN_TEXT = "0.000000100000000000"


def n_of_r(reflection_rate: float) -> float:
    return 4.0 * math.pi / (1.0 - reflection_rate) ** 2


def r_of_n(capacity: float) -> float:
    return 1.0 - math.sqrt(4.0 * math.pi / capacity)


def precision_bounds(
    points: list[dict[str, float | str]],
) -> tuple[Decimal, Decimal, int]:
    r_min = Decimal(f"{min(float(point['R']) for point in points):.18f}")
    r_max = Decimal(f"{max(float(point['R']) for point in points):.18f}")
    margin = Decimal(PRECISION_MARGIN_TEXT)
    step = Decimal(PRECISION_DELTA_TEXT)
    start = r_min - margin
    stop = r_max + margin
    count = int((stop - start) // step) + 1
    return start, stop, count


def reference_points() -> list[dict[str, float | str]]:
    n_geometry = (137.0 + math.sqrt(137.0**2 + 2.0 * math.pi**2)) / 2.0
    r_resonance = math.cos(23.0 * math.pi / 124.0) ** 2
    points: list[dict[str, float | str]] = [
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
    r_obs = float(points[0]["R"])
    for point in points:
        point["delta_R_from_R_obs"] = float(point["R"]) - r_obs
    return points


def read_source(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"source CSV is empty: {path}")
    return rows


def find_grid_step(rows: list[dict[str, str]]) -> float:
    values = sorted({float(row["R"]) for row in rows})
    local_gaps = [
        right - left
        for left, right in zip(values, values[1:])
        if 0.0 < right - left < 1.0e-5
    ]
    return median(local_gaps)


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def prepare_rows(
    source_rows: list[dict[str, str]],
    points: list[dict[str, float | str]],
    half_window: float,
    grid_step: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], float]:
    r_obs = float(points[0]["R"])
    r_min = min(float(point["R"]) for point in points) - half_window
    r_max = max(float(point["R"]) for point in points) + half_window
    global_max_depth = max(
        float(row["best_prefix_gray_depth_no_phase"]) for row in source_rows
    )

    zoom_rows: list[dict[str, object]] = []
    for source in source_rows:
        reflection_rate = float(source["R"])
        if not r_min <= reflection_rate <= r_max:
            continue
        nearest = min(points, key=lambda point: abs(reflection_rate - float(point["R"])))
        row: dict[str, object] = dict(source)
        row.update(
            {
                "offset_from_R_obs": f"{reflection_rate - r_obs:+.18e}",
                "normalized_depth_global": (
                    float(source["best_prefix_gray_depth_no_phase"])
                    / global_max_depth
                ),
                "nearest_reference": nearest["label"],
                "delta_R_to_nearest_reference": (
                    f"{reflection_rate - float(nearest['R']):+.18e}"
                ),
            }
        )
        zoom_rows.append(row)
    if not zoom_rows:
        raise ValueError("no source rows in zoom window")

    nearest_rows: list[dict[str, object]] = []
    for point in points:
        target_r = float(point["R"])
        source = min(source_rows, key=lambda row: abs(float(row["R"]) - target_r))
        source_r = float(source["R"])
        nearest_rows.append(
            {
                "reference": point["label"],
                "definition": point["definition"],
                "target_R": f"{target_r:.18f}",
                "target_N": (
                    "137.035999177"
                    if point["label"] == "R_obs"
                    else f"{float(point['N']):.15f}"
                ),
                "delta_R_from_R_obs": (
                    f"{float(point['delta_R_from_R_obs']):+.18e}"
                ),
                "nearest_source_R": source["R_input_text"],
                "source_minus_target_R": f"{source_r - target_r:+.18e}",
                "absolute_offset_in_old_grid_steps": (
                    f"{abs(source_r - target_r) / grid_step:.9f}"
                ),
                "best_prefix_gray_depth_no_phase": source[
                    "best_prefix_gray_depth_no_phase"
                ],
                "normalized_depth_global": (
                    f"{float(source['best_prefix_gray_depth_no_phase']) / global_max_depth:.15f}"
                ),
                "best_prefix_gray_error_no_phase": source[
                    "best_prefix_gray_error_no_phase"
                ],
                "best_step": source["best_step"],
                "condition": source["best_condition_id"],
            }
        )
    return zoom_rows, nearest_rows, global_max_depth


def make_plot(
    zoom_rows: list[dict[str, object]],
    points: list[dict[str, float | str]],
    grid_step: float,
    png_path: Path,
    svg_path: Path,
) -> None:
    r_obs = float(points[0]["R"])
    r_values = np.array([float(row["R"]) for row in zoom_rows])
    depths = np.array([float(row["normalized_depth_global"]) for row in zoom_rows])
    colors = {"R_obs": "#111111", "R_g": "#2ca02c", "R_124_23": "#d627a5"}
    labels = {
        "R_obs": r"$R_{\rm obs}$",
        "R_g": r"$R_g$",
        "R_124_23": r"$R_{124,23}$",
    }

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2))
    fig.suptitle(
        "Previous full-sweep control near the low-energy triplet",
        fontsize=16,
        fontweight="bold",
    )

    axes[0].plot(
        (r_values - r_obs) * 1.0e6,
        depths,
        color="#2b6cb0",
        linewidth=1.8,
        marker="o",
        markersize=2.5,
        label=r"saved sweep, $\Delta R=10^{-7}$",
    )
    for point in points:
        axes[0].axvline(
            (float(point["R"]) - r_obs) * 1.0e6,
            color=colors[str(point["label"])],
            linestyle="--",
            linewidth=1.7,
            label=labels[str(point["label"])],
        )
    axes[0].set_xlabel(r"$(R-R_{\rm obs})\times10^6$")
    axes[0].set_ylabel("normalized candidate depth (global max = 1)")
    axes[0].set_title("Local peak retained in the 2026-07-15 sweep")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=9)

    mask = np.abs(r_values - r_obs) <= MICRO_HALF_WINDOW
    axes[1].plot(
        (r_values[mask] - r_obs) * 1.0e9,
        depths[mask],
        color="#2b6cb0",
        linewidth=1.8,
        marker="o",
        markersize=6,
        label="saved grid points",
    )
    for point in points:
        offset = (float(point["R"]) - r_obs) * 1.0e9
        axes[1].axvline(
            offset,
            color=colors[str(point["label"])],
            linestyle="--",
            linewidth=2.0,
            label=labels[str(point["label"])],
        )
        axes[1].annotate(
            f"{offset:+.2f}",
            xy=(offset, 0.985),
            xycoords=("data", "axes fraction"),
            xytext=(3, 0),
            textcoords="offset points",
            rotation=90,
            va="top",
            fontsize=8,
            color=colors[str(point["label"])],
        )
    axes[1].set_xlim(-MICRO_HALF_WINDOW * 1.0e9, MICRO_HALF_WINDOW * 1.0e9)
    axes[1].set_xlabel(r"$(R-R_{\rm obs})\times10^9$")
    axes[1].set_ylabel("normalized candidate depth (global max = 1)")
    axes[1].set_title(
        "Three references lie inside one old-grid interval\n"
        + rf"(old $\Delta R={grid_step:.1e}$)"
    )
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=9, loc="lower right")

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def make_markdown(
    path: Path,
    source_path: Path,
    points: list[dict[str, float | str]],
    nearest_rows: list[dict[str, object]],
    grid_step: float,
    zoom_count: int,
) -> None:
    triplet_span = float(points[-1]["R"]) - float(points[0]["R"])
    r_min, r_max, point_count = precision_bounds(points)
    same_row = len({str(row["nearest_source_R"]) for row in nearest_rows}) == 1

    lines = [
        "# 低エネルギー3基準点：前回全域スイープ対照データ",
        "",
        "## 1. 対象",
        "",
        r"\[",
        r"R_{\mathrm{obs}}=R(137.035999177),\qquad",
        r"R_g=R\!\left(\frac{137+\sqrt{137^2+2\pi^2}}{2}\right),\qquad",
        r"R_{124,23}=\cos^2\!\left(\frac{23\pi}{124}\right).",
        r"\]",
        "",
        "入力CSV:",
        "",
        source_path.as_posix(),
        "",
        rf"前回全域スイープの刻みは \(\Delta R_{{\rm old}}={grid_step:.1e}\)。",
        f"局所抽出行数は {zoom_count} 行。",
        "",
        "## 2. 3基準点と前回CSVの最寄り点",
        "",
        "| 基準点 | target R | target N | R−R_obs | 前回CSV最寄りR | CSV−target | depth | error | step |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in nearest_rows:
        lines.append(
            "| {reference} | {target_R} | {target_N} | {delta_R_from_R_obs} | "
            "{nearest_source_R} | {source_minus_target_R} | "
            "{best_prefix_gray_depth_no_phase} | "
            "{best_prefix_gray_error_no_phase} | {best_step} |".format(**row)
        )

    lines.extend(["", "## 3. 対照データの判定", ""])
    if same_row:
        lines.extend(
            [
                "3基準点は、前回CSVでは**すべて同一の格子点**へ写る。",
                "",
                r"\[",
                rf"R_{{124,23}}-R_{{\rm obs}}={triplet_span:.18e}"
                rf"<\Delta R_{{\rm old}}.",
                r"\]",
                "",
                "前回全域スイープは低エネルギー候補帯と局所最大の存在を示すが、",
                "3点間の順位、谷、分裂を判別する解像度を持たない。",
                "これは \(10^{-10}\) スイープの必要性を直接示す対照結果である。",
            ]
        )
    else:
        lines.append("3基準点は前回CSVでも異なる格子点へ対応した。")

    lines.extend(
        [
            "",
            "## 4. 拡大図",
            "",
            "![Previous full-sweep control near the low-energy triplet](previous_full_sweep_three_point_zoom_v1.png)",
            "",
            "左図は保存済み候補帯の局所形状、右図は3基準点と旧格子点の位置関係を示す。",
            "",
            "## 5. 次段階の精密スイープ案",
            "",
            "| 項目 | 値 |",
            "|---|---:|",
            f"| min_R | {r_min} |",
            f"| max_R | {r_max} |",
            f"| delta_R | {PRECISION_DELTA_R:.1e} |",
            f"| 予定点数 | {point_count} |",
            "| 直接評価 | 3基準点をグリッドとは別に厳密値で評価 |",
            "",
            "一つの連続区間で3点全体を覆い、全評価点をCSVへ保存する。",
            "これにより候補選別規則による欠落を避け、ピーク位置と非対称性を比較する。",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--half-window", type=float, default=LOCAL_HALF_WINDOW)
    args = parser.parse_args()

    source_rows = read_source(args.source_csv)
    points = reference_points()
    grid_step = find_grid_step(source_rows)
    try:
        source_display = args.source_csv.resolve().relative_to(
            HERE.parent.parent
        )
    except ValueError:
        source_display = args.source_csv
    zoom_rows, nearest_rows, _ = prepare_rows(
        source_rows, points, args.half_window, grid_step
    )

    zoom_fields = list(source_rows[0].keys()) + [
        "offset_from_R_obs",
        "normalized_depth_global",
        "nearest_reference",
        "delta_R_to_nearest_reference",
    ]
    write_csv(
        args.output_dir / "previous_full_sweep_three_point_zoom_v1.csv",
        zoom_rows,
        zoom_fields,
    )
    write_csv(
        args.output_dir / "previous_full_sweep_three_point_nearest_v1.csv",
        nearest_rows,
        list(nearest_rows[0].keys()),
    )
    make_plot(
        zoom_rows,
        points,
        grid_step,
        args.output_dir / "previous_full_sweep_three_point_zoom_v1.png",
        args.output_dir / "previous_full_sweep_three_point_zoom_v1.svg",
    )
    make_markdown(
        args.output_dir / "previous_full_sweep_three_point_control_v1.md",
        source_display,
        points,
        nearest_rows,
        grid_step,
        len(zoom_rows),
    )

    r_min, r_max, point_count = precision_bounds(points)
    plan = {
        "source_csv": source_display.as_posix(),
        "source_grid_step": grid_step,
        "reference_points": points,
        "zoom_half_window": args.half_window,
        "zoom_row_count": len(zoom_rows),
        "all_references_share_one_nearest_source_row": (
            len({str(row["nearest_source_R"]) for row in nearest_rows}) == 1
        ),
        "precision_sweep_plan": {
            "min_R": str(r_min),
            "max_R": str(r_max),
            "delta_R": f"{PRECISION_DELTA_R:.1e}",
            "point_count": point_count,
            "evaluate_reference_points_exactly": True,
            "save_all_evaluated_rows": True,
        },
    }
    (args.output_dir / "three_point_precision_sweep_plan_v1.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(plan, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
