#!/usr/bin/env python3
"""Create Stage A1 candidate figures without selecting any event."""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PACKAGE_ROOT / "config_candidates.json"
PROCESSED_DIR = PACKAGE_ROOT / "processed"
FIGURE_DIR = PACKAGE_ROOT / "figures"
LOG_DIR = PACKAGE_ROOT / "logs"
GROWTH_MANIFEST_PATH = LOG_DIR / "growth_analysis_manifest.json"
RANK_MANIFEST_PATH = LOG_DIR / "rank_analysis_manifest.json"
FIGURE_MANIFEST_PATH = LOG_DIR / "figures_manifest.json"
TEXT_LOG_PATH = LOG_DIR / "make_calibration_figures.log"

os.environ.setdefault("MPLCONFIGDIR", str(LOG_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(LOG_DIR / "cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


FIGURE_NAMES = [
    "figure01_f_and_log10",
    "figure02_regression_slopes",
    "figure03_regression_r2",
    "figure04_growth_intervals_by_window",
    "figure05_growth_end_candidate_distribution",
    "figure06_q1_q4_and_rank_q",
    "figure07_q_ratios",
    "figure08_rank4_onset_heatmap",
    "figure09_growth_end_vs_rank4_difference",
    "figure10_crossing_1167_zoom",
]


def resolve(relative: str) -> Path:
    return (PACKAGE_ROOT / relative).resolve()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def save(fig, stem: str) -> None:
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=150)
    fig.savefig(FIGURE_DIR / f"{stem}.svg")
    plt.close(fig)


def main() -> int:
    started = time.perf_counter()
    planned = [
        FIGURE_DIR / f"{name}.{extension}"
        for name in FIGURE_NAMES
        for extension in ("png", "svg")
    ] + [FIGURE_MANIFEST_PATH, TEXT_LOG_PATH]
    if any(path.exists() for path in planned):
        raise RuntimeError("校正図出力の上書きを避けて停止")
    if not GROWTH_MANIFEST_PATH.is_file() or not RANK_MANIFEST_PATH.is_file():
        raise RuntimeError("成長候補またはrank候補解析の成功記録がない")
    growth_manifest = json.loads(GROWTH_MANIFEST_PATH.read_text(encoding="utf-8"))
    rank_manifest = json.loads(RANK_MANIFEST_PATH.read_text(encoding="utf-8"))
    if growth_manifest.get("success") is not True or rank_manifest.get("success") is not True:
        raise RuntimeError("先行候補解析が成功していない")

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    crossing = int(config["existing_crossing"])
    windows = [int(value) for value in config["regression_windows"]]
    r2_thresholds = [float(value) for value in config["r2_thresholds"]]
    rank_thresholds = [float(value) for value in config["rank_relative_thresholds"]]
    rank_persistences = [int(value) for value in config["rank_persistence_records"]]

    f_rows = load_csv(resolve(config["inputs"]["fcurve"]["path"]))
    paper_rows = load_csv(resolve(config["inputs"]["paper7_long"]["path"]))
    growth_metrics = load_csv(PROCESSED_DIR / "f_growth_metrics.csv")
    intervals = load_csv(PROCESSED_DIR / "growth_intervals_all_candidates.csv")
    growth_ends = load_csv(PROCESSED_DIR / "growth_end_all_candidates.csv")
    rank_metrics = load_csv(PROCESSED_DIR / "q_rank_candidate_metrics.csv")
    rank_onsets = load_csv(PROCESSED_DIR / "rank4_onset_all_candidates.csv")

    f_step = np.array([int(row["tau"]) for row in f_rows])
    f_value = np.array([float(row["f"]) for row in f_rows])
    step = np.array([int(row["step"]) for row in growth_metrics])
    log_f = np.array([as_float(row["log_f"]) for row in growth_metrics])
    log10_f = np.array([as_float(row["log10_f"]) for row in growth_metrics])
    slopes = {
        window: np.array(
            [as_float(row[f"slope_w{window:03d}"]) for row in growth_metrics]
        )
        for window in windows
    }
    r2_values = {
        window: np.array(
            [as_float(row[f"r2_w{window:03d}"]) for row in growth_metrics]
        )
        for window in windows
    }

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    (LOG_DIR / "matplotlib").mkdir(parents=True, exist_ok=True)
    (LOG_DIR / "cache").mkdir(parents=True, exist_ok=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(windows)))

    # Figure 1
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    axes[0].plot(f_step, f_value, color="black", lw=0.8)
    axes[0].axvline(crossing, color="tab:red", ls="--", lw=1.0)
    axes[0].set_ylabel("f")
    axes[0].set_title("Figure 1a: reproduced f(t), full available range")
    axes[1].plot(step, log10_f, color="tab:blue", lw=0.8)
    axes[1].axvline(crossing, color="tab:red", ls="--", lw=1.0, label="existing crossing=1167")
    axes[1].set_xlabel("absolute step")
    axes[1].set_ylabel("log10(f)")
    axes[1].legend()
    axes[1].set_title("Figure 1b: log10(f); f<=0 would remain NaN")
    fig.tight_layout()
    save(fig, FIGURE_NAMES[0])

    # Figure 2
    fig, ax = plt.subplots(figsize=(13, 6))
    for color, window in zip(colors, windows):
        ax.plot(step, slopes[window], lw=0.65, color=color, label=f"window={window}")
    ax.axhline(0, color="black", lw=0.7)
    ax.axvline(crossing, color="tab:red", ls="--", lw=0.9)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("OLS slope of log(f)")
    ax.set_title("Figure 2: slope comparison for every regression window")
    ax.legend(ncol=3)
    fig.tight_layout()
    save(fig, FIGURE_NAMES[1])

    # Figure 3
    fig, ax = plt.subplots(figsize=(13, 6))
    for color, window in zip(colors, windows):
        ax.plot(step, r2_values[window], lw=0.65, color=color, label=f"window={window}")
    for threshold in r2_thresholds:
        ax.axhline(threshold, color="gray", lw=0.35, alpha=0.45)
    ax.axvline(crossing, color="tab:red", ls="--", lw=0.9)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("R squared")
    ax.set_title("Figure 3: R squared comparison; all candidate thresholds shown equally")
    ax.legend(ncol=3)
    fig.tight_layout()
    save(fig, FIGURE_NAMES[2])

    # Figure 4
    threshold_colors = {
        threshold: plt.cm.plasma(index / max(1, len(r2_thresholds) - 1))
        for index, threshold in enumerate(r2_thresholds)
    }
    fig, axes = plt.subplots(3, 2, figsize=(15, 14), sharex=True, sharey=True)
    for ax, window in zip(axes.ravel(), windows):
        ax.plot(step, log_f, color="black", lw=0.55)
        subset = [row for row in intervals if int(row["window"]) == window]
        for row in subset:
            threshold = float(row["r2_threshold"])
            ax.axvspan(
                int(row["interval_start"]),
                int(row["interval_end"]),
                color=threshold_colors[threshold],
                alpha=0.012,
                lw=0,
            )
        ax.axvline(crossing, color="tab:red", ls="--", lw=0.7)
        ax.set_title(f"window={window}; candidate rows={len(subset)}")
        ax.set_ylabel("log(f)")
    axes[-1, 0].set_xlabel("absolute step")
    axes[-1, 1].set_xlabel("absolute step")
    handles = [
        Line2D([0], [0], color=threshold_colors[value], lw=5, alpha=0.65, label=f"R2 >= {value}")
        for value in r2_thresholds
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3)
    fig.suptitle(
        "Figure 4: every exponential-interval candidate over log(f); no adopted interval",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.98))
    save(fig, FIGURE_NAMES[3])

    # Figure 5
    found_ends = [row for row in growth_ends if row["status"] == "found"]
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    r2_color_map = {
        threshold: plt.cm.viridis(index / max(1, len(r2_thresholds) - 1))
        for index, threshold in enumerate(r2_thresholds)
    }
    duration_values = [int(value) for value in config["growth_minimum_durations"]]
    duration_offsets = {
        duration: (index - (len(duration_values) - 1) / 2) * 0.07
        for index, duration in enumerate(duration_values)
    }
    for ax, condition_name in zip(axes, config["growth_end_conditions"]):
        subset = [row for row in found_ends if row["end_condition"] == condition_name]
        for persistence in config["growth_end_persistence"]:
            current = [
                row for row in subset if int(row["end_persistence"]) == int(persistence)
            ]
            x = [int(row["growth_end_candidate"]) for row in current]
            y = [
                windows.index(int(row["window"]))
                + duration_offsets[int(row["minimum_duration"])]
                for row in current
            ]
            c = [r2_color_map[float(row["r2_threshold"])] for row in current]
            ax.scatter(
                x,
                y,
                c=c,
                s=4 + 0.12 * int(persistence),
                alpha=0.22,
                edgecolors="none",
                rasterized=True,
                label=f"end persistence={persistence}",
            )
        ax.set_yticks(range(len(windows)), [str(value) for value in windows])
        ax.set_ylabel("regression window")
        ax.set_title(f"condition {condition_name}; all found candidates={len(subset)}")
    axes[-1].set_xlabel("growth-end candidate absolute step")
    axes[0].legend(ncol=4, fontsize=8)
    fig.suptitle(
        "Figure 5: growth-end candidate distribution by window, R2 color, durations; no selection"
    )
    fig.tight_layout()
    save(fig, FIGURE_NAMES[4])

    # Figure 6
    q_step = np.array([int(row["step"]) for row in rank_metrics])
    q_values = {
        name: np.array([float(row[name]) for row in rank_metrics])
        for name in ("q1", "q2", "q3", "q4")
    }
    rank_q = np.array(
        [int(row["rank_q_existing_1e-8"]) for row in rank_metrics]
    )
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for name in ("q1", "q2", "q3", "q4"):
        axes[0].plot(q_step, q_values[name], lw=0.8, label=name)
    axes[0].axvline(crossing, color="tab:red", ls="--", lw=0.9)
    axes[0].set_ylabel("q value")
    axes[0].legend(ncol=4)
    axes[0].set_title("Figure 6a: existing q1 to q4")
    axes[1].step(q_step, rank_q, where="post", color="black", lw=0.8)
    axes[1].axvline(crossing, color="tab:red", ls="--", lw=0.9)
    axes[1].set_xlabel("absolute step")
    axes[1].set_ylabel("existing rank_Q")
    axes[1].set_yticks([1, 2, 3, 4])
    axes[1].set_title("Figure 6b: existing rank_Q at relative threshold 1e-8")
    fig.tight_layout()
    save(fig, FIGURE_NAMES[5])

    # Figure 7
    fig, ax = plt.subplots(figsize=(13, 6))
    for column, label in [
        ("q3_over_q1", "q3/q1"),
        ("q4_over_q1", "q4/q1"),
        ("min_q3_q4_over_q1", "min(q3,q4)/q1"),
    ]:
        values = np.array([as_float(row[column]) for row in rank_metrics])
        ax.plot(q_step, values, lw=0.8, label=label)
    ax.axvline(crossing, color="tab:red", ls="--", lw=0.9)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("ratio")
    ax.set_title("Figure 7: q-ratios without event threshold")
    ax.legend()
    fig.tight_layout()
    save(fig, FIGURE_NAMES[6])

    # Figure 8
    heat = np.full((len(rank_thresholds), len(rank_persistences)), np.nan)
    for row in rank_onsets:
        if row["status"] != "found":
            continue
        i = rank_thresholds.index(float(row["relative_threshold"]))
        j = rank_persistences.index(int(row["persistence_records"]))
        heat[i, j] = int(row["rank4_onset_candidate"])
    fig, ax = plt.subplots(figsize=(11, 7))
    image = ax.imshow(heat, aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_xticks(range(len(rank_persistences)), [str(value) for value in rank_persistences])
    ax.set_yticks(
        range(len(rank_thresholds)),
        [f"{value:.0e}" for value in rank_thresholds],
    )
    ax.set_xlabel("persistence in consecutive saved q records")
    ax.set_ylabel("relative rank threshold")
    ax.set_title("Figure 8: first observed rank=4 candidate; no adopted cell")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            text = "NA" if np.isnan(heat[i, j]) else str(int(heat[i, j]))
            ax.text(j, i, text, ha="center", va="center", fontsize=7, color="white")
    fig.colorbar(image, ax=ax, label="onset candidate absolute step")
    fig.tight_layout()
    save(fig, FIGURE_NAMES[7])

    # Figure 9
    onset_steps = np.array(
        [
            int(row["rank4_onset_candidate"])
            for row in rank_onsets
            if row["status"] == "found"
        ],
        dtype=np.int64,
    )
    differences_by_condition = {}
    global_min = None
    global_max = None
    for condition_name in config["growth_end_conditions"]:
        end_steps = np.array(
            [
                int(row["growth_end_candidate"])
                for row in found_ends
                if row["end_condition"] == condition_name
            ],
            dtype=np.int64,
        )
        values = (onset_steps[None, :] - end_steps[:, None]).ravel()
        differences_by_condition[condition_name] = values
        if len(values):
            current_min = int(np.min(values))
            current_max = int(np.max(values))
            global_min = current_min if global_min is None else min(global_min, current_min)
            global_max = current_max if global_max is None else max(global_max, current_max)
    bins = np.linspace(global_min, global_max, 121) if global_min != global_max else 20
    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
    for ax, (condition_name, values) in zip(axes, differences_by_condition.items()):
        ax.hist(values, bins=bins, color="tab:blue", alpha=0.75)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_ylabel("all-pair count")
        ax.set_title(f"end condition {condition_name}; pairs={len(values)}")
    axes[-1].set_xlabel("rank4 onset candidate - growth-end candidate (steps)")
    fig.suptitle(
        "Figure 9: all candidate-pair time differences; no pairing selected"
    )
    fig.tight_layout()
    save(fig, FIGURE_NAMES[8])

    # Figure 10
    zoom_low = crossing + int(config["crossing_zoom_relative_range"][0])
    zoom_high = crossing + int(config["crossing_zoom_relative_range"][1])
    mask = (step >= zoom_low) & (step <= zoom_high)
    paper_step = np.array([int(row["step"]) for row in paper_rows])
    paper_f = np.array([float(row["splitting_fraction"]) for row in paper_rows])
    paper_mask = (paper_step >= zoom_low) & (paper_step <= zoom_high)
    fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    axes[0].plot(step[mask] - crossing, log10_f[mask], color="black", lw=0.9, label="fcurve log10(f)")
    axes[0].scatter(
        paper_step[paper_mask] - crossing,
        np.log10(np.where(paper_f[paper_mask] > 0, paper_f[paper_mask], np.nan)),
        s=8,
        color="tab:orange",
        alpha=0.5,
        label="paper7_long saved splitting_fraction",
    )
    axes[0].axvline(0, color="tab:red", ls="--", lw=1.0)
    axes[0].set_ylabel("log10(f)")
    axes[0].legend()
    axes[0].set_title("Figure 10a: existing crossing=1167 reference zoom")
    for color, window in zip(colors, windows):
        axes[1].plot(
            step[mask] - crossing,
            slopes[window][mask],
            color=color,
            lw=0.8,
            label=f"window={window}",
        )
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].axvline(0, color="tab:red", ls="--", lw=1.0)
    axes[1].set_xlabel("step - existing crossing 1167")
    axes[1].set_ylabel("slope of log(f)")
    axes[1].legend(ncol=3)
    axes[1].set_title("Figure 10b: every window, no adopted growth event")
    fig.tight_layout()
    save(fig, FIGURE_NAMES[9])

    generated = sorted(
        str(path)
        for path in FIGURE_DIR.iterdir()
        if path.is_file() and path.suffix in (".png", ".svg")
    )
    manifest = {
        "stage": "A1",
        "success": len(generated) == 20,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "figure_stems": FIGURE_NAMES,
        "generated_files": generated,
        "figure_count_by_format": {
            "png": len([path for path in generated if path.endswith(".png")]),
            "svg": len([path for path in generated if path.endswith(".svg")]),
        },
        "automatic_event_selection": False,
    }
    FIGURE_MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        f"png_count={manifest['figure_count_by_format']['png']}",
        f"svg_count={manifest['figure_count_by_format']['svg']}",
        f"duration_seconds={manifest['duration_seconds']:.6f}",
        "SUCCESS" if manifest["success"] else "STOP: figure count incomplete",
    ]
    TEXT_LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if manifest["success"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
