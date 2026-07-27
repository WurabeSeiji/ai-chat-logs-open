#!/usr/bin/env python3
"""Create Stage A1b transition-anatomy figures from processed Stage A0 data."""

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
PROCESSED_DIR = PACKAGE_ROOT / "processed"
FIGURE_DIR = PACKAGE_ROOT / "figures"
LOG_DIR = PACKAGE_ROOT / "logs"
ANALYSIS_MANIFEST_PATH = LOG_DIR / "analysis_manifest.json"
FIGURE_MANIFEST_PATH = LOG_DIR / "figures_manifest.json"
TEXT_LOG_PATH = LOG_DIR / "make_transition_figures.log"
CROSSING = 1167
REGRESSION_WINDOWS = [11, 21, 41, 81, 161]
ZOOM_WINDOWS = [
    ("0-500", 0, 500),
    ("500-1000", 500, 1000),
    ("800-1400", 800, 1400),
    ("1000-1800", 1000, 1800),
    ("1400-2500", 1400, 2500),
]
FIGURE_STEMS = [
    "figure01_f_running_max_0_3000",
    "figure02_log10_f_and_slopes_0_3000",
    "figure03_direction_1_to_4_occupation_0_3000",
    "figure04_other_kernel_splitting_0_3000",
    "figure05_q1_to_q4_0_3000",
    "figure06_q_ratios_and_rank_q_0_3000",
    "figure07_first_passage_vs_direction_3_4_actual_records",
    "figure08_first_passage_vs_q3_q4_actual_records",
    "figure09_all_quantities_800_1400",
    "figure10_all_quantities_1000_1800",
    "figure11_first_passage_level_step_intervals",
    "figure12_adjacent_level_mean_exponential_rates",
    "supplement_all_requested_zoom_ranges",
]

os.environ.setdefault("MPLCONFIGDIR", str(LOG_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(LOG_DIR / "cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def add_crossing(ax, start: int = 0, end: int = 3000) -> None:
    if start <= CROSSING <= end:
        ax.axvline(
            CROSSING,
            color="tab:red",
            ls="--",
            lw=0.9,
            label="existing crossing=1167",
        )


def actual_mask(steps: np.ndarray, start: int, end: int) -> np.ndarray:
    return (steps >= start) & (steps <= end)


def draw_all_quantities(
    f_step: np.ndarray,
    log10_f: np.ndarray,
    occ_step: np.ndarray,
    occ_actual: dict[str, np.ndarray],
    display_step: np.ndarray,
    occ_display: dict[str, np.ndarray],
    q_step: np.ndarray,
    q_actual: dict[str, np.ndarray],
    start: int,
    end: int,
    title: str,
):
    fig, axes = plt.subplots(4, 1, figsize=(14, 13), sharex=True)
    fmask = actual_mask(f_step, start, end)
    omask = actual_mask(occ_step, start, end)
    dmask = actual_mask(display_step, start, end)
    qmask = actual_mask(q_step, start, end)

    axes[0].plot(f_step[fmask], log10_f[fmask], color="black", lw=0.9)
    axes[0].set_ylabel("log10(f)")
    axes[0].set_title("f normalization: logarithm only; no min-max scaling")

    for column, label in [
        ("direction_1_occupation", "direction 1"),
        ("direction_2_occupation", "direction 2"),
        ("direction_3_occupation", "direction 3"),
        ("direction_4_occupation", "direction 4"),
    ]:
        axes[1].plot(
            display_step[dmask],
            occ_display[column][dmask],
            lw=0.8,
            label=f"{label} display interpolation",
        )
        axes[1].scatter(
            occ_step[omask],
            occ_actual[column][omask],
            s=8,
            alpha=0.55,
        )
    axes[1].set_ylabel("raw occupation fraction")
    axes[1].set_title(
        "directions: raw fractions; lines are display-only linear interpolation"
    )
    axes[1].legend(ncol=4, fontsize=7)

    for column, label in [
        ("other_rotating_occupation", "other rotating"),
        ("kernel_occupation", "kernel"),
        ("splitting_fraction", "splitting f"),
    ]:
        axes[2].plot(
            display_step[dmask],
            occ_display[column][dmask],
            lw=0.8,
            label=f"{label} display interpolation",
        )
        axes[2].scatter(
            occ_step[omask],
            occ_actual[column][omask],
            s=8,
            alpha=0.55,
        )
    axes[2].set_ylabel("raw fraction")
    axes[2].set_title("other/kernel/splitting: raw fractions")
    axes[2].legend(ncol=3, fontsize=8)

    q1 = q_actual["q1"][qmask]
    for column, label in [
        ("q1", "q1/q1"),
        ("q2", "q2/q1"),
        ("q3", "q3/q1"),
        ("q4", "q4/q1"),
    ]:
        values = np.divide(
            q_actual[column][qmask],
            q1,
            out=np.full_like(q1, np.nan),
            where=q1 != 0,
        )
        axes[3].plot(
            q_step[qmask],
            values,
            marker=".",
            ms=3,
            lw=0.65,
            label=label,
        )
    rank_axis = axes[3].twinx()
    rank_axis.step(
        q_step[qmask],
        q_actual["rank_q"][qmask],
        where="post",
        color="gray",
        lw=0.7,
        alpha=0.7,
        label="rank_q",
    )
    axes[3].set_ylabel("qj/q1 at actual q records")
    rank_axis.set_ylabel("existing rank_q")
    rank_axis.set_yticks([1, 2, 3, 4])
    axes[3].set_title("q normalization: each qj divided by q1; no q interpolation")
    axes[3].legend(ncol=4, fontsize=8, loc="upper left")
    axes[3].set_xlabel("absolute step")

    for ax in axes:
        add_crossing(ax, start, end)
        ax.set_xlim(start, end)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def main() -> int:
    started = time.perf_counter()
    planned = [
        FIGURE_DIR / f"{stem}.{extension}"
        for stem in FIGURE_STEMS
        for extension in ("png", "svg")
    ] + [FIGURE_MANIFEST_PATH, TEXT_LOG_PATH]
    if any(path.exists() for path in planned):
        raise RuntimeError("Stage A1b図出力の上書きを避けて停止")
    if not ANALYSIS_MANIFEST_PATH.is_file():
        raise RuntimeError("analyze_first_transition.py の成功記録がない")
    analysis_manifest = json.loads(
        ANALYSIS_MANIFEST_PATH.read_text(encoding="utf-8")
    )
    if analysis_manifest.get("success") is not True:
        raise RuntimeError("先行解析が成功していない")

    f_rows = load_csv(PROCESSED_DIR / "transition_f_metrics_0_3000.csv")
    occ_rows = load_csv(PROCESSED_DIR / "occupation_actual_records_0_3000.csv")
    display_rows = load_csv(
        PROCESSED_DIR
        / "occupation_display_only_linear_interpolation_0_3000.csv"
    )
    q_rows = load_csv(PROCESSED_DIR / "q_actual_records_0_3000.csv")
    first_passage = load_csv(PROCESSED_DIR / "f_first_passage_levels.csv")
    rates = load_csv(PROCESSED_DIR / "f_decade_growth_rates.csv")
    occ_nearest = load_csv(
        PROCESSED_DIR / "first_passage_nearest_occupation_records.csv"
    )
    q_nearest = load_csv(
        PROCESSED_DIR / "first_passage_nearest_q_records.csv"
    )

    f_step = np.array([int(row["step"]) for row in f_rows])
    f_value = np.array([as_float(row["f"]) for row in f_rows])
    log10_f = np.array([as_float(row["log10_f"]) for row in f_rows])
    running_max = np.array([as_float(row["running_max_f"]) for row in f_rows])
    slopes = {
        window: np.array(
            [as_float(row[f"slope_w{window:03d}"]) for row in f_rows]
        )
        for window in REGRESSION_WINDOWS
    }

    occ_step = np.array([int(row["step"]) for row in occ_rows])
    occ_columns = [
        "direction_1_occupation",
        "direction_2_occupation",
        "direction_3_occupation",
        "direction_4_occupation",
        "other_rotating_occupation",
        "kernel_occupation",
        "splitting_fraction",
    ]
    occ_actual = {
        column: np.array([as_float(row[column]) for row in occ_rows])
        for column in occ_columns
    }
    display_step = np.array([int(row["step"]) for row in display_rows])
    occ_display = {
        column: np.array(
            [
                as_float(row[f"{column}_display_interp"])
                for row in display_rows
            ]
        )
        for column in occ_columns
    }

    q_step = np.array([int(row["step"]) for row in q_rows])
    q_columns = [
        "q1",
        "q2",
        "q3",
        "q4",
        "rank_q",
        "q3_over_q1",
        "q4_over_q1",
        "min_q3_q4_over_q1",
    ]
    q_actual = {
        column: np.array([as_float(row[column]) for row in q_rows])
        for column in q_columns
    }

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    (LOG_DIR / "matplotlib").mkdir(parents=True, exist_ok=True)
    (LOG_DIR / "cache").mkdir(parents=True, exist_ok=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(REGRESSION_WINDOWS)))

    # Figure 1
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.semilogy(f_step, f_value, color="black", lw=0.75, label="f actual every step")
    ax.semilogy(
        f_step,
        running_max,
        color="tab:orange",
        lw=1.1,
        label="running maximum f",
    )
    add_crossing(ax)
    ax.set_xlim(0, 3000)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("f (log scale)")
    ax.set_title("Figure 1: f and running maximum, first-transition observation window")
    ax.legend()
    fig.tight_layout()
    save(fig, FIGURE_STEMS[0])

    # Figure 2
    fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    axes[0].plot(f_step, log10_f, color="black", lw=0.8)
    axes[0].set_ylabel("log10(f)")
    axes[0].set_title("Figure 2a: log10(f)")
    for color, window in zip(colors, REGRESSION_WINDOWS):
        axes[1].plot(
            f_step,
            slopes[window],
            color=color,
            lw=0.75,
            label=f"window={window}",
        )
    axes[1].axhline(0, color="gray", lw=0.6)
    axes[1].set_ylabel("OLS slope of ln(f)")
    axes[1].set_xlabel("absolute step")
    axes[1].set_title("Figure 2b: all requested regression slopes; no interval selected")
    axes[1].legend(ncol=3)
    for ax in axes:
        add_crossing(ax)
        ax.set_xlim(0, 3000)
    fig.tight_layout()
    save(fig, FIGURE_STEMS[1])

    # Figure 3
    fig, ax = plt.subplots(figsize=(13, 6))
    for column in [
        "direction_1_occupation",
        "direction_2_occupation",
        "direction_3_occupation",
        "direction_4_occupation",
    ]:
        label = column.replace("_occupation", "").replace("_", " ")
        ax.plot(
            display_step,
            occ_display[column],
            lw=0.85,
            label=f"{label} display interp",
        )
        ax.scatter(
            occ_step,
            occ_actual[column],
            s=8,
            alpha=0.5,
        )
    add_crossing(ax)
    ax.set_xlim(0, 3000)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("occupation fraction")
    ax.set_title(
        "Figure 3: directions 1-4; markers=actual 25-step records, lines=display-only interpolation"
    )
    ax.legend(ncol=2)
    fig.tight_layout()
    save(fig, FIGURE_STEMS[2])

    # Figure 4
    fig, ax = plt.subplots(figsize=(13, 6))
    for column in [
        "other_rotating_occupation",
        "kernel_occupation",
        "splitting_fraction",
    ]:
        label = column.replace("_occupation", "").replace("_", " ")
        ax.plot(
            display_step,
            occ_display[column],
            lw=0.9,
            label=f"{label} display interp",
        )
        ax.scatter(occ_step, occ_actual[column], s=8, alpha=0.5)
    add_crossing(ax)
    ax.set_xlim(0, 3000)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("raw fraction")
    ax.set_title(
        "Figure 4: other rotating, kernel, splitting; actual markers and display-only interpolation"
    )
    ax.legend()
    fig.tight_layout()
    save(fig, FIGURE_STEMS[3])

    # Figure 5
    fig, ax = plt.subplots(figsize=(13, 6))
    for column in ("q1", "q2", "q3", "q4"):
        ax.plot(
            q_step,
            q_actual[column],
            marker=".",
            ms=3,
            lw=0.7,
            label=f"{column} actual saved",
        )
    add_crossing(ax)
    ax.set_xlim(0, 3000)
    ax.set_xlabel("absolute step")
    ax.set_ylabel("q")
    ax.set_title("Figure 5: q1-q4 at actual saved q records; no q interpolation")
    ax.legend(ncol=4)
    fig.tight_layout()
    save(fig, FIGURE_STEMS[4])

    # Figure 6
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for column, label in [
        ("q3_over_q1", "q3/q1"),
        ("q4_over_q1", "q4/q1"),
        ("min_q3_q4_over_q1", "min(q3,q4)/q1"),
    ]:
        axes[0].plot(
            q_step,
            q_actual[column],
            marker=".",
            ms=3,
            lw=0.7,
            label=label,
        )
    axes[0].set_ylabel("ratio")
    axes[0].set_title("Figure 6a: q ratios at actual saved records")
    axes[0].legend()
    axes[1].step(
        q_step,
        q_actual["rank_q"],
        where="post",
        color="black",
        lw=0.8,
    )
    axes[1].set_ylabel("existing rank_q")
    axes[1].set_yticks([1, 2, 3, 4])
    axes[1].set_xlabel("absolute step")
    axes[1].set_title("Figure 6b: existing rank_q; not equated with q growth or occupation")
    for ax in axes:
        add_crossing(ax)
        ax.set_xlim(0, 3000)
    fig.tight_layout()
    save(fig, FIGURE_STEMS[5])

    # Figure 7
    found_occ = [row for row in occ_nearest if row["status"] == "found"]
    passage_steps = np.array([int(row["first_passage_step"]) for row in found_occ])
    fig, ax = plt.subplots(figsize=(13, 7))
    series = [
        ("before_direction_3_occupation", "d3 before/or-at", "o"),
        ("after_direction_3_occupation", "d3 after/or-at", "x"),
        ("before_direction_4_occupation", "d4 before/or-at", "s"),
        ("after_direction_4_occupation", "d4 after/or-at", "+"),
    ]
    for column, label, marker in series:
        values = np.array([as_float(row[column]) for row in found_occ])
        ax.scatter(passage_steps, values, s=24, marker=marker, label=label)
    add_crossing(ax)
    ax.set_yscale("log")
    ax.set_xlim(0, 1300)
    ax.set_xlabel("f first-passage absolute step")
    ax.set_ylabel("actual saved direction occupation (log)")
    ax.set_title(
        "Figure 7: first-passage coordinates versus bracketing actual direction 3/4 records"
    )
    ax.legend(ncol=2)
    fig.tight_layout()
    save(fig, FIGURE_STEMS[6])

    # Figure 8
    found_q = [row for row in q_nearest if row["status"] == "found"]
    passage_steps_q = np.array([int(row["first_passage_step"]) for row in found_q])
    fig, ax = plt.subplots(figsize=(13, 7))
    for column, label, marker in [
        ("before_q3", "q3 before/or-at", "o"),
        ("after_q3", "q3 after/or-at", "x"),
        ("before_q4", "q4 before/or-at", "s"),
        ("after_q4", "q4 after/or-at", "+"),
    ]:
        values = np.array([as_float(row[column]) for row in found_q])
        ax.scatter(passage_steps_q, values, s=24, marker=marker, label=label)
    add_crossing(ax)
    ax.set_yscale("symlog", linthresh=1e-16)
    ax.set_xlim(0, 1300)
    ax.set_xlabel("f first-passage absolute step")
    ax.set_ylabel("actual saved q value (symlog; linear within +/-1e-16)")
    ax.set_title(
        "Figure 8: first-passage coordinates versus bracketing actual q3/q4 records"
    )
    ax.legend(ncol=2)
    fig.tight_layout()
    save(fig, FIGURE_STEMS[7])

    # Figures 9 and 10
    fig = draw_all_quantities(
        f_step,
        log10_f,
        occ_step,
        occ_actual,
        display_step,
        occ_display,
        q_step,
        q_actual,
        800,
        1400,
        "Figure 9: all existing quantities, 800-1400; axes and normalizations separated",
    )
    save(fig, FIGURE_STEMS[8])
    fig = draw_all_quantities(
        f_step,
        log10_f,
        occ_step,
        occ_actual,
        display_step,
        occ_display,
        q_step,
        q_actual,
        1000,
        1800,
        "Figure 10: all existing quantities, 1000-1800; axes and normalizations separated",
    )
    save(fig, FIGURE_STEMS[9])

    # Figure 11
    found_rates = [
        row
        for row in rates
        if row["status"] in ("found", "zero_step_difference")
    ]
    upper_levels = np.array([as_float(row["upper_level"]) for row in found_rates])
    step_differences = np.array(
        [as_float(row["step_difference"]) for row in found_rates]
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.semilogx(
        upper_levels,
        step_differences,
        marker="o",
        ms=4,
        lw=0.7,
    )
    ax.set_xlabel("upper first-passage level")
    ax.set_ylabel("step difference from adjacent lower level")
    ax.set_title(
        "Figure 11: step intervals between adjacent listed first-passage levels"
    )
    fig.tight_layout()
    save(fig, FIGURE_STEMS[10])

    # Figure 12
    rate_rows = [row for row in rates if row["status"] == "found"]
    geometric_levels = np.sqrt(
        np.array([as_float(row["lower_level"]) for row in rate_rows])
        * np.array([as_float(row["upper_level"]) for row in rate_rows])
    )
    mean_rates = np.array(
        [as_float(row["mean_exponential_rate_per_step"]) for row in rate_rows]
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.semilogx(
        geometric_levels,
        mean_rates,
        marker="o",
        ms=4,
        lw=0.7,
    )
    ax.set_xlabel("geometric mean of adjacent listed levels")
    ax.set_ylabel("ln(level ratio) / step difference")
    ax.set_title(
        "Figure 12: adjacent-level mean exponential rates; no interval adopted"
    )
    fig.tight_layout()
    save(fig, FIGURE_STEMS[11])

    # Supplement: all explicitly requested zoom ranges
    fig, axes = plt.subplots(len(ZOOM_WINDOWS), 3, figsize=(18, 18))
    for row_index, (label, start, end) in enumerate(ZOOM_WINDOWS):
        fmask = actual_mask(f_step, start, end)
        omask = actual_mask(occ_step, start, end)
        dmask = actual_mask(display_step, start, end)
        qmask = actual_mask(q_step, start, end)
        axes[row_index, 0].plot(
            f_step[fmask],
            log10_f[fmask],
            color="black",
            lw=0.8,
        )
        axes[row_index, 0].set_ylabel(f"{label}\nlog10(f)")
        for column in ("direction_3_occupation", "direction_4_occupation"):
            axes[row_index, 1].plot(
                display_step[dmask],
                occ_display[column][dmask],
                lw=0.75,
                label=column.replace("_occupation", ""),
            )
            axes[row_index, 1].scatter(
                occ_step[omask],
                occ_actual[column][omask],
                s=5,
                alpha=0.45,
            )
        axes[row_index, 1].set_ylabel("raw occupation")
        for column, qlabel in [
            ("q3_over_q1", "q3/q1"),
            ("q4_over_q1", "q4/q1"),
        ]:
            axes[row_index, 2].plot(
                q_step[qmask],
                q_actual[column][qmask],
                marker=".",
                ms=2,
                lw=0.6,
                label=qlabel,
            )
        axes[row_index, 2].set_ylabel("q ratio actual")
        for ax in axes[row_index]:
            ax.set_xlim(start, end)
            add_crossing(ax, start, end)
        if row_index == 0:
            axes[row_index, 1].legend(fontsize=7)
            axes[row_index, 2].legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("absolute step")
    fig.suptitle(
        "Supplement: all requested zooms; f=log10, occupation=raw with display interpolation, q=qj/q1 actual only"
    )
    fig.tight_layout()
    save(fig, FIGURE_STEMS[12])

    generated = sorted(
        str(path.resolve())
        for path in FIGURE_DIR.iterdir()
        if path.suffix in (".png", ".svg")
    )
    manifest = {
        "stage": "A1b",
        "success": len(generated) == len(FIGURE_STEMS) * 2,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.perf_counter() - started,
        "figure_stems": FIGURE_STEMS,
        "generated_files": generated,
        "figure_count_by_format": {
            "png": sum(path.endswith(".png") for path in generated),
            "svg": sum(path.endswith(".svg") for path in generated),
        },
        "existing_event_reference_lines": ["crossing=1167"],
        "other_single_event_lines": [],
        "q_interpolation": False,
        "occupation_interpolation_use": "display_only",
    }
    FIGURE_MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        f"png_count={manifest['figure_count_by_format']['png']}",
        f"svg_count={manifest['figure_count_by_format']['svg']}",
        f"duration_seconds={manifest['duration_seconds']:.6f}",
        "SUCCESS" if manifest["success"] else "STOP: figure output incomplete",
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
