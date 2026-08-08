#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""既存の個別CSVから、準安定窓に限定した要約JSONを再生成する。"""

import csv
import glob
import json
import math
import os

from run_relation_amplitude_scaling_v1 import RESULT_DIR, tail_stats


TAIL_KEYS = (
    "f_initial_plane",
    "pr_over_m",
    "entropy_normalized",
    "relation_abs_mean",
    "relation_abs_median",
    "relation_abs_q05",
    "relation_abs_q95",
    "scaled_mean_sqrt_m",
    "scaled_median_sqrt_m",
    "scaled_width_q90_sqrt_m",
    "sigma_ratio_2_1",
)


def refresh(path):
    with open(path, encoding="utf-8") as fh:
        summary = json.load(fh)
    csv_path = os.path.join(
        os.path.dirname(path), summary["trajectory_csv"]
    )
    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    if summary["plateau_tau"] is not None:
        selected = [
            row for row in rows
            if int(row["tau"]) >= int(summary["plateau_tau"])
        ]
        window = "equalization_streak"
    else:
        count = min(int(summary.get("tail_records", 200)), len(rows))
        selected = rows[-count:]
        window = "terminal_window_without_equalization_detection"
    if not selected:
        raise RuntimeError(f"再集計窓が空です: {csv_path}")

    count = len(selected)
    tail = {
        key: tail_stats(selected, key, count)
        for key in TAIL_KEYS
    }
    measured = tail["relation_abs_median"]["mean"]
    prediction = 1.0 / math.sqrt(summary["m"])
    summary["tail_window"] = window
    summary["tail_records"] = count
    summary["tail"] = tail
    summary["duality"] = {
        "measured_relation_abs_median": measured,
        "prediction_1_over_sqrt_m": prediction,
        "measured_times_sqrt_m": measured * math.sqrt(summary["m"]),
        "relative_error": measured / prediction - 1.0,
        "n_times_measured": summary["n"] * measured,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    return summary["n"], summary["seed"], window, count


def main():
    paths = sorted(glob.glob(os.path.join(RESULT_DIR, "summary_N*.json")))
    if not paths:
        raise RuntimeError("再集計対象の要約JSONがありません")
    for path in paths:
        n, seed, window, count = refresh(path)
        print(
            f"N={n} seed={seed}: {window}, "
            f"{count} records"
        )


if __name__ == "__main__":
    main()

