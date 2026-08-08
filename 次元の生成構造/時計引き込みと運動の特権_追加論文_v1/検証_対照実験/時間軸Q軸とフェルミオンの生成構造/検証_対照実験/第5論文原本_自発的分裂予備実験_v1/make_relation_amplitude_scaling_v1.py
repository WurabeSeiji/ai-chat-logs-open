#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""低N・大Nの個別結果を統合し、N=3〜300の振幅則を定量評価する。"""

import argparse
import csv
import glob
import json
import math
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOW_DIR = os.path.join(BASE_DIR, "lowN_metastable_result_v1")
LARGE_DIR = os.path.join(BASE_DIR, "relation_amplitude_scaling_result_v1")
OUTPUT_DIR = LARGE_DIR


def close_delta(value, target):
    return np.isclose(value, target, rtol=0.0, atol=abs(target) * 1e-12)


def load_large_measurement_window(summary, summary_path):
    """既存CSVから、検出済み準安定窓だけを読み直す。"""
    csv_path = os.path.join(
        os.path.dirname(summary_path), summary["trajectory_csv"]
    )
    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if summary["plateau_tau"] is not None:
        selected = [
            row for row in rows
            if int(row["tau"]) >= int(summary["plateau_tau"])
        ]
        window = "detected_equalization_streak"
    else:
        selected = rows[-min(int(summary["tail_records"]), len(rows)):]
        window = "terminal_window_without_equalization_detection"
    if not selected:
        raise RuntimeError(f"準安定窓が空です: {csv_path}")

    def mean(key):
        return float(np.mean([float(row[key]) for row in selected]))

    return {
        "measured": mean("relation_abs_median"),
        "q05": mean("relation_abs_q05"),
        "q95": mean("relation_abs_q95"),
        "pr_over_m": mean("pr_over_m"),
        "scaled_width": mean("scaled_width_q90_sqrt_m"),
        "window": window,
        "window_records": len(selected),
    }


def load_low_measurement_window(summary, summary_path):
    """低N CSVから、既存要約と同じ末尾窓の分布統計を読む。"""
    csv_path = os.path.join(
        os.path.dirname(summary_path), summary["trajectory_csv"]
    )
    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    start_tau = int(summary["tail"]["start_tau"])
    selected = [row for row in rows if int(row["tau"]) >= start_tau]
    if not selected:
        raise RuntimeError(f"低N準安定窓が空です: {csv_path}")

    def mean(key):
        return float(np.mean([float(row[key]) for row in selected]))

    return {
        "measured": mean("relation_abs_median"),
        "q05": mean("relation_abs_q05"),
        "q95": mean("relation_abs_q95"),
        "pr_over_m": mean("pr") / summary["m"],
        "window_records": len(selected),
    }


def load_records(delta):
    records = []
    for path in glob.glob(os.path.join(LOW_DIR, "summary_N*.json")):
        with open(path, encoding="utf-8") as fh:
            summary = json.load(fh)
        if summary["n"] <= 7 and close_delta(summary["delta"], delta):
            window = load_low_measurement_window(summary, path)
            measured = window["measured"]
            q05 = window["q05"]
            q95 = window["q95"]
            records.append({
                "source": "lowN_exact",
                "n": summary["n"],
                "m": summary["m"],
                "seed": summary["seed"],
                "measured": measured,
                "q05": q05,
                "q95": q95,
                "pr_over_m": window["pr_over_m"],
                "scaled_width": (q95 - q05) * math.sqrt(summary["m"]),
                "measurement_window": "terminal_lowN_window",
                "window_records": window["window_records"],
                "closure_deviation": summary["max_closure_deviation"],
                "summary_path": path,
            })

    for path in glob.glob(os.path.join(LARGE_DIR, "summary_N*.json")):
        with open(path, encoding="utf-8") as fh:
            summary = json.load(fh)
        if summary["n"] >= 8 and close_delta(summary["delta"], delta):
            window = load_large_measurement_window(summary, path)
            records.append({
                "source": summary["normalization"],
                "n": summary["n"],
                "m": summary["m"],
                "seed": summary["seed"],
                "measured": window["measured"],
                "q05": window["q05"],
                "q95": window["q95"],
                "pr_over_m": window["pr_over_m"],
                "scaled_width": window["scaled_width"],
                "measurement_window": window["window"],
                "window_records": window["window_records"],
                "closure_deviation": summary["max_closure_deviation"],
                "summary_path": path,
            })
    for record in records:
        record["prediction"] = 1.0 / math.sqrt(record["m"])
        record["scaled_amplitude"] = (
            record["measured"] * math.sqrt(record["m"])
        )
        record["relative_error"] = (
            record["scaled_amplitude"] - 1.0
        )
        record["n_times_amplitude"] = record["n"] * record["measured"]
    records.sort(key=lambda item: (item["n"], item["seed"]))
    return records


def aggregate_by_n(records):
    grouped = defaultdict(list)
    for record in records:
        grouped[record["n"]].append(record)
    result = []
    for n in sorted(grouped):
        rows = grouped[n]
        result.append({
            "n": n,
            "m": rows[0]["m"],
            "n_runs": len(rows),
            "amplitude_median": float(np.median([r["measured"] for r in rows])),
            "amplitude_q16": float(np.quantile([r["measured"] for r in rows], 0.16)),
            "amplitude_q84": float(np.quantile([r["measured"] for r in rows], 0.84)),
            "prediction": rows[0]["prediction"],
            "scaled_amplitude_median": float(
                np.median([r["scaled_amplitude"] for r in rows])
            ),
            "relative_error_median": float(
                np.median([r["relative_error"] for r in rows])
            ),
            "pr_over_m_median": float(
                np.median([r["pr_over_m"] for r in rows])
            ),
            "scaled_width_median": float(
                np.median([r["scaled_width"] for r in rows])
            ),
            "closure_deviation_max": float(
                max(r["closure_deviation"] for r in rows)
            ),
        })
    return result


def log_fit(rows, minimum_n):
    selected = [row for row in rows if row["n"] >= minimum_n]
    x = np.log([row["n"] for row in selected])
    y = np.log([row["amplitude_median"] for row in selected])
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "minimum_n": minimum_n,
        "n_points": len(selected),
        "alpha_in_A_eq_c_N_minus_alpha": float(-slope),
        "c": float(math.exp(intercept)),
        "r_squared": float(1.0 - ss_res / ss_tot) if ss_tot > 0.0 else 1.0,
    }


def write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(rows, output):
    n = np.array([row["n"] for row in rows], dtype=float)
    amplitude = np.array([row["amplitude_median"] for row in rows])
    low = np.array([row["amplitude_q16"] for row in rows])
    high = np.array([row["amplitude_q84"] for row in rows])
    prediction = np.array([row["prediction"] for row in rows])
    scaled_residual = np.array(
        [row["scaled_amplitude_median"] - 1.0 for row in rows]
    )
    pr = np.array([row["pr_over_m_median"] for row in rows])
    width = np.array([row["scaled_width_median"] for row in rows])

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
    ax = axes[0, 0]
    ax.errorbar(
        n, amplitude,
        yerr=np.vstack([amplitude - low, high - amplitude]),
        fmt="o", capsize=3, label="measured median relation amplitude",
    )
    ax.plot(n, prediction, "k--", label=r"$1/\sqrt{M}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$A_{\rm relation}$")
    ax.set_title("Per-relation metastable amplitude")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.semilogx(n, scaled_residual, "o-")
    ax.axhline(0.0, color="k", ls="--")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$A_{\rm relation}\sqrt{M}-1$")
    ax.set_title("Residual from finite-N amplitude law")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1, 0]
    ax.loglog(n, np.maximum(1.0 - pr, 1e-16), "o-")
    ax.set_xlabel("N")
    ax.set_ylabel("1 - PR/M")
    ax.set_title("Equal-amplitude defect")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1, 1]
    ax.loglog(n, np.maximum(width, 1e-16), "o-")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$(Q_{95}-Q_{05})\sqrt{M}$")
    ax.set_title("Collapse width of relation amplitudes")
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delta", type=float, default=1e-15)
    parser.add_argument("--output-dir", default=OUTPUT_DIR)
    args = parser.parse_args()

    records = load_records(args.delta)
    if not records:
        raise RuntimeError("集約対象がありません")
    rows = aggregate_by_n(records)
    os.makedirs(args.output_dir, exist_ok=True)
    tag = f"N00003-00300_delta{args.delta:.0e}"
    csv_path = os.path.join(args.output_dir, f"scaling_{tag}.csv")
    json_path = os.path.join(args.output_dir, f"scaling_{tag}.json")
    figure_path = os.path.join(args.output_dir, f"scaling_{tag}.png")
    write_csv(rows, csv_path)

    scaled_all = np.array([row["scaled_amplitude"] for row in records])
    report = {
        "delta": args.delta,
        "n_runs": len(records),
        "n_values": [row["n"] for row in rows],
        "finite_n_model": "A_relation = 1/sqrt(M) = sqrt(2/[N(N-1)])",
        "scaled_amplitude": {
            "mean": float(np.mean(scaled_all)),
            "std": float(np.std(scaled_all)),
            "max_abs_deviation_from_1": float(np.max(np.abs(scaled_all - 1.0))),
        },
        "log_fits": [
            log_fit(rows, 8),
            log_fit(rows, 20),
        ],
        "max_closure_deviation": float(
            max(record["closure_deviation"] for record in records)
        ),
        "aggregate_csv": os.path.basename(csv_path),
        "figure": os.path.basename(figure_path),
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    make_figure(rows, figure_path)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"wrote {csv_path}")
    print(f"wrote {figure_path}")


if __name__ == "__main__":
    main()
