#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""個別に生成した N=3〜7 のCSV/JSONを集約し、比較図と集約CSVを作る。"""

import argparse
import csv
import glob
import json
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RESULT_DIR = os.path.join(BASE_DIR, "lowN_metastable_result_v1")


def load_summaries(result_dir, delta, n_min, n_max):
    summaries = []
    for path in glob.glob(os.path.join(result_dir, "summary_N*.json")):
        with open(path, encoding="utf-8") as fh:
            summary = json.load(fh)
        if (
            n_min <= summary["n"] <= n_max
            and np.isclose(summary["delta"], delta, rtol=0.0, atol=abs(delta) * 1e-12)
        ):
            summary["_summary_path"] = path
            summaries.append(summary)
    summaries.sort(key=lambda item: (item["n"], item["seed"]))
    return summaries


def load_trajectory(summary, result_dir):
    path = os.path.join(result_dir, summary["trajectory_csv"])
    rows = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    return rows


def tail_value(summary, observable, statistic="mean"):
    value = summary["tail"][observable][statistic]
    return float(value) if value is not None else np.nan


def grouped_values(summaries, observable, statistic="mean"):
    grouped = defaultdict(list)
    for summary in summaries:
        grouped[summary["n"]].append(
            tail_value(summary, observable, statistic)
        )
    return grouped


def median_band(grouped):
    ns = np.array(sorted(grouped), dtype=int)
    median = np.array([np.median(grouped[n]) for n in ns])
    low = np.array([np.quantile(grouped[n], 0.16) for n in ns])
    high = np.array([np.quantile(grouped[n], 0.84) for n in ns])
    return ns, median, low, high


def write_aggregate_csv(summaries, path):
    fields = [
        "n",
        "m",
        "delta",
        "seed",
        "parent_residual",
        "crossing_tau",
        "a_complement_mean",
        "a_complement_oscillation",
        "a_plane2_mean",
        "a_plane2_oscillation",
        "epsilon_half_mean",
        "epsilon_half_oscillation",
        "pr_mean",
        "relation_abs_median_mean",
        "equal_amplitude_prediction_1_over_sqrt_m",
        "sqrt_m_times_relation_abs_median",
        "n_times_relation_abs_median",
        "n_minus_1_times_a_complement",
        "sqrt_n_minus_1_times_a_complement",
        "n_times_a_complement",
        "max_closure_deviation",
        "max_norm_deviation",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            n = summary["n"]
            amplitude = tail_value(summary, "a_complement")
            writer.writerow({
                "n": n,
                "m": summary["m"],
                "delta": summary["delta"],
                "seed": summary["seed"],
                "parent_residual": summary["parent_residual"],
                "crossing_tau": summary["crossing_tau"],
                "a_complement_mean": amplitude,
                "a_complement_oscillation": tail_value(
                    summary, "a_complement", "oscillation_half_q90"
                ),
                "a_plane2_mean": tail_value(summary, "a_plane2"),
                "a_plane2_oscillation": tail_value(
                    summary, "a_plane2", "oscillation_half_q90"
                ),
                "epsilon_half_mean": tail_value(summary, "epsilon_half"),
                "epsilon_half_oscillation": tail_value(
                    summary, "epsilon_half", "oscillation_half_q90"
                ),
                "pr_mean": tail_value(summary, "pr"),
                "relation_abs_median_mean": tail_value(
                    summary, "relation_abs_median"
                ),
                "equal_amplitude_prediction_1_over_sqrt_m":
                    1.0 / np.sqrt(summary["m"]),
                "sqrt_m_times_relation_abs_median":
                    np.sqrt(summary["m"]) * tail_value(
                        summary, "relation_abs_median"
                    ),
                "n_times_relation_abs_median":
                    n * tail_value(summary, "relation_abs_median"),
                "n_minus_1_times_a_complement": (n - 1) * amplitude,
                "sqrt_n_minus_1_times_a_complement": np.sqrt(n - 1) * amplitude,
                "n_times_a_complement": n * amplitude,
                "max_closure_deviation": summary["max_closure_deviation"],
                "max_norm_deviation": summary["max_norm_deviation"],
            })


def representative_by_n(summaries):
    representatives = {}
    for summary in summaries:
        n = summary["n"]
        if n not in representatives or summary["seed"] < representatives[n]["seed"]:
            representatives[n] = summary
    return representatives


def plot_results(summaries, result_dir, delta, output):
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, 5))
    color_for_n = {
        n: colors[i]
        for i, n in enumerate(range(3, 8))
    }
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))

    ax = axes[0, 0]
    for n, summary in representative_by_n(summaries).items():
        rows = load_trajectory(summary, result_dir)
        tau = np.array([float(row["tau"]) for row in rows])
        amplitude = np.array([float(row["a_complement"]) for row in rows])
        ax.plot(tau, amplitude, lw=1.1, color=color_for_n.get(n),
                label=f"N={n}, seed={summary['seed']}")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$A_\perp=\sqrt{f}$")
    ax.set_title("Low-N trajectories: complement amplitude")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    grouped = grouped_values(summaries, "relation_abs_median")
    ns, med, low, high = median_band(grouped)
    ax.errorbar(
        ns,
        med,
        yerr=np.vstack([med - low, high - med]),
        fmt="o-",
        capsize=4,
        label=r"measured median relation amplitude",
    )
    theory_n = np.arange(3, 8)
    theory_m = theory_n * (theory_n - 1) / 2
    ax.plot(
        theory_n,
        1.0 / np.sqrt(theory_m),
        "k--",
        label=r"$1/\sqrt{M}=\sqrt{2/[N(N-1)]}$",
    )
    ax.set_xticks(range(3, 8))
    ax.set_xlabel("N")
    ax.set_ylabel(r"typical $|Z_e|$")
    ax.set_title("Per-relation metastable amplitude versus N")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    grouped_complement = grouped_values(summaries, "a_complement")
    ns_p, med_p, low_p, high_p = median_band(grouped_complement)
    ax.errorbar(
        ns_p,
        med_p,
        yerr=np.vstack([med_p - low_p, high_p - med_p]),
        fmt="o-",
        capsize=4,
        label=r"tail mean of $A_\perp=\sqrt{f}$",
    )
    grouped_osc = grouped_values(
        summaries, "a_complement", "oscillation_half_q90"
    )
    ns_o, med_o, _, _ = median_band(grouped_osc)
    ax.plot(ns_o, med_o, "s--", label="tail oscillation half-width")
    ax.set_xticks(range(3, 8))
    ax.set_xlabel("N")
    ax.set_ylabel("global complement amplitude")
    ax.set_title("Initial-plane complement is a different observable")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    grouped_eps = grouped_values(summaries, "epsilon_half")
    ns_e, med_e, low_e, high_e = median_band(grouped_eps)
    ax.errorbar(
        ns_e,
        med_e,
        yerr=np.vstack([med_e - low_e, high_e - med_e]),
        fmt="o-",
        capsize=4,
        label=r"measured $\varepsilon=1/2-\sigma_2/\sigma_1$",
    )
    theory_n = np.arange(4, 8)
    ax.plot(
        theory_n,
        1.0 / (2.0 * (theory_n - 1)),
        "k--",
        label=r"$1/[2(N-1)]$ reference",
    )
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xticks(range(3, 8))
    ax.set_xlabel("N")
    ax.set_ylabel(r"$\varepsilon$")
    ax.set_title("Spectral deviation kept separate from amplitude")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle(
        rf"Low-$N$ metastable scan, $\delta={delta:.0e}$",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", default=DEFAULT_RESULT_DIR)
    parser.add_argument("--delta", type=float, default=1e-15)
    parser.add_argument("--n-min", type=int, default=3)
    parser.add_argument("--n-max", type=int, default=7)
    args = parser.parse_args()

    summaries = load_summaries(
        args.result_dir, args.delta, args.n_min, args.n_max
    )
    if not summaries:
        raise RuntimeError("集約対象のsummary JSONがありません")
    missing = sorted(
        set(range(args.n_min, args.n_max + 1))
        - {summary["n"] for summary in summaries}
    )
    if missing:
        raise RuntimeError(f"N={missing} の個別結果が不足しています")

    tag = f"N{args.n_min:02d}-{args.n_max:02d}_delta{args.delta:.0e}"
    csv_path = os.path.join(
        args.result_dir, f"aggregate_{tag}.csv"
    )
    json_path = os.path.join(
        args.result_dir, f"aggregate_{tag}.json"
    )
    figure_path = os.path.join(
        args.result_dir, f"metastable_scaling_{tag}.png"
    )
    write_aggregate_csv(summaries, csv_path)
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "delta": args.delta,
                "n_min": args.n_min,
                "n_max": args.n_max,
                "n_runs": len(summaries),
                "runs_per_n": {
                    str(n): sum(item["n"] == n for item in summaries)
                    for n in range(args.n_min, args.n_max + 1)
                },
                "aggregate_csv": os.path.basename(csv_path),
                "figure": os.path.basename(figure_path),
            },
            fh,
            indent=2,
            ensure_ascii=False,
        )
        fh.write("\n")
    plot_results(
        summaries, args.result_dir, args.delta, figure_path
    )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {figure_path}")


if __name__ == "__main__":
    main()
