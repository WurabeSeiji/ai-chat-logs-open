#!/usr/bin/env python3
"""無seed自然軌道を t=110000 まで測定する長時間専用ラッパ。

論文7の実時間発展コードと初期状態構築を変更せず再利用し、観測は
f, q1..q4 と閉鎖診断に限定する。t=55000 と終端で Z, wp を保存する。
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

import run_seedless_natural_figures3_4_v1 as base

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "outputs" / "long_horizon_110000"
XMAX = 110000
FIT_START = 20000
AMPLIFICATION = 20.0
NS = (40, 300)
SAMPLE = {40: 25, 300: 100}
FMT = "%.10e"
COLUMNS = [
    "step",
    "time",
    "N",
    "condition",
    "initial_seed_enabled",
    "metastable_seed_enabled",
    "f_outside_parent",
    "q1",
    "q2",
    "q3",
    "q4",
    "rank_Q",
    "norm_Z",
    "zero_square_abs",
    "crossing_detected",
]


def csv_path(n: int) -> Path:
    return OUTPUT / "raw" / f"N{n:05d}_seedless_f_q3_q4_t110000.csv"


def checkpoint_path(n: int, step: int) -> Path:
    return OUTPUT / "checkpoints" / f"N{n:05d}_state_t{step:06d}.npz"


def fit_one_exp(u: np.ndarray, values: np.ndarray):
    from scipy.optimize import curve_fit

    def model(x, c, a, tau):
        return c + a * np.exp(-x / tau)

    initial = (values[-1], values[0] - values[-1], 5000.0)
    parameters, _ = curve_fit(
        model,
        u,
        values,
        p0=initial,
        bounds=([-np.inf, -np.inf, 100.0], [np.inf, np.inf, 1e7]),
        maxfev=200000,
    )
    return parameters, model(u, *parameters)


def fit_two_exp(u: np.ndarray, values: np.ndarray):
    from scipy.optimize import curve_fit

    def model(x, c, a1, tau1, a2, tau2):
        return c + a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2)

    span = values[0] - values[-1]
    guesses = (
        (values[-1], 0.8 * span, 3000.0, 0.2 * span, 15000.0),
        (values[-1], 0.5 * span, 1500.0, 0.5 * span, 8000.0),
        (values[-1], 1.2 * span, 4000.0, -0.2 * span, 20000.0),
    )
    best = None
    for initial in guesses:
        try:
            parameters, _ = curve_fit(
                model,
                u,
                values,
                p0=initial,
                bounds=(
                    [-np.inf, -np.inf, 100.0, -np.inf, 100.0],
                    [np.inf, np.inf, 1e7, np.inf, 1e7],
                ),
                maxfev=500000,
            )
        except (RuntimeError, ValueError):
            continue
        trend = model(u, *parameters)
        rss = float(np.sum((values - trend) ** 2))
        if best is None or rss < best[0]:
            best = (rss, parameters, trend)
    if best is None:
        raise RuntimeError("two-exponential fit failed")
    _, parameters, trend = best
    return parameters, trend


def run_long(n: int) -> dict:
    if n not in NS:
        raise SystemExit(f"N must be one of {NS}")
    target = csv_path(n)
    meta_path = OUTPUT / "summary" / f"N{n:05d}_long_horizon_meta.json"
    if target.exists() or meta_path.exists():
        raise SystemExit(f"既存の長時間出力を上書きしない: N={n}")

    base.verify_sources()
    base.prepare_import_paths()
    from run_n300_dimension_saturation_v2 import dominant_plane, gram_reduce
    from run_paper7_5color_timeseries import occ

    (
        sys_lr,
        _,
        B_p1,
        _,
        B0,
        p,
        q,
        Z,
        wp,
        parent_residual,
        parent_sigma,
        method,
    ) = base.build_seedless(n)

    target.parent.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "checkpoints").mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    sample_every = SAMPLE[n]
    crossing = None
    max_norm_error = 0.0
    max_zero_square_abs = 0.0
    next_progress = 10000

    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        for t in range(XMAX + 1):
            Z_perp = Z - p * (p @ Z) - q * (q @ Z)
            totZ = float(np.real(np.conj(Z) @ Z))
            f = float(np.real(np.conj(Z_perp) @ Z_perp)) / totZ
            if crossing is None and f > 0.05:
                crossing = t

            if t % sample_every == 0 or t == XMAX:
                parent_plane_occupation = occ(B_p1, Z)
                f_record = 1.0 - parent_plane_occupation / totZ
                gr = gram_reduce(sys_lr, Z)
                _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
                qs = base.qsv4(B0, Bdom)
                rank_q = int(np.sum(qs > base.Q_REL_TAU * qs[0]))
                zero_square_abs = abs(complex(Z @ Z))
                norm_error = abs(totZ - 1.0)
                writer.writerow(
                    {
                        "step": t,
                        "time": t,
                        "N": n,
                        "condition": "A",
                        "initial_seed_enabled": 0,
                        "metastable_seed_enabled": 0,
                        "f_outside_parent": FMT % f_record,
                        "q1": FMT % qs[0],
                        "q2": FMT % qs[1],
                        "q3": FMT % qs[2],
                        "q4": FMT % qs[3],
                        "rank_Q": rank_q,
                        "norm_Z": FMT % np.sqrt(totZ),
                        "zero_square_abs": FMT % zero_square_abs,
                        "crossing_detected": int(crossing is not None),
                    }
                )
                max_norm_error = max(max_norm_error, norm_error)
                max_zero_square_abs = max(max_zero_square_abs, zero_square_abs)

            if t in (55000, XMAX):
                np.savez_compressed(
                    checkpoint_path(n, t),
                    N=np.array(n),
                    step=np.array(t),
                    Z=Z,
                    wp=wp,
                )

            if t >= next_progress:
                elapsed = time.perf_counter() - started
                rate = t / elapsed
                eta = (XMAX - t) / rate if rate > 0 else float("nan")
                print(
                    f"[PROGRESS] N={n} t={t}/{XMAX} "
                    f"elapsed={elapsed:.1f}s ETA={eta:.1f}s",
                    flush=True,
                )
                next_progress += 10000

            if t == XMAX:
                break
            sys_lr.set_theta(np.angle(Z))
            sigma_estimate, wp = sys_lr.sigma_max_power(wp)
            Z = sys_lr.cayley_step(Z, sigma_estimate)

    elapsed = time.perf_counter() - started
    summary = {
        "N": n,
        "M": int(sys_lr.m),
        "condition": "A_seedless_natural",
        "xmax": XMAX,
        "sample_every": sample_every,
        "crossing": crossing,
        "metastable_start": None if crossing is None else crossing + base.GUARD,
        "initial_state_rule": "Z0 = v.copy(); kernel seed is not generated",
        "initial_seed_enabled": False,
        "metastable_seed_enabled": False,
        "benettin_enabled": False,
        "state_feedback_from_observation": False,
        "method_parent_basis": method,
        "parent_residual": float(parent_residual),
        "parent_sigma": [float(value) for value in parent_sigma],
        "max_norm_error": max_norm_error,
        "max_zero_square_abs": max_zero_square_abs,
        "elapsed_seconds": elapsed,
        "csv_sha256": base.sha256(target),
        "checkpoints": [
            str(checkpoint_path(n, 55000)),
            str(checkpoint_path(n, XMAX)),
        ],
    }
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[DONE] N={n} t={XMAX} elapsed={elapsed:.1f}s", flush=True)
    return summary


def reference_csv(n: int) -> Path:
    reproduced = (
        HERE
        / "outputs"
        / "raw"
        / f"N{n:05d}"
        / "paper7_long_timeseries.csv"
    )
    if reproduced.is_file():
        return reproduced
    return base.PAPER8 / "raw" / f"N{n:05d}" / "condition_A_no_seed.csv"


def compare_prefix(n: int) -> dict:
    target = csv_path(n)
    reference = reference_csv(n)
    with target.open(encoding="utf-8") as handle:
        long_rows = {
            int(row["step"]): row
            for row in csv.DictReader(handle)
            if int(row["step"]) <= base.XMAX
        }
    with reference.open(encoding="utf-8") as handle:
        reference_rows = {
            int(row["step"]): row for row in csv.DictReader(handle)
        }
    columns = ("f_outside_parent", "q1", "q2", "q3", "q4")
    first_difference = None
    for step, ref_row in reference_rows.items():
        row = long_rows.get(step)
        if row is None:
            first_difference = {"step": step, "reason": "missing"}
            break
        for column in columns:
            if row[column] != ref_row[column]:
                first_difference = {
                    "step": step,
                    "column": column,
                    "reference": ref_row[column],
                    "long_run": row[column],
                }
                break
        if first_difference:
            break
    result = {
        "N": n,
        "reference": str(reference),
        "long_run": str(target),
        "checked_through": base.XMAX,
        "columns": list(columns),
        "reference_rows": len(reference_rows),
        "long_prefix_rows": len(long_rows),
        "string_equal": first_difference is None,
        "first_difference": first_difference,
    }
    out = OUTPUT / "comparison" / f"N{n:05d}_prefix_through_55000.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"[COMPARE] N={n} prefix_string_equal={result['string_equal']}",
        flush=True,
    )
    return result


def turning_points(t: np.ndarray, values: np.ndarray) -> list[dict]:
    difference = np.diff(values)
    nonzero_indices = np.flatnonzero(difference != 0.0)
    if len(nonzero_indices) < 2:
        return []
    signs = np.sign(difference[nonzero_indices])
    points = []
    for j in range(1, len(signs)):
        if signs[j] == signs[j - 1]:
            continue
        index = int(nonzero_indices[j])
        points.append(
            {
                "step": int(t[index]),
                "value": float(values[index]),
                "kind": "maximum" if signs[j - 1] > 0 else "minimum",
            }
        )
    return points


def fit_metrics(values: np.ndarray, trend: np.ndarray, parameter_count: int):
    residual = values - trend
    rss = float(np.sum(residual**2))
    total_ss = float(np.sum((values - np.mean(values)) ** 2))
    n = len(values)
    return {
        "r_squared": 1.0 - rss / total_ss if total_ss else 1.0,
        "rss": rss,
        "aic": n * np.log(max(rss / n, np.finfo(float).tiny))
        + 2 * parameter_count,
        "residual_std": float(np.std(residual)),
        "residual_peak_to_peak": float(np.ptp(residual)),
    }


def analyze() -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    observables = (
        ("f_outside_parent", "f", "black"),
        ("q3", "q3", "#F58518"),
        ("q4", "q4", "#E45756"),
    )
    report = {
        "fit_window": [FIT_START, XMAX],
        "amplification": AMPLIFICATION,
        "interpretation_rule": (
            "A long-period candidate must survive both one- and "
            "two-exponential baseline subtraction and show repeated extrema."
        ),
        "N": {},
    }
    initial_report = {
        "source": "same long-run CSV prefix, verified against existing t<=55000 data",
        "rule": (
            "Only raw turning points after f>=1e-6 are treated as "
            "resolved dynamical candidates; machine-floor reversals are excluded."
        ),
        "N": {},
    }
    figure_dir = OUTPUT / "figures"
    summary_dir = OUTPUT / "summary"
    figure_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8.4), sharex=True)
    sensitivity_fig, sensitivity_axes = plt.subplots(
        2, 3, figsize=(15, 8.4), sharex=True
    )
    terminal_fig, terminal_axes = plt.subplots(
        2, 3, figsize=(15, 8.4), sharex=True
    )
    initial_fig, initial_axes = plt.subplots(
        2, 3, figsize=(15, 8.4), squeeze=False
    )
    for row_index, n in enumerate(NS):
        with csv_path(n).open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        t_all = np.array([int(row["step"]) for row in rows])
        f_all = np.array([float(row["f_outside_parent"]) for row in rows])
        meta = json.loads(
            (
                summary_dir / f"N{n:05d}_long_horizon_meta.json"
            ).read_text(encoding="utf-8")
        )
        crossing = int(meta["crossing"])
        positive_growth_indices = np.flatnonzero(f_all >= 1e-10)
        initial_start = int(t_all[positive_growth_indices[0]])
        initial_end = crossing + base.GUARD
        initial_mask = (t_all >= initial_start) & (t_all <= initial_end)
        initial_report["N"][str(n)] = {
            "sample_every": SAMPLE[n],
            "crossing": crossing,
            "window": [initial_start, initial_end],
            "observables": {},
        }
        mask = t_all >= FIT_START
        t = t_all[mask]
        u = t - FIT_START
        report["N"][str(n)] = {"observables": {}}

        for column_index, (key, label, color) in enumerate(observables):
            values_all = np.array([float(row[key]) for row in rows])
            values = values_all[mask]
            one_parameters, one_trend = fit_one_exp(u, values)
            two_parameters, two_trend = fit_two_exp(u, values)
            one_residual = values - one_trend
            two_residual = values - two_trend
            raw_points = turning_points(t, values)
            one_points = turning_points(t, one_residual)
            two_points = turning_points(t, two_residual)
            report["N"][str(n)]["observables"][key] = {
                "raw_turning_points": raw_points,
                "one_exponential": {
                    "parameters_c_a_tau": [float(x) for x in one_parameters],
                    **fit_metrics(values, one_trend, 3),
                    "residual_turning_points": one_points,
                },
                "two_exponential": {
                    "parameters_c_a1_tau1_a2_tau2": [
                        float(x) for x in two_parameters
                    ],
                    **fit_metrics(values, two_trend, 5),
                    "residual_turning_points": two_points,
                },
            }

            terminal_mask = t_all >= 55000
            terminal_t = t_all[terminal_mask]
            terminal_values = values_all[terminal_mask]
            terminal_relative = terminal_values - terminal_values[-1]
            terminal_points = turning_points(terminal_t, terminal_values)
            report["N"][str(n)]["observables"][key]["terminal_raw"] = {
                "window": [55000, XMAX],
                "value_at_55000": float(terminal_values[0]),
                "value_at_110000": float(terminal_values[-1]),
                "change": float(terminal_values[-1] - terminal_values[0]),
                "peak_to_peak": float(np.ptp(terminal_values)),
                "turning_points": terminal_points,
            }

            ax = axes[row_index, column_index]
            ax.axhline(0, color="#303030", ls="--", lw=0.9)
            ax.plot(
                t,
                AMPLIFICATION * one_residual,
                color=color,
                lw=0.9,
            )
            ax.set_title(f"N={n}  {label}: one-exp residual ×20")
            ax.grid(alpha=0.15)
            if column_index == 0:
                ax.set_ylabel("20 × residual")

            sax = sensitivity_axes[row_index, column_index]
            sax.axhline(0, color="#303030", ls="--", lw=0.9)
            sax.plot(
                t,
                AMPLIFICATION * one_residual,
                color="#B0B0B0",
                lw=0.8,
                label="one-exp",
            )
            sax.plot(
                t,
                AMPLIFICATION * two_residual,
                color=color,
                lw=0.9,
                label="two-exp",
            )
            sax.set_title(f"N={n}  {label}: baseline sensitivity")
            sax.grid(alpha=0.15)
            if column_index == 0:
                sax.set_ylabel("20 × residual")

            tax = terminal_axes[row_index, column_index]
            tax.axhline(0, color="#303030", ls="--", lw=0.9)
            tax.plot(terminal_t, terminal_relative, color=color, lw=1.0)
            tax.set_title(f"N={n}  {label}: raw − terminal value")
            tax.grid(alpha=0.15)
            if column_index == 0:
                tax.set_ylabel("raw data − x(110000)")

            iax = initial_axes[row_index, column_index]
            initial_t = t_all[initial_mask]
            initial_values = values_all[initial_mask]
            resolved_mask = initial_mask & (f_all >= 1e-6)
            resolved_t = t_all[resolved_mask]
            resolved_values = values_all[resolved_mask]
            resolved_points = turning_points(resolved_t, resolved_values)
            initial_report["N"][str(n)]["observables"][key] = {
                "raw_turning_points_after_f_ge_1e-6": resolved_points,
                "turning_count": len(resolved_points),
                "minimum": float(np.min(initial_values)),
                "maximum": float(np.max(initial_values)),
            }
            if key == "f_outside_parent":
                positive = initial_values > 0.0
                iax.semilogy(
                    initial_t[positive],
                    initial_values[positive],
                    color=color,
                    marker=".",
                    ms=2.5,
                    lw=0.9,
                )
            else:
                iax.plot(
                    initial_t,
                    initial_values,
                    color=color,
                    marker=".",
                    ms=2.5,
                    lw=0.9,
                )
            maxima = [point for point in resolved_points if point["kind"] == "maximum"]
            minima = [point for point in resolved_points if point["kind"] == "minimum"]
            if maxima:
                iax.scatter(
                    [point["step"] for point in maxima],
                    [point["value"] for point in maxima],
                    marker="^",
                    s=25,
                    color="#2A9D8F",
                    zorder=3,
                    label="resolved maximum",
                )
            if minima:
                iax.scatter(
                    [point["step"] for point in minima],
                    [point["value"] for point in minima],
                    marker="v",
                    s=25,
                    color="#7B2CBF",
                    zorder=3,
                    label="resolved minimum",
                )
            iax.axvline(
                crossing,
                color="#666666",
                ls="--",
                lw=0.8,
                label="crossing",
            )
            iax.set_title(f"N={n}  {label}: existing samples")
            iax.set_xlim(initial_start, initial_end)
            iax.grid(alpha=0.15)
            if column_index == 0:
                iax.set_ylabel("raw observable")

    for ax in axes[-1, :]:
        ax.set_xlabel("step (absolute)")
        ax.set_xlim(FIT_START, XMAX)
    for ax in sensitivity_axes[-1, :]:
        ax.set_xlabel("step (absolute)")
        ax.set_xlim(FIT_START, XMAX)
    for ax in terminal_axes[-1, :]:
        ax.set_xlabel("step (absolute)")
        ax.set_xlim(55000, XMAX)
    for ax in initial_axes[-1, :]:
        ax.set_xlabel("step (absolute)")
    sensitivity_axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Seedless natural trajectory to t=110000: "
        "one-exponential asymptotic residual (display ×20)"
    )
    sensitivity_fig.suptitle(
        "Long-horizon residual: one- vs two-exponential baseline "
        "(display ×20)"
    )
    terminal_fig.suptitle(
        "Seedless terminal approach without fitting: "
        "raw observable minus x(110000)"
    )
    initial_fig.suptitle(
        "Initial growth in existing samples "
        "(markers: raw extrema after f ≥ 10⁻⁶)"
    )
    fig.tight_layout()
    sensitivity_fig.tight_layout()
    terminal_fig.tight_layout()
    initial_fig.tight_layout()
    fig.savefig(
        figure_dir / "figure_long_horizon_one_exp_residual_x20.png", dpi=150
    )
    fig.savefig(
        figure_dir / "figure_long_horizon_one_exp_residual_x20.svg"
    )
    sensitivity_fig.savefig(
        figure_dir / "figure_long_horizon_baseline_sensitivity_x20.png",
        dpi=150,
    )
    sensitivity_fig.savefig(
        figure_dir / "figure_long_horizon_baseline_sensitivity_x20.svg"
    )
    terminal_fig.savefig(
        figure_dir / "figure_long_horizon_terminal_raw_relative.png",
        dpi=150,
    )
    terminal_fig.savefig(
        figure_dir / "figure_long_horizon_terminal_raw_relative.svg"
    )
    initial_fig.savefig(
        figure_dir / "figure_initial_growth_existing_samples.png",
        dpi=150,
    )
    initial_fig.savefig(
        figure_dir / "figure_initial_growth_existing_samples.svg"
    )
    plt.close(fig)
    plt.close(sensitivity_fig)
    plt.close(terminal_fig)
    plt.close(initial_fig)
    report_path = summary_dir / "long_horizon_asymptotic_analysis.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (
        summary_dir / "initial_growth_existing_sample_assessment.json"
    ).write_text(
        json.dumps(initial_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[ANALYZE] {report_path}", flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("N", type=int)
    compare_parser = subparsers.add_parser("compare-prefix")
    compare_parser.add_argument("N", type=int)
    subparsers.add_parser("analyze")
    args = parser.parse_args()
    if args.command == "run":
        run_long(args.N)
    elif args.command == "compare-prefix":
        compare_prefix(args.N)
    elif args.command == "analyze":
        analyze()


if __name__ == "__main__":
    main()
