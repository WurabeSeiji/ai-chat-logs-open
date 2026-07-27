#!/usr/bin/env python3
"""t=0から幾何級数的成長域までを、数値床と成長残差に分けて図示する。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = HERE / "outputs" / "long_horizon_110000" / "raw"
OUTPUT = HERE / "outputs" / "initial_geometric_growth"
NS = (40, 300)
SAMPLE = {40: 25, 300: 100}
F_START = 1e-10
F_END = 1e-2
AMPLIFICATION = 20.0
OBSERVABLES = (
    ("f_outside_parent", "f", "black"),
    ("q3", "q3", "#F58518"),
    ("q4", "q4", "#E45756"),
)


def source_csv(n: int) -> Path:
    return INPUT / f"N{n:05d}_seedless_f_q3_q4_t110000.csv"


def turning_points(t: np.ndarray, values: np.ndarray) -> list[dict]:
    difference = np.diff(values)
    nonzero_indices = np.flatnonzero(difference != 0.0)
    if len(nonzero_indices) < 2:
        return []
    signs = np.sign(difference[nonzero_indices])
    points = []
    for index in range(1, len(signs)):
        if signs[index] == signs[index - 1]:
            continue
        value_index = int(nonzero_indices[index])
        points.append(
            {
                "step": int(t[value_index]),
                "value": float(values[value_index]),
                "kind": (
                    "maximum" if signs[index - 1] > 0 else "minimum"
                ),
            }
        )
    return points


def fit_report(
    t: np.ndarray,
    log_values: np.ndarray,
    degree: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    coefficients = np.polyfit(t, log_values, degree)
    trend = np.polyval(coefficients, t)
    residual = log_values - trend
    residual_ss = float(np.sum(residual**2))
    total_ss = float(np.sum((log_values - np.mean(log_values)) ** 2))
    return coefficients, residual, {
        "degree": degree,
        "coefficients_high_to_low": [float(x) for x in coefficients],
        "r_squared": (
            1.0 - residual_ss / total_ss if total_ss else 1.0
        ),
        "log10_residual_std": float(np.std(residual)),
        "log10_residual_peak_to_peak": float(np.ptp(residual)),
        "displayed_peak_to_peak_x20": float(
            AMPLIFICATION * np.ptp(residual)
        ),
    }


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = OUTPUT / "figures"
    summary_dir = OUTPUT / "summary"
    figure_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    report = {
        "definition": {
            "growth_start": f"first sampled f >= {F_START}",
            "growth_end": f"last sampled f <= {F_END}",
            "geometric_model": "log10(x) = intercept + slope * t",
            "display": "20 * (log10(data) - fitted log10 growth)",
            "baseline_sensitivity": (
                "linear geometric model versus quadratic log-growth model"
            ),
        },
        "N": {},
    }
    raw_fig, raw_axes = plt.subplots(2, 3, figsize=(15, 8.4))
    residual_fig, residual_axes = plt.subplots(2, 3, figsize=(15, 8.4))
    sensitivity_fig, sensitivity_axes = plt.subplots(
        2, 3, figsize=(15, 8.4)
    )

    for row_index, n in enumerate(NS):
        with source_csv(n).open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        t_all = np.array([int(row["step"]) for row in rows])
        f_all = np.array([float(row["f_outside_parent"]) for row in rows])
        growth_indices = np.flatnonzero(
            (f_all >= F_START) & (f_all <= F_END)
        )
        if not len(growth_indices):
            raise RuntimeError(f"N={n}: growth interval not found")
        first_index = int(growth_indices[0])
        last_index = int(growth_indices[-1])
        growth_slice = slice(first_index, last_index + 1)
        t_growth = t_all[growth_slice]
        growth_start = int(t_growth[0])
        growth_end = int(t_growth[-1])
        raw_mask = t_all <= growth_end
        t_raw = t_all[raw_mask]
        report["N"][str(n)] = {
            "source_csv": str(source_csv(n)),
            "sample_every": SAMPLE[n],
            "growth_window": [growth_start, growth_end],
            "growth_sample_count": int(len(t_growth)),
            "pre_growth_window": [0, growth_start],
            "observables": {},
        }

        for column_index, (key, label, color) in enumerate(OBSERVABLES):
            values_all = np.array([float(row[key]) for row in rows])
            values_growth = values_all[growth_slice]
            if np.any(values_growth <= 0):
                raise RuntimeError(f"N={n} {key}: non-positive growth value")
            log_values = np.log10(values_growth)
            linear_coefficients, linear_residual, linear_metrics = fit_report(
                t_growth, log_values, 1
            )
            _, quadratic_residual, quadratic_metrics = fit_report(
                t_growth, log_values, 2
            )
            raw_points = turning_points(t_growth, values_growth)
            report["N"][str(n)]["observables"][key] = {
                "minimum_in_growth_window": float(np.min(values_growth)),
                "maximum_in_growth_window": float(np.max(values_growth)),
                "raw_turning_points": raw_points,
                "raw_turning_count": len(raw_points),
                "linear_log_growth": linear_metrics,
                "quadratic_log_growth": quadratic_metrics,
            }

            raw_ax = raw_axes[row_index, column_index]
            raw_values = values_all[raw_mask]
            plotted_values = np.maximum(np.abs(raw_values), 1e-18)
            raw_ax.semilogy(
                t_raw,
                plotted_values,
                color=color,
                marker=".",
                ms=2.5,
                lw=0.8,
                label=f"|{label}|",
            )
            fitted_growth = 10.0 ** np.polyval(
                linear_coefficients, t_growth
            )
            raw_ax.semilogy(
                t_growth,
                fitted_growth,
                color="#277DA1",
                ls="--",
                lw=1.1,
                label="geometric fit",
            )
            if key == "f_outside_parent":
                raw_ax.axhspan(
                    1e-18,
                    1e-13,
                    color="#999999",
                    alpha=0.12,
                    label="numerical-floor region",
                )
            else:
                raw_ax.axhspan(
                    1e-10,
                    5e-8,
                    color="#999999",
                    alpha=0.12,
                    label="numerical-floor region",
                )
            raw_ax.axvline(growth_start, color="#666666", ls=":", lw=0.8)
            raw_ax.axvline(growth_end, color="#666666", ls=":", lw=0.8)
            raw_ax.set_title(f"N={n}  {label}: t=0 to growth end")
            raw_ax.set_xlim(0, growth_end)
            raw_ax.grid(alpha=0.15)
            if column_index == 0:
                raw_ax.set_ylabel("absolute raw value (log scale)")

            residual_ax = residual_axes[row_index, column_index]
            residual_ax.axhline(0, color="#303030", ls="--", lw=0.9)
            residual_ax.plot(
                t_growth,
                AMPLIFICATION * linear_residual,
                color=color,
                marker=".",
                ms=3,
                lw=0.9,
            )
            for point in raw_points:
                residual_index = int(
                    np.flatnonzero(t_growth == point["step"])[0]
                )
                residual_ax.scatter(
                    point["step"],
                    AMPLIFICATION * linear_residual[residual_index],
                    marker="^" if point["kind"] == "maximum" else "v",
                    s=28,
                    color=(
                        "#2A9D8F"
                        if point["kind"] == "maximum"
                        else "#7B2CBF"
                    ),
                    zorder=3,
                )
            residual_ax.set_title(
                f"N={n}  {label}: geometric residual ×20"
            )
            residual_ax.set_xlim(growth_start, growth_end)
            residual_ax.grid(alpha=0.15)
            if column_index == 0:
                residual_ax.set_ylabel(
                    "20 × log10(data / geometric fit)"
                )

            sensitivity_ax = sensitivity_axes[row_index, column_index]
            sensitivity_ax.axhline(0, color="#303030", ls="--", lw=0.9)
            sensitivity_ax.plot(
                t_growth,
                AMPLIFICATION * linear_residual,
                color="#B0B0B0",
                lw=0.9,
                label="linear log-growth",
            )
            sensitivity_ax.plot(
                t_growth,
                AMPLIFICATION * quadratic_residual,
                color=color,
                lw=0.9,
                label="quadratic log-growth",
            )
            sensitivity_ax.set_title(
                f"N={n}  {label}: growth-baseline sensitivity"
            )
            sensitivity_ax.set_xlim(growth_start, growth_end)
            sensitivity_ax.grid(alpha=0.15)
            if column_index == 0:
                sensitivity_ax.set_ylabel("20 × log residual")

    for axes in (raw_axes, residual_axes, sensitivity_axes):
        for ax in axes[-1, :]:
            ax.set_xlabel("step (absolute)")
    raw_axes[0, 0].legend(fontsize=7, loc="best")
    sensitivity_axes[0, 0].legend(fontsize=7, loc="best")
    raw_fig.suptitle(
        "Seedless initial interval: numerical floor and geometric growth"
    )
    residual_fig.suptitle(
        "Residual around geometric-growth envelope "
        "(display ×20 in log10 ratio)"
    )
    sensitivity_fig.suptitle(
        "Initial growth residual: baseline sensitivity "
        "(linear versus quadratic log growth)"
    )
    raw_fig.tight_layout()
    residual_fig.tight_layout()
    sensitivity_fig.tight_layout()
    for figure, name in (
        (raw_fig, "figure_initial_raw_to_geometric_growth"),
        (residual_fig, "figure_initial_geometric_growth_residual_x20"),
        (
            sensitivity_fig,
            "figure_initial_growth_baseline_sensitivity_x20",
        ),
    ):
        figure.savefig(figure_dir / f"{name}.png", dpi=150)
        figure.savefig(figure_dir / f"{name}.svg")
        plt.close(figure)
    report_path = summary_dir / "initial_geometric_growth_analysis.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[DONE] {report_path}", flush=True)


if __name__ == "__main__":
    main()
