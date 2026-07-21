"""N体ランク飽和予備実験の論文用図生成 v1

nbody_rank_saturation_preliminary_result_v1.json を読み、論文用PNGを生成する。
"""

from __future__ import annotations

import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "nbody_rank_saturation_preliminary_result_v1"
MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SERIES_BLUE = "#2a78d6"
SERIES_GREEN = "#008300"
SERIES_MAGENTA = "#e87ba4"
NEUTRAL_GRAY = "#52514e"


def load_results() -> dict:
    with (OUT_DIR / "nbody_rank_saturation_preliminary_result_v1.json").open(
        encoding="utf-8"
    ) as handle:
        return json.load(handle)


def figure_three_layer_scaling(results: dict) -> None:
    rows = results["experiment_b_rank_law"]
    body_counts = [row["body_count"] for row in rows]
    relation_counts = [row["relation_count"] for row in rows]
    observed_ranks = [row["expected_rank"] for row in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        body_counts, relation_counts, marker="o", color=SERIES_BLUE,
        label=r"pair-relation waves  $M_N=\binom{N}{2}$  ($O(N^2)$)",
    )
    ax.plot(
        body_counts, observed_ranks, marker="s", color=SERIES_GREEN,
        label=r"generator rank  $2\min(N,\lfloor M/2\rfloor)$  ($O(N)$)",
    )
    ax.axhline(
        3, color=NEUTRAL_GRAY, linestyle="--", linewidth=1.2,
        label="uniquely readable spatial directions  (3)",
    )
    ax.set_xticks(body_counts)
    ax.set_xlabel("body count N")
    ax.set_ylabel("count")
    ax.set_title("Three-layer scaling: relation waves, rotation modes, spatial directions")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "three_layer_scaling_v1.png", dpi=180)
    plt.close(fig)


def figure_conservation_errors(results: dict) -> None:
    rows = results["experiment_a_summaries"]
    body_counts = [row["body_count"] for row in rows]
    closure = [row["max_closure_error"] for row in rows]
    amplitude = [row["max_amplitude_drift"] for row in rows]
    covariance = [row["max_trajectory_covariance_error"] for row in rows]
    tolerance = results["parameters"]["invariant_tol"]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(body_counts, closure, marker="o", color=SERIES_BLUE,
                label="max squared-closure error")
    ax.semilogy(body_counts, amplitude, marker="s", color=SERIES_GREEN,
                label="max absolute-square sum drift")
    ax.semilogy(body_counts, covariance, marker="^", color=SERIES_MAGENTA,
                label="max relabeling-covariance error")
    ax.axhline(tolerance, color=NEUTRAL_GRAY, linestyle="--", linewidth=1.2,
               label=r"tolerance $10^{-10}$")
    ax.set_xticks(body_counts)
    ax.set_xlabel("body count N")
    ax.set_ylabel("maximum error over 32 trials x 720 steps")
    ax.set_title("Conservation and covariance errors under the fixed generator")
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="center right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "nbody_conservation_errors_v1.png", dpi=180)
    plt.close(fig)


def figure_kernel_dimension(results: dict) -> None:
    rows = results["experiment_a_summaries"]
    body_counts = [row["body_count"] for row in rows]
    nullities = [row["nullity_values"][0] for row in rows]
    bound_counts = [n for n in body_counts if n >= 5]
    bounds = [n * (n - 5) // 2 for n in bound_counts]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(body_counts, nullities, marker="o", color=SERIES_BLUE,
            label=r"observed  $\dim\ker K$  (all 32 trials identical)")
    ax.plot(bound_counts, bounds, linestyle="--", marker="s", color=SERIES_GREEN,
            label=r"lower bound  $M_N-2N=\frac{N(N-5)}{2}$  ($N\geq5$)")
    ax.set_xticks(body_counts)
    ax.set_xlabel("body count N")
    ax.set_ylabel("kernel dimension")
    ax.set_title("Residual subspace growth: the kernel reappears and dominates")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "kernel_dimension_growth_v1.png", dpi=180)
    plt.close(fig)


def figure_normal_gauge_freedom(results: dict) -> None:
    rows = results["experiment_c_normal_uniqueness"]
    dims = [row["display_dimension"] for row in rows]
    angles = [row["max_candidate_line_angle_deg"] for row in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar([str(d) for d in dims], angles, width=0.5, color=SERIES_BLUE)
    for bar, angle in zip(bars, angles):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
                f"{angle:.4f}" + "°", ha="center", va="bottom",
                color=NEUTRAL_GRAY)
    ax.set_xlabel("display dimension d")
    ax.set_ylabel("max angle between gauge-equivalent normals (deg)")
    ax.set_title(
        "Normal uniqueness: only d = 3 pins the normal "
        "(projector invariant in all cases)"
    )
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "normal_gauge_freedom_v1.png", dpi=180)
    plt.close(fig)


def main() -> None:
    results = load_results()
    figure_three_layer_scaling(results)
    figure_conservation_errors(results)
    figure_kernel_dimension(results)
    figure_normal_gauge_freedom(results)
    for name in (
        "three_layer_scaling_v1.png",
        "nbody_conservation_errors_v1.png",
        "kernel_dimension_growth_v1.png",
        "normal_gauge_freedom_v1.png",
    ):
        print(OUT_DIR / name)


if __name__ == "__main__":
    main()
