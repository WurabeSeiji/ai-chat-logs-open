from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from run_gray_cat_ab_metastable_interface_preliminary_v1 import (
    apply_exchange,
    apply_stability_gain,
    normalize_pair,
    state_from_s_phi,
)
from run_gray_cat_c_readout_window_preliminary_v1 import apply_c_backaction, c_visibility
from run_gray_cat_d_observation_response_preliminary_v1 import apply_d_backaction, d_visibility


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "gray_cat_observed_value_transition_figures_v1"


COLORS = {
    "ink": "#1f2933",
    "muted": "#667085",
    "grid": "#d8dee9",
    "A": "#2676b8",
    "B": "#c47b2c",
    "S": "#475467",
    "C": "#4b8b5b",
    "D": "#b94d4d",
    "AB": "#eef4ff",
    "ABC": "#f2f7ee",
    "ABCD": "#fff1f1",
}


@dataclass(frozen=True)
class Scenario:
    key: str
    title: str
    epsilon: float
    phi: float
    s0: float
    stability_gain: float
    c_g: float
    c_backaction_scale: float
    d_g: float
    d_backaction_scale: float
    ab_steps: int = 240
    abc_steps: int = 240
    abcd_steps: int = 900


def setup_font() -> None:
    plt.rcParams["font.family"] = [
        "Hiragino Sans",
        "YuGothic",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["svg.fonttype"] = "path"


def measure(a: complex, b: complex) -> Tuple[float, float, float, float]:
    p_a = abs(a) ** 2
    p_b = abs(b) ** 2
    q = p_a + p_b
    if q <= 0.0:
        raise ValueError("zero AB norm")
    p_a /= q
    p_b /= q
    return p_a, p_b, p_a - p_b, q


def record(
    rows: List[Dict[str, float | str]],
    scenario: Scenario,
    step: int,
    stage: str,
    a: complex,
    b: complex,
    c_vis: float,
    d_vis: float,
) -> None:
    p_a, p_b, s_value, q = measure(a, b)
    c_s = c_vis * s_value if stage in ("ABC", "ABCD") else math.nan
    d_s = d_vis * s_value if stage == "ABCD" else math.nan
    rows.append(
        {
            "scenario": scenario.key,
            "stage": stage,
            "step": step,
            "p_A": p_a,
            "p_B": p_b,
            "S": s_value,
            "Q": q,
            "C_A": 0.5 * (1.0 + c_s) if not math.isnan(c_s) else math.nan,
            "C_B": 0.5 * (1.0 - c_s) if not math.isnan(c_s) else math.nan,
            "D_A": 0.5 * (1.0 + d_s) if not math.isnan(d_s) else math.nan,
            "D_B": 0.5 * (1.0 - d_s) if not math.isnan(d_s) else math.nan,
        }
    )


def run_scenario(scenario: Scenario) -> List[Dict[str, float | str]]:
    a, b = state_from_s_phi(scenario.s0, scenario.phi)
    rows: List[Dict[str, float | str]] = []
    c_vis = c_visibility(scenario.c_g, 0.02)
    d_vis = d_visibility(scenario.d_g, 0.02)
    c_gain = scenario.c_g * scenario.c_backaction_scale
    d_gain = scenario.d_g * scenario.d_backaction_scale
    step = 0

    for _ in range(scenario.ab_steps):
        record(rows, scenario, step, "AB", a, b, c_vis, d_vis)
        a, b = apply_exchange(a, b, scenario.epsilon)
        a, b = apply_stability_gain(a, b, scenario.stability_gain)
        step += 1

    for _ in range(scenario.abc_steps):
        p_a, p_b, s_value, _ = measure(a, b)
        s_c = c_vis * s_value
        record(rows, scenario, step, "ABC", a, b, c_vis, d_vis)
        a, b = apply_exchange(a, b, scenario.epsilon)
        a, b = apply_stability_gain(a, b, scenario.stability_gain)
        a, b = apply_c_backaction(a, b, c_gain, s_c)
        step += 1

    for _ in range(scenario.abcd_steps):
        p_a, p_b, s_value, _ = measure(a, b)
        s_d = d_vis * s_value
        record(rows, scenario, step, "ABCD", a, b, c_vis, d_vis)
        a, b = apply_exchange(a, b, scenario.epsilon)
        a, b = apply_stability_gain(a, b, scenario.stability_gain)
        a, b = apply_d_backaction(a, b, d_gain, s_d)
        a, b = normalize_pair(a, b)
        step += 1

    record(rows, scenario, step, "ABCD", a, b, c_vis, d_vis)
    return rows


def write_csv(path: Path, rows: List[Dict[str, float | str]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def stage_bounds(rows: List[Dict[str, float | str]]) -> Dict[str, Tuple[float, float]]:
    out: Dict[str, Tuple[float, float]] = {}
    for stage in ("AB", "ABC", "ABCD"):
        xs = [float(row["step"]) for row in rows if row["stage"] == stage]
        out[stage] = (min(xs), max(xs))
    return out


def plot_rows(scenarios: List[Scenario], rows_by_key: Dict[str, List[Dict[str, float | str]]]) -> None:
    setup_font()
    fig, axes = plt.subplots(len(scenarios), 1, figsize=(13.5, 8.8), sharex=False)
    if len(scenarios) == 1:
        axes = [axes]

    for ax, scenario in zip(axes, scenarios):
        rows = rows_by_key[scenario.key]
        bounds = stage_bounds(rows)
        for stage, color in (("AB", "AB"), ("ABC", "ABC"), ("ABCD", "ABCD")):
            left, right = bounds[stage]
            ax.axvspan(left, right, color=COLORS[color], alpha=0.85, zorder=0)
            ax.text((left + right) / 2.0, 1.045, stage, ha="center", va="center", fontsize=10, color=COLORS["muted"])
        for boundary in (bounds["ABC"][0], bounds["ABCD"][0]):
            ax.axvline(boundary, color=COLORS["grid"], lw=1.2)

        x = [float(row["step"]) for row in rows]
        p_a = [float(row["p_A"]) for row in rows]
        p_b = [float(row["p_B"]) for row in rows]
        s_scaled = [0.5 * (float(row["S"]) + 1.0) for row in rows]
        read_a = [
            float(row["C_A"])
            if row["stage"] == "ABC"
            else float(row["D_A"])
            if row["stage"] == "ABCD"
            else math.nan
            for row in rows
        ]
        read_b = [
            float(row["C_B"])
            if row["stage"] == "ABC"
            else float(row["D_B"])
            if row["stage"] == "ABCD"
            else math.nan
            for row in rows
        ]

        ax.plot(x, p_a, color=COLORS["A"], lw=1.9, label="p_A")
        ax.plot(x, p_b, color=COLORS["B"], lw=1.9, label="p_B")
        ax.plot(x, s_scaled, color=COLORS["S"], lw=1.15, ls="--", label="(S+1)/2")
        ax.plot(x, read_a, color=COLORS["A"], lw=1.0, ls=":", label="C_A or D_A")
        ax.plot(x, read_b, color=COLORS["B"], lw=1.0, ls=":", label="C_B or D_B")
        ax.set_ylim(-0.05, 1.08)
        ax.set_ylabel("probability / scaled S")
        ax.set_title(scenario.title, loc="left", fontsize=12, fontweight="bold")
        ax.grid(True, axis="y", color=COLORS["grid"], lw=0.8)

    axes[-1].set_xlabel("internal update step")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.985))
    fig.suptitle("白猫・黒猫・灰色猫: 実測値による AB -> ABC -> ABCD 遷移", fontsize=17, fontweight="bold", y=1.025)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(OUT_DIR / "gray_cat_ab_abc_abcd_observed_values_three_scenarios_v1.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT_DIR / "gray_cat_ab_abc_abcd_observed_values_three_scenarios_v1.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scenarios = [
        Scenario(
            key="metastable_d_selection",
            title="準安定混在: C弱読出しでは保持し、D強観測で片側へ選択",
            epsilon=0.01,
            phi=0.0,
            s0=0.01,
            stability_gain=0.0,
            c_g=1.0,
            c_backaction_scale=1.0e-5,
            d_g=0.1,
            d_backaction_scale=1.0,
        ),
        Scenario(
            key="gray_eigen_kept",
            title="灰色固有相: C弱読出しでもD強観測でも灰色を保持",
            epsilon=0.0,
            phi=0.0,
            s0=0.0,
            stability_gain=-0.01,
            c_g=1.0,
            c_backaction_scale=1.0e-5,
            d_g=1.0,
            d_backaction_scale=1.0,
        ),
        Scenario(
            key="strong_c_selection",
            title="Cが強すぎる場合: ABC段階で白側へ選択",
            epsilon=0.001,
            phi=0.0,
            s0=0.01,
            stability_gain=0.0,
            c_g=1.0,
            c_backaction_scale=0.01,
            d_g=1.0,
            d_backaction_scale=1.0,
        ),
    ]
    rows_by_key = {scenario.key: run_scenario(scenario) for scenario in scenarios}
    all_rows: List[Dict[str, float | str]] = []
    for scenario in scenarios:
        all_rows.extend(rows_by_key[scenario.key])
    write_csv(OUT_DIR / "gray_cat_ab_abc_abcd_observed_values_timeseries_v1.csv", all_rows)
    plot_rows(scenarios, rows_by_key)


if __name__ == "__main__":
    main()
