from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from draw_gray_cat_ab_abc_abcd_state_diagrams_v1 import (
    COLORS as STATE_COLORS,
    add_label,
    draw_arrow,
    draw_card,
    draw_cat,
    draw_observer,
    draw_or_pair,
)
from draw_gray_cat_observed_value_transition_figures_v1 import (
    COLORS as VALUE_COLORS,
    measure,
    record,
    stage_bounds,
    write_csv,
)
from run_gray_cat_ab_metastable_interface_preliminary_v1 import (
    apply_exchange,
    apply_stability_gain,
    normalize_pair,
    state_from_s_phi,
)
from run_gray_cat_c_readout_window_preliminary_v1 import apply_c_backaction, c_visibility
from run_gray_cat_d_observation_response_preliminary_v1 import apply_d_backaction, d_visibility


BASE_DIR = Path(__file__).resolve().parent
STATE_OUT_DIR = BASE_DIR / "gray_cat_state_transition_figures_v1"
VALUE_OUT_DIR = BASE_DIR / "gray_cat_observed_value_transition_figures_v1"


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
    plt.rcParams["font.family"] = ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["svg.fonttype"] = "path"


def draw_scenario_en(
    ax,
    y: float,
    row_title: str,
    states: tuple[str, str, str],
    subtitles: tuple[str, str, str],
    row_color: str,
) -> None:
    x_positions = (2.1, 5.2, 8.3)
    for x, col_title, subtitle in zip(x_positions, ("AB", "ABC", "ABCD"), subtitles):
        draw_card(ax, x, y, 2.25, 1.85, col_title, subtitle, edge=row_color)

    for x, state in zip(x_positions, states):
        if state == "or":
            draw_or_pair(ax, x, y - 0.03, 0.48)
        else:
            draw_cat(ax, x, y - 0.03, 0.43, state)

    draw_observer(ax, x_positions[1] + 0.92, y + 0.52, "C", "weak", "blue")
    draw_observer(ax, x_positions[2] + 0.92, y + 0.52, "D", "strong", "red")
    draw_arrow(ax, 3.35, y, 3.95, y)
    draw_arrow(ax, 6.45, y, 7.05, y)
    add_label(ax, 0.82, y + 0.35, row_title, size=10, weight="bold", color=row_color, ha="right")


def draw_state_overview_en() -> None:
    setup_font()
    fig, ax = plt.subplots(figsize=(13.5, 8.8))
    ax.set_xlim(0, 10.0)
    ax.set_ylim(0, 7.4)
    ax.axis("off")
    add_label(ax, 5.0, 7.05, "White Cat, Black Cat, Gray Cat: AB -> ABC -> ABCD", size=18, weight="bold")
    add_label(
        ax,
        5.0,
        6.65,
        "C reads weakly. D observes strongly. Gray metastable and gray eigen phases are separated.",
        size=11,
        color="muted",
    )

    draw_scenario_en(
        ax,
        5.55,
        "Metastable\nmix",
        ("mix", "mix", "or"),
        ("make white + black", "read without breaking", "select white or black"),
        "green",
    )
    draw_scenario_en(
        ax,
        3.45,
        "Gray eigen\nphase",
        ("gray", "gray", "gray"),
        ("gray cat is established", "read as gray", "remains gray"),
        "amber",
    )
    draw_scenario_en(
        ax,
        1.35,
        "C too\nstrong",
        ("mix", "or", "or"),
        ("make white + black", "falls to white or black by C", "read after selection"),
        "red",
    )

    fig.savefig(STATE_OUT_DIR / "gray_cat_ab_abc_abcd_three_scenarios_en_v1.png", dpi=220, bbox_inches="tight")
    fig.savefig(STATE_OUT_DIR / "gray_cat_ab_abc_abcd_three_scenarios_en_v1.svg", bbox_inches="tight")
    plt.close(fig)


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
        _, _, s_value, _ = measure(a, b)
        s_c = c_vis * s_value
        record(rows, scenario, step, "ABC", a, b, c_vis, d_vis)
        a, b = apply_exchange(a, b, scenario.epsilon)
        a, b = apply_stability_gain(a, b, scenario.stability_gain)
        a, b = apply_c_backaction(a, b, c_gain, s_c)
        step += 1

    for _ in range(scenario.abcd_steps):
        _, _, s_value, _ = measure(a, b)
        s_d = d_vis * s_value
        record(rows, scenario, step, "ABCD", a, b, c_vis, d_vis)
        a, b = apply_exchange(a, b, scenario.epsilon)
        a, b = apply_stability_gain(a, b, scenario.stability_gain)
        a, b = apply_d_backaction(a, b, d_gain, s_d)
        a, b = normalize_pair(a, b)
        step += 1

    record(rows, scenario, step, "ABCD", a, b, c_vis, d_vis)
    return rows


def plot_observed_values_en(scenarios: List[Scenario], rows_by_key: Dict[str, List[Dict[str, float | str]]]) -> None:
    setup_font()
    fig, axes = plt.subplots(len(scenarios), 1, figsize=(13.5, 8.8), sharex=False)
    if len(scenarios) == 1:
        axes = [axes]

    for ax, scenario in zip(axes, scenarios):
        rows = rows_by_key[scenario.key]
        bounds = stage_bounds(rows)
        for stage, color in (("AB", "AB"), ("ABC", "ABC"), ("ABCD", "ABCD")):
            left, right = bounds[stage]
            ax.axvspan(left, right, color=VALUE_COLORS[color], alpha=0.85, zorder=0)
            ax.text((left + right) / 2.0, 1.045, stage, ha="center", va="center", fontsize=10, color=VALUE_COLORS["muted"])
        for boundary in (bounds["ABC"][0], bounds["ABCD"][0]):
            ax.axvline(boundary, color=VALUE_COLORS["grid"], lw=1.2)

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

        ax.plot(x, p_a, color=VALUE_COLORS["A"], lw=1.9, label="p_A")
        ax.plot(x, p_b, color=VALUE_COLORS["B"], lw=1.9, label="p_B")
        ax.plot(x, s_scaled, color=VALUE_COLORS["S"], lw=1.15, ls="--", label="(S+1)/2")
        ax.plot(x, read_a, color=VALUE_COLORS["A"], lw=1.0, ls=":", label="C_A or D_A")
        ax.plot(x, read_b, color=VALUE_COLORS["B"], lw=1.0, ls=":", label="C_B or D_B")
        ax.set_ylim(-0.05, 1.08)
        ax.set_ylabel("probability / scaled S")
        ax.set_title(scenario.title, loc="left", fontsize=12, fontweight="bold")
        ax.grid(True, axis="y", color=VALUE_COLORS["grid"], lw=0.8)

    axes[-1].set_xlabel("internal update step")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.985))
    fig.suptitle("White Cat, Black Cat, Gray Cat: Observed AB -> ABC -> ABCD Transitions", fontsize=17, fontweight="bold", y=1.025)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(VALUE_OUT_DIR / "gray_cat_ab_abc_abcd_observed_values_three_scenarios_en_v1.png", dpi=220, bbox_inches="tight")
    fig.savefig(VALUE_OUT_DIR / "gray_cat_ab_abc_abcd_observed_values_three_scenarios_en_v1.svg", bbox_inches="tight")
    plt.close(fig)


def draw_observed_values_en() -> None:
    scenarios = [
        Scenario(
            key="metastable_d_selection",
            title="Metastable mix: weak C readout preserves it, strong D selects one side",
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
            title="Gray eigen phase: gray is kept under weak C and strong D",
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
            title="C is too strong: selection to the white side occurs at ABC",
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
    write_csv(VALUE_OUT_DIR / "gray_cat_ab_abc_abcd_observed_values_timeseries_en_v1.csv", all_rows)
    plot_observed_values_en(scenarios, rows_by_key)


def main() -> None:
    STATE_OUT_DIR.mkdir(parents=True, exist_ok=True)
    VALUE_OUT_DIR.mkdir(parents=True, exist_ok=True)
    draw_state_overview_en()
    draw_observed_values_en()


if __name__ == "__main__":
    main()
