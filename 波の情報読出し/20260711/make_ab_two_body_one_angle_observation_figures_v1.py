from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
PRELIM_DIR = BASE_DIR / "ab_two_body_one_angle_harmonic_readout_preliminary_result_v1"
SWEEP_DIR = BASE_DIR / "ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1"
OUT_DIR = BASE_DIR / "ab_two_body_one_angle_harmonic_readout_observation_figures_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "figure.dpi": 140,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
    }
)


COLORS = {
    "black": "#1f2328",
    "gray": "#6e7781",
    "light_gray": "#d0d7de",
    "blue": "#0969da",
    "green": "#1a7f37",
    "red": "#cf222e",
    "teal": "#007d8a",
}


def read_csv(path: Path) -> List[Dict[str, Any]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rows_for(
    rows: List[Dict[str, Any]],
    *,
    case_id: str,
    protocol: str,
    readout_mode: str,
) -> List[Dict[str, Any]]:
    selected = [
        row
        for row in rows
        if row["case_id"] == case_id and row["protocol"] == protocol and row["readout_mode"] == readout_mode
    ]
    selected.sort(key=lambda row: int(row["step"]))
    return selected


def save_dual(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.svg")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=220)
    plt.close(fig)


def float_values(rows: List[Dict[str, Any]], key: str) -> List[float]:
    return [float(row[key]) for row in rows]


def int_values(rows: List[Dict[str, Any]], key: str) -> List[int]:
    return [int(row[key]) for row in rows]


def make_harmonic_state_figure(series_rows: List[Dict[str, Any]]) -> None:
    off = rows_for(series_rows, case_id="near_pi_05deg", protocol="Protocol_B", readout_mode="readout_off")
    normal = rows_for(series_rows, case_id="near_pi_05deg", protocol="Protocol_B", readout_mode="readout_normal")
    strong = rows_for(series_rows, case_id="near_pi_05deg", protocol="Protocol_B", readout_mode="readout_strong")

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.2))
    fig.suptitle("AB one-angle readout: harmonic state", fontsize=15, fontweight="bold")

    ax = axes[0, 0]
    x = float_values(off, "internal_signed_deviation_rad")
    y = float_values(off, "closure_complement")
    ax.plot(x, y, color=COLORS["blue"], linewidth=1.4)
    ax.scatter([x[0]], [y[0]], s=28, color=COLORS["black"], zorder=3, label="start")
    ax.scatter([x[-1]], [y[-1]], s=28, facecolors="white", edgecolors=COLORS["blue"], zorder=3, label="end")
    ax.axhline(0.0, color=COLORS["light_gray"], linewidth=0.8)
    ax.axvline(0.0, color=COLORS["light_gray"], linewidth=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("read deviation component")
    ax.set_ylabel("closure complement")
    ax.set_title("closure auxiliary plane")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[0, 1]
    steps = int_values(off, "step")
    ax.plot(steps, float_values(off, "V_AB"), color=COLORS["black"], linewidth=1.2, label="readout_off")
    ax.plot(steps, float_values(normal, "V_AB"), color=COLORS["blue"], linewidth=1.0, label="readout_normal")
    ax.plot(steps, float_values(strong, "V_AB"), color=COLORS["red"], linewidth=1.0, label="readout_strong")
    ax.set_xlabel("step")
    ax.set_ylabel("V_AB")
    ax.set_title("label-free harmonic readout")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[1, 0]
    ax.plot(steps, float_values(off, "D_AB_near_deg"), color=COLORS["blue"], linewidth=1.1, label="near arc")
    ax.plot(steps, float_values(off, "D_AB_far_deg"), color=COLORS["green"], linewidth=1.1, label="far arc")
    ax.axhline(180.0, color=COLORS["black"], linestyle="--", linewidth=0.8)
    ax.set_xlabel("step")
    ax.set_ylabel("degrees")
    ax.set_title("two unlabeled arcs around pi")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[1, 1]
    ax.plot(steps, float_values(off, "f_AB_center"), color=COLORS["black"], linewidth=1.2, label="center display")
    ax.plot(steps, float_values(off, "f_AB_circle"), color=COLORS["teal"], linewidth=1.0, linestyle="--", label="circle display")
    ax.set_xlabel("step")
    ax.set_ylabel("f_AB readout")
    ax.set_title("f_AB display consistency, readout_off")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    fig.text(
        0.5,
        0.01,
        "Axes are readout coordinates for this figure only; they do not assert absolute background coordinates.",
        ha="center",
        fontsize=8,
        color=COLORS["gray"],
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    save_dual(fig, "ab_two_body_one_angle_harmonic_state_v1")


def make_protocol_degeneracy_figure(series_rows: List[Dict[str, Any]], protocol_rows: List[Dict[str, Any]]) -> None:
    f_rows = rows_for(series_rows, case_id="near_pi_10deg", protocol="Protocol_F", readout_mode="readout_off")
    b_rows = rows_for(series_rows, case_id="near_pi_10deg", protocol="Protocol_B", readout_mode="readout_off")
    steps = int_values(f_rows, "step")

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.4))
    fig.suptitle("Protocol F/B: internal display differs, label-free readout collapses", fontsize=14, fontweight="bold")

    ax = axes[0]
    ax.plot(steps, float_values(f_rows, "protocol_display_deviation_rad"), color=COLORS["black"], label="Protocol F display")
    ax.plot(steps, float_values(b_rows, "protocol_display_deviation_rad"), color=COLORS["red"], label="Protocol B display")
    ax.set_xlabel("step")
    ax.set_ylabel("display deviation [rad]")
    ax.set_title("signed/unsigned display")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[1]
    ax.plot(steps, float_values(f_rows, "V_AB"), color=COLORS["black"], linewidth=1.2, label="Protocol F")
    ax.plot(steps, float_values(b_rows, "V_AB"), color=COLORS["blue"], linestyle="--", linewidth=1.1, label="Protocol B")
    ax.set_xlabel("step")
    ax.set_ylabel("V_AB")
    ax.set_title("same label-free V_AB")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    selected = [row for row in protocol_rows if row["case_id"] == "near_pi_10deg"]
    labels = [row["readout_mode"].replace("readout_", "") for row in selected]
    x = np.arange(len(labels))
    d_diff = [float(row["max_D_AB_near_protocol_diff"]) for row in selected]
    v_diff = [float(row["max_V_AB_protocol_diff"]) for row in selected]
    display_diff = [float(row["max_protocol_display_deviation_diff"]) for row in selected]
    width = 0.26
    ax = axes[2]
    ax.bar(x - width, d_diff, width=width, color=COLORS["blue"], label="max D diff")
    ax.bar(x, v_diff, width=width, color=COLORS["green"], label="max V diff")
    ax.bar(x + width, display_diff, width=width, color=COLORS["red"], label="display diff")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_yscale("symlog", linthresh=1.0e-16)
    ax.set_title("degeneracy check")
    ax.set_ylabel("max difference")
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    save_dual(fig, "ab_two_body_one_angle_protocol_degeneracy_v1")


def make_readout_leak_response_figure(leak_rows: List[Dict[str, Any]], period_rows: List[Dict[str, Any]]) -> None:
    leak = [float(row["per_step_leak"]) for row in leak_rows]
    max_f = [float(row["max_normalized_f_AB_projection_error"]) for row in leak_rows]
    mean_f = [float(row["mean_normalized_f_AB_projection_error"]) for row in leak_rows]
    max_decay = [float(row["max_abs_decay_rate_V_AB"]) for row in leak_rows]
    min_envelope = [float(row["min_envelope_ratio_final_over_initial"]) for row in leak_rows]

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.0))
    fig.suptitle("Readout wave response: damping and f_AB perturbation", fontsize=15, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(leak, max_decay, marker="o", color=COLORS["red"], label="max |decay|")
    ax.set_xscale("symlog", linthresh=1.0e-7)
    ax.set_yscale("symlog", linthresh=1.0e-12)
    ax.set_xlabel("per-step readout leak")
    ax.set_ylabel("|decay_rate_V_AB|")
    ax.set_title("envelope decay")
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax = axes[0, 1]
    ax.plot(leak, max_f, marker="o", color=COLORS["black"], label="max")
    ax.plot(leak, mean_f, marker="o", color=COLORS["blue"], label="mean")
    ax.set_xscale("symlog", linthresh=1.0e-7)
    ax.set_yscale("symlog", linthresh=1.0e-12)
    ax.set_xlabel("per-step readout leak")
    ax.set_ylabel("normalized f_AB error")
    ax.set_title("f_AB projection perturbation")
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax = axes[1, 0]
    ax.plot(leak, min_envelope, marker="o", color=COLORS["green"])
    ax.set_xscale("symlog", linthresh=1.0e-7)
    ax.set_xlabel("per-step readout leak")
    ax.set_ylabel("min final/initial envelope")
    ax.set_title("remaining oscillation envelope")

    ax = axes[1, 1]
    ax.plot(
        [int(row["period_steps"]) for row in period_rows],
        [float(row["max_normalized_f_AB_projection_error"]) for row in period_rows],
        marker="o",
        color=COLORS["teal"],
        label="max normalized f error",
    )
    ax.plot(
        [int(row["period_steps"]) for row in period_rows],
        [float(row["mean_normalized_f_AB_projection_error"]) for row in period_rows],
        marker="o",
        color=COLORS["gray"],
        label="mean normalized f error",
    )
    ax.set_xlabel("period steps")
    ax.set_ylabel("normalized f_AB error")
    ax.set_title("period dependence")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    save_dual(fig, "ab_two_body_one_angle_readout_leak_response_v1")


def make_observation_summary_figure(prelim: Dict[str, Any], sweep: Dict[str, Any]) -> None:
    prelim_v = prelim["aggregate_verdict"]
    sweep_v = sweep["aggregate_verdict"]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    fig.suptitle("AB one-angle preliminary observations", fontsize=15, fontweight="bold")

    ax = axes[0]
    labels = [
        "max Q_closed",
        "max protocol D diff",
        "max protocol V diff",
        "readout_off decay",
        "readout_strong decay",
        "non-strong f error",
    ]
    values = [
        float(prelim_v["max_Q_closed_abs"]),
        float(prelim_v["max_D_AB_near_protocol_diff"]),
        float(prelim_v["max_V_AB_protocol_diff"]),
        float(prelim_v["readout_off_decay_max_abs"]),
        float(prelim_v["readout_strong_decay_min_abs"]),
        float(prelim_v["max_f_AB_projection_consistency_error_nonstrong"]),
    ]
    display_values = [max(value, 1.0e-18) for value in values]
    ax.barh(
        labels,
        display_values,
        color=[COLORS["gray"], COLORS["blue"], COLORS["blue"], COLORS["green"], COLORS["red"], COLORS["teal"]],
    )
    for idx, value in enumerate(values):
        ax.text(max(display_values[idx] * 1.4, 2.0e-18), idx, f"{value:.2e}", va="center", fontsize=8)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_title("primary verdict scales")
    ax.set_xlabel("absolute scale; exact zeros shown at floor")

    ax = axes[1]
    ax.axis("off")
    checklist = [
        ("sweep cases", f"{int(sweep_v['case_summary_count'])} cases", COLORS["black"]),
        ("F/B label-free collapse", "PASS" if sweep_v["label_free_protocol_degenerate_all_cases"] else "FAIL", COLORS["green"]),
        ("oscillation detected", "PASS" if sweep_v["oscillation_detected_all_cases"] else "FAIL", COLORS["green"]),
        ("decay monotonic by leak", "PASS" if sweep_v["decay_abs_monotonic_by_leak_all_grids"] else "FAIL", COLORS["green"]),
        (
            "f_AB error monotonic by leak",
            "PASS" if sweep_v["normalized_f_error_monotonic_by_leak_all_grids"] else "FAIL",
            COLORS["green"],
        ),
        (
            "strong readout perturbation",
            "CONTROL DETECTED" if sweep_v["strong_leak_perturbs_projection_all_cases"] else "NOT DETECTED",
            COLORS["red"],
        ),
        ("readout_off max decay", f"{float(sweep_v['readout_off_decay_max_abs']):.2e}", COLORS["green"]),
        ("max normalized f_AB error", f"{float(sweep_v['max_normalized_f_AB_projection_error']):.2e}", COLORS["red"]),
    ]
    ax.set_title("sweep verdict dashboard", pad=10)
    y0 = 0.90
    dy = 0.10
    for idx, (name, value, color) in enumerate(checklist):
        y = y0 - idx * dy
        ax.scatter([0.05], [y], s=70, color=color, transform=ax.transAxes, clip_on=False)
        ax.text(0.10, y, name, transform=ax.transAxes, va="center", fontsize=10, color=COLORS["black"])
        ax.text(0.72, y, value, transform=ax.transAxes, va="center", fontsize=10, color=color, fontweight="bold")

    fig.tight_layout(rect=(0, 0.02, 1, 0.92))
    save_dual(fig, "ab_two_body_one_angle_observation_summary_v1")


def log_fit_slope(xs: List[float], ys: List[float]) -> float:
    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    mask = (x > 0.0) & (y > 0.0)
    slope, _ = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return float(slope)


def make_phase_difference_scaling_figure(case_rows: List[Dict[str, Any]]) -> None:
    selected = [
        row
        for row in case_rows
        if row["protocol"] == "Protocol_B"
        and row["readout_mode"] == "readout_off"
        and int(row["period_steps"]) == 96
    ]
    selected.sort(key=lambda row: float(row["initial_deviation_deg"]))
    deviations_deg = [float(row["initial_deviation_deg"]) for row in selected]
    deviations_rad = [math.radians(value) for value in deviations_deg]
    v_max = [float(row["V_AB_max"]) for row in selected]
    omega_step = float(selected[0]["omega_step"])
    omega_discrete_sq = 4.0 * math.sin(omega_step / 2.0) ** 2
    f_expected = [omega_discrete_sq * value for value in deviations_rad]
    v_slope = log_fit_slope(deviations_rad, v_max)
    f_slope = log_fit_slope(deviations_rad, f_expected)

    leak_selected = [
        row
        for row in case_rows
        if row["protocol"] == "Protocol_B"
        and row["readout_mode"] == "leak_5e-5"
        and int(row["period_steps"]) == 96
    ]
    leak_selected.sort(key=lambda row: float(row["initial_deviation_deg"]))
    leak_decay = [abs(float(row["decay_rate_V_AB"])) for row in leak_selected]
    leak_f_error = [float(row["normalized_f_AB_projection_error"]) for row in leak_selected]

    ref_dev = deviations_rad[3]
    ref_f = f_expected[3]
    inv_ref = [ref_f * (ref_dev / value) for value in deviations_rad]
    inv2_ref = [ref_f * (ref_dev / value) ** 2 for value in deviations_rad]

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    fig.suptitle("Phase-difference scaling in the one-angle AB readout", fontsize=15, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(deviations_deg, v_max, marker="o", color=COLORS["blue"], label=f"V_AB max, slope {v_slope:.2f}")
    ax.plot(
        deviations_deg,
        [v_max[3] * (value / ref_dev) ** 2 for value in deviations_rad],
        linestyle="--",
        color=COLORS["gray"],
        label="delta^2 reference",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("initial deviation from pi [deg]")
    ax.set_ylabel("V_AB max")
    ax.set_title("symmetric deviation readout")
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax = axes[0, 1]
    ax.plot(deviations_deg, f_expected, marker="o", color=COLORS["black"], label=f"f_AB expected, slope {f_slope:.2f}")
    ax.plot(deviations_deg, inv_ref, linestyle="--", color=COLORS["red"], label="1/delta reference")
    ax.plot(deviations_deg, inv2_ref, linestyle=":", color=COLORS["red"], label="1/delta^2 reference")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("initial deviation from pi [deg]")
    ax.set_ylabel("f_AB scale")
    ax.set_title("no inverse-power behavior in S1")
    ax.legend(frameon=False, fontsize=8, loc="best")

    ax = axes[1, 0]
    ax.plot(deviations_deg, leak_decay, marker="o", color=COLORS["red"], label="|decay|, leak_5e-5")
    ax.plot(deviations_deg, leak_f_error, marker="o", color=COLORS["teal"], label="normalized f error, leak_5e-5")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("initial deviation from pi [deg]")
    ax.set_ylabel("scale")
    ax.set_title("readout perturbation is nearly deviation-independent")
    ax.legend(frameon=False, fontsize=8, loc="best")

    ax = axes[1, 1]
    ax.axis("off")
    lines = [
        ("S1 conclusion", COLORS["black"], 0.88, 12, "bold"),
        ("V_AB follows delta^2 by definition/readout.", COLORS["blue"], 0.74, 10, "normal"),
        ("f_AB scale follows delta^1 in this one-angle map.", COLORS["black"], 0.62, 10, "normal"),
        ("1/delta and 1/delta^2 references run the wrong way.", COLORS["red"], 0.50, 10, "normal"),
        ("Readout leakage affects decay and f_AB consistency,", COLORS["teal"], 0.38, 10, "normal"),
        ("but not as an inverse-power phase-distance law.", COLORS["teal"], 0.30, 10, "normal"),
        ("Therefore inverse-square remains an S2/S3 test,", COLORS["black"], 0.18, 10, "bold"),
        ("not a result of this one-angle AB experiment.", COLORS["black"], 0.10, 10, "bold"),
    ]
    for text, color, y, size, weight in lines:
        ax.text(0.04, y, text, transform=ax.transAxes, color=color, fontsize=size, fontweight=weight)

    fig.text(
        0.5,
        0.01,
        "delta is the unlabeled deviation from pi in the one-angle readout, not an assumed external distance.",
        ha="center",
        fontsize=8,
        color=COLORS["gray"],
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    save_dual(fig, "ab_two_body_one_angle_phase_difference_scaling_v1")


def write_index() -> None:
    lines = [
        "# AB二体一角度円周位相調和読出し 主要観測図 v1",
        "",
        "## 図一覧",
        "",
        "| 図 | SVG | PNG | 内容 |",
        "|---|---|---|---|",
        "| 調和状態 | [ab_two_body_one_angle_harmonic_state_v1.svg](ab_two_body_one_angle_harmonic_state_v1.svg) | [ab_two_body_one_angle_harmonic_state_v1.png](ab_two_body_one_angle_harmonic_state_v1.png) | 閉鎖補助平面、`V_AB` 脈動、二弧読出し、`f_AB` 整合 |",
        "| Protocol 縮退 | [ab_two_body_one_angle_protocol_degeneracy_v1.svg](ab_two_body_one_angle_protocol_degeneracy_v1.svg) | [ab_two_body_one_angle_protocol_degeneracy_v1.png](ab_two_body_one_angle_protocol_degeneracy_v1.png) | `Protocol F/B` の内部表示差とラベルなし読出し縮退 |",
        "| 読出し波応答 | [ab_two_body_one_angle_readout_leak_response_v1.svg](ab_two_body_one_angle_readout_leak_response_v1.svg) | [ab_two_body_one_angle_readout_leak_response_v1.png](ab_two_body_one_angle_readout_leak_response_v1.png) | 読出し漏れによる減衰と `f_AB` 射影不整合 |",
        "| 位相差スケーリング | [ab_two_body_one_angle_phase_difference_scaling_v1.svg](ab_two_body_one_angle_phase_difference_scaling_v1.svg) | [ab_two_body_one_angle_phase_difference_scaling_v1.png](ab_two_body_one_angle_phase_difference_scaling_v1.png) | 一角度系での偏差依存性と逆冪候補の否定対照 |",
        "| 統合観測サマリー | [ab_two_body_one_angle_observation_summary_v1.svg](ab_two_body_one_angle_observation_summary_v1.svg) | [ab_two_body_one_angle_observation_summary_v1.png](ab_two_body_one_angle_observation_summary_v1.png) | 主要判定値とスイープ判定フラグ |",
        "",
        "図中の軸は読出し補助表示であり、絶対背景座標の存在を仮定しない。",
        "",
    ]
    (OUT_DIR / "ab_two_body_one_angle_observation_figures_index_v1.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    series_rows = read_csv(PRELIM_DIR / "ab_two_body_one_angle_harmonic_readout_series_v1.csv")
    protocol_rows = read_csv(PRELIM_DIR / "ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.csv")
    leak_rows = read_csv(SWEEP_DIR / "ab_two_body_one_angle_parameter_sweep_leak_summary_v1.csv")
    period_rows = read_csv(SWEEP_DIR / "ab_two_body_one_angle_parameter_sweep_period_summary_v1.csv")
    case_rows = read_csv(SWEEP_DIR / "ab_two_body_one_angle_parameter_sweep_case_summary_v1.csv")
    prelim = read_json(PRELIM_DIR / "ab_two_body_one_angle_harmonic_readout_preliminary_result_v1.json")
    sweep = read_json(SWEEP_DIR / "ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1.json")

    make_harmonic_state_figure(series_rows)
    make_protocol_degeneracy_figure(series_rows, protocol_rows)
    make_readout_leak_response_figure(leak_rows, period_rows)
    make_phase_difference_scaling_figure(case_rows)
    make_observation_summary_figure(prelim, sweep)
    write_index()

    print(json.dumps({"figure_dir": str(OUT_DIR), "figure_count": 10}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
