from __future__ import annotations

"""公開済み AB 加速度実験 v4 の保存CSVから論文用図表を再構成する。

このスクリプトは既存CSVを読み出して二階差分を再表示するだけであり、
運動更新、散乱計算、パラメータ掃引は実行しない。
"""

import csv
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = BASE_DIR.parents[1]
SOURCE_DIR = (
    REPO_ROOT
    / "波の情報読出し"
    / "20260711"
    / "ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4"
)
SERIES_PATH = SOURCE_DIR / "ab_two_body_fermionic_reflection_harmonic_series_v4.csv"
SUMMARY_PATH = SOURCE_DIR / "ab_two_body_fermionic_reflection_harmonic_case_summary_v4.csv"

FIGURE_DIR = BASE_DIR / "figures"
TABLE_DIR = BASE_DIR / "tables"
FIGURE_DIR.mkdir(exist_ok=True)
TABLE_DIR.mkdir(exist_ok=True)

MPL_DIR = BASE_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


PERIOD_STEPS = 96
OMEGA_STEP = 2.0 * math.pi / PERIOD_STEPS
OMEGA_DISCRETE_SQ = 4.0 * math.sin(OMEGA_STEP / 2.0) ** 2
SELECTED_PROTOCOLS = ["pass_through", "fermionic_reflection_pi"]
SELECTED_READOUT_MODE = "readout_off"
REPRESENTATIVE_CASE = "near_pi_05deg"
PLOT_LAST_STEP = 192

PROTOCOL_LABELS = {
    "pass_through": "pass-through",
    "fermionic_reflection_pi": "fermionic reflection",
}
PROTOCOL_STYLES = {
    "pass_through": {"color": "#1f2933", "linestyle": "-", "marker": "o"},
    "fermionic_reflection_pi": {
        "color": "#087f8c",
        "linestyle": "--",
        "marker": "x",
    },
}
DEVIATION_COLORS = {
    2.0: "#315c8a",
    5.0: "#087f8c",
    10.0: "#b77812",
    20.0: "#b54a3a",
}


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def group_selected_series(
    rows: Sequence[Dict[str, str]],
) -> Dict[Tuple[str, str], List[Dict[str, str]]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = {}
    for row in rows:
        protocol = row["scattering_protocol"]
        if protocol not in SELECTED_PROTOCOLS:
            continue
        if row["readout_mode"] != SELECTED_READOUT_MODE:
            continue
        key = (row["case_id"], protocol)
        groups.setdefault(key, []).append(row)
    for selected in groups.values():
        selected.sort(key=lambda row: int(row["step"]))
    return groups


def selected_summary_map(
    rows: Sequence[Dict[str, str]],
) -> Dict[Tuple[str, str], Dict[str, str]]:
    result: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in rows:
        protocol = row["scattering_protocol"]
        if protocol not in SELECTED_PROTOCOLS:
            continue
        if row["readout_mode"] != SELECTED_READOUT_MODE:
            continue
        result[(row["case_id"], protocol)] = row
    return result


def arrays(rows: Sequence[Dict[str, str]]) -> Dict[str, np.ndarray]:
    return {
        "step": np.array([int(row["step"]) for row in rows], dtype=int),
        "chi": np.array([float(row["chi_read"]) for row in rows], dtype=float),
        "f_center": np.array([float(row["f_AB_center"]) for row in rows], dtype=float),
        "f_circle": np.array([float(row["f_AB_circle"]) for row in rows], dtype=float),
        "q_closed": np.array([float(row["Q_closed"]) for row in rows], dtype=float),
        "unitarity_error": np.array(
            [float(row["scattering_unitarity_error"]) for row in rows], dtype=float
        ),
    }


def signed_second_difference(values: np.ndarray) -> np.ndarray:
    result = np.full_like(values, np.nan, dtype=float)
    if len(values) >= 3:
        result[1:-1] = values[2:] - 2.0 * values[1:-1] + values[:-2]
    return result


def linear_fit_with_intercept(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    prediction = slope * x + intercept
    residual = y - prediction
    total = y - np.mean(y)
    ss_res = float(np.dot(residual, residual))
    ss_tot = float(np.dot(total, total))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return float(slope), float(intercept), float(r_squared)


def build_table(
    groups: Dict[Tuple[str, str], List[Dict[str, str]]],
    summary_map: Dict[Tuple[str, str], Dict[str, str]],
) -> List[Dict[str, Any]]:
    output: List[Dict[str, Any]] = []
    ordered_keys = sorted(
        groups,
        key=lambda key: (
            float(groups[key][0]["initial_deviation_deg"]),
            SELECTED_PROTOCOLS.index(key[1]),
        ),
    )
    for key in ordered_keys:
        data = arrays(groups[key])
        summary = summary_map[key]
        second = signed_second_difference(data["chi"])
        valid = np.isfinite(second)
        slope, intercept, r_squared = linear_fit_with_intercept(
            data["chi"][valid], second[valid]
        )
        expected_slope = -OMEGA_DISCRETE_SQ
        slope_relative_error = abs(slope - expected_slope) / abs(expected_slope)
        source_projection_error = np.abs(data["f_center"] - data["f_circle"])
        output.append(
            {
                "case_id": key[0],
                "initial_deviation_deg": float(groups[key][0]["initial_deviation_deg"]),
                "scattering_protocol": key[1],
                "readout_mode": SELECTED_READOUT_MODE,
                "step_count": int(summary["step_count"]),
                "period_steps": PERIOD_STEPS,
                "sign_change_count": int(summary["sign_change_count_chi_read"]),
                "oscillation_detected": summary["oscillation_detected"],
                "max_abs_second_difference": float(np.nanmax(np.abs(second))),
                "fitted_second_difference_slope": slope,
                "fitted_intercept": intercept,
                "expected_slope_minus_omega_discrete_sq": expected_slope,
                "slope_relative_error": slope_relative_error,
                "linear_fit_R2": r_squared,
                "max_f_AB_center_circle_abs_error": float(
                    np.max(source_projection_error[1:-1])
                ),
                "max_Q_closed_abs": float(np.max(np.abs(data["q_closed"]))),
                "max_scattering_unitarity_error": float(
                    np.max(np.abs(data["unitarity_error"]))
                ),
                "envelope_ratio_final_over_initial": float(
                    summary["envelope_ratio_final_over_initial"]
                ),
                "decay_rate_V_AB": float(summary["decay_rate_V_AB"]),
            }
        )
    return output


def make_figure(
    groups: Dict[Tuple[str, str], List[Dict[str, str]]],
    table_rows: Sequence[Dict[str, Any]],
) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Hiragino Sans",
                "Yu Gothic",
                "Noto Sans CJK JP",
                "DejaVu Sans",
            ],
            "axes.unicode_minus": False,
            "font.size": 11,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 10.2))
    fig.suptitle(
        "Published AB acceleration experiment: nonzero second difference and the harmonic acceleration relation",
        fontsize=14,
        fontweight="bold",
        y=0.985,
    )

    # A/B: representative time series.
    for protocol in SELECTED_PROTOCOLS:
        rows = groups[(REPRESENTATIVE_CASE, protocol)]
        data = arrays(rows)
        mask = data["step"] <= PLOT_LAST_STEP
        style = PROTOCOL_STYLES[protocol]
        axes[0, 0].plot(
            data["step"][mask],
            data["chi"][mask],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=PROTOCOL_LABELS[protocol],
        )
        second = signed_second_difference(data["chi"])
        axes[0, 1].plot(
            data["step"][mask],
            second[mask],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=PROTOCOL_LABELS[protocol],
        )

    axes[0, 0].axhline(0.0, color="#9aa6b2", linewidth=0.8)
    axes[0, 0].set_title("A  Periodic motion of the relative position phase χ (initial deviation 5°)")
    axes[0, 0].set_xlabel("step")
    axes[0, 0].set_ylabel(r"$\chi_s$ [rad]")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].axhline(0.0, color="#9aa6b2", linewidth=0.8)
    axes[0, 1].set_title("B  Signed second difference computed directly from the stored series")
    axes[0, 1].set_xlabel("step")
    axes[0, 1].set_ylabel(r"$\Delta^2\chi_s$")
    axes[0, 1].legend(frameon=False)

    # C: all cases collapse onto the same harmonic acceleration law.
    x_all: List[float] = []
    y_all: List[float] = []
    for key, rows in sorted(groups.items()):
        data = arrays(rows)
        second = signed_second_difference(data["chi"])
        valid = np.isfinite(second)
        deviation = float(rows[0]["initial_deviation_deg"])
        protocol = key[1]
        axes[1, 0].scatter(
            data["chi"][valid][::4],
            second[valid][::4],
            s=13,
            alpha=0.46,
            color=DEVIATION_COLORS[deviation],
            marker=PROTOCOL_STYLES[protocol]["marker"],
            linewidths=0.8,
        )
        x_all.extend(data["chi"][valid].tolist())
        y_all.extend(second[valid].tolist())

    x_line = np.linspace(min(x_all), max(x_all), 400)
    axes[1, 0].plot(
        x_line,
        -OMEGA_DISCRETE_SQ * x_line,
        color="#111827",
        linewidth=2.0,
        linestyle="--",
        label=r"$\Delta^2\chi=-\omega_d^2\chi$",
    )
    axes[1, 0].axhline(0.0, color="#9aa6b2", linewidth=0.8)
    axes[1, 0].axvline(0.0, color="#9aa6b2", linewidth=0.8)
    axes[1, 0].set_title("C  4 initial deviations and 2 protocols collapse onto the same relation")
    axes[1, 0].set_xlabel(r"$\chi_s$ [rad]")
    axes[1, 0].set_ylabel(r"$\Delta^2\chi_s$")
    deviation_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor=color,
            label=f"{deviation:g}°",
        )
        for deviation, color in DEVIATION_COLORS.items()
    ]
    protocol_handles = [
        Line2D(
            [0],
            [0],
            marker=PROTOCOL_STYLES[protocol]["marker"],
            linestyle="none",
            color="#4b5563",
            label=PROTOCOL_LABELS[protocol],
        )
        for protocol in SELECTED_PROTOCOLS
    ]
    law_handle = Line2D(
        [0],
        [0],
        color="#111827",
        linestyle="--",
        linewidth=2.0,
        label=r"$-\omega_d^2\chi$",
    )
    axes[1, 0].legend(
        handles=deviation_handles + protocol_handles + [law_handle],
        ncol=2,
        frameon=False,
        fontsize=9,
    )

    # D: two independent displays of the acceleration magnitude.
    representative = arrays(groups[(REPRESENTATIVE_CASE, "pass_through")])
    mask = representative["step"] <= PLOT_LAST_STEP
    axes[1, 1].plot(
        representative["step"][mask],
        representative["f_center"][mask],
        color="#1f2933",
        linewidth=2.0,
        label=r"center display $\omega_d^2|\chi|$",
    )
    axes[1, 1].plot(
        representative["step"][mask],
        representative["f_circle"][mask],
        color="#087f8c",
        linestyle="--",
        linewidth=1.7,
        label=r"circumferential display $|\Delta^2\chi|$",
    )
    representative_table = next(
        row
        for row in table_rows
        if row["case_id"] == REPRESENTATIVE_CASE
        and row["scattering_protocol"] == "pass_through"
    )
    axes[1, 1].text(
        0.03,
        0.93,
        "max display difference = "
        f"{float(representative_table['max_f_AB_center_circle_abs_error']):.3e}",
        transform=axes[1, 1].transAxes,
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#c8d3dd"},
    )
    axes[1, 1].set_title("D  Center display agrees with the circumferential second difference (5°)")
    axes[1, 1].set_xlabel("step")
    axes[1, 1].set_ylabel(r"$f_{AB}$ / acceleration-like readout")
    axes[1, 1].legend(frameon=False)

    for ax in axes.flat:
        ax.grid(alpha=0.22)

    fig.text(
        0.5,
        0.012,
        "Source: published v4 harmonic series CSV. No new motion update, scattering calculation, or parameter sweep.",
        ha="center",
        fontsize=10,
        color="#596879",
    )
    fig.tight_layout(rect=[0.02, 0.035, 0.98, 0.955])
    stem = FIGURE_DIR / "existing_ab_acceleration_second_difference_en_v1"
    fig.savefig(stem.with_suffix(".png"), dpi=210)
    fig.savefig(stem.with_suffix(".svg"))
    plt.close(fig)


def write_markdown_table(rows: Sequence[Dict[str, Any]]) -> None:
    lines = [
        "# 公開済みAB加速度実験・二階差分集計 v1",
        "",
        "**データ源:** `波の情報読出し/20260711/ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4/`",
        "",
        "既存の harmonic series CSV を再集計した。新しい運動更新、散乱計算、パラメータ掃引は行っていない。",
        "",
        "離散加速度関係は、端点を除いて次式で評価した。",
        "",
        "$$",
        "\\Delta^2\\chi_s=\\chi_{s+1}-2\\chi_s+\\chi_{s-1}",
        "$$",
        "",
        "$$",
        "\\Delta^2\\chi_s=-\\omega_d^2\\chi_s,",
        "\\qquad",
        "\\omega_d^2=4\\sin^2\\left(\\frac{\\pi}{96}\\right)",
        "$$",
        "",
        "| 初期偏差 | プロトコル | 最大 $|\\Delta^2\\chi|$ | 回帰傾き | 理論傾き | $R^2$ | 最大表示差 | $\\max|Q_{closed}|$ |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{float(row['initial_deviation_deg']):.0f}° | "
            f"{PROTOCOL_LABELS[str(row['scattering_protocol'])]} | "
            f"{float(row['max_abs_second_difference']):.9e} | "
            f"{float(row['fitted_second_difference_slope']):.12e} | "
            f"{float(row['expected_slope_minus_omega_discrete_sq']):.12e} | "
            f"{float(row['linear_fit_R2']):.12f} | "
            f"{float(row['max_f_AB_center_circle_abs_error']):.3e} | "
            f"{float(row['max_Q_closed_abs']):.3e} |"
        )
    lines.extend(
        [
            "",
            "## 判定",
            "",
            "- 全8条件で符号付き二階差分時系列は恒等的ゼロではなく、周期的な非零振幅を持つ。",
            "- 全8条件で $\\Delta^2\\chi_s=-\\omega_d^2\\chi_s$ の回帰が成立する。",
            "- 通過型とフェルミオン反射型は同じ傾きへ一致する。",
            "- 全条件で $Q_{closed}=0$ が保存される。",
            "- この表が示すのは加速度様二階構造の存在であり、逆二乗則は調和位相条件との代数的接続として別に示す。",
            "",
        ]
    )
    (TABLE_DIR / "既存AB加速度発生集計_v1.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    series_rows = read_csv(SERIES_PATH)
    summary_rows = read_csv(SUMMARY_PATH)
    groups = group_selected_series(series_rows)
    summary_map = selected_summary_map(summary_rows)
    expected_group_count = 4 * len(SELECTED_PROTOCOLS)
    if len(groups) != expected_group_count:
        raise RuntimeError(f"unexpected selected group count: {len(groups)}")
    if set(groups) != set(summary_map):
        raise RuntimeError("series/summary key mismatch")

    table_rows = build_table(groups, summary_map)
    # 英語版は図のみ生成する。日本語の表（tables/）は上書きしない。
    make_figure(groups, table_rows)

    max_slope_relative_error = max(float(row["slope_relative_error"]) for row in table_rows)
    min_r_squared = min(float(row["linear_fit_R2"]) for row in table_rows)
    max_q_closed = max(float(row["max_Q_closed_abs"]) for row in table_rows)
    print(
        {
            "source_series": str(SERIES_PATH.relative_to(REPO_ROOT)),
            "source_summary": str(SUMMARY_PATH.relative_to(REPO_ROOT)),
            "new_simulation_performed": False,
            "table_case_count": len(table_rows),
            "expected_omega_discrete_sq": OMEGA_DISCRETE_SQ,
            "max_slope_relative_error": max_slope_relative_error,
            "min_linear_fit_R2": min_r_squared,
            "max_Q_closed_abs": max_q_closed,
        }
    )


if __name__ == "__main__":
    main()
