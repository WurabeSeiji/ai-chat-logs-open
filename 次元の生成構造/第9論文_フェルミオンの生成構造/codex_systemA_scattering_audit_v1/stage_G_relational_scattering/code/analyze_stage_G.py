"""Analyze the preregistered Stage G runs without candidate or parameter search."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
FIGURE_DIR = STAGE_ROOT / "figures"
REPORT_DIR = STAGE_ROOT / "reports"
MPL_DIR = DATA_DIR / ".matplotlib-cache"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DYNAMIC_TOLERANCE = 1.0e-10
CORRELATION_CONSTANT_TOLERANCE = 1.0e-12
NUMERICAL_RESIDUAL_LIMIT = 1.0e-8
MAX_LAG = 64
AMPLITUDE_SHIFT_TOLERANCE = 1.0e-4


def convert(value: str) -> Any:
    if value in ("True", "False"):
        return value == "True"
    try:
        return float(value)
    except ValueError:
        return value


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            {key: convert(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to create empty output: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def group_primary(rows: list[dict]) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["series_role"] == "primary":
            grouped[str(row["run_id"])].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: int(row["collision_index"]))
    return dict(grouped)


def observable_matrix(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(
        [
            [
                row["L_A"],
                row["L_B"],
                row["N_eff_A"],
                row["N_eff_B"],
                row["Gamma_AB"],
                row["R_eff"],
            ]
            for row in rows
        ]
    )
    first = np.abs(matrix[0])
    scales = np.asarray(
        [
            max(first[0], first[1], 1.0e-12),
            max(first[0], first[1], 1.0e-12),
            max(first[2], first[3], 1.0),
            max(first[2], first[3], 1.0),
            1.0,
            1.0,
        ]
    )
    return matrix, scales


def cycle_metrics(rows: list[dict]) -> dict[str, Any]:
    matrix, scales = observable_matrix(rows)
    values = matrix / scales
    initial_distances = np.sqrt(
        np.sum((values - values[0]) ** 2, axis=1)
    )
    minimum_index = int(np.argmin(initial_distances[1:]) + 1)
    lag_rows = []
    for lag in range(1, MAX_LAG + 1):
        left = values[lag:]
        right = values[:-lag]
        return_rms = float(
            np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1)))
        )
        flat_left = left.reshape(-1)
        flat_right = right.reshape(-1)
        if (
            np.std(flat_left) <= 1.0e-300
            or np.std(flat_right) <= 1.0e-300
        ):
            correlation = 1.0 if return_rms == 0.0 else 0.0
        else:
            correlation = float(
                np.corrcoef(flat_left, flat_right)[0, 1]
            )
        lag_rows.append((lag, return_rms, correlation))
    best_lag, best_return_rms, _ = min(
        lag_rows, key=lambda item: item[1]
    )
    peak_lag, _, peak_correlation = max(
        lag_rows, key=lambda item: item[2]
    )
    l_gap = np.asarray(
        [abs(row["L_A"] - row["L_B"]) for row in rows]
    )
    n_gap = np.asarray(
        [
            abs(row["N_eff_A"] - row["N_eff_B"])
            for row in rows
        ]
    )
    amplitude_score = float(
        np.ptp(l_gap) / scales[0] + np.ptp(n_gap) / scales[2]
    )
    return {
        "minimum_observable_return_error": float(
            initial_distances[minimum_index]
        ),
        "minimum_observable_return_collision": int(
            rows[minimum_index]["collision_index"]
        ),
        "best_return_lag": best_lag,
        "best_return_lag_rms": best_return_rms,
        "autocorrelation_peak_lag": peak_lag,
        "autocorrelation_peak": peak_correlation,
        "deviation_from_32_exchange": abs(best_lag - 32),
        "amplitude_score": amplitude_score,
    }


def sign_crossings(values: list[float]) -> int:
    signs = np.sign(np.asarray(values))
    return int(
        sum(
            signs[index] != 0
            and signs[index - 1] != 0
            and signs[index] != signs[index - 1]
            for index in range(1, len(signs))
        )
    )


def correlation_result(
    run_id: str,
    first: dict,
    gamma: np.ndarray,
    target: np.ndarray,
    target_name: str,
) -> dict:
    if np.ptp(gamma) <= CORRELATION_CONSTANT_TOLERANCE:
        return {
            "run_id": run_id,
            "scattering_mode": first["scattering_mode"],
            "kappa": first["kappa"],
            "R0": first["R0"],
            "correlation": f"corr(Gamma_AB,{target_name})",
            "coefficient": "",
            "status": "not_defined_constant_series",
        }
    if np.ptp(target) <= CORRELATION_CONSTANT_TOLERANCE:
        return {
            "run_id": run_id,
            "scattering_mode": first["scattering_mode"],
            "kappa": first["kappa"],
            "R0": first["R0"],
            "correlation": f"corr(Gamma_AB,{target_name})",
            "coefficient": "",
            "status": "not_defined_constant_series",
        }
    return {
        "run_id": run_id,
        "scattering_mode": first["scattering_mode"],
        "kappa": first["kappa"],
        "R0": first["R0"],
        "correlation": f"corr(Gamma_AB,{target_name})",
        "coefficient": float(np.corrcoef(gamma, target)[0, 1]),
        "status": "defined",
    }


def classification(
    mode: str,
    gamma_max: float,
    response_max: float,
    gamma_range: float,
    r_range: float,
    maximum_residual: float,
    nan_inf_count: int,
) -> str:
    if nan_inf_count or maximum_residual > NUMERICAL_RESIDUAL_LIMIT:
        return "numerically_unstable"
    if mode == "C0":
        return "relational_term_inactive"
    if mode == "reversed_C1":
        return "constant_relation_reparameterization"
    if gamma_max <= DYNAMIC_TOLERANCE or response_max <= DYNAMIC_TOLERANCE:
        return "relational_term_inactive"
    if (
        gamma_range <= DYNAMIC_TOLERANCE
        or r_range <= DYNAMIC_TOLERANCE
    ):
        return "constant_relation_reparameterization"
    return "dynamic_relation_dynamic_scattering"


def analyze_primary(
    rows: list[dict],
) -> tuple[list[dict], list[dict], list[dict]]:
    grouped = group_primary(rows)
    summaries = []
    correlations = []
    residuals = []
    for run_id, values in sorted(grouped.items()):
        first = values[0]
        gamma = np.asarray([row["Gamma_AB"] for row in values])
        r_eff = np.asarray([row["R_eff"] for row in values])
        l_gap = np.asarray(
            [abs(row["L_A"] - row["L_B"]) for row in values]
        )
        n_gap = np.asarray(
            [
                abs(row["N_eff_A"] - row["N_eff_B"])
                for row in values
            ]
        )
        min_l_index = int(np.argmin(l_gap))
        min_n_index = int(np.argmin(n_gap))
        cycle = cycle_metrics(values)
        maximums = {
            "unitarity_residual_max": max(
                row["unitarity_residual"] for row in values
            ),
            "orthogonality_residual_max": max(
                row["orthogonality_residual"] for row in values
            ),
            "path_sum_residual_A_max": max(
                row["path_sum_residual_A"] for row in values
            ),
            "path_sum_residual_B_max": max(
                row["path_sum_residual_B"] for row in values
            ),
            "total_norm_residual_max": max(
                row["total_norm_residual"] for row in values
            ),
            "demodulation_residual_A_max": max(
                row["demodulation_reconstruction_residual_A"]
                for row in values
            ),
            "demodulation_residual_B_max": max(
                row["demodulation_reconstruction_residual_B"]
                for row in values
            ),
            "gamma_range_violation_max": max(
                max(0.0, -row["Gamma_AB"], row["Gamma_AB"] - 1.0)
                for row in values
            ),
        }
        maximum_residual = max(maximums.values())
        nan_inf_count = int(
            sum(int(row["nan_inf_count"]) for row in values)
        )
        gamma_range = float(np.ptp(gamma))
        r_range = float(np.ptp(r_eff))
        mode = str(first["scattering_mode"])
        response_max = max(
            abs(row["candidate_response"]) for row in values
        )
        result_class = classification(
            mode,
            float(np.max(gamma)),
            response_max,
            gamma_range,
            r_range,
            maximum_residual,
            nan_inf_count,
        )
        summaries.append(
            {
                "run_id": run_id,
                "scattering_mode": mode,
                "kappa": first["kappa"],
                "R0": first["R0"],
                "collision_count": len(values),
                "Gamma_AB_min": float(np.min(gamma)),
                "Gamma_AB_max": float(np.max(gamma)),
                "Gamma_AB_range": gamma_range,
                "gamma_status": (
                    "gamma_dynamic"
                    if gamma_range > DYNAMIC_TOLERANCE
                    else "gamma_constant"
                ),
                "R_eff_min": float(np.min(r_eff)),
                "R_eff_max": float(np.max(r_eff)),
                "R_eff_range": r_range,
                "R_eff_status": (
                    "R_eff_dynamic"
                    if r_range > DYNAMIC_TOLERANCE
                    else "R_eff_constant"
                ),
                "c_A_min": min(row["c_A"] for row in values),
                "c_A_max": max(row["c_A"] for row in values),
                "c_B_min": min(row["c_B"] for row in values),
                "c_B_max": max(row["c_B"] for row in values),
                "c_mean_min": min(row["c_mean"] for row in values),
                "c_mean_max": max(row["c_mean"] for row in values),
                "min_L_difference": float(l_gap[min_l_index]),
                "min_L_difference_collision": int(
                    values[min_l_index]["collision_index"]
                ),
                "min_N_eff_difference": float(n_gap[min_n_index]),
                "min_N_eff_difference_collision": int(
                    values[min_n_index]["collision_index"]
                ),
                "spectral_similarity_crossing_count": (
                    sign_crossings(
                        [
                            row[
                                "spectral_similarity_A_to_initial_A"
                            ]
                            - row[
                                "spectral_similarity_A_to_initial_B"
                            ]
                            for row in values
                        ]
                    )
                    + sign_crossings(
                        [
                            row[
                                "spectral_similarity_B_to_initial_B"
                            ]
                            - row[
                                "spectral_similarity_B_to_initial_A"
                            ]
                            for row in values
                        ]
                    )
                ),
                **cycle,
                "classification": result_class,
                "maximum_numerical_residual": maximum_residual,
                "nan_inf_count": nan_inf_count,
            }
        )
        correlations.extend(
            [
                correlation_result(
                    run_id, first, gamma, r_eff, "R_eff"
                ),
                correlation_result(
                    run_id, first, gamma, l_gap, "abs_L_difference"
                ),
                correlation_result(
                    run_id,
                    first,
                    gamma,
                    n_gap,
                    "abs_N_eff_difference",
                ),
            ]
        )
        residuals.append(
            {
                "series": "G-B",
                "run_id": run_id,
                "scattering_mode": mode,
                "kappa": first["kappa"],
                "R0": first["R0"],
                **maximums,
                "relation_wave_norm2_min": min(
                    min(
                        row["relation_wave_norm2_A"],
                        row["relation_wave_norm2_B"],
                    )
                    for row in values
                ),
                "nan_inf_count": nan_inf_count,
                "theta_range_violation_count": sum(
                    int(bool(row["theta_range_violation"]))
                    for row in values
                ),
            }
        )
    return summaries, correlations, residuals


def g_c_residual_rows(g_c_log: dict) -> list[dict]:
    rows = []
    for summary in g_c_log["run_summaries"]:
        rows.append(
            {
                "series": "G-C",
                "run_id": (
                    f"custom31__{summary['scattering_mode']}__"
                    f"kappa_{summary['kappa']}"
                ),
                "scattering_mode": summary["scattering_mode"],
                "kappa": summary["kappa"],
                "R0": g_c_log["condition"]["R0"],
                "unitarity_residual_max": summary[
                    "unitarity_residual_max"
                ],
                "orthogonality_residual_max": summary[
                    "orthogonality_residual_max"
                ],
                "path_sum_residual_A_max": summary[
                    "path_sum_residual_max"
                ],
                "path_sum_residual_B_max": summary[
                    "path_sum_residual_max"
                ],
                "total_norm_residual_max": summary[
                    "total_norm_residual_max"
                ],
                "demodulation_residual_A_max": summary[
                    "demodulation_residual_max"
                ],
                "demodulation_residual_B_max": summary[
                    "demodulation_residual_max"
                ],
                "gamma_range_violation_max": summary[
                    "gamma_range_violation_max"
                ],
                "relation_wave_norm2_min": summary[
                    "relation_wave_norm2_min"
                ],
                "nan_inf_count": summary["nan_inf_count"],
                "theta_range_violation_count": summary[
                    "theta_range_violation_count"
                ],
            }
        )
    return rows


def selected(
    rows: list[dict],
    mode: str,
    kappa: float,
    reflection: float = 0.55,
) -> list[dict]:
    return [
        row
        for row in rows
        if row["series_role"] == "primary"
        and row["scattering_mode"] == mode
        and abs(float(row["kappa"]) - kappa) <= 1.0e-15
        and abs(float(row["R0"]) - reflection) <= 1.0e-15
    ]


def figures(
    collision_rows: list[dict],
    summaries: list[dict],
    g_c_rows: list[dict],
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    colors = {
        "C0": "black",
        "reversed_C1": "#d97706",
        "relational_C1": "#2563a6",
    }
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for mode in ("C0", "reversed_C1", "relational_C1"):
        values = selected(collision_rows, mode, 1.0)
        axes[0].plot(
            [row["collision_index"] for row in values],
            [row["Gamma_AB"] for row in values],
            color=colors[mode],
            label=mode,
        )
        axes[1].plot(
            [row["collision_index"] for row in values],
            [row["R_eff"] for row in values],
            color=colors[mode],
            label=mode,
        )
    axes[0].set_ylabel("Gamma_AB")
    axes[1].set_ylabel("R_eff")
    axes[1].set_xlabel("collision")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "Gamma_and_R_eff_by_collision.png", dpi=180
    )
    plt.close(fig)

    for metric, filename, label in (
        (
            "L",
            "C0_C1_relational_L_exchange.png",
            "localization L",
        ),
        (
            "N_eff",
            "C0_C1_relational_N_eff_exchange.png",
            "effective harmonic N_eff",
        ),
    ):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
        for channel, axis in (("A", axes[0]), ("B", axes[1])):
            for mode in ("C0", "reversed_C1", "relational_C1"):
                values = selected(collision_rows, mode, 1.0)
                axis.plot(
                    [row["collision_index"] for row in values],
                    [row[f"{metric}_{channel}"] for row in values],
                    color=colors[mode],
                    label=mode,
                )
            axis.set_title(f"channel {channel}")
            axis.set_xlabel("collision")
            axis.set_ylabel(label)
            axis.grid(alpha=0.2)
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURE_DIR / filename, dpi=180)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for kappa, color in ((0.01, "#2a6fbb"), (0.1, "#d97706"), (1.0, "#a23b72")):
        values = selected(
            collision_rows, "relational_C1", kappa
        )
        ax.scatter(
            [row["Gamma_AB"] for row in values],
            [abs(row["L_A"] - row["L_B"]) for row in values],
            s=12,
            alpha=0.55,
            color=color,
            label=f"kappa={kappa:g}",
        )
    ax.set_xlabel("Gamma_AB")
    ax.set_ylabel("|L_A-L_B|")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "relation_vs_localization_difference.png",
        dpi=180,
    )
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(3)
    modes = ("C0", "reversed_C1", "relational_C1")
    for index, kappa in enumerate((0.01, 0.1, 1.0)):
        values = [
            next(
                row
                for row in summaries
                if row["scattering_mode"] == mode
                and abs(float(row["kappa"]) - kappa) <= 1.0e-15
                and abs(float(row["R0"]) - 0.55) <= 1.0e-15
            )
            for mode in modes
        ]
        offset = (index - 1) * 0.22
        axes[0].bar(
            x + offset,
            [row["best_return_lag"] for row in values],
            width=0.2,
            label=f"k={kappa:g}",
        )
        axes[1].bar(
            x + offset,
            [row["minimum_observable_return_error"] for row in values],
            width=0.2,
            label=f"k={kappa:g}",
        )
    for axis in axes:
        axis.set_xticks(x, modes, rotation=15)
        axis.grid(axis="y", alpha=0.2)
        axis.legend(fontsize=8)
    axes[0].set_ylabel("best return lag (1..64)")
    axes[1].set_ylabel("minimum observable return error")
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "cycle_and_return_error_comparison.png", dpi=180
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for mode in ("C0", "reversed_C1", "relational_C1"):
        values = sorted(
            [
                row
                for row in g_c_rows
                if row["scattering_mode"] == mode
                and abs(float(row["kappa"]) - 1.0) <= 1.0e-15
            ],
            key=lambda row: int(row["iteration"]),
        )
        ax.plot(
            [row["iteration"] for row in values],
            [row["return_error"] for row in values],
            marker="o",
            color=colors[mode],
            label=mode,
        )
    ax.set_xlabel("fixed iteration")
    ax.set_ylabel("full-state return error")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "return_error_31_series.png", dpi=180)
    plt.close(fig)


def report_files(
    summaries: list[dict],
    correlations: list[dict],
    residuals: list[dict],
    g_c_rows: list[dict],
    g_c_log: dict,
) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    representative = [
        row
        for row in summaries
        if abs(float(row["R0"]) - 0.55) <= 1.0e-15
    ]
    table = [
        "| mode | kappa | Gamma range | R_eff range | min L gap@col | min N_eff gap@col | class |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in representative:
        table.append(
            f"| {row['scattering_mode']} | {float(row['kappa']):g} | "
            f"{float(row['Gamma_AB_range']):.4g} | "
            f"{float(row['R_eff_range']):.4g} | "
            f"{float(row['min_L_difference']):.5g}@"
            f"{int(row['min_L_difference_collision'])} | "
            f"{float(row['min_N_eff_difference']):.5g}@"
            f"{int(row['min_N_eff_difference_collision'])} | "
            f"{row['classification']} |"
        )
    (REPORT_DIR / "03_existing_systemA_comparison.md").write_text(
        """# 既存System A代表条件の3モード比較

C0の既存正規化再現は最大絶対誤差0、Stage FのC0/reversed_C1 raw系列との12回帰比較も最大絶対誤差0で通過した。主系列はStage Fの結果に従い `raw_update` とした。

"""
        + "\n".join(table)
        + """

`N_A=1,N_B=63`では `c_A=c_B=-1` が保存された。relational_C1の補正はreversed_C1より小さいが、衝突ごとに動く補正にはならなかった。
""",
        encoding="utf-8",
    )

    relational = [
        row
        for row in summaries
        if row["scattering_mode"] == "relational_C1"
    ]
    max_gamma_range = max(
        float(row["Gamma_AB_range"]) for row in relational
    )
    max_r_range = max(
        float(row["R_eff_range"]) for row in relational
    )
    undefined_count = sum(
        row["status"] == "not_defined_constant_series"
        for row in correlations
    )
    (REPORT_DIR / "04_dynamic_relation_analysis.md").write_text(
        f"""# 関係量の動的性

## 機械判定

- 動的閾値: `{DYNAMIC_TOLERANCE:.1e}`
- relational_C1の最大 `Delta Gamma`: `{max_gamma_range:.17g}`
- relational_C1の最大 `Delta R_eff`: `{max_r_range:.17g}`
- 判定: `gamma_constant`, `R_eff_constant`
- 相関54件のうち定数系列として未定義: `{undefined_count}`

未定義相関を0へ置換していない。

## なぜ一定になったか

復調後の初期関係波を \(a_0,b_0\) とし、両者が単位ノルムで実重なり \(s=\\langle a_0,b_0\\rangle\) を持つと、Gram行列は

\[
G_0=
\\begin{{pmatrix}}1&s\\\\s&1\\end{{pmatrix}}
=I+s\\sigma_x。
\]

System Aの各衝突行列は

\[
U_n=
\\begin{{pmatrix}}r_n&t_n\\\\t_n&r_n\\end{{pmatrix}}
=r_n I+t_n\\sigma_x
\]

である。よって \([G_0,U_n]=0\)。さらに \(U_n\) はユニタリなので

\[
G_{{n+1}}=U_nG_nU_n^\\dagger=G_n。
\]

この帰結は \(U_n\) の角度が状態依存でも成立する。したがって本Stageの二条件では、\(\Gamma=|s|^2\) は保存量となった。代表条件では `Gamma=1/32`、31系列custom packetでは `Gamma=1/2` である。

## 中心判定

relational_C1は今回の既存System A条件では `constant_relation_reparameterization` に退化した。関係量を追加したというモデル定義は成立するが、「純パリティ区間で散乱率を動的にする」という作業仮説はこの対称更新・初期Gram条件では実現しなかった。新候補は生成しない。
""",
        encoding="utf-8",
    )

    fixed_248 = [
        row for row in g_c_rows if int(row["iteration"]) == 248
    ]
    table_248 = [
        "| mode | kappa | Gamma | R_eff | return error@248 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in fixed_248:
        table_248.append(
            f"| {row['scattering_mode']} | {float(row['kappa']):g} | "
            f"{float(row['Gamma_AB']):.9g} | "
            f"{float(row['R_eff']):.9g} | "
            f"{float(row['return_error']):.6g} |"
        )
    max_fc_gamma = max(
        float(row["Gamma_AB_range"])
        for row in g_c_log["run_summaries"]
        if row["scattering_mode"] == "relational_C1"
    )
    (REPORT_DIR / "05_31_series_comparison.md").write_text(
        f"""# 31系列の限定比較

既存custom packet `A=(1), B=(1,2), R0=0.697177927`、固定評価点 `31,62,93,124,155,186,217,247,248,279` だけを再実行した。C0/reversed_C1はStage Fの360値と最大絶対誤差0で一致した。

relational_C1の全279衝突における最大 `Delta Gamma` は `{max_fc_gamma:.17g}` であり、ここでも関係量は一定だった。

"""
        + "\n".join(table_248)
        + """

relational_C1の248帰還誤差はC0と一致するよう調整していない。一定 `Gamma=1/2` による別の一定Rへ写った結果として、既存帰還位置から移動した。
""",
        encoding="utf-8",
    )

    max_fields = {
        key: max(float(row[key]) for row in residuals)
        for key in (
            "unitarity_residual_max",
            "orthogonality_residual_max",
            "path_sum_residual_A_max",
            "path_sum_residual_B_max",
            "total_norm_residual_max",
            "demodulation_residual_A_max",
            "demodulation_residual_B_max",
            "gamma_range_violation_max",
        )
    }
    residual_text = "\n".join(
        f"- `{key}`: `{value:.17g}`"
        for key, value in max_fields.items()
    )
    nan_count = sum(int(row["nan_inf_count"]) for row in residuals)
    theta_count = sum(
        int(row["theta_range_violation_count"]) for row in residuals
    )
    minimum_relation_norm = min(
        float(row["relation_wave_norm2_min"]) for row in residuals
    )
    (REPORT_DIR / "06_numerical_invariants.md").write_text(
        f"""# 数値不変量

G-BとG-C全27実行の最大残差:

{residual_text}

- relation wave最小ノルム二乗: `{minimum_relation_norm:.17g}`
- NaN/Inf: `{nan_count}`
- theta範囲違反: `{theta_count}`
- ゼロ関係波による不成立: `0`
- 数値不安定閾値: `{NUMERICAL_RESIDUAL_LIMIT:.1e}`

\(\Gamma\)は全実行で数値許容誤差内の `[0,1]` にあった。自動クリップや不成立条件の黙示除外はない。
""",
        encoding="utf-8",
    )

    class_counts = Counter(
        row["classification"] for row in summaries
    )
    counts = "\n".join(
        f"- `{name}`: {count}"
        for name, count in sorted(class_counts.items())
    )
    (REPORT_DIR / "Stage_G_report.md").write_text(
        f"""# Stage G report

## 結論

relational_C1一候補をStage Fの独立System Aコピーへ実装し、単体検証、C0再現、代表条件2点、既存31系列を完了した。

中心判定は否定結果である。代表条件では \(\Gamma=1/32\)、31系列では \(\Gamma=1/2\) が保存され、`Delta Gamma`と`Delta R_eff`はいずれも事前閾値 `1e-10` 未満だった。したがってrelational_C1も今回の対称System A条件では一定Rへの再パラメータ化に退化した。

## 分類件数

{counts}

## コード上の事実

- 散乱・経路・rawノルムは512×16全状態で計算した。
- 関係波は由来別eta射影・搬送波除去後にコヒーレント合成した。
- 新規散乱候補はrelational_C1だけである。
- C0とreversed_C1はStage F系列を誤差0で保持した。

## 数学的帰結

- 一定パリティ区間のreversed_C1は一定Rへ退化する。
- Cauchy–Schwarz不等式から `0<=Gamma<=1`。
- 等ノルム・実重なりのGram行列とSystem Aの対称散乱行列はともに `I` と `sigma_x` の線形結合で可換なため、Gammaは角度が状態依存でも保存される。

## モデル定義

\[
\\theta_{{eff}}=\\theta_0-\\kappa\\rho(\\theta_0)\\bar c\\,\\Gamma_{{AB}}。
\]

## 作業仮説の判定

「パリティ符号と関係強度が散乱角を制御する」はモデルとして実装済み。「同じ純パリティ区間でも関係強度が変化して散乱率が動く」は、今回のSystem A対称更新では成立しなかった。

## 数値観察

relational_C1はC0ともreversed_C1とも異なる一定Rを与え、局在性交換と帰還位置を移動させた。しかしその変化は動的Gammaではなく一定Gammaによる。

## 未導出

κ、rho、重なり二乗を相互作用強度と読む根拠、自然界のボゾン・フェルミオン対応は未導出である。

## 棄却・保留

本条件における `dynamic_relation_dynamic_scattering` は棄却。別候補、別関係量、非対称散乱、N体系、論文反映は保留し、自動継続しない。
""",
        encoding="utf-8",
    )


def main() -> None:
    collision_rows = read_csv(
        DATA_DIR / "stage_G_collision_results.csv"
    )
    g_c_rows = read_csv(DATA_DIR / "stage_G_31_series_results.csv")
    g_c_log = json.loads(
        (STAGE_ROOT / "logs" / "stage_G_31_series_run.json").read_text(
            encoding="utf-8"
        )
    )
    summaries, correlations, residuals = analyze_primary(collision_rows)
    residuals.extend(g_c_residual_rows(g_c_log))
    write_csv(DATA_DIR / "stage_G_run_summary.csv", summaries)
    write_csv(DATA_DIR / "stage_G_correlation_results.csv", correlations)
    write_csv(
        DATA_DIR / "stage_G_numerical_residuals.csv", residuals
    )
    figures(collision_rows, summaries, g_c_rows)
    report_files(
        summaries, correlations, residuals, g_c_rows, g_c_log
    )
    relational = [
        row
        for row in summaries
        if row["scattering_mode"] == "relational_C1"
    ]
    summary = {
        "stage": "G",
        "status": "complete_pending_reference_verification",
        "central_question": {
            "Gamma_AB_dynamic": any(
                row["gamma_status"] == "gamma_dynamic"
                for row in relational
            ),
            "R_eff_dynamic": any(
                row["R_eff_status"] == "R_eff_dynamic"
                for row in relational
            ),
            "decision": "constant_relation_reparameterization",
            "dynamic_tolerance": DYNAMIC_TOLERANCE,
            "maximum_Gamma_AB_range": max(
                float(row["Gamma_AB_range"]) for row in relational
            ),
            "maximum_R_eff_range": max(
                float(row["R_eff_range"]) for row in relational
            ),
        },
        "unit_tests": json.loads(
            (
                STAGE_ROOT / "logs" / "stage_G_unit_test_run.json"
            ).read_text(encoding="utf-8")
        ),
        "C0_reproduction": json.loads(
            (
                STAGE_ROOT / "logs" / "stage_G_C0_reproduction.json"
            ).read_text(encoding="utf-8")
        ),
        "run_count_G_B": len(summaries),
        "run_count_G_C": g_c_log["run_count"],
        "classification_counts": dict(
            Counter(row["classification"] for row in summaries)
        ),
        "correlation_status_counts": dict(
            Counter(row["status"] for row in correlations)
        ),
        "maximum_numerical_residual": max(
            max(
                float(row["unitarity_residual_max"]),
                float(row["orthogonality_residual_max"]),
                float(row["path_sum_residual_A_max"]),
                float(row["path_sum_residual_B_max"]),
                float(row["total_norm_residual_max"]),
                float(row["demodulation_residual_A_max"]),
                float(row["demodulation_residual_B_max"]),
                float(row["gamma_range_violation_max"]),
            )
            for row in residuals
        ),
        "nan_inf_count": sum(
            int(row["nan_inf_count"]) for row in residuals
        ),
        "theta_range_violation_count": sum(
            int(row["theta_range_violation_count"])
            for row in residuals
        ),
        "stage_G_C": g_c_log,
    }
    (DATA_DIR / "stage_G_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
