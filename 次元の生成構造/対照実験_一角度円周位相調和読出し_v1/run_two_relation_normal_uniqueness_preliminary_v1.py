from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_ab_two_body_one_angle_harmonic_readout_preliminary_v1 import (
    Params as OneAngleParams,
    ReadoutMode,
    closure_rotation_series,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "two_relation_normal_uniqueness_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi


@dataclass(frozen=True)
class ExperimentParams:
    ambient_dimensions: Tuple[int, ...] = (2, 3, 4, 5, 6)
    phase_offset_degs: Tuple[float, ...] = (5.0, 15.0, 30.0, 60.0)
    trial_count: int = 64
    candidate_phase_count: int = 144
    random_seed: int = 20260720
    rank_tol: float = 1.0e-12
    invariant_tol: float = 1.0e-11
    distinct_candidate_angle_deg: float = 80.0
    linear_variation_min: float = 1.5


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def random_orthogonal(rng: np.random.Generator, dimension: int) -> np.ndarray:
    matrix = rng.normal(size=(dimension, dimension))
    q, r = np.linalg.qr(matrix)
    signs = np.where(np.diag(r) >= 0.0, 1.0, -1.0)
    return q @ np.diag(signs)


def phase_seed_pair(dimension: int, phase_offset_rad: float) -> Tuple[np.ndarray, np.ndarray]:
    if dimension < 2:
        raise ValueError("dimension must be at least 2")
    relation_1 = np.zeros(dimension, dtype=float)
    relation_2 = np.zeros(dimension, dtype=float)
    relation_1[0] = 1.0
    relation_2[0] = math.cos(phase_offset_rad)
    relation_2[1] = math.sin(phase_offset_rad)
    return relation_1, relation_2


def ordered_sweep_basis(relation_1: np.ndarray, relation_2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    e1 = relation_1 / np.linalg.norm(relation_1)
    residual = relation_2 - float(np.dot(relation_2, e1)) * e1
    residual_norm = float(np.linalg.norm(residual))
    if residual_norm <= 1.0e-14:
        raise ValueError("the two relations must have a nonzero phase offset")
    e2 = residual / residual_norm
    return e1, e2


def normal_basis(
    relation_1: np.ndarray,
    relation_2: np.ndarray,
    rank_tol: float,
) -> Tuple[int, np.ndarray, np.ndarray]:
    relation_matrix = np.stack([relation_1, relation_2], axis=0)
    _, singular_values, vh = np.linalg.svd(relation_matrix, full_matrices=True)
    rank = int(np.sum(singular_values > rank_tol))
    basis = vh[rank:].T.copy()
    projector = basis @ basis.T if basis.shape[1] else np.zeros(
        (relation_1.size, relation_1.size), dtype=float
    )
    return rank, basis, projector


def projective_angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    cosine = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    cosine = max(-1.0, min(1.0, abs(cosine)))
    return math.degrees(math.acos(cosine))


def candidate_family_diagnostics(
    relation_1: np.ndarray,
    relation_2: np.ndarray,
    basis: np.ndarray,
    candidate_phase_count: int,
) -> Dict[str, float]:
    nullity = int(basis.shape[1])
    if nullity == 0:
        return {
            "candidate_sample_count": 0,
            "candidate_relation_signature_spread": 0.0,
            "candidate_projective_angle_max_deg": 0.0,
            "candidate_linear_probe_range": 0.0,
            "candidate_quadratic_norm_range": 0.0,
        }

    if nullity == 1:
        candidates = np.array([basis[:, 0]], dtype=float)
    else:
        phases = np.linspace(0.0, TAU, candidate_phase_count, endpoint=False)
        candidates = np.array(
            [math.cos(phase) * basis[:, 0] + math.sin(phase) * basis[:, 1] for phase in phases],
            dtype=float,
        )

    signatures = np.column_stack(
        [
            candidates @ relation_1,
            candidates @ relation_2,
            np.sum(candidates * candidates, axis=1),
        ]
    )
    signature_spread = float(np.max(np.ptp(signatures, axis=0))) if len(candidates) > 1 else 0.0
    reference = candidates[0]
    projective_angles = [projective_angle_deg(reference, candidate) for candidate in candidates]

    diagnostic_probe = basis[:, 0]
    linear_components = candidates @ diagnostic_probe
    quadratic_norms = np.sum(candidates * candidates, axis=1)
    return {
        "candidate_sample_count": int(len(candidates)),
        "candidate_relation_signature_spread": signature_spread,
        "candidate_projective_angle_max_deg": float(max(projective_angles)),
        "candidate_linear_probe_range": float(np.ptp(linear_components)),
        "candidate_quadratic_norm_range": float(np.ptp(quadratic_norms)),
    }


def one_angle_control_series() -> Tuple[List[complex], Dict[str, float]]:
    params = OneAngleParams()
    readout_off = ReadoutMode("readout_off", 0.0, False)
    initial_deviation_rad = math.radians(5.0)
    series = closure_rotation_series(initial_deviation_rad, readout_off, params)
    r2_values = np.array([abs(value) ** 2 for value in series], dtype=float)
    cycle_repeat_errors = [
        abs(series[index + 96] - series[index])
        for index in range(len(series) - 96)
    ]
    diagnostics = {
        "step_count": params.step_count,
        "period_steps": 96,
        "omega_step": params.omega_step,
        "initial_deviation_deg": 5.0,
        "R2_initial": float(r2_values[0]),
        "max_R2_drift": float(np.max(np.abs(r2_values - r2_values[0]))),
        "max_one_period_repeat_error": float(max(cycle_repeat_errors)),
    }
    return series, diagnostics


def trial_row(
    dimension: int,
    phase_offset_deg: float,
    trial_index: int,
    q: np.ndarray,
    one_angle_series: List[complex],
    params: ExperimentParams,
) -> Dict[str, Any]:
    phase_offset_rad = math.radians(phase_offset_deg)
    base_1, base_2 = phase_seed_pair(dimension, phase_offset_rad)
    relation_1 = q @ base_1
    relation_2 = q @ base_2

    seed_rank, basis, projector = normal_basis(relation_1, relation_2, params.rank_tol)
    _, base_basis, base_projector = normal_basis(base_1, base_2, params.rank_tol)
    _, _, swap_projector = normal_basis(relation_2, relation_1, params.rank_tol)
    expected_projector = q @ base_projector @ q.T

    e1, e2 = ordered_sweep_basis(relation_1, relation_2)
    mapped_series = np.array(
        [value.real * e1 + value.imag * e2 for value in one_angle_series],
        dtype=float,
    )
    mapped_r2 = np.sum(mapped_series * mapped_series, axis=1)
    expected_r2 = abs(one_angle_series[0]) ** 2

    nullity = int(basis.shape[1])
    if nullity:
        normal = basis[:, 0]
        xyz_rank = int(np.linalg.matrix_rank(np.stack([relation_1, relation_2, normal]), tol=params.rank_tol))
        xz_abs = abs(float(np.dot(relation_1, normal)))
        yz_abs = abs(float(np.dot(relation_2, normal)))
    else:
        xyz_rank = seed_rank
        xz_abs = 0.0
        yz_abs = 0.0

    family = candidate_family_diagnostics(
        relation_1,
        relation_2,
        basis,
        params.candidate_phase_count,
    )
    relation_angle = math.degrees(
        math.acos(
            max(
                -1.0,
                min(1.0, float(np.dot(relation_1, relation_2))),
            )
        )
    )
    return {
        "ambient_dimension": dimension,
        "phase_offset_deg": phase_offset_deg,
        "trial_index": trial_index,
        "seed_relation_rank": seed_rank,
        "normal_nullity": nullity,
        "unique_normal_line": bool(nullity == 1),
        "unique_readout_direction_count": 2 + int(nullity == 1),
        "XYZ_rank_using_one_normal_candidate": xyz_rank,
        "measured_relation_angle_deg": relation_angle,
        "relation_angle_error_deg": abs(relation_angle - phase_offset_deg),
        "relation_1_dot_normal_abs": xz_abs,
        "relation_2_dot_normal_abs": yz_abs,
        "normal_projector_trace": float(np.trace(projector)),
        "normal_projector_idempotence_error": float(np.linalg.norm(projector @ projector - projector)),
        "normal_projector_covariance_error": float(np.linalg.norm(projector - expected_projector)),
        "label_swap_projector_error": float(np.linalg.norm(projector - swap_projector)),
        "mapped_one_angle_R2_max_drift": float(np.max(np.abs(mapped_r2 - expected_r2))),
        "base_normal_nullity": int(base_basis.shape[1]),
        "absolute_background_axis_used_for_selection": False,
        "external_direction_name_used_for_selection": False,
        "diagnostic_probe_used_for_selection": False,
        **family,
    }


def summarize_dimensions(rows: List[Dict[str, Any]], params: ExperimentParams) -> List[Dict[str, Any]]:
    summaries: List[Dict[str, Any]] = []
    for dimension in params.ambient_dimensions:
        selected = [row for row in rows if int(row["ambient_dimension"]) == dimension]
        nullities = [int(row["normal_nullity"]) for row in selected]
        unique_counts = [int(row["unique_readout_direction_count"]) for row in selected]
        xyz_ranks = [int(row["XYZ_rank_using_one_normal_candidate"]) for row in selected]
        summary = {
            "ambient_dimension": dimension,
            "case_count": len(selected),
            "normal_nullity_min": min(nullities),
            "normal_nullity_max": max(nullities),
            "unique_normal_line_all": bool_all(bool(row["unique_normal_line"]) for row in selected),
            "unique_readout_direction_count_min": min(unique_counts),
            "unique_readout_direction_count_max": max(unique_counts),
            "XYZ_rank_min": min(xyz_ranks),
            "XYZ_rank_max": max(xyz_ranks),
            "max_relation_angle_error_deg": max(float(row["relation_angle_error_deg"]) for row in selected),
            "max_relation_dot_normal_abs": max(
                max(float(row["relation_1_dot_normal_abs"]), float(row["relation_2_dot_normal_abs"]))
                for row in selected
            ),
            "max_projector_idempotence_error": max(
                float(row["normal_projector_idempotence_error"]) for row in selected
            ),
            "max_projector_covariance_error": max(
                float(row["normal_projector_covariance_error"]) for row in selected
            ),
            "max_label_swap_projector_error": max(float(row["label_swap_projector_error"]) for row in selected),
            "max_mapped_one_angle_R2_drift": max(
                float(row["mapped_one_angle_R2_max_drift"]) for row in selected
            ),
            "max_candidate_relation_signature_spread": max(
                float(row["candidate_relation_signature_spread"]) for row in selected
            ),
            "min_candidate_projective_angle_max_deg": min(
                float(row["candidate_projective_angle_max_deg"]) for row in selected
            ),
            "min_candidate_linear_probe_range": min(
                float(row["candidate_linear_probe_range"]) for row in selected
            ),
            "max_candidate_quadratic_norm_range": max(
                float(row["candidate_quadratic_norm_range"]) for row in selected
            ),
        }
        summary["dimension_classification_pass"] = bool(
            (dimension == 2 and summary["normal_nullity_max"] == 0)
            or (
                dimension == 3
                and summary["normal_nullity_min"] == 1
                and summary["normal_nullity_max"] == 1
                and summary["unique_normal_line_all"]
                and summary["XYZ_rank_min"] == 3
            )
            or (
                dimension >= 4
                and summary["normal_nullity_min"] >= 2
                and not summary["unique_normal_line_all"]
                and summary["min_candidate_projective_angle_max_deg"]
                >= params.distinct_candidate_angle_deg
                and summary["min_candidate_linear_probe_range"] >= params.linear_variation_min
                and summary["max_candidate_relation_signature_spread"] <= params.invariant_tol
                and summary["max_candidate_quadratic_norm_range"] <= params.invariant_tol
            )
        )
        summaries.append(summary)
    return summaries


def aggregate_verdict(
    rows: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
    one_angle: Dict[str, float],
    params: ExperimentParams,
) -> Dict[str, Any]:
    by_dimension = {int(row["ambient_dimension"]): row for row in summaries}
    d3 = by_dimension[3]
    high_dimensions = [row for row in summaries if int(row["ambient_dimension"]) >= 4]
    max_covariance_error = max(float(row["normal_projector_covariance_error"]) for row in rows)
    max_swap_error = max(float(row["label_swap_projector_error"]) for row in rows)
    max_r2_drift = max(float(row["mapped_one_angle_R2_max_drift"]) for row in rows)
    d4plus_nonunique = bool_all(
        int(row["normal_nullity_min"]) >= 2 and not bool(row["unique_normal_line_all"])
        for row in high_dimensions
    )
    d4plus_pair_signatures_degenerate = bool_all(
        float(row["max_candidate_relation_signature_spread"]) <= params.invariant_tol
        for row in high_dimensions
    )
    d4plus_quadratic_invariant = bool_all(
        float(row["max_candidate_quadratic_norm_range"]) <= params.invariant_tol
        for row in high_dimensions
    )
    d4plus_linear_components_vary = bool_all(
        float(row["min_candidate_linear_probe_range"]) >= params.linear_variation_min
        for row in high_dimensions
    )
    return {
        "experiment": "two_relation_normal_uniqueness_preliminary_v1",
        "trial_case_count": len(rows),
        "one_angle_control_kernel_used": True,
        "one_angle_control_max_R2_drift": one_angle["max_R2_drift"],
        "one_angle_control_max_one_period_repeat_error": one_angle["max_one_period_repeat_error"],
        "two_relation_phase_offset_includes_5deg": 5.0 in params.phase_offset_degs,
        "dimension_3_unique_third_direction_all": bool(
            d3["unique_normal_line_all"] and d3["XYZ_rank_min"] == 3 and d3["XYZ_rank_max"] == 3
        ),
        "dimension_4plus_nonunique_normal_family_all": d4plus_nonunique,
        "dimension_4plus_pair_signatures_degenerate_all": d4plus_pair_signatures_degenerate,
        "dimension_4plus_linear_components_vary_all": d4plus_linear_components_vary,
        "dimension_4plus_quadratic_norm_invariant_all": d4plus_quadratic_invariant,
        "max_projector_covariance_error": max_covariance_error,
        "max_label_swap_projector_error": max_swap_error,
        "max_mapped_one_angle_R2_drift": max_r2_drift,
        "absolute_background_axis_used_for_selection": False,
        "external_direction_name_used_for_selection": False,
        "ambient_dimension_used_as_test_parameter": True,
        "physical_dimension_selection_derived": False,
        "imaginary_axis_identification_tested": False,
        "preliminary_experiment_valid": bool(
            one_angle["max_R2_drift"] <= params.invariant_tol
            and one_angle["max_one_period_repeat_error"] <= params.invariant_tol
            and bool_all(bool(row["dimension_classification_pass"]) for row in summaries)
            and max_covariance_error <= params.invariant_tol
            and max_swap_error <= params.invariant_tol
            and max_r2_drift <= params.invariant_tol
            and d4plus_pair_signatures_degenerate
            and d4plus_linear_components_vary
            and d4plus_quadratic_invariant
        ),
    }


def make_plots(
    one_angle_series: List[complex],
    one_angle: Dict[str, float],
    summaries: List[Dict[str, Any]],
    params: ExperimentParams,
) -> None:
    dimensions = [int(row["ambient_dimension"]) for row in summaries]
    nullities = [int(row["normal_nullity_min"]) for row in summaries]
    readable_counts = [int(row["unique_readout_direction_count_min"]) for row in summaries]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(dimensions, nullities, marker="o", label="normal-candidate nullity")
    ax.plot(dimensions, readable_counts, marker="s", label="uniquely readable direction count")
    ax.axvline(3, color="black", linestyle="--", linewidth=1.0, alpha=0.5)
    ax.set_xticks(dimensions)
    ax.set_xlabel("test representation dimension d")
    ax.set_ylabel("count")
    ax.set_title("Two-relation readout: unique third direction only at d=3")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "two_relation_dimension_uniqueness_v1.png", dpi=180)
    plt.close(fig)

    base_1, base_2 = phase_seed_pair(4, math.radians(5.0))
    _, basis, _ = normal_basis(base_1, base_2, params.rank_tol)
    phases = np.linspace(0.0, TAU, params.candidate_phase_count, endpoint=False)
    candidates = np.array(
        [math.cos(phase) * basis[:, 0] + math.sin(phase) * basis[:, 1] for phase in phases],
        dtype=float,
    )
    component_1 = candidates @ basis[:, 0]
    component_2 = candidates @ basis[:, 1]
    quadratic = component_1 * component_1 + component_2 * component_2
    fig, ax = plt.subplots(figsize=(9, 5))
    phase_degs = np.degrees(phases)
    ax.plot(phase_degs, component_1, label="hidden linear component 1")
    ax.plot(phase_degs, component_2, label="hidden linear component 2")
    ax.plot(phase_degs, quadratic, label="quadratic sum", linewidth=2.2)
    ax.set_xlabel("unresolved candidate phase (deg)")
    ax.set_ylabel("diagnostic readout")
    ax.set_title("d=4 unresolved normal family: linear variation and quadratic conservation")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dimension4_hidden_linear_and_quadratic_readout_v1.png", dpi=180)
    plt.close(fig)

    steps = np.arange(len(one_angle_series))
    r2_values = np.array([abs(value) ** 2 for value in one_angle_series], dtype=float)
    drift = r2_values - float(one_angle["R2_initial"])
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(steps, [value.real for value in one_angle_series], label="real phase component")
    axes[0].plot(steps, [value.imag for value in one_angle_series], label="closure complement")
    axes[1].plot(steps, drift, label="R^2 drift")
    axes[0].set_ylabel("one-angle state")
    axes[1].set_ylabel("R^2 - R^2(0)")
    axes[1].set_xlabel("step")
    axes[0].set_title("Copied one-angle kernel used by the preliminary experiment")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "copied_one_angle_kernel_R2_control_v1.png", dpi=180)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    one_angle = result["one_angle_control"]
    summaries = result["dimension_summaries"]
    lines: List[str] = [
        "# 二体位相関係からの第三方向一意性予備実験検証メモ v1",
        "",
        "## 目的",
        "",
        "無改変対照テストを通過した一角度円周位相調和読出しをスイープ核として使い、二つの無名な位相関係から第三方向を追加の特権軸なしに一意に再構成できる条件を検査した。",
        "",
        "## 検査方法",
        "",
        "二つの関係ベクトルを `u`, `v` とし、両者の位相差を `5, 15, 30, 60 deg` とした。各ケースをランダム直交変換で基底交換し、二体関係が張る行列",
        "",
        "```math",
        "A = [u^T; v^T]",
        "```",
        "",
        "の零空間を法線候補集合として計算した。候補集合への射影は、",
        "",
        "```math",
        "P_perp = I - A^T (A A^T)^{-1} A",
        "```",
        "",
        "である。選択に外部軸名は使用していない。",
        "",
        "## 一角度核の保存検算",
        "",
        f"- `step_count`: `{one_angle['step_count']}`",
        f"- `period_steps`: `{one_angle['period_steps']}`",
        f"- `max_R2_drift`: `{one_angle['max_R2_drift']:.16e}`",
        f"- `max_one_period_repeat_error`: `{one_angle['max_one_period_repeat_error']:.16e}`",
        "",
        "## 次元別結果",
        "",
        "| test d | normal nullity | unique third line | unique readable directions | XYZ rank | candidate angle min | linear range min | quadratic range max | pass |",
        "|---:|---:|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summaries:
        lines.append(
            f"| {row['ambient_dimension']} | {row['normal_nullity_min']} | {row['unique_normal_line_all']} | "
            f"{row['unique_readout_direction_count_min']} | {row['XYZ_rank_min']} | "
            f"{row['min_candidate_projective_angle_max_deg']:.8f} | "
            f"{row['min_candidate_linear_probe_range']:.8f} | "
            f"{row['max_candidate_quadratic_norm_range']:.8e} | "
            f"{row['dimension_classification_pass']} |"
        )

    lines.extend(
        [
            "",
            "## 統合判定",
            "",
            "| 量 | 値 |",
            "|---|---:|",
        ]
    )
    for key, value in verdict.items():
        lines.append(f"| `{key}` | `{value}` |")

    lines.extend(
        [
            "",
            "## 観測事実",
            "",
            "- `d=3` では二関係平面の法線候補空間が1次元となり、符号を同一視した第三方向は全ケースで一意だった。二関係と法線候補を並べた行列のランクは3だった。",
            "- `d>=4` では法線候補空間が2次元以上となり、同じ二体関係読出しを持つ互いに異なる候補が残った。",
            "- `d>=4` の候補族では、任意の診断基底に沿う一次成分は変化したが、候補成分の二乗和は数値精度内で一定だった。",
            "- 法線候補射影はランダム基底交換に対して共変であり、二関係のラベル交換でも変化しなかった。",
            "",
            "## 分類",
            "",
            "- `d=3` における第三方向の一意性: 本予備実験の数値結果。",
            "- `d>=4` における二体読出しだけからの法線選択不能: 本予備実験の数値結果。",
            "- 一次候補を選べない場合にも二乗和が読めること: 本予備実験の数値結果。",
            "- なぜ完全系が `d=3` の表示を選ぶか: 本実験では未導出。`d` は比較用パラメータとして与えた。",
            "- 選択されない候補を虚数軸と同定すること: 本実験では未検査。",
            "",
            "## 出力",
            "",
            "- `two_relation_normal_uniqueness_preliminary_result_v1.json`",
            "- `two_relation_normal_uniqueness_trials_v1.csv`",
            "- `two_relation_normal_uniqueness_dimension_summary_v1.csv`",
            "- `two_relation_dimension_uniqueness_v1.png`",
            "- `dimension4_hidden_linear_and_quadratic_readout_v1.png`",
            "- `copied_one_angle_kernel_R2_control_v1.png`",
        ]
    )
    (OUT_DIR / "two_relation_normal_uniqueness_preliminary_report_v1.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def run() -> Dict[str, Any]:
    params = ExperimentParams()
    rng = np.random.default_rng(params.random_seed)
    one_angle_series, one_angle = one_angle_control_series()
    rows: List[Dict[str, Any]] = []
    for dimension in params.ambient_dimensions:
        for phase_offset_deg in params.phase_offset_degs:
            for trial_index in range(params.trial_count):
                q = random_orthogonal(rng, dimension)
                rows.append(
                    trial_row(
                        dimension,
                        phase_offset_deg,
                        trial_index,
                        q,
                        one_angle_series,
                        params,
                    )
                )

    summaries = summarize_dimensions(rows, params)
    verdict = aggregate_verdict(rows, summaries, one_angle, params)
    result = {
        "experiment": "two_relation_normal_uniqueness_preliminary_v1",
        "params": asdict(params),
        "one_angle_control": one_angle,
        "dimension_summaries": summaries,
        "aggregate_verdict": verdict,
        "note": (
            "The ambient dimension is a comparison parameter, not a derived output. "
            "The experiment tests when two relational directions determine a unique normal line "
            "and whether unresolved first-order candidates retain a quadratic invariant."
        ),
    }
    (OUT_DIR / "two_relation_normal_uniqueness_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(OUT_DIR / "two_relation_normal_uniqueness_trials_v1.csv", rows)
    write_csv(OUT_DIR / "two_relation_normal_uniqueness_dimension_summary_v1.csv", summaries)
    make_plots(one_angle_series, one_angle, summaries, params)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
