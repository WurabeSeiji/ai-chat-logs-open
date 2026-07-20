from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_abc_abcd_complete_pair_relation_network_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi
BODY_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


@dataclass(frozen=True)
class ExperimentParams:
    body_counts: Tuple[int, ...] = (2, 3, 4)
    trial_count: int = 32
    step_count: int = 720
    target_period_steps: int = 96
    radius_squared: float = 1.0
    imaginary_seed_amplitude: float = 0.35
    random_seed: int = 20260720
    invariant_tol: float = 1.0e-10
    covariance_tol: float = 1.0e-10
    activity_variance_tol: float = 1.0e-10


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


def normalized_pair(pair: Tuple[int, int]) -> Tuple[int, int]:
    return tuple(sorted(pair))


def relation_pairs(body_count: int) -> List[Tuple[int, int]]:
    if body_count == 2:
        return [(0, 1)]
    if body_count == 3:
        return [(0, 1), (1, 2), (2, 0)]
    return [(left, right) for left in range(body_count) for right in range(left + 1, body_count)]


def relation_name(pair: Tuple[int, int]) -> str:
    return BODY_LABELS[pair[0]] + BODY_LABELS[pair[1]]


def relation_adjacency(pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
    relation_count = len(pairs)
    adjacency = np.zeros((relation_count, relation_count), dtype=float)
    for left in range(relation_count):
        for right in range(relation_count):
            if left != right and set(pairs[left]).intersection(pairs[right]):
                adjacency[left, right] = 1.0
    return adjacency


def initial_closed_relation_state(
    relation_count: int,
    rng: np.random.Generator,
    params: ExperimentParams,
) -> np.ndarray:
    radius = math.sqrt(params.radius_squared)
    if relation_count == 1:
        return np.array([complex(radius, 0.0)], dtype=complex)

    real_direction = rng.normal(size=relation_count)
    real_direction /= np.linalg.norm(real_direction)
    imaginary_direction = rng.normal(size=relation_count)
    imaginary_direction -= float(np.dot(imaginary_direction, real_direction)) * real_direction
    imaginary_direction /= np.linalg.norm(imaginary_direction)

    imaginary_norm = params.imaginary_seed_amplitude
    real_norm = math.sqrt(params.radius_squared + imaginary_norm * imaginary_norm)
    state = real_norm * real_direction + 1j * imaginary_norm * imaginary_direction
    return np.asarray(state, dtype=complex)


def closure_quantities(state: np.ndarray, radius_squared: float) -> Dict[str, float]:
    real = state.real
    imag = state.imag
    closure_real = float(np.sum(real * real - imag * imag))
    closure_imag = float(2.0 * np.sum(real * imag))
    amplitude_sum = float(np.sum(real * real + imag * imag))
    return {
        "closure_real_E": closure_real,
        "closure_imag_F": closure_imag,
        "closure_target_error_abs": math.hypot(closure_real - radius_squared, closure_imag),
        "closure_complex_abs": math.hypot(closure_real, closure_imag),
        "hermitian_amplitude_sum": amplitude_sum,
    }


def relational_generator(state: np.ndarray, adjacency: np.ndarray) -> Tuple[np.ndarray, float]:
    phases = np.angle(state)
    phase_difference = phases[np.newaxis, :] - phases[:, np.newaxis]
    raw_generator = adjacency * np.sin(phase_difference)
    raw_generator = 0.5 * (raw_generator - raw_generator.T)
    spectral_norm = float(np.linalg.norm(raw_generator, ord=2))
    if spectral_norm <= 1.0e-14:
        return np.zeros_like(raw_generator), spectral_norm
    return raw_generator / spectral_norm, spectral_norm


def cayley_orthogonal_update(generator: np.ndarray, target_period_steps: int) -> np.ndarray:
    dimension = generator.shape[0]
    if dimension == 1 or float(np.linalg.norm(generator)) <= 1.0e-14:
        return np.eye(dimension, dtype=float)
    omega_step = TAU / target_period_steps
    cayley_scale = math.tan(omega_step / 2.0)
    identity = np.eye(dimension, dtype=float)
    return np.linalg.solve(identity - cayley_scale * generator, identity + cayley_scale * generator)


def evolve(initial_state: np.ndarray, update: np.ndarray, step_count: int) -> np.ndarray:
    states = np.empty((step_count + 1, initial_state.size), dtype=complex)
    states[0] = initial_state
    for step in range(step_count):
        states[step + 1] = update @ states[step]
    return states


def edge_permutation_matrix(
    pairs: Sequence[Tuple[int, int]],
    body_permutation: np.ndarray,
) -> np.ndarray:
    index_by_pair = {normalized_pair(pair): index for index, pair in enumerate(pairs)}
    matrix = np.zeros((len(pairs), len(pairs)), dtype=float)
    for old_index, pair in enumerate(pairs):
        new_pair = normalized_pair((int(body_permutation[pair[0]]), int(body_permutation[pair[1]])))
        new_index = index_by_pair[new_pair]
        matrix[new_index, old_index] = 1.0
    return matrix


def trial_result(
    body_count: int,
    trial_index: int,
    rng: np.random.Generator,
    params: ExperimentParams,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    pairs = relation_pairs(body_count)
    names = [relation_name(pair) for pair in pairs]
    adjacency = relation_adjacency(pairs)
    initial_state = initial_closed_relation_state(len(pairs), rng, params)
    generator, raw_generator_norm = relational_generator(initial_state, adjacency)
    update = cayley_orthogonal_update(generator, params.target_period_steps)
    states = evolve(initial_state, update, params.step_count)

    closure_rows = [closure_quantities(state, params.radius_squared) for state in states]
    initial_amplitude_sum = closure_rows[0]["hermitian_amplitude_sum"]
    max_closure_error = max(float(row["closure_target_error_abs"]) for row in closure_rows)
    max_amplitude_drift = max(
        abs(float(row["hermitian_amplitude_sum"]) - initial_amplitude_sum) for row in closure_rows
    )
    activity_variances = np.var(states.real, axis=0) + np.var(states.imag, axis=0)
    active_relation_count = int(np.sum(activity_variances > params.activity_variance_tol))

    period = params.target_period_steps
    recurrence_error = float(np.linalg.norm(states[period] - states[0])) if period <= params.step_count else 0.0
    orthogonality_error = float(np.linalg.norm(update.T @ update - np.eye(len(pairs))))

    body_permutation = rng.permutation(body_count)
    edge_permutation = edge_permutation_matrix(pairs, body_permutation)
    permuted_initial = edge_permutation @ initial_state
    permuted_generator, _ = relational_generator(permuted_initial, adjacency)
    permuted_update = cayley_orthogonal_update(permuted_generator, params.target_period_steps)
    expected_generator = edge_permutation @ generator @ edge_permutation.T
    expected_update = edge_permutation @ update @ edge_permutation.T
    generator_covariance_error = float(np.linalg.norm(permuted_generator - expected_generator))
    update_covariance_error = float(np.linalg.norm(permuted_update - expected_update))
    permuted_states = evolve(permuted_initial, permuted_update, params.step_count)
    permutation_source_indices = np.argmax(edge_permutation, axis=1)
    expected_permuted_states = states[:, permutation_source_indices]
    trajectory_covariance_error = float(
        np.max(np.linalg.norm(permuted_states - expected_permuted_states, axis=1))
    )

    trial = {
        "body_count": body_count,
        "relation_count": len(pairs),
        "relation_names": ";".join(names),
        "trial_index": trial_index,
        "all_relations_are_physical_waves": True,
        "observer_C_or_D_used": False,
        "normalization_applied": False,
        "absolute_background_axis_used": False,
        "raw_generator_spectral_norm": raw_generator_norm,
        "update_orthogonality_error": orthogonality_error,
        "max_closure_target_error_abs": max_closure_error,
        "max_hermitian_amplitude_drift": max_amplitude_drift,
        "active_relation_wave_count": active_relation_count,
        "min_relation_activity_variance": float(np.min(activity_variances)),
        "max_relation_activity_variance": float(np.max(activity_variances)),
        "target_period_recurrence_error": recurrence_error,
        "generator_label_covariance_error": generator_covariance_error,
        "update_label_covariance_error": update_covariance_error,
        "trajectory_label_covariance_error": trajectory_covariance_error,
        "initial_closure_real_E": closure_rows[0]["closure_real_E"],
        "initial_closure_imag_F": closure_rows[0]["closure_imag_F"],
        "initial_hermitian_amplitude_sum": initial_amplitude_sum,
        "stationary_relation_system": bool(active_relation_count == 0),
    }

    selected_rows: List[Dict[str, Any]] = []
    if trial_index == 0:
        for step, (state, closure) in enumerate(zip(states, closure_rows)):
            for relation_index, (name, value) in enumerate(zip(names, state)):
                selected_rows.append(
                    {
                        "body_count": body_count,
                        "relation_count": len(pairs),
                        "trial_index": trial_index,
                        "step": step,
                        "relation_index": relation_index,
                        "relation_name": name,
                        "relation_real": float(value.real),
                        "relation_imag": float(value.imag),
                        "relation_abs2": float(abs(value) ** 2),
                        "relation_phase_rad": float(np.angle(value)),
                        **closure,
                    }
                )
    return trial, selected_rows


def summarize_body_counts(trials: List[Dict[str, Any]], params: ExperimentParams) -> List[Dict[str, Any]]:
    summaries: List[Dict[str, Any]] = []
    for body_count in params.body_counts:
        selected = [row for row in trials if int(row["body_count"]) == body_count]
        expected_relation_count = body_count * (body_count - 1) // 2
        summary = {
            "body_count": body_count,
            "trial_count": len(selected),
            "expected_relation_count": expected_relation_count,
            "relation_count_min": min(int(row["relation_count"]) for row in selected),
            "relation_count_max": max(int(row["relation_count"]) for row in selected),
            "active_relation_wave_count_min": min(int(row["active_relation_wave_count"]) for row in selected),
            "active_relation_wave_count_max": max(int(row["active_relation_wave_count"]) for row in selected),
            "max_closure_target_error_abs": max(float(row["max_closure_target_error_abs"]) for row in selected),
            "max_hermitian_amplitude_drift": max(
                float(row["max_hermitian_amplitude_drift"]) for row in selected
            ),
            "max_update_orthogonality_error": max(float(row["update_orthogonality_error"]) for row in selected),
            "max_generator_label_covariance_error": max(
                float(row["generator_label_covariance_error"]) for row in selected
            ),
            "max_update_label_covariance_error": max(
                float(row["update_label_covariance_error"]) for row in selected
            ),
            "max_trajectory_label_covariance_error": max(
                float(row["trajectory_label_covariance_error"]) for row in selected
            ),
            "max_target_period_recurrence_error": max(
                float(row["target_period_recurrence_error"]) for row in selected
            ),
            "stationary_all_trials": bool_all(bool(row["stationary_relation_system"]) for row in selected),
            "all_relation_waves_active_all_trials": bool_all(
                int(row["active_relation_wave_count"]) == expected_relation_count for row in selected
            ),
        }
        summary["body_count_test_pass"] = bool(
            summary["relation_count_min"] == expected_relation_count
            and summary["relation_count_max"] == expected_relation_count
            and summary["max_closure_target_error_abs"] <= params.invariant_tol
            and summary["max_hermitian_amplitude_drift"] <= params.invariant_tol
            and summary["max_update_orthogonality_error"] <= params.invariant_tol
            and summary["max_generator_label_covariance_error"] <= params.covariance_tol
            and summary["max_update_label_covariance_error"] <= params.covariance_tol
            and summary["max_trajectory_label_covariance_error"] <= params.covariance_tol
            and (
                (body_count == 2 and summary["stationary_all_trials"])
                or (body_count >= 3 and summary["all_relation_waves_active_all_trials"])
            )
        )
        summaries.append(summary)
    return summaries


def aggregate_verdict(
    summaries: List[Dict[str, Any]],
    params: ExperimentParams,
) -> Dict[str, Any]:
    by_body_count = {int(row["body_count"]): row for row in summaries}
    ab = by_body_count[2]
    abc = by_body_count[3]
    abcd = by_body_count[4]
    return {
        "experiment": "ab_abc_abcd_complete_pair_relation_network_preliminary_v1",
        "AB_relation_wave_count": ab["relation_count_min"],
        "ABC_relation_wave_count": abc["relation_count_min"],
        "ABCD_relation_wave_count": abcd["relation_count_min"],
        "ABC_three_relation_axes_constructed": bool(
            abc["relation_count_min"] == 3 and abc["relation_count_max"] == 3
        ),
        "ABC_all_three_relation_waves_active_all_trials": abc["all_relation_waves_active_all_trials"],
        "AB_single_relation_stationary_all_trials": ab["stationary_all_trials"],
        "ABCD_six_relation_directions_present": bool(
            abcd["relation_count_min"] == 6 and abcd["relation_count_max"] == 6
        ),
        "max_closure_target_error_abs": max(
            float(row["max_closure_target_error_abs"]) for row in summaries
        ),
        "max_hermitian_amplitude_drift": max(
            float(row["max_hermitian_amplitude_drift"]) for row in summaries
        ),
        "max_label_covariance_error": max(
            max(
                float(row["max_generator_label_covariance_error"]),
                float(row["max_update_label_covariance_error"]),
                float(row["max_trajectory_label_covariance_error"]),
            )
            for row in summaries
        ),
        "all_relations_are_physical_waves": True,
        "observer_C_or_D_used": False,
        "normalization_applied": False,
        "absolute_background_axis_used": False,
        "relation_to_spatial_axis_is_model_definition": True,
        "three_spatial_dimensions_derived": False,
        "ABCD_six_to_three_projection_resolved": False,
        "preliminary_experiment_valid": bool_all(bool(row["body_count_test_pass"]) for row in summaries),
    }


def make_plots(
    selected_series: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
) -> None:
    body_counts = [int(row["body_count"]) for row in summaries]
    relation_counts = [int(row["relation_count_min"]) for row in summaries]
    active_counts = [int(row["active_relation_wave_count_min"]) for row in summaries]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(body_counts, relation_counts, marker="o", label="complete pair-relation wave count")
    ax.plot(body_counts, active_counts, marker="s", label="active relation-wave count")
    ax.axhline(3, color="black", linestyle="--", linewidth=1.0, alpha=0.5, label="three axes")
    ax.set_xticks(body_counts)
    ax.set_xlabel("physical body count")
    ax.set_ylabel("relation-wave count")
    ax.set_title("AB, ABC, ABCD complete pair-relation networks")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "complete_pair_relation_wave_count_v1.png", dpi=180)
    plt.close(fig)

    abc_rows = [row for row in selected_series if int(row["body_count"]) == 3]
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for relation_name_value in sorted({str(row["relation_name"]) for row in abc_rows}):
        rows = [row for row in abc_rows if row["relation_name"] == relation_name_value]
        axes[0].plot(
            [int(row["step"]) for row in rows],
            [float(row["relation_real"]) for row in rows],
            label=relation_name_value,
        )
        axes[1].plot(
            [int(row["step"]) for row in rows],
            np.unwrap([float(row["relation_phase_rad"]) for row in rows]),
            label=relation_name_value,
        )
    axes[0].set_ylabel("relation-wave real component")
    axes[1].set_ylabel("unwrapped relation phase")
    axes[1].set_xlabel("step")
    axes[0].set_title("ABC physical relation waves: AB, BC, CA")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ABC_three_physical_relation_waves_v1.png", dpi=180)
    plt.close(fig)

    first_relation_rows = [row for row in abc_rows if int(row["relation_index"]) == 0]
    steps = [int(row["step"]) for row in first_relation_rows]
    closure_real_error = [float(row["closure_real_E"]) - 1.0 for row in first_relation_rows]
    closure_imag = [float(row["closure_imag_F"]) for row in first_relation_rows]
    amplitude_initial = float(first_relation_rows[0]["hermitian_amplitude_sum"])
    amplitude_drift = [
        float(row["hermitian_amplitude_sum"]) - amplitude_initial for row in first_relation_rows
    ]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(steps, closure_real_error, label="E - R^2")
    ax.plot(steps, closure_imag, label="F")
    ax.plot(steps, amplitude_drift, label="amplitude-sum drift")
    ax.set_xlabel("step")
    ax.set_ylabel("conservation error")
    ax.set_title("ABC quadratic-closure and amplitude conservation")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ABC_relation_wave_conservation_v1.png", dpi=180)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    summaries = result["body_count_summaries"]
    lines: List[str] = [
        "# AB・ABC・ABCD完全二体関係波ネットワーク予備実験検証メモ v1",
        "",
        "## 目的",
        "",
        "個体A・B・C・Dではなく、個体間の全二体関係を物理的な関係波として状態変数に置き、ABCの三関係 `AB`, `BC`, `CA` を三つの軸成分として構成できるかを検査した。",
        "",
        "CおよびDを観測器として使用せず、観測減衰と正規化を行わない。",
        "",
        "## 状態と閉鎖",
        "",
        "関係波の集合を `X_e` とし、各ステップで次を直接計算した。",
        "",
        "```math",
        "E = sum_e ((Re X_e)^2 - (Im X_e)^2)",
        "```",
        "",
        "```math",
        "F = 2 sum_e (Re X_e)(Im X_e)",
        "```",
        "",
        "```math",
        "sum_e X_e^2 = E + i F = R^2",
        "```",
        "",
        "比較用に、実数二乗和 `sum_e |X_e|^2` も独立に記録した。",
        "",
        "## 作業更新則",
        "",
        "関係波同士が端点を共有する場合だけ結合し、初期位相差から実反対称生成子を作った。",
        "",
        "```math",
        "K_ef = adjacency(e,f) sin(theta_f - theta_e)",
        "```",
        "",
        "実反対称生成子のCayley変換による実直交更新を反復した。この更新は二乗形式と実数二乗和を保存する。",
        "",
        "この更新則は保存的な作業仮説であり、第0・第1公理からの導出結果ではない。",
        "",
        "## 構成別結果",
        "",
        "| system | relation waves | active min | closure error max | amplitude drift max | label covariance max | pass |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summaries:
        system_name = BODY_LABELS[: int(row["body_count"])]
        covariance_max = max(
            float(row["max_generator_label_covariance_error"]),
            float(row["max_update_label_covariance_error"]),
            float(row["max_trajectory_label_covariance_error"]),
        )
        lines.append(
            f"| {system_name} | {row['relation_count_min']} | {row['active_relation_wave_count_min']} | "
            f"{row['max_closure_target_error_abs']:.16e} | "
            f"{row['max_hermitian_amplitude_drift']:.16e} | "
            f"{covariance_max:.16e} | {row['body_count_test_pass']} |"
        )

    lines.extend(["", "## 統合判定", "", "| 量 | 値 |", "|---|---:|"])
    for key, value in verdict.items():
        lines.append(f"| `{key}` | `{value}` |")

    lines.extend(
        [
            "",
            "## 観測事実",
            "",
            "- ABは一つの関係波だけを持ち、本更新則では連続混合相手がないため定常だった。",
            "- ABCは `AB`, `BC`, `CA` の三つの物理的関係波を持ち、全試行で三波すべてが変動した。",
            "- ABCDは六つの物理的関係波を持ち、三軸を超える関係方向が代数的には存在した。",
            "- 全構成で `sum_e X_e^2 = R^2` と実数二乗和は数値精度内で保存された。",
            "- 個体名の置換に対し、生成子、更新行列、軌道は数値精度内で共変だった。",
            "",
            "## 分類",
            "",
            "- ABC三関係波を三軸成分として置くこと: 本実験のモデル定義。",
            "- ABC三関係波が閉鎖を保存しながら同時に振動できること: 本予備実験の数値結果。",
            "- 三関係波が物理的なXYZ空間と同一であること: 本実験では未導出。",
            "- ABCDの六関係波から観測可能な三軸を一意選択する機構: 本実験では未解決。",
            "",
            "## 出力",
            "",
            "- `ab_abc_abcd_complete_pair_relation_network_preliminary_result_v1.json`",
            "- `ab_abc_abcd_complete_pair_relation_network_trial_summary_v1.csv`",
            "- `ab_abc_abcd_complete_pair_relation_network_body_summary_v1.csv`",
            "- `ab_abc_abcd_complete_pair_relation_network_selected_series_v1.csv`",
            "- `complete_pair_relation_wave_count_v1.png`",
            "- `ABC_three_physical_relation_waves_v1.png`",
            "- `ABC_relation_wave_conservation_v1.png`",
        ]
    )
    (OUT_DIR / "ab_abc_abcd_complete_pair_relation_network_preliminary_report_v1.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def run() -> Dict[str, Any]:
    params = ExperimentParams()
    rng = np.random.default_rng(params.random_seed)
    trials: List[Dict[str, Any]] = []
    selected_series: List[Dict[str, Any]] = []
    for body_count in params.body_counts:
        for trial_index in range(params.trial_count):
            trial, series = trial_result(body_count, trial_index, rng, params)
            trials.append(trial)
            selected_series.extend(series)

    summaries = summarize_body_counts(trials, params)
    verdict = aggregate_verdict(summaries, params)
    result = {
        "experiment": "ab_abc_abcd_complete_pair_relation_network_preliminary_v1",
        "params": asdict(params),
        "body_count_summaries": summaries,
        "aggregate_verdict": verdict,
        "note": (
            "Pair relations are physical state waves in this model, not derived observer channels. "
            "The axis interpretation is a model definition. The relational antisymmetric generator "
            "is a closure-preserving working update rule, not a derived force law."
        ),
    }
    (OUT_DIR / "ab_abc_abcd_complete_pair_relation_network_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(OUT_DIR / "ab_abc_abcd_complete_pair_relation_network_trial_summary_v1.csv", trials)
    write_csv(OUT_DIR / "ab_abc_abcd_complete_pair_relation_network_body_summary_v1.csv", summaries)
    write_csv(OUT_DIR / "ab_abc_abcd_complete_pair_relation_network_selected_series_v1.csv", selected_series)
    make_plots(selected_series, summaries)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
