"""N体完全二体関係波・生成子ランク線形上界と三方向飽和の予備実験 v1

固定生成子本線（初期位相から一回構成、以後固定）。
実験A: N=3..9 保存量・共変性・ランク・零空間射影のN拡張検査
実験B: 一般位置ランク則 rank K = 2*min(N, floor(M/2)) の検証と縮退探索
実験C: 二関係面の法線非一意性（d=3のみ一意、d>=4でO(d-2)ゲージ族）
実験D: N=5のABCD部分系：局所生成子と主生成子部分行列の一致検査
実験E: N>=6の多次元核：射影子一意・内部基底非一意の検査
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, asdict
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

TAU = 2.0 * math.pi
BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "nbody_rank_saturation_preliminary_result_v1"


@dataclass
class ExperimentParams:
    body_counts_a: Tuple[int, ...] = (3, 4, 5, 6, 7, 8, 9)
    trial_count_a: int = 32
    step_count: int = 720
    target_period_steps: int = 144
    radius_squared: float = 1.0
    imaginary_seed_amplitude: float = 0.35
    random_seed: int = 20260721
    invariant_tol: float = 1.0e-10
    covariance_tol: float = 1.0e-10
    rank_tol: float = 1.0e-10
    body_counts_b: Tuple[int, ...] = (3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    trial_count_b: int = 256
    normal_dims_c: Tuple[int, ...] = (3, 4, 5)
    gauge_sample_count_c: int = 64
    trial_count_c: int = 32
    subsystem_trial_count_d: int = 8
    kernel_trial_count_e: int = 8


def relation_pairs(body_count: int) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(body_count) for j in range(i + 1, body_count)]


def incidence_matrix(body_count: int, pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
    matrix = np.zeros((body_count, len(pairs)), dtype=float)
    for column, (i, j) in enumerate(pairs):
        matrix[i, column] = 1.0
        matrix[j, column] = 1.0
    return matrix


def relation_adjacency(pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
    count = len(pairs)
    matrix = np.zeros((count, count), dtype=float)
    for a in range(count):
        for b in range(count):
            if a != b and set(pairs[a]) & set(pairs[b]):
                matrix[a, b] = 1.0
    return matrix


def initial_closed_relation_state(
    relation_count: int, rng: np.random.Generator, params: ExperimentParams
) -> np.ndarray:
    real_direction = rng.normal(size=relation_count)
    real_direction /= np.linalg.norm(real_direction)
    imaginary_direction = rng.normal(size=relation_count)
    imaginary_direction -= float(np.dot(imaginary_direction, real_direction)) * real_direction
    imaginary_direction /= np.linalg.norm(imaginary_direction)
    s = params.imaginary_seed_amplitude
    real_norm = math.sqrt(params.radius_squared + s * s)
    return np.asarray(real_norm * real_direction + 1j * s * imaginary_direction, dtype=complex)


def raw_generator_from_phases(phases: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    phase_difference = phases[np.newaxis, :] - phases[:, np.newaxis]
    raw = adjacency * np.sin(phase_difference)
    return 0.5 * (raw - raw.T)


def normalized_generator(raw: np.ndarray) -> Tuple[np.ndarray, float]:
    spectral_norm = float(np.linalg.norm(raw, ord=2))
    if spectral_norm <= 1.0e-14:
        return np.zeros_like(raw), spectral_norm
    return raw / spectral_norm, spectral_norm


def cayley_orthogonal_update(generator: np.ndarray, target_period_steps: int) -> np.ndarray:
    dimension = generator.shape[0]
    if dimension == 1 or float(np.linalg.norm(generator)) <= 1.0e-14:
        return np.eye(dimension, dtype=float)
    cayley_scale = math.tan((TAU / target_period_steps) / 2.0)
    identity = np.eye(dimension, dtype=float)
    return np.linalg.solve(identity - cayley_scale * generator, identity + cayley_scale * generator)


def evolve(initial_state: np.ndarray, update: np.ndarray, step_count: int) -> np.ndarray:
    states = np.empty((step_count + 1, initial_state.size), dtype=complex)
    states[0] = initial_state
    for step in range(step_count):
        states[step + 1] = update @ states[step]
    return states


def closure_error_series(states: np.ndarray, radius_squared: float) -> Tuple[float, float]:
    real = states.real
    imag = states.imag
    closure_real = np.sum(real * real - imag * imag, axis=1)
    closure_imag = 2.0 * np.sum(real * imag, axis=1)
    closure_err = np.hypot(closure_real - radius_squared, closure_imag)
    amplitude = np.sum(real * real + imag * imag, axis=1)
    return float(np.max(closure_err)), float(np.max(np.abs(amplitude - amplitude[0])))


def spectral_structure(generator: np.ndarray, rank_tol: float) -> Tuple[int, int, np.ndarray]:
    _, singular_values, right_vectors_h = np.linalg.svd(generator, full_matrices=True)
    scale = max(1.0, float(singular_values[0]) if singular_values.size else 0.0)
    rank = int(np.sum(singular_values > rank_tol * scale))
    null_basis = right_vectors_h[rank:].T
    return rank, generator.shape[0] - rank, null_basis


def edge_permutation_matrix(
    pairs: Sequence[Tuple[int, int]], body_permutation: np.ndarray
) -> np.ndarray:
    index_by_pair = {tuple(sorted(pair)): index for index, pair in enumerate(pairs)}
    matrix = np.zeros((len(pairs), len(pairs)), dtype=float)
    for old_index, pair in enumerate(pairs):
        new_pair = tuple(sorted((int(body_permutation[pair[0]]), int(body_permutation[pair[1]]))))
        matrix[index_by_pair[new_pair], old_index] = 1.0
    return matrix


def expected_rank(body_count: int, relation_count: int) -> int:
    return 2 * min(body_count, relation_count // 2)


def vertex_decomposition_residual(
    raw_generator: np.ndarray, incidence: np.ndarray, phases: np.ndarray
) -> float:
    cosines = np.cos(phases)
    sines = np.sin(phases)
    total = np.zeros_like(raw_generator)
    for k in range(incidence.shape[0]):
        c_k = incidence[k, :] * cosines
        s_k = incidence[k, :] * sines
        total += np.outer(c_k, s_k) - np.outer(s_k, c_k)
    return float(np.max(np.abs(total - raw_generator)))


# ---------------------------------------------------------------- 実験A
def experiment_a(params: ExperimentParams) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    rng = np.random.default_rng(params.random_seed)
    trial_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for body_count in params.body_counts_a:
        pairs = relation_pairs(body_count)
        relation_count = len(pairs)
        adjacency = relation_adjacency(pairs)
        incidence = incidence_matrix(body_count, pairs)
        per_trial: List[Dict[str, Any]] = []
        for trial_index in range(params.trial_count_a):
            state0 = initial_closed_relation_state(relation_count, rng, params)
            phases = np.angle(state0)
            raw = raw_generator_from_phases(phases, adjacency)
            generator, _ = normalized_generator(raw)
            update = cayley_orthogonal_update(generator, params.target_period_steps)
            orthogonality_error = float(
                np.max(np.abs(update.T @ update - np.eye(relation_count)))
            )
            states = evolve(state0, update, params.step_count)
            closure_err, amplitude_drift = closure_error_series(states, params.radius_squared)

            rank, nullity, null_basis = spectral_structure(generator, params.rank_tol)
            decomposition_residual = vertex_decomposition_residual(raw, incidence, phases)

            kernel_projection_drift = 0.0
            if nullity > 0:
                projections = states @ null_basis
                kernel_projection_drift = float(np.max(np.abs(projections - projections[0])))

            body_permutation = rng.permutation(body_count)
            edge_perm = edge_permutation_matrix(pairs, body_permutation)
            permuted_state0 = edge_perm @ state0
            raw_p = raw_generator_from_phases(np.angle(permuted_state0), adjacency)
            generator_p, _ = normalized_generator(raw_p)
            generator_covariance_error = float(
                np.max(np.abs(generator_p - edge_perm @ generator @ edge_perm.T))
            )
            update_p = cayley_orthogonal_update(generator_p, params.target_period_steps)
            states_p = evolve(permuted_state0, update_p, params.step_count)
            trajectory_covariance_error = float(
                np.max(np.abs(states_p - states @ edge_perm.T))
            )

            per_trial.append(
                {
                    "body_count": body_count,
                    "trial_index": trial_index,
                    "relation_count": relation_count,
                    "generator_rank": rank,
                    "generator_nullity": nullity,
                    "rotation_plane_count": rank // 2,
                    "expected_rank": expected_rank(body_count, relation_count),
                    "rank_law_holds": bool(rank == expected_rank(body_count, relation_count)),
                    "max_closure_error": closure_err,
                    "max_amplitude_drift": amplitude_drift,
                    "max_orthogonality_error": orthogonality_error,
                    "vertex_decomposition_residual": decomposition_residual,
                    "max_kernel_projection_drift": kernel_projection_drift,
                    "generator_covariance_error": generator_covariance_error,
                    "trajectory_covariance_error": trajectory_covariance_error,
                }
            )
        trial_rows.extend(per_trial)
        summary = {
            "body_count": body_count,
            "relation_count": relation_count,
            "trial_count": params.trial_count_a,
            "rank_values": sorted({row["generator_rank"] for row in per_trial}),
            "nullity_values": sorted({row["generator_nullity"] for row in per_trial}),
            "rank_law_all_trials": all(row["rank_law_holds"] for row in per_trial),
            "max_closure_error": max(row["max_closure_error"] for row in per_trial),
            "max_amplitude_drift": max(row["max_amplitude_drift"] for row in per_trial),
            "max_orthogonality_error": max(row["max_orthogonality_error"] for row in per_trial),
            "max_vertex_decomposition_residual": max(
                row["vertex_decomposition_residual"] for row in per_trial
            ),
            "max_kernel_projection_drift": max(
                row["max_kernel_projection_drift"] for row in per_trial
            ),
            "max_generator_covariance_error": max(
                row["generator_covariance_error"] for row in per_trial
            ),
            "max_trajectory_covariance_error": max(
                row["trajectory_covariance_error"] for row in per_trial
            ),
        }
        summary["passed"] = bool(
            summary["rank_law_all_trials"]
            and summary["max_closure_error"] <= params.invariant_tol
            and summary["max_amplitude_drift"] <= params.invariant_tol
            and summary["max_orthogonality_error"] <= params.invariant_tol
            and summary["max_vertex_decomposition_residual"] <= params.invariant_tol
            and summary["max_kernel_projection_drift"] <= params.invariant_tol
            and summary["max_generator_covariance_error"] <= params.covariance_tol
            and summary["max_trajectory_covariance_error"] <= params.covariance_tol
        )
        summaries.append(summary)
    return trial_rows, summaries


# ---------------------------------------------------------------- 実験B
def experiment_b(params: ExperimentParams) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(params.random_seed + 1)
    rows: List[Dict[str, Any]] = []
    for body_count in params.body_counts_b:
        pairs = relation_pairs(body_count)
        relation_count = len(pairs)
        adjacency = relation_adjacency(pairs)
        generic_ranks: List[int] = []
        for _ in range(params.trial_count_b):
            phases = rng.uniform(0.0, TAU, size=relation_count)
            raw = raw_generator_from_phases(phases, adjacency)
            generator, _ = normalized_generator(raw)
            rank, _, _ = spectral_structure(generator, params.rank_tol)
            generic_ranks.append(rank)

        equal_phases = np.zeros(relation_count)
        raw_equal = raw_generator_from_phases(equal_phases, adjacency)
        rank_equal = int(np.linalg.matrix_rank(raw_equal, tol=1.0e-12))

        binary_phases = np.where(np.arange(relation_count) % 2 == 0, 0.0, math.pi)
        raw_binary = raw_generator_from_phases(binary_phases, adjacency)
        rank_binary = int(np.linalg.matrix_rank(raw_binary, tol=1.0e-12))

        two_value_phases = np.where(np.arange(relation_count) % 2 == 0, 0.0, 0.5)
        raw_two = raw_generator_from_phases(two_value_phases, adjacency)
        generator_two, _ = normalized_generator(raw_two)
        rank_two, _, _ = spectral_structure(generator_two, params.rank_tol)

        rows.append(
            {
                "body_count": body_count,
                "relation_count": relation_count,
                "expected_rank": expected_rank(body_count, relation_count),
                "generic_rank_values": sorted(set(generic_ranks)),
                "generic_rank_law_all": all(
                    rank == expected_rank(body_count, relation_count) for rank in generic_ranks
                ),
                "degenerate_equal_phase_rank": rank_equal,
                "degenerate_zero_pi_phase_rank": rank_binary,
                "degenerate_two_value_rank": rank_two,
            }
        )
    return rows


# ---------------------------------------------------------------- 実験C
def experiment_c(params: ExperimentParams) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(params.random_seed + 2)
    rows: List[Dict[str, Any]] = []
    for dimension in params.normal_dims_c:
        max_projector_drift = 0.0
        max_candidate_line_angle = 0.0
        for _ in range(params.trial_count_c):
            u = rng.normal(size=dimension)
            v = rng.normal(size=dimension)
            u /= np.linalg.norm(u)
            v -= float(np.dot(v, u)) * u
            v /= np.linalg.norm(v)
            plane_matrix = np.vstack([u, v])
            projector = np.eye(dimension) - plane_matrix.T @ np.linalg.solve(
                plane_matrix @ plane_matrix.T, plane_matrix
            )
            _, _, vh = np.linalg.svd(plane_matrix, full_matrices=True)
            null_basis = vh[2:].T
            candidates: List[np.ndarray] = []
            for _ in range(params.gauge_sample_count_c):
                gauge = np.linalg.qr(rng.normal(size=(dimension - 2, dimension - 2)))[0]
                rotated = null_basis @ gauge
                candidates.append(rotated[:, 0])
                drift = float(np.max(np.abs(rotated @ rotated.T - projector)))
                max_projector_drift = max(max_projector_drift, drift)
            for a in range(len(candidates)):
                for b in range(a + 1, len(candidates)):
                    cosine = min(1.0, abs(float(np.dot(candidates[a], candidates[b]))))
                    angle = math.degrees(math.acos(cosine))
                    max_candidate_line_angle = max(max_candidate_line_angle, angle)
        rows.append(
            {
                "display_dimension": dimension,
                "normal_candidate_dimension": dimension - 2,
                "max_projector_gauge_drift": max_projector_drift,
                "max_candidate_line_angle_deg": max_candidate_line_angle,
                "normal_unique_up_to_sign": bool(max_candidate_line_angle <= 1.0e-4),
            }
        )
    return rows


# ---------------------------------------------------------------- 実験D
def experiment_d(params: ExperimentParams) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(params.random_seed + 3)
    body_count = 5
    pairs = relation_pairs(body_count)
    adjacency = relation_adjacency(pairs)
    rows: List[Dict[str, Any]] = []
    for trial_index in range(params.subsystem_trial_count_d):
        state0 = initial_closed_relation_state(len(pairs), rng, params)
        phases = np.angle(state0)
        raw_global = raw_generator_from_phases(phases, adjacency)
        for subset in combinations(range(body_count), 4):
            edge_indices = [
                index for index, pair in enumerate(pairs) if set(pair) <= set(subset)
            ]
            local_pairs = [pairs[index] for index in edge_indices]
            local_adjacency = relation_adjacency(local_pairs)
            local_raw = raw_generator_from_phases(phases[edge_indices], local_adjacency)
            submatrix = raw_global[np.ix_(edge_indices, edge_indices)]
            submatrix_match_error = float(np.max(np.abs(local_raw - submatrix)))
            local_generator, _ = normalized_generator(local_raw)
            local_rank, local_nullity, _ = spectral_structure(local_generator, params.rank_tol)
            rows.append(
                {
                    "trial_index": trial_index,
                    "subsystem_bodies": "".join(chr(ord("A") + body) for body in subset),
                    "local_generator_rank": local_rank,
                    "local_rotation_plane_count": local_rank // 2,
                    "local_nullity": local_nullity,
                    "submatrix_match_error": submatrix_match_error,
                }
            )
    return rows


# ---------------------------------------------------------------- 実験E
def experiment_e(params: ExperimentParams) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(params.random_seed + 4)
    rows: List[Dict[str, Any]] = []
    for body_count in (6, 7, 8, 9):
        pairs = relation_pairs(body_count)
        relation_count = len(pairs)
        adjacency = relation_adjacency(pairs)
        for trial_index in range(params.kernel_trial_count_e):
            state0 = initial_closed_relation_state(relation_count, rng, params)
            raw = raw_generator_from_phases(np.angle(state0), adjacency)
            generator, _ = normalized_generator(raw)
            rank, nullity, null_basis = spectral_structure(generator, params.rank_tol)
            projector = null_basis @ null_basis.T
            idempotence_error = float(np.max(np.abs(projector @ projector - projector)))
            annihilation_error = float(np.max(np.abs(generator @ projector)))
            gauge_one = np.linalg.qr(rng.normal(size=(nullity, nullity)))[0]
            gauge_two = np.linalg.qr(rng.normal(size=(nullity, nullity)))[0]
            basis_one = null_basis @ gauge_one
            basis_two = null_basis @ gauge_two
            projector_gauge_error = float(
                np.max(np.abs(basis_one @ basis_one.T - basis_two @ basis_two.T))
            )
            first_direction_angle = math.degrees(
                math.acos(
                    min(1.0, abs(float(np.dot(basis_one[:, 0], basis_two[:, 0]))))
                )
            )
            rows.append(
                {
                    "body_count": body_count,
                    "trial_index": trial_index,
                    "relation_count": relation_count,
                    "generator_rank": rank,
                    "kernel_dimension": nullity,
                    "kernel_dimension_lower_bound": relation_count - 2 * body_count,
                    "kernel_projector_idempotence_error": idempotence_error,
                    "kernel_generator_annihilation_error": annihilation_error,
                    "kernel_projector_gauge_error": projector_gauge_error,
                    "kernel_basis_first_direction_angle_deg": first_direction_angle,
                }
            )
    return rows


# ---------------------------------------------------------------- 出力
def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def assert_all_finite(rows: List[Dict[str, Any]]) -> None:
    for row in rows:
        for key, value in row.items():
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f"non-finite value: {key}={value} in {row}")


def main() -> None:
    # macOS Accelerate BLAS は matmul で疑似的な浮動小数警告を出すため、
    # 警告抑制の代わりに全記録値の有限性を assert_all_finite で明示検査する。
    np.seterr(all="ignore")
    params = ExperimentParams()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    trial_rows_a, summaries_a = experiment_a(params)
    rows_b = experiment_b(params)
    rows_c = experiment_c(params)
    rows_d = experiment_d(params)
    rows_e = experiment_e(params)

    assert_all_finite(trial_rows_a)
    assert_all_finite(summaries_a)
    assert_all_finite(rows_c)
    assert_all_finite(rows_d)
    assert_all_finite(rows_e)

    write_csv(RESULT_DIR / "experiment_a_trials_v1.csv", trial_rows_a)
    write_csv(RESULT_DIR / "experiment_a_summary_v1.csv",
              [{**s, "rank_values": ";".join(map(str, s["rank_values"])),
                "nullity_values": ";".join(map(str, s["nullity_values"]))} for s in summaries_a])
    write_csv(RESULT_DIR / "experiment_b_rank_law_v1.csv",
              [{**r, "generic_rank_values": ";".join(map(str, r["generic_rank_values"]))} for r in rows_b])
    write_csv(RESULT_DIR / "experiment_c_normal_uniqueness_v1.csv", rows_c)
    write_csv(RESULT_DIR / "experiment_d_subsystems_v1.csv", rows_d)
    write_csv(RESULT_DIR / "experiment_e_kernel_v1.csv", rows_e)

    payload = {
        "parameters": asdict(params),
        "experiment_a_summaries": summaries_a,
        "experiment_b_rank_law": rows_b,
        "experiment_c_normal_uniqueness": rows_c,
        "experiment_d_subsystem_row_count": len(rows_d),
        "experiment_d_max_submatrix_match_error": max(
            row["submatrix_match_error"] for row in rows_d
        ),
        "experiment_d_local_rank_values": sorted(
            {row["local_generator_rank"] for row in rows_d}
        ),
        "experiment_e_rows": rows_e,
        "all_experiment_a_passed": all(summary["passed"] for summary in summaries_a),
        "all_rank_law_generic": all(row["generic_rank_law_all"] for row in rows_b),
    }
    with (RESULT_DIR / "nbody_rank_saturation_preliminary_result_v1.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    print(json.dumps(
        {
            "all_experiment_a_passed": payload["all_experiment_a_passed"],
            "all_rank_law_generic": payload["all_rank_law_generic"],
            "experiment_a": [
                {
                    "N": s["body_count"], "M": s["relation_count"],
                    "rank": s["rank_values"], "nullity": s["nullity_values"],
                    "passed": s["passed"],
                    "max_closure": s["max_closure_error"],
                    "max_H_drift": s["max_amplitude_drift"],
                    "max_traj_cov": s["max_trajectory_covariance_error"],
                    "max_vtx_decomp": s["max_vertex_decomposition_residual"],
                    "max_ker_drift": s["max_kernel_projection_drift"],
                }
                for s in summaries_a
            ],
            "experiment_b": [
                {"N": r["body_count"], "M": r["relation_count"],
                 "expected": r["expected_rank"], "generic": r["generic_rank_values"],
                 "law": r["generic_rank_law_all"],
                 "deg_equal": r["degenerate_equal_phase_rank"],
                 "deg_0pi": r["degenerate_zero_pi_phase_rank"],
                 "deg_2val": r["degenerate_two_value_rank"]}
                for r in rows_b
            ],
            "experiment_c": rows_c,
            "experiment_d": {
                "max_submatrix_match_error": payload["experiment_d_max_submatrix_match_error"],
                "local_rank_values": payload["experiment_d_local_rank_values"],
            },
            "experiment_e_max_gauge_error": max(
                r["kernel_projector_gauge_error"] for r in rows_e
            ),
            "experiment_e_min_basis_angle": min(
                r["kernel_basis_first_direction_angle_deg"] for r in rows_e
            ),
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
