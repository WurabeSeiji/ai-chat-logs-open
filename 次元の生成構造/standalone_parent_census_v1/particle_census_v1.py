#!/usr/bin/env python3
"""Standalone structural census for ``closed_wave_trajectory_v1``.

This reader does not import or execute the generator.  It reads only the
versioned file contract, resolves temporal harmonics on the system-resolution
clock, and keeps algebraic/internal readouts separate from physical labels.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA_VERSION = "closed_wave_trajectory_v1"
CENSUS_VERSION = "closed_wave_particle_census_v1"
RANDOM_REFERENCE_ALPHA = 0.01  # reporting convention, not a model parameter


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read a closed-wave trajectory and make a standalone structural census."
    )
    parser.add_argument("--input", type=Path, required=True, help="generator run directory")
    parser.add_argument("--output", type=Path, required=True, help="new census directory")
    parser.add_argument(
        "--max-print", type=int, default=24, help="maximum high-power rows printed to stdout"
    )
    return parser.parse_args()


def load_contract(input_dir: Path) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray]:
    manifest_path = input_dir / "manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(f"missing manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA_VERSION:
        raise SystemExit(
            f"unsupported schema {manifest.get('schema')!r}; expected {SCHEMA_VERSION!r}"
        )

    edges = np.load(input_dir / manifest["arrays"]["edges"]["file"], mmap_mode="r")
    steps = np.load(input_dir / manifest["arrays"]["steps"]["file"], mmap_mode="r")
    trajectory = np.load(
        input_dir / manifest["arrays"]["trajectory"]["file"], mmap_mode="r"
    )
    expected_shape = tuple(manifest["arrays"]["trajectory"]["shape"])
    if trajectory.shape != expected_shape:
        raise SystemExit(
            f"trajectory shape mismatch: file={trajectory.shape}, manifest={expected_shape}"
        )
    if trajectory.ndim != 2 or edges.shape != (trajectory.shape[1], 2):
        raise SystemExit("invalid trajectory/edge topology in input contract")
    if steps.shape != (trajectory.shape[0],):
        raise SystemExit("steps array length does not match trajectory")
    return manifest, edges, steps, trajectory


def lag_fidelity(values: np.ndarray, lag: int) -> float | None:
    overlap = lag_overlap(values, lag)
    return float(abs(overlap)) if overlap is not None else None


def lag_overlap(values: np.ndarray, lag: int) -> complex | None:
    if lag <= 0 or len(values) <= lag:
        return None
    left = values[:-lag]
    right = values[lag:]
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom == 0.0:
        return None
    return complex(np.vdot(left, right) / denom)


def correlation_lifetime(values: np.ndarray, recurrence_order: int) -> tuple[int | None, str]:
    if recurrence_order < 1:
        return None, "unavailable"
    largest_tested = 0
    for lag in range(recurrence_order, len(values), recurrence_order):
        fidelity = lag_fidelity(values, lag)
        if fidelity is None:
            continue
        largest_tested = lag
        if fidelity < math.exp(-1.0):
            return lag, "measured_first_e_fold_crossing"
    if largest_tested:
        return largest_tested, "lower_bound_no_e_fold_crossing"
    return None, "unavailable"


def random_vector_coherence_reference(
    complex_dimension: int, comparisons: int, alpha: float
) -> float:
    """Bonferroni reference for |<u,v>| of independent complex unit vectors.

    This is only a transparent random-vector reference, not proof that
    consecutive nonlinear trajectory windows are statistically independent.
    """

    if complex_dimension <= 1:
        return 1.0
    tail_probability = alpha / max(comparisons, 1)
    squared = 1.0 - tail_probability ** (1.0 / (complex_dimension - 1))
    return float(math.sqrt(max(0.0, min(1.0, squared))))


def representative_mode(
    vectors: np.ndarray,
) -> tuple[np.ndarray, dict[str, float | None]]:
    """Phase-align one spatial vector per window and return their representative."""

    norms = np.linalg.norm(vectors, axis=1)
    reference_index = int(np.argmax(norms))
    reference = vectors[reference_index]
    aligned = np.zeros_like(vectors)
    for index, vector in enumerate(vectors):
        if norms[index] == 0.0 or norms[reference_index] == 0.0:
            continue
        overlap = np.vdot(reference, vector)
        phase = float(np.angle(overlap)) if abs(overlap) else 0.0
        aligned[index] = vector * np.exp(-1j * phase)
    representative = np.mean(aligned, axis=0)
    representative_norm = float(np.linalg.norm(representative))
    if representative_norm:
        representative /= representative_norm

    coherences: list[float] = []
    phase_steps: list[complex] = []
    for left, right, norm_left, norm_right in zip(
        vectors[:-1], vectors[1:], norms[:-1], norms[1:]
    ):
        denom = float(norm_left * norm_right)
        if denom == 0.0:
            continue
        overlap = np.vdot(left, right)
        coherences.append(float(abs(overlap) / denom))
        if abs(overlap):
            phase_steps.append(overlap / abs(overlap))

    gram = np.einsum("im,jm->ij", vectors, np.conj(vectors), optimize=False)
    eigenvalues = np.linalg.eigvalsh(gram).real
    eigenvalues[eigenvalues < 0.0] = 0.0
    gram_trace = float(np.sum(eigenvalues))
    rank1_fraction = float(eigenvalues[-1] / gram_trace) if gram_trace else 0.0
    phase_step = float(np.degrees(np.angle(np.sum(phase_steps)))) if phase_steps else None
    return representative, {
        "coherence_mean": float(np.mean(coherences)) if coherences else 0.0,
        "coherence_min": float(np.min(coherences)) if coherences else 0.0,
        "rank1_fraction": rank1_fraction,
        "phase_step_per_window_deg": phase_step,
    }


def numerical_rank_from_rows(rows: np.ndarray) -> tuple[int, list[float], float]:
    """Rank through the small row Gram matrix, with the standard dtype tolerance."""

    if rows.size == 0:
        return 0, [], 0.0
    gram = np.einsum("im,jm->ij", rows, np.conj(rows), optimize=False)
    eigenvalues = np.linalg.eigvalsh(gram).real
    eigenvalues[eigenvalues < 0.0] = 0.0
    singular_values = np.sqrt(eigenvalues[::-1])
    if not len(singular_values) or singular_values[0] == 0.0:
        return 0, singular_values.tolist(), 0.0
    real_dtype = np.empty((), dtype=rows.dtype).real.dtype
    # The Gram eigenvalues are squared singular values.  Roundoff is O(eps)
    # at the eigenvalue level; taking sqrt first would inflate it to O(sqrt(eps))
    # and can create impossible ranks larger than the column dimension.
    eigenvalue_tolerance = (
        max(rows.shape) * np.finfo(real_dtype).eps * float(eigenvalues[-1])
    )
    tolerance = math.sqrt(eigenvalue_tolerance)
    rank = int(np.sum(eigenvalues > eigenvalue_tolerance))
    rank = min(rank, rows.shape[0], rows.shape[1])
    return rank, singular_values.tolist(), tolerance


def project_trajectory(
    trajectory: np.ndarray,
    representative: np.ndarray,
    used_samples: int,
    chunk_rows: int = 4096,
) -> np.ndarray:
    """Compute <representative|Z_t> without relying on platform BLAS quirks."""

    coefficient = np.empty(used_samples, dtype=np.complex128)
    conjugate = np.conj(np.asarray(representative, dtype=np.complex128))
    for start in range(0, used_samples, chunk_rows):
        stop = min(used_samples, start + chunk_rows)
        block = np.asarray(trajectory[start:stop], dtype=np.complex128)
        coefficient[start:stop] = np.einsum(
            "tm,m->t", block, conjugate, optimize=False
        )
    return coefficient


def fmt_float(value: Any, digits: int = 5) -> str:
    if value is None:
        return "—"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return "—"
    return f"{number:.{digits}g}"


STATUS_JA = {
    "exact_finite_order_closed_mode": "厳密有限位数・単独零閉鎖",
    "persistent_closed_candidate": "持続する単独零閉鎖候補",
    "persistent_needs_bundle_closure": "持続するが倍音束での閉鎖が必要",
    "sea_or_unresolved": "海または未解決",
    "pending_recurrence": "回帰判定待ち",
}

SPIN_STRUCTURE_JA = {
    "half_integer_internal_layer": "半整数型の内部二重被覆",
    "integer_internal_layer": "整数型の内部層",
    "unresolved_non_null_single_mode": "未確定（単独では零閉鎖しない）",
}

PAIR_BRANCH_JA = {
    "+m/n": "正方向の共役枝",
    "-m/n": "負方向の共役枝",
    "neutral": "中立",
}

SOURCE_JA = {
    "white_null": "白色雑音から作った零閉鎖初期状態",
    "single_mode": "単一波の対照初期状態",
}


def human_readable_row(row: dict[str, Any]) -> dict[str, str]:
    if row["lifetime_steps"] is None:
        lifetime = "測定不能"
    elif row["lifetime_status"] == "lower_bound_no_e_fold_crossing":
        lifetime = f"{row['lifetime_steps']}以上"
    else:
        lifetime = str(row["lifetime_steps"])
    bf_display = {"B": "B型", "F": "F型"}.get(row["BF_type"], "未定義")
    return {
        "モード": f"k={int(row['signed_k']):+d}",
        "波数k": str(row["signed_k"]),
        "有限位数住所m/n": row["address"],
        "回帰位数n": str(row["state_order"]),
        "全体に占める強度比": fmt_float(row["power_fraction"], 6),
        "窓間コヒーレンス平均": fmt_float(row["coherence_mean"], 5),
        "同一空間方向への集中度": fmt_float(row["rank1_fraction"], 5),
        "単独零閉鎖誤差": fmt_float(row["single_mode_closure_residual"], 5),
        "B/F住所予測": bf_display,
        "粒子・反粒子の共役枝": PAIR_BRANCH_JA.get(row["pair_branch"], row["pair_branch"]),
        "住所から読むモデル電荷量": fmt_float(row["model_charge_magnitude"], 7),
        "正負共役対のグラム行列式detΓ": fmt_float(row["conjugate_pair_det_gamma"], 5),
        "状態と二次量の被覆比": f"{row['cover_ratio']}:1",
        "内部スピン構造": SPIN_STRUCTURE_JA.get(
            row["spin_structure"], row["spin_structure"]
        ),
        "相関寿命（ステップ）": lifetime,
        "総合判定": STATUS_JA.get(row["status"], row["status"]),
    }


def markdown_subtable(human_rows: list[dict[str, str]], headers: list[str]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for human_row in human_rows:
        lines.append("| " + " | ".join(human_row[header] for header in headers) + " |")
    return "\n".join(lines)


def markdown_table(rows: list[dict[str, Any]]) -> str:
    human_rows = [human_readable_row(row) for row in rows]
    wave_headers = [
        "モード",
        "有限位数住所m/n",
        "回帰位数n",
        "全体に占める強度比",
        "窓間コヒーレンス平均",
        "同一空間方向への集中度",
        "単独零閉鎖誤差",
        "相関寿命（ステップ）",
        "総合判定",
    ]
    property_headers = [
        "モード",
        "B/F住所予測",
        "粒子・反粒子の共役枝",
        "住所から読むモデル電荷量",
        "正負共役対のグラム行列式detΓ",
        "状態と二次量の被覆比",
        "内部スピン構造",
    ]
    return "\n".join(
        [
            "### A. 波・閉鎖・寿命",
            "",
            markdown_subtable(human_rows, wave_headers),
            "",
            "### B. 粒子属性の読出し候補",
            "",
            markdown_subtable(human_rows, property_headers),
        ]
    )


def main() -> None:
    args = parse_args()
    started = time.time()
    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    if output_dir.exists():
        raise SystemExit(f"output already exists; refusing to overwrite: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)

    manifest, edges, steps, trajectory = load_contract(input_dir)
    theory = manifest["theory_inputs"]
    resolution = int(theory["resolution_mathcal_N"])
    n_body = int(theory["n_body"])
    m_relations = int(theory["m_relations"])
    stride = int(manifest["numerical_observation"]["record_stride"])
    if stride != 1:
        raise SystemExit(
            "census v1 requires record_stride=1 to avoid harmonic-address aliasing; "
            f"input stride is {stride}"
        )
    if len(steps) < 2 or np.any(np.diff(steps) != 1):
        raise SystemExit("saved steps are not consecutive")

    n_windows = len(trajectory) // resolution
    if n_windows < 2:
        raise SystemExit(
            f"at least two complete resolution windows are required; "
            f"have {len(trajectory)} states for resolution {resolution}"
        )
    used_samples = n_windows * resolution
    total_power_by_bin = np.zeros(resolution, dtype=float)
    storage_real_dtype = np.empty((), dtype=trajectory.dtype).real.dtype
    machine_epsilon = float(np.finfo(storage_real_dtype).eps)

    audit = manifest.get("invariant_audit", {})
    numerical_error_envelope = float(
        math.sqrt(machine_epsilon)
        + abs(float(audit.get("max_norm_drift", 0.0)))
        + abs(float(audit.get("max_closure_abs", 0.0)))
    )
    comparison_count = max(1, (n_windows - 1) * resolution)
    random_coherence_reference = random_vector_coherence_reference(
        m_relations, comparison_count, RANDOM_REFERENCE_ALPHA
    )

    with tempfile.TemporaryDirectory(prefix="census_work_", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        spectra = np.lib.format.open_memmap(
            temp_dir / "spectra.npy",
            mode="w+",
            dtype=np.complex128,
            shape=(n_windows, resolution, m_relations),
        )
        for window in range(n_windows):
            start = window * resolution
            segment = np.asarray(trajectory[start : start + resolution], dtype=np.complex128)
            transformed = np.fft.fft(segment, axis=0) / math.sqrt(resolution)
            spectra[window] = transformed
            total_power_by_bin += np.sum(np.abs(transformed) ** 2, axis=1)
        spectra.flush()
        total_power_by_bin /= n_windows
        total_spectral_power = float(np.sum(total_power_by_bin))
        largest_bin_power = float(np.max(total_power_by_bin))
        power_floor = (
            max(resolution, m_relations)
            * machine_epsilon
            * max(largest_bin_power, np.finfo(float).tiny)
        )
        occupied_bins = np.where(total_power_by_bin > power_floor)[0]

        representatives = np.lib.format.open_memmap(
            temp_dir / "representatives.npy",
            mode="w+",
            dtype=np.complex128,
            shape=(resolution, m_relations),
        )
        representative_metrics: list[dict[str, float | None] | None] = [None] * resolution
        for bin_index in occupied_bins:
            vectors = np.asarray(spectra[:, int(bin_index), :])
            representative, metrics = representative_mode(vectors)
            # Canonical global phase: the largest component is real nonnegative.
            pivot = int(np.argmax(np.abs(representative)))
            if abs(representative[pivot]):
                representative *= np.exp(-1j * np.angle(representative[pivot]))
            representatives[int(bin_index)] = representative
            representative_metrics[int(bin_index)] = metrics
        representatives.flush()

        occupied_representatives = np.asarray(representatives[occupied_bins])
        wave_rank, singular_values, rank_tolerance = numerical_rank_from_rows(
            occupied_representatives
        )
        if wave_rank > m_relations:
            raise RuntimeError(
                f"internal rank error: W={wave_rank} exceeds relation dimension M={m_relations}"
            )

        preliminary_rows: list[dict[str, Any]] = []
        recurrence_bins: list[int] = []
        for bin_index_raw in occupied_bins:
            bin_index = int(bin_index_raw)
            signed_k = bin_index if bin_index <= resolution // 2 else bin_index - resolution
            if signed_k == 0:
                divisor = resolution
                address_m = 0
                state_order = 1
            else:
                divisor = math.gcd(abs(signed_k), resolution)
                address_m = signed_k // divisor
                state_order = resolution // divisor
            quadratic_order = state_order // math.gcd(state_order, 2)
            cover_ratio = state_order // quadratic_order
            representative = np.asarray(representatives[bin_index])
            norm2 = float(np.real(np.vdot(representative, representative)))
            closure_residual = (
                abs(complex(representative @ representative)) / norm2 if norm2 else math.inf
            )
            real_imag_rank = int(
                np.linalg.matrix_rank(
                    np.column_stack([representative.real, representative.imag])
                )
            )
            metrics = representative_metrics[bin_index]
            assert metrics is not None
            coherence_significant = metrics["coherence_min"] >= random_coherence_reference
            cartan_ready = (
                real_imag_rank == 2 and closure_residual <= numerical_error_envelope
            )
            if coherence_significant:
                recurrence_bins.append(bin_index)

            if signed_k > 0:
                pair_branch = "+m/n"
            elif signed_k < 0:
                pair_branch = "-m/n"
            else:
                pair_branch = "neutral"

            conjugate_bin = (-bin_index) % resolution
            pair_det_gamma: float | None = None
            if conjugate_bin in occupied_bins:
                a = representative
                b = np.asarray(representatives[conjugate_bin])
                norm_a = float(np.real(np.vdot(a, a)))
                norm_b = float(np.real(np.vdot(b, b)))
                overlap = np.vdot(a, b)
                pair_det_gamma = max(0.0, norm_a * norm_b - float(abs(overlap) ** 2))

            if cartan_ready and cover_ratio == 2:
                spin_structure = "half_integer_internal_layer"
            elif cartan_ready:
                spin_structure = "integer_internal_layer"
            else:
                spin_structure = "unresolved_non_null_single_mode"

            bf_prediction = (
                ("F" if abs(signed_k) % 2 else "B") if resolution % 2 == 0 else None
            )
            preliminary_rows.append(
                {
                    "mode_id": f"k{signed_k:+d}",
                    "bin_index": bin_index,
                    "signed_k": signed_k,
                    "address_m": address_m,
                    "address_n": state_order,
                    "address": f"{address_m}/{state_order}",
                    "state_order": state_order,
                    "quadratic_order": quadratic_order,
                    "cover_ratio": cover_ratio,
                    "power": float(total_power_by_bin[bin_index]),
                    "power_fraction": float(
                        total_power_by_bin[bin_index] / total_spectral_power
                    ),
                    **metrics,
                    "single_mode_closure_residual": closure_residual,
                    "real_imag_rank": real_imag_rank,
                    "cartan_ready": bool(cartan_ready),
                    "BF_address_prediction": bf_prediction,
                    "BF_measured": None,
                    "BF_type": bf_prediction if bf_prediction is not None else "—",
                    "BF_status": (
                        "register_address_prediction_only"
                        if bf_prediction is not None
                        else "undefined_without_half_register_shift"
                    ),
                    "half_cycle_overlap_real": None,
                    "half_cycle_overlap_imag": None,
                    "half_cycle_fidelity": None,
                    "BF_behavior_residual_to_address_prediction": None,
                    "pair_branch": pair_branch,
                    "conjugate_address": f"{-address_m}/{state_order}",
                    "model_charge_magnitude": float(
                        math.sin(math.pi * address_m / state_order) ** 2
                    ),
                    "conjugate_pair_det_gamma": pair_det_gamma,
                    "coherence_random_reference": random_coherence_reference,
                    "coherence_significant_vs_random_reference": bool(coherence_significant),
                    "state_recurrence_fidelity": None,
                    "quadratic_recurrence_fidelity": None,
                    "spin_structure": spin_structure,
                    "physical_spin": "unresolved_parent_cartan_map",
                    "lifetime_steps": None,
                    "lifetime_status": "unavailable",
                    "status": "pending_recurrence" if coherence_significant else "sea_or_unresolved",
                }
            )

        rows_by_bin = {int(row["bin_index"]): row for row in preliminary_rows}
        for bin_index in recurrence_bins:
            row = rows_by_bin[bin_index]
            representative = np.asarray(representatives[bin_index])
            coefficient = project_trajectory(trajectory, representative, used_samples)
            state_fidelity = lag_fidelity(coefficient, int(row["state_order"]))
            squared_fidelity = lag_fidelity(
                coefficient * coefficient, int(row["quadratic_order"])
            )
            half_cycle_overlap = (
                lag_overlap(coefficient, resolution // 2) if resolution % 2 == 0 else None
            )
            lifetime_steps, lifetime_status = correlation_lifetime(
                coefficient, int(row["state_order"])
            )
            row["state_recurrence_fidelity"] = state_fidelity
            row["quadratic_recurrence_fidelity"] = squared_fidelity
            if half_cycle_overlap is not None:
                row["half_cycle_overlap_real"] = float(half_cycle_overlap.real)
                row["half_cycle_overlap_imag"] = float(half_cycle_overlap.imag)
                row["half_cycle_fidelity"] = float(abs(half_cycle_overlap))
                row["BF_measured"] = "F" if half_cycle_overlap.real < 0.0 else "B"
                prediction = row["BF_address_prediction"]
                target = -1.0 if prediction == "F" else 1.0
                behavior_residual = float(abs(half_cycle_overlap - target))
                row["BF_behavior_residual_to_address_prediction"] = behavior_residual
                row["BF_status"] = (
                    "confirmed_projected_half_cycle"
                    if behavior_residual <= numerical_error_envelope
                    else "address_prediction_not_confirmed_by_half_cycle"
                )
            row["lifetime_steps"] = lifetime_steps
            row["lifetime_status"] = lifetime_status

            recurrence_exact = (
                state_fidelity is not None
                and squared_fidelity is not None
                and 1.0 - state_fidelity <= numerical_error_envelope
                and 1.0 - squared_fidelity <= numerical_error_envelope
            )
            if recurrence_exact and row["cartan_ready"]:
                row["status"] = "exact_finite_order_closed_mode"
            elif row["cartan_ready"]:
                row["status"] = "persistent_closed_candidate"
            else:
                row["status"] = "persistent_needs_bundle_closure"

        del representatives
        del spectra

    rows = sorted(preliminary_rows, key=lambda row: row["power"], reverse=True)
    bound_statuses = {"exact_finite_order_closed_mode", "persistent_closed_candidate"}
    bound_power = float(sum(row["power"] for row in rows if row["status"] in bound_statuses))
    exact_count = sum(row["status"] == "exact_finite_order_closed_mode" for row in rows)
    persistent_closed_count = sum(row["status"] in bound_statuses for row in rows)

    summary = {
        "census_version": CENSUS_VERSION,
        "input_schema": SCHEMA_VERSION,
        "input_directory": str(input_dir),
        "source_condition": manifest["source"],
        "generator_sha256": manifest.get("generator_sha256"),
        "generator_seed": manifest.get("reproducibility", {}).get("seed"),
        "generator_steps": manifest.get("numerical_observation", {}).get("steps"),
        "input_invariant_audit": manifest.get("invariant_audit", {}),
        "resolution_mathcal_N": resolution,
        "n_body": n_body,
        "m_relations": m_relations,
        "complete_windows": n_windows,
        "used_samples": used_samples,
        "occupied_harmonic_bins": int(len(occupied_bins)),
        "independent_wave_directions_W": wave_rank,
        "W_le_M": bool(wave_rank <= m_relations),
        "mode_singular_values": singular_values,
        "numerical_rank_tolerance": rank_tolerance,
        "power_floor": power_floor,
        "random_vector_reference_alpha": RANDOM_REFERENCE_ALPHA,
        "random_vector_coherence_reference": random_coherence_reference,
        "numerical_error_envelope": numerical_error_envelope,
        "exact_finite_order_closed_modes": exact_count,
        "persistent_closed_candidates": persistent_closed_count,
        "bound_candidate_power_fraction": (
            bound_power / total_spectral_power if total_spectral_power else 0.0
        ),
        "sea_or_unresolved_power_fraction": (
            1.0 - bound_power / total_spectral_power if total_spectral_power else 1.0
        ),
        "spin_claim_boundary": (
            "internal 2:1 cover and Cartan-ready null plane are measured; "
            "physical spin remains unresolved until the parent-direction/Gram-Cartan map is verified"
        ),
        "runtime_seconds": time.time() - started,
    }

    script_path = Path(__file__).resolve()
    payload = {
        "summary": summary,
        "reader": script_path.name,
        "reader_sha256": sha256_file(script_path),
        "columns": {
            "address": "DFT register address only; persistence is a separate column",
            "model_charge_magnitude": "sin^2(pi*m/n), model-internal address readout",
            "conjugate_pair_det_gamma": "Gram determinant of +/-k representatives, not physical mass",
            "physical_spin": "not assigned in v1",
        },
        "rows": rows,
    }
    (output_dir / "census.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    fieldnames = list(rows[0].keys()) if rows else []
    with (output_dir / "census.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)

    human_rows = [human_readable_row(row) for row in rows]
    human_fieldnames = list(human_rows[0].keys()) if human_rows else []
    with (output_dir / "census_ja.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=human_fieldnames)
        if human_fieldnames:
            writer.writeheader()
            writer.writerows(human_rows)

    markdown = [
        "# 独立粒子構造一覧表 v1",
        "",
        f"- 入力条件: {SOURCE_JA.get(manifest['source'], manifest['source'])}",
        f"- $N={n_body}$, $\\mathcal N={resolution}$, $M={m_relations}$",
        f"- 独立波方向数: $W={wave_rank}\\le M$",
        f"- 完全分解能窓: {n_windows}",
        f"- 厳密有限位数・単独零閉鎖モード: {exact_count}",
        f"- 持続する単独零閉鎖候補: {persistent_closed_count}",
        f"- 候補割当強度: {summary['bound_candidate_power_fraction']:.8f}",
        f"- 海または未解決強度: {summary['sea_or_unresolved_power_fraction']:.8f}",
        "",
        "離散フーリエ変換（DFT）で読んだ住所は有限窓の数え上げである。粒子候補かどうかは、窓間コヒーレンス、",
        "同一方向への集中度、実回帰、単独零閉鎖誤差を別々に見て判定する。",
        "",
        "- **強度比**: 全時間周波数成分のうち、そのモードが占める割合",
        "- **窓間コヒーレンス**: 隣接する分解能窓で空間パターンが似ている度合い",
        f"- **同一方向への集中度**: 全{n_windows}窓を通じて一つの空間方向を保つ度合い",
        "- **単独零閉鎖誤差**: 0に近いほど、そのモード単独で $\\sum w_n^2=0$",
        "- **束閉鎖が必要**: 単独では閉じず、複数倍音を束ねる必要がある状態",
        "- **B/F住所予測**: ボース型（B型）／フェルミ型（F型）を住所の偶奇から読む予測。動力学的確定とは別",
        "- **正負共役対のグラム量**: 正負の波数を対にした内部量。物理的質量とはまだ同定しない",
        "- **内部スピン構造**: 内部二重被覆まで。物理スピン値はまだ未同定",
        "",
        markdown_table(rows),
        "",
    ]
    (output_dir / "census.md").write_text("\n".join(markdown), encoding="utf-8")

    print(
        f"N={n_body} resolution={resolution} M={m_relations} "
        f"W={wave_rank} exact={exact_count} persistent_closed={persistent_closed_count}"
    )
    print(markdown_table(rows[: max(0, args.max_print)]))
    print(f"saved: {output_dir}")


if __name__ == "__main__":
    main()
