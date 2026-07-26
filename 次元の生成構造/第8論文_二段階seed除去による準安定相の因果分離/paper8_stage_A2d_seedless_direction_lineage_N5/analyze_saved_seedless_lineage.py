#!/usr/bin/env python3
"""保存済み完全無seed Bdom標本からD34/P34方向系譜を再構成する。"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/paper8_a2d_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
A2A = PAPER8 / "paper8_stage_A2a_seedless_N5"
A2C = PAPER8 / "paper8_stage_A2c_direction_lineage_N5"
CONFIG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
PROCESSED = HERE / "processed"
FIGURES = HERE / "figures"
REPORTS = HERE / "reports"
LOGS = HERE / "logs"
FMT = ".17e"

INPUT_FILES = (
    "dominant_plane_steps.npy",
    "dominant_plane_values.npy",
    "f_timeseries.csv",
    "q_timeseries.csv",
    "occupation_timeseries.csv",
    "first_passage_measurements.csv",
    "run_summary.json",
)
BITWISE_INPUT_FILES = tuple(name for name in INPUT_FILES if name != "run_summary.json")


def f17(value: float) -> str:
    return format(float(value), FMT)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"empty rows: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig: plt.Figure, path: Path, dpi: int | None = None) -> None:
    fig.savefig(path, dpi=dpi)
    if path.suffix == ".svg":
        source = path.read_text(encoding="utf-8")
        cleaned = "\n".join(line.rstrip() for line in source.splitlines()) + "\n"
        path.write_text(cleaned, encoding="utf-8")


def ensure_output_dirs(allow_overwrite: bool) -> None:
    for directory in (PROCESSED, FIGURES, REPORTS, LOGS):
        if directory.exists() and any(directory.iterdir()) and not allow_overwrite:
            raise SystemExit(f"EXECUTION_FAILED: 出力先が空ではない: {directory}")
        directory.mkdir(parents=True, exist_ok=True)


def band_id(value: float) -> str:
    for band in CONFIG["q_resolution_bands"]:
        lower = band["lower"]
        upper = band["upper"]
        if (lower is None or value >= float(lower)) and (upper is None or value < float(upper)):
            return band["id"]
    raise ValueError(value)


def s4_new_dirs(B0: np.ndarray, Bdom: np.ndarray) -> np.ndarray:
    """論文7原本s4_new_dirsと同じ二行の定義。"""
    residual = Bdom - B0 @ (B0.T @ Bdom)
    basis, _ = np.linalg.qr(residual)
    return basis[:, :2]


def subspace_metrics(
    A: np.ndarray,
    PA: np.ndarray,
    B: np.ndarray,
    PB: np.ndarray,
) -> dict[str, float]:
    singular = np.clip(np.linalg.svd(A.T @ B, compute_uv=False), -1.0, 1.0)
    angles = np.arccos(singular)
    return {
        "theta_1_rad": float(angles[0]),
        "theta_2_rad": float(angles[1]),
        "maximum_principal_angle_rad": float(np.max(angles)),
        "overlap": float(np.trace(PA @ PB) / 2.0),
        "projector_distance": float(np.linalg.norm(PA - PB) / math.sqrt(2.0)),
        "minimum_singular_value": float(np.min(singular)),
    }


def input_audit() -> tuple[Path, dict]:
    run_dirs = [A2A / "raw" / run_id for run_id in CONFIG["input_runs"]]
    for run_dir in run_dirs:
        if not run_dir.is_dir():
            raise SystemExit(f"INPUT_MISSING: {run_dir}")
        for name in INPUT_FILES:
            if not (run_dir / name).is_file():
                raise SystemExit(f"INPUT_MISSING: {run_dir / name}")

    comparisons = {}
    hashes = {}
    for name in INPUT_FILES:
        path1 = run_dirs[0] / name
        path2 = run_dirs[1] / name
        hash1 = sha256(path1)
        hash2 = sha256(path2)
        hashes[name] = {
            "exec1": hash1,
            "exec2": hash2,
            "byte_identical": hash1 == hash2,
        }
        if name in BITWISE_INPUT_FILES:
            comparisons[name] = hash1 == hash2
    if not all(comparisons.values()):
        failed = [name for name, passed in comparisons.items() if not passed]
        raise SystemExit(f"INPUT_MISMATCH: A2a exec1/exec2 {failed}")

    summary = json.loads((run_dirs[0] / "run_summary.json").read_text(encoding="utf-8"))
    summary_exec2 = json.loads(
        (run_dirs[1] / "run_summary.json").read_text(encoding="utf-8")
    )
    summary_without_run_id = {key: value for key, value in summary.items() if key != "run_id"}
    summary_exec2_without_run_id = {
        key: value for key, value in summary_exec2.items() if key != "run_id"
    }
    if summary_without_run_id != summary_exec2_without_run_id:
        raise SystemExit("INPUT_MISMATCH: A2a run summaries differ beyond run_id")
    required_summary = {
        "status": "COMPLETED",
        "n": CONFIG["n"],
        "initial_state_rule": "Z0 = v.copy()",
    }
    if summary.get("status") != required_summary["status"] or summary.get("n") != required_summary["n"]:
        raise SystemExit("INPUT_MISMATCH: A2a summary status/N")
    if summary.get("explicit_state_seed_added") is not False:
        raise SystemExit("INPUT_MISMATCH: A2a explicit seed is not OFF")
    if summary.get("zero_closure_kernel_seed_called") is not False:
        raise SystemExit("INPUT_MISMATCH: kernel seed function was called")

    manifest = {
        "stage": CONFIG["stage"],
        "status": "INPUTS_VERIFIED",
        "input_runs": [str(path.relative_to(PAPER8)) for path in run_dirs],
        "all_numerical_input_files_byte_identical_between_execs": all(
            comparisons.values()
        ),
        "run_summaries_identical_except_run_id": True,
        "input_hashes": hashes,
        "input_summary_checks": {
            "status": summary.get("status"),
            "n": summary.get("n"),
            "initial_state_rule": summary.get("initial_state_rule"),
            "explicit_state_seed_added": summary.get("explicit_state_seed_added"),
            "zero_closure_kernel_seed_called": summary.get("zero_closure_kernel_seed_called"),
        },
    }
    return run_dirs[0], manifest


def analyze(run_dir: Path, manifest: dict) -> dict:
    steps = np.load(run_dir / "dominant_plane_steps.npy")
    Bdom = np.load(run_dir / "dominant_plane_values.npy")
    if steps.ndim != 1 or Bdom.shape != (len(steps), 10, 2):
        raise SystemExit(f"INPUT_MISMATCH: steps/Bdom shape {steps.shape} {Bdom.shape}")
    if len(np.unique(steps)) != len(steps) or not np.all(np.diff(steps) > 0):
        raise SystemExit("INPUT_MISMATCH: dominant-plane steps are not strictly increasing")
    if int(steps[0]) != 0 or int(steps[-1]) != CONFIG["max_step"]:
        raise SystemExit("INPUT_MISMATCH: dominant-plane time range")
    if not np.all(np.isfinite(Bdom)):
        raise SystemExit("INPUT_MISMATCH: nonfinite Bdom")

    f_rows = load_csv(run_dir / "f_timeseries.csv")
    q_rows = load_csv(run_dir / "q_timeseries.csv")
    passage_rows = load_csv(run_dir / "first_passage_measurements.csv")
    f_by_step = {int(row["step"]): float(row["f"]) for row in f_rows}
    q_by_step = {int(row["step"]): row for row in q_rows}
    if list(sorted(f_by_step)) != list(range(CONFIG["max_step"] + 1)):
        raise SystemExit("INPUT_MISMATCH: f does not cover all integer steps")
    if list(q_by_step) != [int(step) for step in steps]:
        raise SystemExit("INPUT_MISMATCH: q steps and Bdom steps differ")

    crossing = next(step for step in range(CONFIG["max_step"] + 1) if f_by_step[step] > 0.05)
    B0_index = int(np.flatnonzero(steps == 0)[0])
    B0 = Bdom[B0_index].copy()

    D34 = np.asarray([s4_new_dirs(B0, basis) for basis in Bdom])
    P34 = np.einsum("nmi,nki->nmk", D34, D34)
    bdom_orth = np.linalg.norm(np.swapaxes(Bdom, 1, 2) @ Bdom - np.eye(2), axis=(1, 2))
    d34_orth = np.linalg.norm(np.swapaxes(D34, 1, 2) @ D34 - np.eye(2), axis=(1, 2))
    p34_idem = np.linalg.norm(P34 @ P34 - P34, axis=(1, 2))
    p34_sym = np.linalg.norm(P34 - np.swapaxes(P34, 1, 2), axis=(1, 2))

    late_cfg = CONFIG["late_reference"]
    late_mask = (steps >= late_cfg["start_step"]) & (steps <= late_cfg["end_step"])
    if not np.any(late_mask):
        raise SystemExit("ANALYSIS_FAILED: no late-reference points")
    Pmean = np.mean(P34[late_mask], axis=0)
    evals, evecs = np.linalg.eigh(0.5 * (Pmean + Pmean.T))
    order = np.argsort(evals)[::-1]
    late_eigenvalues = evals[order]
    Dlate = evecs[:, order[:2]]
    Plate = Dlate @ Dlate.T

    seeded_path = A2C / "processed" / "early_vs_late_direction_overlap.csv"
    seeded_by_step = {}
    if seeded_path.is_file():
        seeded_by_step = {int(row["step"]): row for row in load_csv(seeded_path)}
    seeded_D_path = A2C / "raw" / "D34_basis_all_steps.npy"
    seeded_P_path = A2C / "raw" / "P34_all_steps.npy"
    seeded_D_at_samples = None
    seeded_P_at_samples = None
    if seeded_D_path.is_file() and seeded_P_path.is_file():
        seeded_D_all = np.load(seeded_D_path)
        seeded_P_all = np.load(seeded_P_path)
        expected_D_shape = (CONFIG["max_step"] + 1, 10, 2)
        expected_P_shape = (CONFIG["max_step"] + 1, 10, 10)
        if seeded_D_all.shape != expected_D_shape or seeded_P_all.shape != expected_P_shape:
            raise SystemExit("INPUT_MISMATCH: seeded A2c D34/P34 shape")
        seeded_D_at_samples = seeded_D_all[steps]
        seeded_P_at_samples = seeded_P_all[steps]
        manifest["seeded_A2c_comparison_inputs"] = {
            "D34_basis_all_steps": {
                "path": str(seeded_D_path.relative_to(PAPER8)),
                "sha256": sha256(seeded_D_path),
            },
            "P34_all_steps": {
                "path": str(seeded_P_path.relative_to(PAPER8)),
                "sha256": sha256(seeded_P_path),
            },
        }

    time_rows = []
    numeric_late = []
    ratios = []
    for index, step_value in enumerate(steps):
        step = int(step_value)
        qrow = q_by_step[step]
        ratio = min(float(qrow["q3_over_q1"]), float(qrow["q4_over_q1"]))
        ratios.append(ratio)
        metrics = subspace_metrics(D34[index], P34[index], Dlate, Plate)
        numeric_late.append(metrics)
        seeded = seeded_by_step.get(step)
        cross_trajectory = (
            subspace_metrics(
                D34[index],
                P34[index],
                seeded_D_at_samples[index],
                seeded_P_at_samples[index],
            )
            if seeded_D_at_samples is not None and seeded_P_at_samples is not None
            else None
        )
        time_rows.append({
            "sample_index": index,
            "step": step,
            "f": f17(f_by_step[step]),
            "q3_over_q1": qrow["q3_over_q1"],
            "q4_over_q1": qrow["q4_over_q1"],
            "min_q3_q4_over_q1": f17(ratio),
            "q_resolution_band": band_id(ratio),
            "theta_1_rad": f17(metrics["theta_1_rad"]),
            "theta_2_rad": f17(metrics["theta_2_rad"]),
            "maximum_principal_angle_rad": f17(metrics["maximum_principal_angle_rad"]),
            "overlap": f17(metrics["overlap"]),
            "projector_distance": f17(metrics["projector_distance"]),
            "minimum_singular_value": f17(metrics["minimum_singular_value"]),
            "Bdom_orthogonality_error": f17(bdom_orth[index]),
            "D34_orthogonality_error": f17(d34_orth[index]),
            "P34_idempotence_error": f17(p34_idem[index]),
            "P34_symmetry_error": f17(p34_sym[index]),
            "seeded_A2c_overlap_at_same_step": (
                seeded["overlap"] if seeded is not None else ""
            ),
            "seedless_minus_seeded_overlap": (
                f17(metrics["overlap"] - float(seeded["overlap"]))
                if seeded is not None else ""
            ),
            "seedless_vs_seeded_same_step_overlap": (
                f17(cross_trajectory["overlap"])
                if cross_trajectory is not None else ""
            ),
            "seedless_vs_seeded_same_step_max_angle_rad": (
                f17(cross_trajectory["maximum_principal_angle_rad"])
                if cross_trajectory is not None else ""
            ),
            "seedless_vs_seeded_same_step_projector_distance": (
                f17(cross_trajectory["projector_distance"])
                if cross_trajectory is not None else ""
            ),
        })

    ratios_array = np.asarray(ratios)
    continuity_rows = []
    continuity_numeric = []
    for index in range(len(steps) - 1):
        step_from = int(steps[index])
        step_to = int(steps[index + 1])
        delta_step = step_to - step_from
        metrics = subspace_metrics(
            D34[index], P34[index], D34[index + 1], P34[index + 1]
        )
        continuity_numeric.append((step_from, step_to, delta_step, metrics))
        continuity_rows.append({
            "sample_index_from": index,
            "sample_index_to": index + 1,
            "step_from": step_from,
            "step_to": step_to,
            "delta_step": delta_step,
            "q_resolution_band_to": band_id(ratios_array[index + 1]),
            "theta_1_rad": f17(metrics["theta_1_rad"]),
            "theta_2_rad": f17(metrics["theta_2_rad"]),
            "maximum_principal_angle_rad": f17(metrics["maximum_principal_angle_rad"]),
            "maximum_principal_angle_per_step_descriptive": f17(
                metrics["maximum_principal_angle_rad"] / delta_step
            ),
            "overlap": f17(metrics["overlap"]),
            "projector_distance": f17(metrics["projector_distance"]),
            "projector_distance_per_step_descriptive": f17(
                metrics["projector_distance"] / delta_step
            ),
        })

    passage_out = []
    index_by_step = {int(step): index for index, step in enumerate(steps)}
    for passage in passage_rows:
        if passage["status"] != "found":
            continue
        step = int(passage["first_passage_step"])
        if step not in index_by_step:
            raise SystemExit(f"ANALYSIS_FAILED: first-passage step not sampled: {step}")
        index = index_by_step[step]
        metrics = numeric_late[index]
        passage_out.append({
            "level_index": passage["level_index"],
            "level": passage["level"],
            "level_label": passage["level_label"],
            "step": step,
            "f_at_step": passage["f_at_first_passage"],
            "min_q3_q4_over_q1": f17(ratios_array[index]),
            "q_resolution_band": band_id(ratios_array[index]),
            "D34_vs_late_overlap": f17(metrics["overlap"]),
            "D34_vs_late_max_angle_rad": f17(metrics["maximum_principal_angle_rad"]),
            "D34_vs_late_projector_distance": f17(metrics["projector_distance"]),
        })

    band_rows = []
    continuity_by_to = {row[1]: row[3] for row in continuity_numeric}
    for band in CONFIG["q_resolution_bands"]:
        indices = np.flatnonzero([
            band_id(value) == band["id"] for value in ratios_array
        ])
        if len(indices) == 0:
            band_rows.append({
                "q_resolution_band": band["id"],
                "point_count": 0,
                "first_step": "",
                "last_step": "",
                "median_D34_vs_late_overlap": "",
                "median_D34_vs_late_max_angle_rad": "",
                "maximum_sample_to_sample_angle_rad": "",
            })
            continue
        band_steps = [int(steps[index]) for index in indices]
        cont = [continuity_by_to[step] for step in band_steps if step in continuity_by_to]
        band_rows.append({
            "q_resolution_band": band["id"],
            "point_count": len(indices),
            "first_step": min(band_steps),
            "last_step": max(band_steps),
            "median_D34_vs_late_overlap": f17(np.median([
                numeric_late[index]["overlap"] for index in indices
            ])),
            "median_D34_vs_late_max_angle_rad": f17(np.median([
                numeric_late[index]["maximum_principal_angle_rad"] for index in indices
            ])),
            "maximum_sample_to_sample_angle_rad": (
                f17(max(item["maximum_principal_angle_rad"] for item in cont))
                if cont else ""
            ),
        })

    resolved = np.flatnonzero(ratios_array >= CONFIG["q_resolution_threshold"])
    pre_resolved = resolved[steps[resolved] < crossing]
    if len(resolved) == 0:
        classification = CONFIG["classification"]["insufficient"]
        evaluation_start_step = None
        evaluation_overlap_median = None
        evaluation_angle_median = None
    else:
        evaluation_start_step = int(steps[resolved[0]])
        evaluation_indices = np.arange(resolved[0], len(steps))
        evaluation_overlap_median = float(np.median([
            numeric_late[index]["overlap"] for index in evaluation_indices
        ]))
        evaluation_angle_median = float(np.median([
            numeric_late[index]["maximum_principal_angle_rad"]
            for index in evaluation_indices
        ]))
        cls = CONFIG["classification"]
        if (
            evaluation_overlap_median >= cls["continuous_late_overlap_median_min"]
            and evaluation_angle_median <= cls["continuous_late_max_angle_median_max_rad"]
        ):
            classification = "SAMPLED_FIRST_DIRECTIONS_CONTINUOUS"
        else:
            classification = cls["fallback_if_resolved"]

    rapid_start, rapid_end = CONFIG["display_windows"]["crossing_zoom"]
    rapid = [
        item for item in continuity_numeric
        if rapid_start <= item[0] and item[1] <= rapid_end
    ]
    seeded_pre_overlaps = [
        float(seeded_by_step[int(steps[index])]["overlap"])
        for index in pre_resolved
        if int(steps[index]) in seeded_by_step
    ]
    seedless_pre_overlaps = [
        numeric_late[index]["overlap"] for index in pre_resolved
    ]
    cross_trajectory_overlap = np.asarray([
        float(row["seedless_vs_seeded_same_step_overlap"])
        if row["seedless_vs_seeded_same_step_overlap"] else np.nan
        for row in time_rows
    ])
    cross_trajectory_angle = np.asarray([
        float(row["seedless_vs_seeded_same_step_max_angle_rad"])
        if row["seedless_vs_seeded_same_step_max_angle_rad"] else np.nan
        for row in time_rows
    ])

    def cross_trajectory_window(mask: np.ndarray) -> dict | None:
        valid = mask & np.isfinite(cross_trajectory_overlap)
        if not np.any(valid):
            return None
        return {
            "sample_count": int(np.sum(valid)),
            "overlap_median": float(np.median(cross_trajectory_overlap[valid])),
            "overlap_minimum": float(np.min(cross_trajectory_overlap[valid])),
            "max_angle_median_rad": float(np.median(cross_trajectory_angle[valid])),
            "max_angle_maximum_rad": float(np.max(cross_trajectory_angle[valid])),
        }

    cross_windows = {
        "pre_crossing_resolved": cross_trajectory_window(
            (steps < crossing) & (ratios_array >= CONFIG["q_resolution_threshold"])
        ),
        "transition_crossing_to_1799": cross_trajectory_window(
            (steps >= crossing) & (steps < 1800)
        ),
        "late_reference_1800_to_2500": cross_trajectory_window(
            (steps >= 1800) & (steps <= 2500)
        ),
        "post_2500_to_5000": cross_trajectory_window(
            (steps > 2500) & (steps <= 5000)
        ),
    }
    resolved_cross_indices = np.flatnonzero(
        (steps >= evaluation_start_step)
        & np.isfinite(cross_trajectory_overlap)
    ) if evaluation_start_step is not None else np.asarray([], dtype=int)
    below_095 = resolved_cross_indices[
        cross_trajectory_overlap[resolved_cross_indices] < 0.95
    ]
    below_05 = resolved_cross_indices[
        cross_trajectory_overlap[resolved_cross_indices] < 0.5
    ]

    summary = {
        "stage": CONFIG["stage"],
        "status": "ANALYSIS_COMPLETE",
        "trajectory": "N=5 float64 Z0=v; both explicit seeds OFF",
        "new_dynamics_executed": False,
        "input_execs_bitwise_identical": True,
        "sample_count": int(len(steps)),
        "step_range": [int(steps[0]), int(steps[-1])],
        "step_gap_values": sorted(int(value) for value in np.unique(np.diff(steps))),
        "step_gap_counts": {
            str(int(value)): int(np.sum(np.diff(steps) == value))
            for value in np.unique(np.diff(steps))
        },
        "crossing_step": crossing,
        "B0_source": "saved Bdom at step 0",
        "late_reference_window": [
            late_cfg["start_step"], late_cfg["end_step"]
        ],
        "late_reference_sample_count": int(np.sum(late_mask)),
        "late_mean_projector_eigenvalues": [
            float(value) for value in late_eigenvalues
        ],
        "late_window_overlap_median": float(np.median([
            numeric_late[index]["overlap"] for index in np.flatnonzero(late_mask)
        ])),
        "late_window_overlap_minimum": float(np.min([
            numeric_late[index]["overlap"] for index in np.flatnonzero(late_mask)
        ])),
        "resolved_threshold": CONFIG["q_resolution_threshold"],
        "resolved_sample_count": int(len(resolved)),
        "evaluation_start_step": evaluation_start_step,
        "evaluation_overlap_median": evaluation_overlap_median,
        "evaluation_max_angle_median_rad": evaluation_angle_median,
        "pre_crossing_resolved_sample_count": int(len(pre_resolved)),
        "pre_crossing_resolved_step_range": (
            [int(steps[pre_resolved[0]]), int(steps[pre_resolved[-1]])]
            if len(pre_resolved) else None
        ),
        "pre_crossing_resolved_late_overlap_median": (
            float(np.median(seedless_pre_overlaps))
            if seedless_pre_overlaps else None
        ),
        "pre_crossing_resolved_late_overlap_range": (
            [float(np.min(seedless_pre_overlaps)), float(np.max(seedless_pre_overlaps))]
            if seedless_pre_overlaps else None
        ),
        "seeded_A2c_same_sample_pre_crossing_overlap_median": (
            float(np.median(seeded_pre_overlaps))
            if seeded_pre_overlaps else None
        ),
        "seedless_minus_seeded_same_sample_pre_crossing_overlap_median": (
            float(np.median(np.asarray(seedless_pre_overlaps) - np.asarray(seeded_pre_overlaps)))
            if len(seedless_pre_overlaps) == len(seeded_pre_overlaps)
            and seedless_pre_overlaps else None
        ),
        "seedless_vs_seeded_same_step_windows": cross_windows,
        "seedless_vs_seeded_first_resolved_step_overlap_below_0_95": (
            int(steps[below_095[0]]) if len(below_095) else None
        ),
        "seedless_vs_seeded_first_resolved_step_overlap_below_0_5": (
            int(steps[below_05[0]]) if len(below_05) else None
        ),
        "crossing_window_sample_pair_count": len(rapid),
        "crossing_window_max_sample_to_sample_angle_rad": (
            max(item[3]["maximum_principal_angle_rad"] for item in rapid)
            if rapid else None
        ),
        "crossing_window_max_sample_to_sample_projector_distance": (
            max(item[3]["projector_distance"] for item in rapid)
            if rapid else None
        ),
        "crossing_window_max_angle_per_step_descriptive": (
            max(item[3]["maximum_principal_angle_rad"] / item[2] for item in rapid)
            if rapid else None
        ),
        "classification": classification,
        "classification_scope": (
            "sampled early-vs-late lineage only; no Tperp and no one-step maximum"
        ),
        "numerical_health": {
            "max_Bdom_orthogonality_error": float(np.max(bdom_orth)),
            "max_D34_orthogonality_error": float(np.max(d34_orth)),
            "max_P34_idempotence_error": float(np.max(p34_idem)),
            "max_P34_symmetry_error": float(np.max(p34_sym)),
            "all_finite": bool(
                np.all(np.isfinite(D34))
                and np.all(np.isfinite(P34))
            ),
        },
        "limitations": [
            "dominant planes are sampled mainly every 5 steps",
            "sample-to-sample maxima are not one-step maxima",
            "Tperp comparison at exact seedless t0 is not performed",
            "no single physical direction-establishment step is selected",
            "H1/H2/H0 is not judged"
        ],
    }

    write_csv(PROCESSED / "seedless_sampled_lineage_timeseries.csv", time_rows)
    write_csv(PROCESSED / "seedless_sampled_continuity.csv", continuity_rows)
    write_csv(PROCESSED / "seedless_lineage_by_f_level.csv", passage_out)
    write_csv(PROCESSED / "seedless_lineage_by_q_band.csv", band_rows)
    (PROCESSED / "seedless_lineage_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (LOGS / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    make_figures(steps, time_rows, continuity_rows, crossing)
    make_report(summary)
    return summary


def make_figures(
    steps: np.ndarray,
    time_rows: list[dict],
    continuity_rows: list[dict],
    crossing: int,
) -> None:
    overlap = np.asarray([float(row["overlap"]) for row in time_rows])
    max_angle = np.asarray([
        float(row["maximum_principal_angle_rad"]) for row in time_rows
    ])
    ratios = np.asarray([
        float(row["min_q3_q4_over_q1"]) for row in time_rows
    ])
    seeded = np.asarray([
        float(row["seeded_A2c_overlap_at_same_step"])
        if row["seeded_A2c_overlap_at_same_step"] else np.nan
        for row in time_rows
    ])
    cross_trajectory = np.asarray([
        float(row["seedless_vs_seeded_same_step_overlap"])
        if row["seedless_vs_seeded_same_step_overlap"] else np.nan
        for row in time_rows
    ])

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    axes[0].plot(steps, overlap, color="#1f77b4", lw=1.4, label="seedless sampled")
    if np.any(np.isfinite(seeded)):
        axes[0].plot(steps, seeded, color="#ff7f0e", lw=0.9, alpha=0.8,
                     label="seeded A2c at same steps")
    if np.any(np.isfinite(cross_trajectory)):
        axes[0].plot(
            steps,
            cross_trajectory,
            color="#111111",
            lw=0.9,
            ls="--",
            alpha=0.9,
            label="seedless vs seeded at same step",
        )
    axes[0].axhline(0.95, color="0.5", ls="--", lw=0.8)
    axes[0].set_ylabel("overlap with late D34")
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].legend(loc="lower right")

    axes[1].plot(steps, max_angle, color="#d62728", lw=1.1)
    axes[1].set_ylabel("max principal angle [rad]")

    positive = ratios > 0.0
    axes[2].semilogy(steps[positive], ratios[positive], color="#2ca02c", lw=1.0)
    axes[2].axhline(CONFIG["q_resolution_threshold"], color="0.3", ls="--", lw=0.9)
    axes[2].set_ylabel("min(q3,q4)/q1")
    axes[2].set_xlabel("step")

    for axis in axes:
        axis.axvline(crossing, color="black", ls=":", lw=1.0)
        axis.axvspan(
            CONFIG["late_reference"]["start_step"],
            CONFIG["late_reference"]["end_step"],
            color="0.8", alpha=0.25,
        )
        axis.grid(alpha=0.2)
    fig.suptitle("Complete-seedless N=5: sampled D34 lineage")
    fig.tight_layout()
    for suffix in ("png", "svg"):
        save_figure(
            fig,
            FIGURES / f"figure01_seedless_early_vs_late_lineage.{suffix}",
            dpi=180 if suffix == "png" else None,
        )
    plt.close(fig)

    step_to = np.asarray([int(row["step_to"]) for row in continuity_rows])
    gaps = np.asarray([int(row["delta_step"]) for row in continuity_rows])
    angles = np.asarray([
        float(row["maximum_principal_angle_rad"]) for row in continuity_rows
    ])
    angles_per_step = np.asarray([
        float(row["maximum_principal_angle_per_step_descriptive"])
        for row in continuity_rows
    ])
    distances = np.asarray([
        float(row["projector_distance"]) for row in continuity_rows
    ])

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    axes[0].plot(step_to, angles, color="#d62728", lw=0.9)
    axes[0].set_ylabel("sample-pair angle [rad]")
    axes[1].plot(step_to, angles_per_step, color="#9467bd", lw=0.9)
    axes[1].set_ylabel("angle / delta step")
    axes[2].plot(step_to, distances, color="#17becf", lw=0.9)
    axes[2].set_ylabel("projector distance")
    axes[2].set_xlabel("step to")
    for axis in axes:
        axis.axvline(crossing, color="black", ls=":", lw=1.0)
        axis.grid(alpha=0.2)
    twin = axes[0].twinx()
    twin.step(step_to, gaps, color="0.4", lw=0.5, alpha=0.5)
    twin.set_ylabel("delta step", color="0.4")
    fig.suptitle("Complete-seedless N=5: sampled subspace rotation")
    fig.tight_layout()
    for suffix in ("png", "svg"):
        save_figure(
            fig,
            FIGURES / f"figure02_seedless_sampled_rotation.{suffix}",
            dpi=180 if suffix == "png" else None,
        )
    plt.close(fig)


def make_report(summary: dict) -> None:
    pre_range = summary["pre_crossing_resolved_step_range"]
    health = summary["numerical_health"]
    cross = summary["seedless_vs_seeded_same_step_windows"]
    cross_pre = cross["pre_crossing_resolved"]
    cross_transition = cross["transition_crossing_to_1799"]
    cross_late = cross["late_reference_1800_to_2500"]
    cross_post = cross["post_2500_to_5000"]
    limitations = "\n".join(f"- {item}" for item in summary["limitations"])
    report = f"""# Stage A2d 保存済み完全無seed N=5方向系譜 報告書

## 実行状態

**{summary["classification"]}**

本解析はStage A2aで保存済みの完全無seed軌道だけを使用した。
新しい時間発展、状態更新、seed注入、横摂動、Benettin再投入は行っていない。

## 入力と標本

- 軌道: `{summary["trajectory"]}`
- A2a exec 1/2 bitwise一致: `{summary["input_execs_bitwise_identical"]}`
- 保存点数: `{summary["sample_count"]}`
- step範囲: `{summary["step_range"][0]}..{summary["step_range"][1]}`
- step間隔: `{summary["step_gap_counts"]}`
- crossing: `{summary["crossing_step"]}`
- B0: `{summary["B0_source"]}`

## 後期代表D34

- 後期窓: `{summary["late_reference_window"][0]}..{summary["late_reference_window"][1]}`
- 使用標本数: `{summary["late_reference_sample_count"]}`
- 後期窓内overlap中央値: `{summary["late_window_overlap_median"]:.12f}`
- 後期窓内overlap最小値: `{summary["late_window_overlap_minimum"]:.12f}`

## 急拡大前の解像可能D34

- 解像基準: `min(q3,q4)/q1 >= {summary["resolved_threshold"]:.0e}`
- crossing前の標本数: `{summary["pre_crossing_resolved_sample_count"]}`
- step範囲: `{pre_range}`
- 無seedD34対late overlap中央値: `{summary["pre_crossing_resolved_late_overlap_median"]:.12f}`
- 無seedD34対late overlap範囲: `{summary["pre_crossing_resolved_late_overlap_range"]}`
- 同じ保存stepでのseedありA2c overlap中央値:
  `{summary["seeded_A2c_same_sample_pre_crossing_overlap_median"]:.12f}`
- 無seed−seedあり overlap差の中央値:
  `{summary["seedless_minus_seeded_same_sample_pre_crossing_overlap_median"]:.12e}`

急拡大前の解像可能部分空間は、急拡大後の後期代表部分空間と
単純に同一とは分類されない。完全無seed条件でも、既存seedありA2cと
同程度の低い早期対後期overlapが得られた。

## 無seed軌道とseedあり軌道の同一step直接比較

これは各軌道の後期代表方向へのoverlap比較ではなく、同一stepにおける
二つの `P34(t)` 自体のoverlapである。

| 区間 | 標本数 | overlap中央値 | overlap最小値 | 最大主角中央値(rad) | 最大主角最大値(rad) |
|---|---:|---:|---:|---:|---:|
| crossing前・解像可能 | {cross_pre["sample_count"]} | {cross_pre["overlap_median"]:.12f} | {cross_pre["overlap_minimum"]:.12f} | {cross_pre["max_angle_median_rad"]:.12e} | {cross_pre["max_angle_maximum_rad"]:.12e} |
| crossing〜1799 | {cross_transition["sample_count"]} | {cross_transition["overlap_median"]:.12f} | {cross_transition["overlap_minimum"]:.12f} | {cross_transition["max_angle_median_rad"]:.12e} | {cross_transition["max_angle_maximum_rad"]:.12e} |
| 1800〜2500 | {cross_late["sample_count"]} | {cross_late["overlap_median"]:.12f} | {cross_late["overlap_minimum"]:.12f} | {cross_late["max_angle_median_rad"]:.12e} | {cross_late["max_angle_maximum_rad"]:.12e} |
| 2500より後 | {cross_post["sample_count"]} | {cross_post["overlap_median"]:.12f} | {cross_post["overlap_minimum"]:.12f} | {cross_post["max_angle_median_rad"]:.12e} | {cross_post["max_angle_maximum_rad"]:.12e} |

- overlapが0.95を初めて下回る解像後step:
  `{summary["seedless_vs_seeded_first_resolved_step_overlap_below_0_95"]}`
- overlapが0.5を初めて下回る解像後step:
  `{summary["seedless_vs_seeded_first_resolved_step_overlap_below_0_5"]}`

急拡大前から三方向閉包の成立直後まで、無seedとseedありのD34部分空間は
ほぼ同一である。したがって、明示的初期seedは最初に生成される方向部分空間を
選択していない。step 2500以後の分岐は、生成方向の初期選択とは分離して
長期軌道差として扱う。

## 標本間回転

- crossing窓の標本対数: `{summary["crossing_window_sample_pair_count"]}`
- 最大標本間主角: `{summary["crossing_window_max_sample_to_sample_angle_rad"]:.12f}` rad
- 最大標本間射影距離:
  `{summary["crossing_window_max_sample_to_sample_projector_distance"]:.12f}`
- `angle/delta_step` の最大記述値:
  `{summary["crossing_window_max_angle_per_step_descriptive"]:.12f}` rad/step

これらは不等間隔標本間の記述量であり、未保存stepを含む1-step最大値ではない。

## 数値健全性

- Bdom直交誤差最大: `{health["max_Bdom_orthogonality_error"]:.3e}`
- D34直交誤差最大: `{health["max_D34_orthogonality_error"]:.3e}`
- P34冪等誤差最大: `{health["max_P34_idempotence_error"]:.3e}`
- P34対称誤差最大: `{health["max_P34_symmetry_error"]:.3e}`
- 全配列finite: `{health["all_finite"]}`

## 判定

分類: **{summary["classification"]}**

分類範囲: `{summary["classification_scope"]}`

保存済み完全無seed軌道においても、急拡大前の解像可能D34は後期D34の
単純な微小振幅版ではない。観測された低overlapと標本間回転は、
方向部分空間が急拡大過程で回転・混合を伴って再編される記述を支持する。

## データだけでは言えないこと

{limitations}
"""
    (REPORTS / "stage_A2d_seedless_direction_lineage_N5_report.md").write_text(
        report, encoding="utf-8"
    )


def main() -> None:
    allow_overwrite = "--overwrite" in sys.argv[1:]
    ensure_output_dirs(allow_overwrite)
    run_dir, manifest = input_audit()
    summary = analyze(run_dir, manifest)
    print(
        "ANALYSIS_COMPLETE "
        f"classification={summary['classification']} "
        f"pre_overlap={summary['pre_crossing_resolved_late_overlap_median']:.12f} "
        f"samples={summary['sample_count']}"
    )


if __name__ == "__main__":
    main()
