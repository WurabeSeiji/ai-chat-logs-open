#!/usr/bin/env python3
"""既存D34/P34/Tperpを主角・射影量で比較し、固定分類規則を適用する。"""

from __future__ import annotations

import csv
import hashlib
import importlib
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = HERE / "raw"
PROCESSED = HERE / "processed"
LOGS = HERE / "logs"
CFG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
EXPECTED = json.loads((HERE / "expected_hashes.json").read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv_md(stem: str, fields: list[str], rows: list[dict]) -> None:
    csv_path = PROCESSED / f"{stem}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    md_path = PROCESSED / f"{stem}.md"
    with md_path.open("w", encoding="utf-8") as fh:
        fh.write("| " + " | ".join(fields) + " |\n")
        fh.write("|" + "|".join(["---"] * len(fields)) + "|\n")
        for row in rows:
            fh.write("| " + " | ".join(str(row.get(c, "")).replace("|", "\\|") for c in fields) + " |\n")


def f17(x: float) -> str:
    return format(float(x), ".17e")


def import_existing_functions():
    src = EXPECTED["sources"]
    for item in src.values():
        path = REPO / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise SystemExit(f"SOURCE_MISMATCH: {path}")
    engine = (REPO / src["run_n_scaling_lowrank_v1.py"]["path"]).parent
    code = (REPO / src["run_n300_dimension_saturation_v2.py"]["path"]).parent
    p7code = (REPO / src["run_paper7_5color_timeseries.py"]["path"]).parent
    for p in (engine, code, p7code):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    sat = importlib.import_module("run_n300_dimension_saturation_v2")
    color = importlib.import_module("run_paper7_5color_timeseries")
    return sat.principal_angles, color.align_2d


def band_id(value: float) -> str:
    for b in CFG["q_resolution_bands"]:
        lower_ok = True if b["lower"] is None else value >= float(b["lower"])
        upper_ok = True if b["upper"] is None else value < float(b["upper"])
        if lower_ok and upper_ok:
            return b["id"]
    raise ValueError(value)


def compare_subspaces(A: np.ndarray, PA: np.ndarray, B: np.ndarray, PB: np.ndarray, principal_angles) -> dict:
    theta_deg = principal_angles(A, B)
    theta_rad = np.radians(theta_deg)
    singular = np.linalg.svd(A.T @ B, compute_uv=False)
    return {
        "theta_1_rad": float(theta_rad[0]),
        "theta_2_rad": float(theta_rad[1]),
        "theta_1_deg": float(theta_deg[0]),
        "theta_2_deg": float(theta_deg[1]),
        "overlap": float(np.trace(PA @ PB) / 2.0),
        "projector_distance": float(np.linalg.norm(PA - PB) / math.sqrt(2.0)),
        "minimum_singular_value": float(np.min(singular)),
        "maximum_principal_angle_rad": float(np.max(theta_rad)),
    }


def basis_column_relation(A: np.ndarray, B: np.ndarray, align_2d) -> dict:
    overlap = A.T @ B
    diagonal_score = abs(overlap[0, 0]) + abs(overlap[1, 1])
    swapped_score = abs(overlap[0, 1]) + abs(overlap[1, 0])
    swapped = swapped_score > diagonal_score
    matches = [(0, 1), (1, 0)] if swapped else [(0, 0), (1, 1)]
    sign_flip_count = sum(overlap[i, j] < 0 for i, j in matches)
    U, _, Vt = np.linalg.svd(overlap)
    polar = U @ Vt
    signed_permutations = []
    for swap in (False, True):
        P = np.array([[0.0, 1.0], [1.0, 0.0]]) if swap else np.eye(2)
        for s0 in (-1.0, 1.0):
            for s1 in (-1.0, 1.0):
                signed_permutations.append(P @ np.diag([s0, s1]))
    closest = min(signed_permutations, key=lambda S: np.linalg.norm(polar - S))
    cosine = float(np.clip(np.trace(closest.T @ polar) / 2.0, -1.0, 1.0))
    inside_rotation = float(math.acos(cosine))
    aligned = align_2d(A, B)
    return {
        "basis_column_swap": swapped,
        "sign_flip_count": int(sign_flip_count),
        "rotation_inside_D34_rad": inside_rotation,
        "raw_column_1_inner": float(overlap[0, 0]),
        "raw_column_2_inner": float(overlap[1, 1]),
        "aligned_column_1_abs_inner": float(abs(A[:, 0] @ aligned[:, 0])),
        "aligned_column_2_abs_inner": float(abs(A[:, 1] @ aligned[:, 1])),
    }


def main() -> None:
    replay = json.loads((LOGS / "replay_summary.json").read_text(encoding="utf-8"))
    trans = json.loads((LOGS / "transverse_reconstruction_summary.json").read_text(encoding="utf-8"))
    if replay.get("status") != "REPLAY_COMPLETE" or not replay.get("trajectory_f_bitwise_matches_stage_a0"):
        raise SystemExit("TRAJECTORY_MISMATCH: replay gate")
    if trans.get("status") != "TRANSVERSE_RECONSTRUCTED" or not trans.get("stage_a0_first_record_strings_all_match"):
        raise SystemExit("TRANSVERSE_RECONSTRUCTION_MISMATCH")
    if PROCESSED.exists() and any(PROCESSED.iterdir()):
        raise SystemExit("EXECUTION_FAILED: processed/が空ではないため上書きを拒否")
    PROCESSED.mkdir(parents=True, exist_ok=True)
    principal_angles, align_2d = import_existing_functions()

    D = np.load(RAW / "D34_basis_all_steps.npy")
    P = np.load(RAW / "P34_all_steps.npy")
    S4 = np.load(RAW / "S4_basis_all_steps.npy")
    PS4 = np.load(RAW / "PS4_all_steps.npy")
    Tbasis = np.load(RAW / "Tperp_real_basis.npy")
    Tproj = np.load(RAW / "Tperp_projectors.npy")
    metrics = load_csv(RAW / "trajectory_direction_metrics.csv")
    snapshots = load_csv(RAW / "direction_basis_snapshots.csv")
    passages = load_csv(RAW / "f_first_passage_levels.csv")
    trans_rows = load_csv(RAW / "transverse_direction_reconstruction.csv")
    max_step = CFG["trajectory"]["max_step"]
    if D.shape != (max_step + 1, 10, 2) or P.shape != (max_step + 1, 10, 10):
        raise SystemExit("EXECUTION_FAILED: basis array shape mismatch")

    late_start = CFG["late_reference"]["start_step"]
    late_end = CFG["late_reference"]["end_step"]
    Pmean = np.mean(P[late_start:late_end + 1], axis=0)
    evals, evecs = np.linalg.eigh(0.5 * (Pmean + Pmean.T))
    order = np.argsort(evals)[::-1]
    late_eigenvalues = evals[order]
    Dlate = evecs[:, order[:2]]
    Plate = Dlate @ Dlate.T
    np.save(RAW / "D34_late_basis.npy", Dlate)
    np.save(RAW / "D34_late_projector.npy", Plate)
    np.save(RAW / "D34_late_mean_projector.npy", Pmean)
    np.save(RAW / "D34_late_mean_projector_eigenvalues.npy", late_eigenvalues)

    # 1. snapshots
    write_csv_md("direction_basis_snapshots", list(snapshots[0]), snapshots)

    # 2. projector quality
    quality_fields = [
        "step", "f", "min_q3_q4_over_q1", "q_resolution_band",
        "Bdom_orthogonality_error", "D34_orthogonality_error",
        "S4_orthogonality_error", "P34_idempotence_error",
        "P34_symmetry_error", "S4_direct_equivalence_error",
        "dominant_eigenvalue_gap", "dominant_eigenvalue_gap_relative",
        "dominant_lifted_residual", "norm_error", "closure_error",
    ]
    quality_rows = []
    for r in metrics:
        quality_rows.append({
            **{k: r[k] for k in quality_fields if k not in {"q_resolution_band"}},
            "q_resolution_band": band_id(float(r["min_q3_q4_over_q1"])),
        })
    write_csv_md("direction_projector_quality", quality_fields, quality_rows)

    # 3. consecutive continuity and basis-column bookkeeping
    continuity_rows = []
    for step in range(max_step):
        comp = compare_subspaces(D[step], P[step], D[step + 1], P[step + 1], principal_angles)
        columns = basis_column_relation(D[step], D[step + 1], align_2d)
        continuity_rows.append({
            "step_from": step,
            "step_to": step + 1,
            "delta_step": 1,
            "f_from": metrics[step]["f"],
            "f_to": metrics[step + 1]["f"],
            "q_resolution_band_to": band_id(float(metrics[step + 1]["min_q3_q4_over_q1"])),
            **{k: f17(v) for k, v in comp.items()},
            "basis_column_swap": str(columns["basis_column_swap"]),
            "sign_flip_count": columns["sign_flip_count"],
            **{k: f17(v) for k, v in columns.items() if k not in {"basis_column_swap", "sign_flip_count"}},
        })
    continuity_fields = list(continuity_rows[0])
    write_csv_md("consecutive_subspace_continuity", continuity_fields, continuity_rows)

    # 4. all D34(t) vs D34_late
    late_rows = []
    late_comp_numeric = []
    for step in range(max_step + 1):
        comp = compare_subspaces(D[step], P[step], Dlate, Plate, principal_angles)
        late_comp_numeric.append(comp)
        late_rows.append({
            "step": step,
            "f": metrics[step]["f"],
            "min_q3_q4_over_q1": metrics[step]["min_q3_q4_over_q1"],
            "q_resolution_band": band_id(float(metrics[step]["min_q3_q4_over_q1"])),
            **{k: f17(v) for k, v in comp.items()},
        })
    late_fields = list(late_rows[0])
    write_csv_md("early_vs_late_direction_overlap", late_fields, late_rows)

    # 5. all D34(t) vs three unique original Tperp spans
    trans_overlap = np.empty((max_step + 1, Tbasis.shape[0]))
    early_trans_rows = []
    for step in range(max_step + 1):
        for seed in range(Tbasis.shape[0]):
            comp = compare_subspaces(D[step], P[step], Tbasis[seed], Tproj[seed], principal_angles)
            trans_overlap[step, seed] = comp["overlap"]
            early_trans_rows.append({
                "step": step,
                "seed": seed,
                "epsilons_sharing_this_direction": ";".join(f"{x:.0e}" for x in CFG["transverse"]["epsilons"]),
                "t0": CFG["transverse"]["t0"],
                "f": metrics[step]["f"],
                "min_q3_q4_over_q1": metrics[step]["min_q3_q4_over_q1"],
                "q_resolution_band": band_id(float(metrics[step]["min_q3_q4_over_q1"])),
                **{k: f17(v) for k, v in comp.items()},
            })
    early_trans_fields = list(early_trans_rows[0])
    write_csv_md("early_vs_transverse_overlap", early_trans_fields, early_trans_rows)

    # 6. D34_late vs each unique Tperp
    late_trans_rows = []
    late_trans_numeric = []
    for seed in range(Tbasis.shape[0]):
        comp = compare_subspaces(Dlate, Plate, Tbasis[seed], Tproj[seed], principal_angles)
        late_trans_numeric.append(comp)
        s4_error = float(np.linalg.norm(S4[CFG["transverse"]["t0"]].T @ Tbasis[seed]))
        late_trans_rows.append({
            "seed": seed,
            "epsilons_sharing_this_direction": ";".join(f"{x:.0e}" for x in CFG["transverse"]["epsilons"]),
            "t0": CFG["transverse"]["t0"],
            **{k: f17(v) for k, v in comp.items()},
            "same_t0_S4_vs_Tperp_orthogonality_error": f17(s4_error),
        })
    late_trans_fields = list(late_trans_rows[0])
    write_csv_md("late_direction_vs_transverse_overlap", late_trans_fields, late_trans_rows)

    # 7. f-level lineage
    flevel_rows = []
    for psg in passages:
        if psg["status"] != "found":
            flevel_rows.append({
                "level_index": psg["level_index"], "level": psg["level"], "level_label": psg["level_label"],
                "status": "not_found", "step": "",
            })
            continue
        step = int(psg["first_passage_step"])
        trans_values = trans_overlap[step]
        best_seed = int(np.argmax(trans_values))
        row = {
            "level_index": psg["level_index"],
            "level": psg["level"],
            "level_label": psg["level_label"],
            "status": "found",
            "step": step,
            "f_at_step": metrics[step]["f"],
            "min_q3_q4_over_q1": metrics[step]["min_q3_q4_over_q1"],
            "q_resolution_band": band_id(float(metrics[step]["min_q3_q4_over_q1"])),
            "D34_vs_late_overlap": f17(late_comp_numeric[step]["overlap"]),
            "D34_vs_late_max_angle_rad": f17(late_comp_numeric[step]["maximum_principal_angle_rad"]),
            "D34_vs_late_projector_distance": f17(late_comp_numeric[step]["projector_distance"]),
            "D34_vs_Tperp_seed0_overlap": f17(trans_values[0]),
            "D34_vs_Tperp_seed1_overlap": f17(trans_values[1]),
            "D34_vs_Tperp_seed2_overlap": f17(trans_values[2]),
            "D34_vs_Tperp_max_overlap": f17(trans_values[best_seed]),
            "D34_vs_Tperp_best_seed": best_seed,
        }
        flevel_rows.append(row)
    flevel_fields = list(next(r for r in flevel_rows if r["status"] == "found"))
    # not_found rows receive the complete schema
    flevel_rows = [{c: r.get(c, "") for c in flevel_fields} for r in flevel_rows]
    write_csv_md("lineage_by_f_level", flevel_fields, flevel_rows)

    # 8. q-resolution bands
    continuity_by_to = {int(r["step_to"]): r for r in continuity_rows}
    qband_rows = []
    for band in CFG["q_resolution_bands"]:
        indices = [i for i, r in enumerate(metrics) if band_id(float(r["min_q3_q4_over_q1"])) == band["id"]]
        if not indices:
            qband_rows.append({"q_resolution_band": band["id"], "point_count": 0})
            continue
        late_ov = np.asarray([late_comp_numeric[i]["overlap"] for i in indices])
        late_ang = np.asarray([late_comp_numeric[i]["maximum_principal_angle_rad"] for i in indices])
        max_trans = np.max(trans_overlap[indices], axis=1)
        cont = [continuity_by_to[i] for i in indices if i in continuity_by_to]
        qband_rows.append({
            "q_resolution_band": band["id"],
            "point_count": len(indices),
            "first_step": min(indices),
            "last_step": max(indices),
            "minimum_q_ratio": f17(min(float(metrics[i]["min_q3_q4_over_q1"]) for i in indices)),
            "maximum_q_ratio": f17(max(float(metrics[i]["min_q3_q4_over_q1"]) for i in indices)),
            "median_D34_vs_late_overlap": f17(np.median(late_ov)),
            "median_D34_vs_late_max_angle_rad": f17(np.median(late_ang)),
            "median_max_D34_vs_Tperp_overlap": f17(np.median(max_trans)),
            "median_consecutive_max_angle_rad": f17(np.median([
                float(r["maximum_principal_angle_rad"]) for r in cont
            ])) if cont else "",
            "maximum_consecutive_max_angle_rad": f17(max([
                float(r["maximum_principal_angle_rad"]) for r in cont
            ])) if cont else "",
            "basis_column_swap_count": sum(r["basis_column_swap"] == "True" for r in cont),
            "sign_flip_total": sum(int(r["sign_flip_count"]) for r in cont),
            "median_direction_3_occupation": f17(np.median([
                float(metrics[i]["direction_3_occupation"]) for i in indices
            ])),
            "median_direction_4_occupation": f17(np.median([
                float(metrics[i]["direction_4_occupation"]) for i in indices
            ])),
        })
    qband_fields = list(max(qband_rows, key=len))
    qband_rows = [{c: r.get(c, "") for c in qband_fields} for r in qband_rows]
    write_csv_md("lineage_by_q_resolution_band", qband_fields, qband_rows)

    # 9. exact transverse reconstruction table
    write_csv_md("transverse_direction_reconstruction", list(trans_rows[0]), trans_rows)

    # fixed descriptive classification
    ratios = np.asarray([float(r["min_q3_q4_over_q1"]) for r in metrics])
    resolved = np.flatnonzero(ratios >= 1e-6)
    if len(resolved) == 0:
        classification = "INSUFFICIENT_RESOLUTION"
        evaluation_start = None
        late_overlap_median = None
        late_angle_median = None
        per_seed_transverse_medians = []
        transverse_aggregate = None
    else:
        evaluation_start = int(resolved[0])
        evaluation_indices = np.arange(evaluation_start, max_step + 1)
        late_overlap_median = float(np.median([
            late_comp_numeric[i]["overlap"] for i in evaluation_indices
        ]))
        late_angle_median = float(np.median([
            late_comp_numeric[i]["maximum_principal_angle_rad"] for i in evaluation_indices
        ]))
        per_seed_transverse_medians = [
            float(np.median(trans_overlap[evaluation_indices, seed]))
            for seed in range(Tbasis.shape[0])
        ]
        transverse_aggregate = max(per_seed_transverse_medians)
        c = CFG["classification"]
        if (
            late_overlap_median >= c["late_overlap_median_min"]
            and late_angle_median <= c["late_max_principal_angle_median_max_rad"]
            and transverse_aggregate <= c["transverse_overlap_median_max_for_continuous"]
        ):
            classification = "FIRST_DIRECTIONS_CONTINUOUS"
        elif (
            transverse_aggregate >= c["transverse_overlap_median_min_for_match"]
            and late_overlap_median <= c["late_overlap_median_max_for_transverse_match"]
        ):
            classification = "MATCHES_LATE_TRANSVERSE_GERM"
        else:
            classification = "ROTATING_OR_MIXED_LINEAGE"

    # 10. health
    health_rows = []
    arrays = {
        "D34_basis_all_steps": D,
        "P34_all_steps": P,
        "S4_basis_all_steps": S4,
        "PS4_all_steps": PS4,
        "Tperp_real_basis": Tbasis,
        "Tperp_projectors": Tproj,
        "D34_late_basis": Dlate,
        "D34_late_projector": Plate,
    }
    health_ok = True
    for name, arr in arrays.items():
        nonfinite = int(np.sum(~np.isfinite(arr)))
        health_ok = health_ok and nonfinite == 0
        health_rows.append({
            "artifact": name,
            "shape": "x".join(str(x) for x in arr.shape),
            "nonfinite_count": nonfinite,
            "minimum": f17(np.min(arr)),
            "maximum": f17(np.max(arr)),
        })
    maxima = {
        "max_D34_orthogonality_error": max(float(r["D34_orthogonality_error"]) for r in metrics),
        "max_S4_orthogonality_error": max(float(r["S4_orthogonality_error"]) for r in metrics),
        "max_P34_idempotence_error": max(float(r["P34_idempotence_error"]) for r in metrics),
        "max_P34_symmetry_error": max(float(r["P34_symmetry_error"]) for r in metrics),
        "max_norm_error": max(float(r["norm_error"]) for r in metrics),
        "max_closure_error": max(float(r["closure_error"]) for r in metrics),
    }
    for name, value in maxima.items():
        health_rows.append({
            "artifact": "trajectory_metrics", "shape": str(max_step + 1),
            "nonfinite_count": 0, "minimum": "", "maximum": f17(value),
            "metric": name,
        })
    health_fields = ["artifact", "shape", "nonfinite_count", "minimum", "maximum", "metric"]
    health_rows = [{c: r.get(c, "") for c in health_fields} for r in health_rows]
    write_csv_md("numerical_health", health_fields, health_rows)

    rapid_start, rapid_end = CFG["display_windows"]["crossing_zoom"]
    rapid_cont = [
        r for r in continuity_rows
        if rapid_start < int(r["step_to"]) <= rapid_end
    ]
    pre_indices = np.arange(0, replay["crossing"])
    pre_resolved = pre_indices[ratios[pre_indices] >= 1e-6]
    summary = {
        "stage": "A2c",
        "status": "ANALYSIS_COMPLETE" if health_ok else "ANALYSIS_INCOMPLETE",
        "classification": classification,
        "classification_evaluation_start_step": evaluation_start,
        "classification_evaluation_end_step": max_step,
        "classification_late_overlap_median": late_overlap_median,
        "classification_late_max_angle_median_rad": late_angle_median,
        "classification_D34_Tperp_overlap_median_by_seed": per_seed_transverse_medians,
        "classification_transverse_aggregate_max_seed_median": transverse_aggregate,
        "D34_late_window": [late_start, late_end],
        "D34_late_member_count": late_end - late_start + 1,
        "D34_late_mean_projector_eigenvalues": [float(x) for x in late_eigenvalues],
        "late_vs_transverse": late_trans_numeric,
        "t0_S4_vs_Tperp_max_orthogonality_error": max(
            float(r["same_t0_S4_vs_Tperp_orthogonality_error"]) for r in late_trans_rows
        ),
        "pre_crossing_resolved_point_count": int(len(pre_resolved)),
        "pre_crossing_resolved_late_overlap_median": (
            float(np.median([late_comp_numeric[i]["overlap"] for i in pre_resolved]))
            if len(pre_resolved) else None
        ),
        "pre_crossing_resolved_max_transverse_overlap_median": (
            float(np.median(np.max(trans_overlap[pre_resolved], axis=1)))
            if len(pre_resolved) else None
        ),
        "crossing_window_max_consecutive_angle_rad": max(
            float(r["maximum_principal_angle_rad"]) for r in rapid_cont
        ),
        "crossing_window_max_projector_distance": max(
            float(r["projector_distance"]) for r in rapid_cont
        ),
        "crossing_window_basis_swap_count": sum(r["basis_column_swap"] == "True" for r in rapid_cont),
        "crossing_window_sign_flip_total": sum(int(r["sign_flip_count"]) for r in rapid_cont),
        "trajectory_f_bitwise_matches_stage_a0": replay["trajectory_f_bitwise_matches_stage_a0"],
        "q_diagnostics_bitwise_match_stage_a0": replay["q_diagnostics_bitwise_match_stage_a0"],
        "numerical_health_passed": health_ok,
        "single_direction_establishment_step_selected": False,
        "H1_H2_H0_judged": False,
    }
    (PROCESSED / "lineage_analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        f"ANALYSIS_COMPLETE: classification={classification}, "
        f"evaluation_start={evaluation_start}, health={health_ok}"
    )


if __name__ == "__main__":
    main()
