#!/usr/bin/env python3
"""Stage A0 N=5軌道を再実行し、既存方向基底と射影行列を全step保存する。"""

from __future__ import annotations

import csv
import hashlib
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = HERE / "raw"
LOGS = HERE / "logs"
CFG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
EXPECTED = json.loads((HERE / "expected_hashes.json").read_text(encoding="utf-8"))
FMT = ".17e"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def f17(x: float) -> str:
    return format(float(x), FMT)


def bitwise_equal_float(a: float, b: float) -> bool:
    return np.float64(a).view(np.uint64) == np.float64(b).view(np.uint64)


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def gate_and_import():
    gate = LOGS / "source_verification.json"
    if not gate.is_file() or json.loads(gate.read_text(encoding="utf-8")).get("status") != "VERIFIED":
        raise SystemExit("SOURCE_MISMATCH: verify_sources.py gate")
    for group in ("sources", "dependencies", "stage_a0"):
        for item in EXPECTED[group].values():
            path = REPO / item["path"]
            if not path.is_file() or sha256(path) != item["sha256"]:
                raise SystemExit(f"SOURCE_MISMATCH: {path}")
    src = EXPECTED["sources"]
    engine = (REPO / src["run_n_scaling_lowrank_v1.py"]["path"]).parent
    code = (REPO / src["run_n300_dimension_saturation_v2.py"]["path"]).parent
    p7code = (REPO / src["run_paper7_5color_timeseries.py"]["path"]).parent
    for p in (engine, code, p7code):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    low = importlib.import_module("run_n_scaling_lowrank_v1")
    exact = importlib.import_module("run_plane_flow_exact_v1")
    sat = importlib.import_module("run_n300_dimension_saturation_v2")
    color = importlib.import_module("run_paper7_5color_timeseries")
    transverse = importlib.import_module("run_paper7_transverse")
    return low, exact, sat, color, transverse


def main() -> None:
    if RAW.exists() and any(RAW.iterdir()):
        raise SystemExit("EXECUTION_FAILED: raw/が空ではないため上書きを拒否")
    RAW.mkdir(parents=True, exist_ok=True)
    low, exact, sat, color, transverse = gate_and_import()

    n = CFG["n"]
    max_step = CFG["trajectory"]["max_step"]
    M = n * (n - 1) // 2
    rng = np.random.default_rng(CFG["trajectory"]["parent_prng_seed"])
    sys_lr = low.LowRankSystem(n)
    v, parent_residual, parent_sigma = low.make_parent(
        sys_lr, rng, iters=CFG["trajectory"]["parent_iters"],
        tol=CFG["trajectory"]["parent_tolerance"]
    )
    _, B_p1, B_rot, _ = exact.parent_plane_split_exact(sys_lr, v)
    gr0 = sat.gram_reduce(sys_lr, v)
    _, B0, _, _, _ = sat.dominant_plane(sys_lr, gr0)
    g = low.zero_closure_kernel_seed(sys_lr, rng)
    Z = v + CFG["trajectory"]["delta"] * g
    Z = Z / np.linalg.norm(Z)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=M)

    f_reference = load_csv(REPO / EXPECTED["stage_a0"]["fcurve"]["path"])
    f_ref = {int(r["tau"]): float(r["f"]) for r in f_reference if int(r["tau"]) <= max_step}
    q_reference = load_csv(REPO / EXPECTED["stage_a0"]["q_svd"]["path"])
    q_ref = {int(r["step"]): r for r in q_reference if int(r["step"]) <= max_step}
    if sorted(f_ref) != list(range(max_step + 1)):
        raise SystemExit("EXECUTION_FAILED: Stage A0 f時間軸が0..5000を覆わない")

    shape2 = (max_step + 1, M, 2)
    shape4 = (max_step + 1, M, 4)
    shapeP = (max_step + 1, M, M)
    states = np.empty((max_step + 1, M), dtype=np.complex128)
    warmstarts = np.empty((max_step + 1, M), dtype=np.float64)
    Bdom_all = np.empty(shape2, dtype=np.float64)
    D34_all = np.empty(shape2, dtype=np.float64)
    display_all = np.empty(shape2, dtype=np.float64)
    P34_all = np.empty(shapeP, dtype=np.float64)
    S4_all = np.empty(shape4, dtype=np.float64)
    PS4_all = np.empty(shapeP, dtype=np.float64)
    metric_rows: list[dict] = []
    passages: list[dict] = []
    pending = set(range(len(CFG["f_levels"])))
    previous_display = None
    crossing = None
    q_comparison_count = 0
    q_bitwise_all = True
    maximum_q_abs_error = 0.0
    maximum_f_abs_error = 0.0
    started = time.perf_counter()

    def f_value(state: np.ndarray) -> float:
        zp = state - p * (p @ state) - q * (q @ state)
        return float(np.real(np.vdot(zp, zp))) / float(np.real(np.vdot(state, state)))

    for step in range(max_step + 1):
        states[step] = Z
        warmstarts[step] = wp
        f = f_value(Z)
        ref_f = f_ref[step]
        maximum_f_abs_error = max(maximum_f_abs_error, abs(f - ref_f))
        if not bitwise_equal_float(f, ref_f):
            raise SystemExit(f"TRAJECTORY_MISMATCH: f bitwise mismatch at step {step}")
        if crossing is None and f > 0.05:
            crossing = step

        gr = sat.gram_reduce(sys_lr, Z)
        dominant_eigenvalue, Bdom, _, lifted_residual, Bdom_orth = sat.dominant_plane(sys_lr, gr)
        q4, _ = sat.qsv4(B0, Bdom)
        rank_q = int(np.sum(q4 > 1e-8 * q4[0]))
        D34 = color.s4_new_dirs(B0, Bdom)
        P34 = D34 @ D34.T
        S4 = transverse.s4_basis(sys_lr, B0, Z)
        PS4 = S4 @ S4.T
        direct_s4, _ = np.linalg.qr(np.column_stack([B0, Bdom]))
        direct_PS4 = direct_s4[:, :4] @ direct_s4[:, :4].T
        s4_equivalence_error = float(np.linalg.norm(PS4 - direct_PS4))

        projected = B_rot @ (B_rot.T @ D34)
        display_basis, _ = np.linalg.qr(projected)
        display_basis = color.align_2d(previous_display, display_basis[:, :2])
        previous_display = display_basis

        total = float(np.real(np.vdot(Z, Z)))
        E_p1 = color.occ(B_p1, Z)
        E_other = color.occ(B_rot, Z)
        E_kernel = total - E_p1 - E_other
        E_d3 = color.occ(display_basis[:, [0]], Z)
        E_d4 = color.occ(display_basis[:, [1]], Z)
        direction_3 = E_d3 / total
        direction_4 = E_d4 / total
        positive_mu = np.sort(gr["mu"][gr["mu"] > 0.0])[::-1]
        if len(positive_mu) < 2:
            raise SystemExit(f"EXECUTION_FAILED: dominant eigenvalue gap undefined at step {step}")
        eigen_gap = float(positive_mu[0] - positive_mu[1])
        eigen_gap_relative = eigen_gap / float(positive_mu[0])
        q3_over_q1 = float(q4[2] / q4[0])
        q4_over_q1 = float(q4[3] / q4[0])
        qmin = min(q3_over_q1, q4_over_q1)
        D_orth = float(np.linalg.norm(D34.T @ D34 - np.eye(2)))
        S4_orth = float(np.linalg.norm(S4.T @ S4 - np.eye(4)))
        P_idem = float(np.linalg.norm(P34 @ P34 - P34))
        P_sym = float(np.linalg.norm(P34 - P34.T))

        if step in q_ref:
            rr = q_ref[step]
            q_comparison_count += 1
            for key, value in (
                ("dominant_eigenvalue", dominant_eigenvalue),
                ("q1", q4[0]), ("q2", q4[1]), ("q3", q4[2]), ("q4", q4[3]),
            ):
                reference_value = float(rr[key])
                maximum_q_abs_error = max(maximum_q_abs_error, abs(float(value) - reference_value))
                q_bitwise_all = q_bitwise_all and bitwise_equal_float(float(value), reference_value)
            q_bitwise_all = q_bitwise_all and int(rr["rank_q"]) == rank_q
            q_bitwise_all = q_bitwise_all and int(rr["gram_rank"]) == int(gr["diag"]["r_G"])

        newly_crossed = [i for i in sorted(pending) if f >= float(CFG["f_levels"][i])]
        for i in newly_crossed:
            pending.remove(i)
            passages.append({
                "level_index": i,
                "level": f17(CFG["f_levels"][i]),
                "level_label": f"{float(CFG['f_levels'][i]):.0e}" if float(CFG["f_levels"][i]) < 0.05 else str(CFG["f_levels"][i]),
                "status": "found",
                "first_passage_step": step,
                "f_at_first_passage": f17(f),
            })

        Bdom_all[step] = Bdom
        D34_all[step] = D34
        display_all[step] = display_basis
        P34_all[step] = P34
        S4_all[step] = S4
        PS4_all[step] = PS4
        metric_rows.append({
            "step": step,
            "f": f17(f),
            "log10_f": f17(math.log10(f)) if f > 0 else "nan",
            "q1": f17(q4[0]), "q2": f17(q4[1]), "q3": f17(q4[2]), "q4": f17(q4[3]),
            "q3_over_q1": f17(q3_over_q1),
            "q4_over_q1": f17(q4_over_q1),
            "min_q3_q4_over_q1": f17(qmin),
            "rank_q": rank_q,
            "gram_rank": int(gr["diag"]["r_G"]),
            "dominant_eigenvalue": f17(dominant_eigenvalue),
            "dominant_eigenvalue_gap": f17(eigen_gap),
            "dominant_eigenvalue_gap_relative": f17(eigen_gap_relative),
            "direction_3_occupation": f17(direction_3),
            "direction_4_occupation": f17(direction_4),
            "kernel_occupation": f17(E_kernel / total),
            "Bdom_orthogonality_error": f17(Bdom_orth),
            "D34_orthogonality_error": f17(D_orth),
            "S4_orthogonality_error": f17(S4_orth),
            "P34_idempotence_error": f17(P_idem),
            "P34_symmetry_error": f17(P_sym),
            "S4_direct_equivalence_error": f17(s4_equivalence_error),
            "dominant_lifted_residual": f17(lifted_residual),
            "norm_error": f17(abs(total - 1.0)),
            "closure_error": f17(abs(complex(Z @ Z))),
        })

        if step == max_step:
            break
        sys_lr.set_theta(np.angle(Z))
        sigma_estimate, wp = sys_lr.sigma_max_power(wp)
        if not np.isfinite(sigma_estimate) or sigma_estimate <= 0:
            raise SystemExit(f"EXECUTION_FAILED: invalid sigma at step {step}")
        try:
            Z = sys_lr.cayley_step(Z, sigma_estimate)
        except np.linalg.LinAlgError as exc:
            raise SystemExit(f"EXECUTION_FAILED: linear algebra at step {step}: {exc}")
        if not np.all(np.isfinite(Z.real)) or not np.all(np.isfinite(Z.imag)):
            raise SystemExit(f"EXECUTION_FAILED: nonfinite state after step {step}")

    if crossing != 1167:
        raise SystemExit(f"TRAJECTORY_MISMATCH: crossing={crossing}, expected=1167")
    if pending:
        for i in sorted(pending):
            passages.append({
                "level_index": i,
                "level": f17(CFG["f_levels"][i]),
                "level_label": f"{float(CFG['f_levels'][i]):.0e}" if float(CFG["f_levels"][i]) < 0.05 else str(CFG["f_levels"][i]),
                "status": "not_found",
                "first_passage_step": "",
                "f_at_first_passage": "",
            })
    passages.sort(key=lambda r: int(r["level_index"]))
    if not q_bitwise_all:
        raise SystemExit(f"TRAJECTORY_MISMATCH: Stage A0 q diagnostics differ; max_abs={maximum_q_abs_error}")

    np.save(RAW / "trajectory_states.npy", states)
    np.save(RAW / "trajectory_warmstarts.npy", warmstarts)
    np.save(RAW / "B0.npy", B0)
    np.save(RAW / "Bdom_all_steps.npy", Bdom_all)
    np.save(RAW / "D34_basis_all_steps.npy", D34_all)
    np.save(RAW / "D34_display_aligned_all_steps.npy", display_all)
    np.save(RAW / "P34_all_steps.npy", P34_all)
    np.save(RAW / "S4_basis_all_steps.npy", S4_all)
    np.save(RAW / "PS4_all_steps.npy", PS4_all)
    metric_fields = list(metric_rows[0])
    write_csv(RAW / "trajectory_direction_metrics.csv", metric_fields, metric_rows)
    write_csv(RAW / "f_first_passage_levels.csv", list(passages[0]), passages)

    fixed = set(int(x) for x in CFG["fixed_steps"])
    passage_map: dict[int, list[str]] = {}
    for row in passages:
        if row["status"] == "found":
            passage_map.setdefault(int(row["first_passage_step"]), []).append(row["level_label"])
    late_start = CFG["late_reference"]["start_step"]
    late_end = CFG["late_reference"]["end_step"]
    selected = sorted(fixed | set(passage_map) | set(range(late_start, late_end + 1)))
    snapshots = []
    for step in selected:
        r = metric_rows[step]
        snapshots.append({
            "snapshot_index": len(snapshots),
            "step": step,
            "is_fixed_step": str(step in fixed),
            "f_level_labels": ";".join(passage_map.get(step, [])),
            "is_late_reference_member": str(late_start <= step <= late_end),
            "basis_array_index": step,
            "B0_file": "raw/B0.npy",
            "Bdom_file": "raw/Bdom_all_steps.npy",
            "D34_basis_file": "raw/D34_basis_all_steps.npy",
            "P34_file": "raw/P34_all_steps.npy",
            "S4_basis_file": "raw/S4_basis_all_steps.npy",
            "PS4_file": "raw/PS4_all_steps.npy",
            "f": r["f"],
            "q3_over_q1": r["q3_over_q1"],
            "q4_over_q1": r["q4_over_q1"],
            "direction_3_occupation": r["direction_3_occupation"],
            "direction_4_occupation": r["direction_4_occupation"],
        })
    write_csv(RAW / "direction_basis_snapshots.csv", list(snapshots[0]), snapshots)

    summary = {
        "stage": "A2c",
        "status": "REPLAY_COMPLETE",
        "n": n,
        "m": M,
        "max_step": max_step,
        "crossing": crossing,
        "parent_residual": float(parent_residual),
        "parent_sigma": [float(x) for x in parent_sigma],
        "trajectory_f_bitwise_matches_stage_a0": True,
        "trajectory_f_comparison_count": len(f_ref),
        "trajectory_f_max_absolute_error": maximum_f_abs_error,
        "q_diagnostics_bitwise_match_stage_a0": q_bitwise_all,
        "q_comparison_count": q_comparison_count,
        "q_max_absolute_error": maximum_q_abs_error,
        "first_passage_found_count": sum(r["status"] == "found" for r in passages),
        "snapshot_count": len(snapshots),
        "all_step_basis_count": max_step + 1,
        "elapsed_seconds": time.perf_counter() - started,
    }
    (LOGS / "replay_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        f"REPLAY_COMPLETE: f bitwise={summary['trajectory_f_bitwise_matches_stage_a0']}, "
        f"q bitwise={q_bitwise_all}, crossing={crossing}, bases={max_step + 1}"
    )


if __name__ == "__main__":
    main()
