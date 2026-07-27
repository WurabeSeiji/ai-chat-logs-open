#!/usr/bin/env python3
"""N=5 float64、Z0=v、明示seedなしの固定5000-step軌道を同条件で2回生成する。"""

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
CONFIG_PATH = HERE / "config_locked.json"
HASHES_PATH = HERE / "expected_hashes.json"
VERIFY_LOG = LOGS / "source_verification.json"
FLOAT_FMT = ".17e"

F_COLUMNS = ["step", "f", "log10_f", "norm_error", "closure_error", "conservation_error"]
Q_COLUMNS = [
    "step", "q1", "q2", "q3", "q4", "q3_over_q1", "q4_over_q1",
    "rank_q", "gram_rank", "dominant_eigenvalue",
]
OCC_COLUMNS = [
    "step", "direction_1_occupation", "direction_2_occupation",
    "direction_3_occupation", "direction_4_occupation",
    "other_rotating_occupation", "kernel_occupation", "occupation_sum",
    "splitting_fraction",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fmt(value: float) -> str:
    return format(float(value), FLOAT_FMT)


def load_fixed_inputs() -> tuple[dict, dict]:
    cfg = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    expected = json.loads(HASHES_PATH.read_text(encoding="utf-8"))
    if not VERIFY_LOG.is_file():
        raise SystemExit("EXECUTION_FAILED: verify_sources.pyの成功記録がない")
    verified = json.loads(VERIFY_LOG.read_text(encoding="utf-8"))
    if verified.get("status") != "VERIFIED":
        raise SystemExit("SOURCE_MISMATCH: source verification gate")
    for group in ("sources", "dependencies", "stage_a0_inputs", "stage_a1b_inputs"):
        for item in expected[group].values():
            path = REPO / item["path"]
            if not path.is_file() or sha256(path) != item["sha256"]:
                raise SystemExit(f"SOURCE_MISMATCH: {path}")
    if cfg["n"] != 5 or cfg["dtype"] != "float64" or cfg["max_step"] != 5000:
        raise SystemExit("EXECUTION_FAILED: locked N/dtype/max_step mismatch")
    return cfg, expected


def import_originals(expected: dict):
    src = expected["sources"]
    engine_dir = (REPO / src["run_n_scaling_lowrank_v1.py"]["path"]).parent
    v2_code = (REPO / src["run_n300_dimension_saturation_v2.py"]["path"]).parent
    p7_code = (REPO / src["run_paper7_5color_timeseries.py"]["path"]).parent
    for path in (str(engine_dir), str(v2_code), str(p7_code)):
        if path not in sys.path:
            sys.path.insert(0, path)
    lowrank = importlib.import_module("run_n_scaling_lowrank_v1")
    exact = importlib.import_module("run_plane_flow_exact_v1")
    saturation = importlib.import_module("run_n300_dimension_saturation_v2")
    color = importlib.import_module("run_paper7_5color_timeseries")
    required = {
        "LowRankSystem": lowrank.LowRankSystem,
        "make_parent": lowrank.make_parent,
        "parent_plane_split_exact": exact.parent_plane_split_exact,
        "gram_reduce": saturation.gram_reduce,
        "dominant_plane": saturation.dominant_plane,
        "qsv4": saturation.qsv4,
        "occ": color.occ,
        "s4_new_dirs": color.s4_new_dirs,
        "align_2d": color.align_2d,
    }
    return required


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def serialize_row(columns: list[str], values: dict, integer_columns: set[str]) -> dict:
    return {
        key: (str(int(values[key])) if key in integer_columns else fmt(values[key]))
        for key in columns
    }


def run_one(run_id: str, cfg: dict, fn: dict) -> dict:
    run_dir = RAW / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        raise SystemExit(f"EXECUTION_FAILED: 既存rawを上書きする可能性: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)

    n = cfg["n"]
    rng = np.random.default_rng(cfg["parent_prng_seed"])
    sys_lr = fn["LowRankSystem"](n)
    v, parent_residual, parent_sigma = fn["make_parent"](
        sys_lr,
        rng,
        iters=cfg["parent_iters"],
        tol=cfg["parent_tolerance"],
    )
    _, B_p1, B_rot, parent_rotation_spectrum = fn["parent_plane_split_exact"](sys_lr, v)
    gr0 = fn["gram_reduce"](sys_lr, v)
    _, B0, _, _, _ = fn["dominant_plane"](sys_lr, gr0)

    # 固定仕様: seed生成関数を呼ばず、親をそのまま初期状態とする。
    Z = v.copy()
    if not np.array_equal(Z, v):
        raise RuntimeError("Z0 is not bitwise equal to v")
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    # 状態には加えない。sigma_max_powerの既存warm-start専用。
    wp = rng.normal(size=sys_lr.m)

    levels = [float(x) for x in cfg["f_levels"]]
    pending = set(range(len(levels)))
    passages: list[dict] = []
    f_rows: list[dict] = []
    q_rows: list[dict] = []
    occ_rows: list[dict] = []
    plane_steps: list[int] = []
    plane_values: list[np.ndarray] = []
    prev_f34 = None
    max_projection_closure = 0.0
    initial_norm = float(np.real(np.vdot(Z, Z)))

    def f_value(state: np.ndarray) -> float:
        transverse = state - p * (p @ state) - q * (q @ state)
        return float(np.real(np.vdot(transverse, transverse))) / float(np.real(np.vdot(state, state)))

    def geometry(state: np.ndarray):
        gr = fn["gram_reduce"](sys_lr, state)
        dominant_eigenvalue, Bdom, _, _, _ = fn["dominant_plane"](sys_lr, gr)
        q4, _ = fn["qsv4"](B0, Bdom)
        if q4.shape[0] != 4:
            raise RuntimeError("qsv4 did not return q1..q4")
        rank_q = int(np.sum(q4 > cfg["rank_q_relative_threshold"] * q4[0]))
        return gr, float(dominant_eigenvalue), Bdom, q4, rank_q

    started = time.perf_counter()
    for step in range(cfg["max_step"] + 1):
        norm_sq = float(np.real(np.vdot(Z, Z)))
        f = f_value(Z)
        log10_f = math.log10(f) if f > 0.0 else float("nan")
        norm_error = abs(norm_sq - 1.0)
        closure_error = abs(complex(Z @ Z))
        conservation_error = abs(norm_sq - initial_norm)
        base_metrics = np.array([f, norm_sq, norm_error, closure_error, conservation_error], dtype=float)
        if not np.all(np.isfinite(base_metrics)):
            raise FloatingPointError(f"NaN/Inf at step {step}")
        f_rows.append(serialize_row(F_COLUMNS, {
            "step": step,
            "f": f,
            "log10_f": log10_f,
            "norm_error": norm_error,
            "closure_error": closure_error,
            "conservation_error": conservation_error,
        }, {"step"}))

        newly_crossed = [idx for idx in sorted(pending) if f >= levels[idx]]
        for idx in newly_crossed:
            pending.remove(idx)
        need_q = (step % cfg["q_every"] == 0) or bool(newly_crossed)
        need_occ = (step % cfg["occupation_every"] == 0) or bool(newly_crossed)
        q_values = None
        occ_values = None
        Bdom = None

        if need_q or need_occ:
            gr, dominant_eigenvalue, Bdom, q4, rank_q = geometry(Z)
            q_values = {
                "step": step,
                "q1": q4[0], "q2": q4[1], "q3": q4[2], "q4": q4[3],
                "q3_over_q1": q4[2] / q4[0] if q4[0] != 0 else float("nan"),
                "q4_over_q1": q4[3] / q4[0] if q4[0] != 0 else float("nan"),
                "rank_q": rank_q,
                "gram_rank": gr["diag"]["r_G"],
                "dominant_eigenvalue": dominant_eigenvalue,
            }
            if not np.all(np.isfinite([
                q_values[k] for k in ("q1", "q2", "q3", "q4", "q3_over_q1",
                                      "q4_over_q1", "dominant_eigenvalue")
            ])):
                raise FloatingPointError(f"nonfinite q measurement at step {step}")

        if need_q:
            q_rows.append(serialize_row(Q_COLUMNS, q_values, {"step", "rank_q", "gram_rank"}))
            plane_steps.append(step)
            plane_values.append(Bdom.copy())

        if need_occ:
            total = float(np.real(np.vdot(Z, Z)))
            E_p1 = fn["occ"](B_p1, Z)
            E_other = fn["occ"](B_rot, Z)
            E_kernel = total - E_p1 - E_other
            e34 = fn["s4_new_dirs"](B0, Bdom)
            projection = B_rot @ (B_rot.T @ e34)
            projected_basis, _ = np.linalg.qr(projection)
            f34 = projected_basis[:, :2]
            f34 = fn["align_2d"](prev_f34, f34)
            prev_f34 = f34
            E_d3 = fn["occ"](f34[:, [0]], Z)
            E_d4 = fn["occ"](f34[:, [1]], Z)
            E_remaining = max(0.0, E_other - E_d3 - E_d4)
            E_a1 = fn["occ"](B_p1[:, [0]], Z)
            E_a2 = fn["occ"](B_p1[:, [1]], Z)
            occupation_sum = (E_p1 + E_d3 + E_d4 + E_remaining + E_kernel) / total
            max_projection_closure = max(max_projection_closure, abs(occupation_sum - 1.0))
            occ_values = {
                "step": step,
                "direction_1_occupation": E_a1 / total,
                "direction_2_occupation": E_a2 / total,
                "direction_3_occupation": E_d3 / total,
                "direction_4_occupation": E_d4 / total,
                "other_rotating_occupation": E_remaining / total,
                "kernel_occupation": E_kernel / total,
                "occupation_sum": occupation_sum,
                "splitting_fraction": 1.0 - E_p1 / total,
            }
            if not np.all(np.isfinite([occ_values[k] for k in OCC_COLUMNS if k != "step"])):
                raise FloatingPointError(f"nonfinite occupation measurement at step {step}")
            occ_rows.append(serialize_row(OCC_COLUMNS, occ_values, {"step"}))

        if newly_crossed:
            if q_values is None or occ_values is None:
                raise RuntimeError("first-passage measurement missing")
            for idx in newly_crossed:
                passages.append({
                    "level_index": idx,
                    "level": fmt(levels[idx]),
                    "level_label": f"{levels[idx]:.0e}" if levels[idx] < 0.05 else str(levels[idx]),
                    "status": "found",
                    "first_passage_step": step,
                    "f_at_first_passage": fmt(f),
                    **{k: (str(int(q_values[k])) if k in {"rank_q", "gram_rank"} else fmt(q_values[k]))
                       for k in Q_COLUMNS if k != "step"},
                    **{k: fmt(occ_values[k]) for k in OCC_COLUMNS if k != "step"},
                })

        if step == cfg["max_step"]:
            break
        try:
            sys_lr.set_theta(np.angle(Z))
            sigma_estimate, wp = sys_lr.sigma_max_power(wp)
            if not np.isfinite(sigma_estimate) or sigma_estimate <= 0.0:
                raise FloatingPointError(f"invalid sigma at step {step}: {sigma_estimate}")
            Z = sys_lr.cayley_step(Z, sigma_estimate)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(f"linear algebra failure after step {step}") from exc
        if not np.all(np.isfinite(Z.real)) or not np.all(np.isfinite(Z.imag)):
            raise FloatingPointError(f"NaN/Inf state after step {step}")

    for idx in sorted(pending):
        passages.append({
            "level_index": idx,
            "level": fmt(levels[idx]),
            "level_label": f"{levels[idx]:.0e}" if levels[idx] < 0.05 else str(levels[idx]),
            "status": "not_found",
            "first_passage_step": "",
            "f_at_first_passage": "",
            **{k: "" for k in Q_COLUMNS if k != "step"},
            **{k: "" for k in OCC_COLUMNS if k != "step"},
        })
    passages.sort(key=lambda row: int(row["level_index"]))

    write_csv(run_dir / "f_timeseries.csv", F_COLUMNS, f_rows)
    write_csv(run_dir / "q_timeseries.csv", Q_COLUMNS, q_rows)
    write_csv(run_dir / "occupation_timeseries.csv", OCC_COLUMNS, occ_rows)
    passage_columns = [
        "level_index", "level", "level_label", "status", "first_passage_step",
        "f_at_first_passage",
        *[x for x in Q_COLUMNS if x != "step"],
        *[x for x in OCC_COLUMNS if x != "step"],
    ]
    write_csv(run_dir / "first_passage_measurements.csv", passage_columns, passages)
    np.save(run_dir / "dominant_plane_steps.npy", np.asarray(plane_steps, dtype=np.int64))
    np.save(run_dir / "dominant_plane_values.npy", np.asarray(plane_values, dtype=np.float64))

    elapsed = time.perf_counter() - started
    summary = {
        "run_id": run_id,
        "stage": "A2a",
        "status": "COMPLETED",
        "n": n,
        "m": sys_lr.m,
        "dtype": str(Z.real.dtype),
        "parent_prng_seed": cfg["parent_prng_seed"],
        "parent_residual": float(parent_residual),
        "parent_sigma_spectrum": [float(x) for x in parent_sigma],
        "parent_rotation_spectrum": [float(x) for x in parent_rotation_spectrum],
        "initial_state_rule": "Z0 = v.copy()",
        "initial_state_bitwise_equal_parent": True,
        "explicit_state_seed_added": False,
        "zero_closure_kernel_seed_called": False,
        "warmstart_rule": "rng.normal after make_parent; not added to state",
        "max_step": cfg["max_step"],
        "final_step": cfg["max_step"],
        "f_row_count": len(f_rows),
        "q_row_count": len(q_rows),
        "occupation_row_count": len(occ_rows),
        "first_passage_found_count": sum(row["status"] == "found" for row in passages),
        "first_passage_total_count": len(passages),
        "initial_f": float(f_rows[0]["f"]),
        "final_f": float(f_rows[-1]["f"]),
        "maximum_norm_error": max(float(row["norm_error"]) for row in f_rows),
        "maximum_closure_error": max(float(row["closure_error"]) for row in f_rows),
        "maximum_conservation_error": max(float(row["conservation_error"]) for row in f_rows),
        "maximum_projection_closure_error": float(max_projection_closure),
    }
    (run_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {"run_id": run_id, "elapsed_seconds": elapsed, "summary": summary}


def main() -> None:
    cfg, expected = load_fixed_inputs()
    RAW.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    if cfg["output_policy"]["fail_if_raw_not_empty"] and any(RAW.iterdir()):
        raise SystemExit("EXECUTION_FAILED: raw/が空ではないため上書きを拒否")
    fn = import_originals(expected)
    runs = []
    try:
        for run_id in cfg["run_ids"]:
            print(f"START {run_id}", flush=True)
            result = run_one(run_id, cfg, fn)
            runs.append(result)
            print(f"COMPLETED {run_id}: {result['elapsed_seconds']:.3f} s", flush=True)
    except Exception as exc:
        failure = {
            "stage": "A2a",
            "status": "EXECUTION_FAILED",
            "completed_runs": [r["run_id"] for r in runs],
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
        (LOGS / "execution_manifest.json").write_text(
            json.dumps(failure, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        raise
    manifest = {
        "stage": "A2a",
        "status": "COMPLETED",
        "execution_order": cfg["run_ids"],
        "runs": [{"run_id": r["run_id"], "elapsed_seconds": r["elapsed_seconds"]} for r in runs],
        "total_elapsed_seconds": sum(r["elapsed_seconds"] for r in runs),
    }
    (LOGS / "execution_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("COMPLETED: two N=5 seedless float64 executions")


if __name__ == "__main__":
    main()
