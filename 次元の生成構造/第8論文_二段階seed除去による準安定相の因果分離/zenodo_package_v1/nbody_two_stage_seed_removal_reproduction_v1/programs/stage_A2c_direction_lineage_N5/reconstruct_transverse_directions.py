#!/usr/bin/env python3
"""第7論文横摂動コードと同じseed/epsilon/t0規則でTperpを正確に再構成する。"""

from __future__ import annotations

import csv
import hashlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = HERE / "raw"
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


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def f17(x: float) -> str:
    return format(float(x), ".17e")


def import_originals():
    src = EXPECTED["sources"]
    for item in list(src.values()) + list(EXPECTED["dependencies"].values()) + list(EXPECTED["stage_a0"].values()):
        path = REPO / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise SystemExit(f"SOURCE_MISMATCH: {path}")
    engine = (REPO / src["run_n_scaling_lowrank_v1.py"]["path"]).parent
    code = (REPO / src["run_n300_dimension_saturation_v2.py"]["path"]).parent
    p7code = (REPO / src["run_paper7_transverse.py"]["path"]).parent
    for p in (engine, code, p7code):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    low = importlib.import_module("run_n_scaling_lowrank_v1")
    transverse = importlib.import_module("run_paper7_transverse")
    return low, transverse


def main() -> None:
    replay_gate = LOGS / "replay_summary.json"
    if not replay_gate.is_file():
        raise SystemExit("EXECUTION_FAILED: replay_and_extract_bases.py gate missing")
    replay = json.loads(replay_gate.read_text(encoding="utf-8"))
    if replay.get("status") != "REPLAY_COMPLETE" or not replay.get("trajectory_f_bitwise_matches_stage_a0"):
        raise SystemExit("TRAJECTORY_MISMATCH: replay gate")
    output_files = [
        RAW / "Tperp_eta_complex.npy",
        RAW / "Tperp_real_basis.npy",
        RAW / "Tperp_projectors.npy",
        RAW / "transverse_direction_reconstruction.csv",
    ]
    if any(p.exists() for p in output_files):
        raise SystemExit("EXECUTION_FAILED: transverse出力を上書きする可能性")

    low, transverse = import_originals()
    n = CFG["n"]
    M = n * (n - 1) // 2
    t0 = CFG["transverse"]["t0"]
    record_step = CFG["transverse"]["verification_record_step"]
    states = np.load(RAW / "trajectory_states.npy")
    warmstarts = np.load(RAW / "trajectory_warmstarts.npy")
    B0 = np.load(RAW / "B0.npy")
    saved_PS4 = np.load(RAW / "PS4_all_steps.npy")
    Z0 = states[t0].copy()
    wp0 = warmstarts[t0].copy()
    sys_lr = low.LowRankSystem(n)
    S4_t0 = transverse.s4_basis(sys_lr, B0, Z0)
    PS4_t0 = S4_t0 @ S4_t0.T
    s4_saved_error = float(np.linalg.norm(PS4_t0 - saved_PS4[t0]))
    if s4_saved_error > 1e-12:
        raise SystemExit(f"EXECUTION_FAILED: s4_basis reconstruction mismatch {s4_saved_error}")

    rng = np.random.default_rng(CFG["transverse"]["direction_prng_seed"])
    eta_all = np.empty((CFG["transverse"]["seeds"], M), dtype=np.complex128)
    basis_all = np.empty((CFG["transverse"]["seeds"], M, 2), dtype=np.float64)
    projector_all = np.empty((CFG["transverse"]["seeds"], M, M), dtype=np.float64)
    expected_rows = load_csv(REPO / EXPECTED["stage_a0"]["transverse_stability_timeseries"]["path"])
    expected_4200 = {
        (int(r["seed"]), float(r["epsilon"])): r
        for r in expected_rows if int(r["step"]) == record_step
    }
    if len(expected_4200) != CFG["transverse"]["seeds"] * len(CFG["transverse"]["epsilons"]):
        raise SystemExit("EXECUTION_FAILED: Stage A0 transverse verification rows missing")

    rows = []
    all_record_strings_match = True
    all_baseline_states_match = True
    for seed in range(CFG["transverse"]["seeds"]):
        eta_r = rng.normal(size=M)
        eta_i = rng.normal(size=M)
        eta_r = eta_r - S4_t0 @ (S4_t0.T @ eta_r)
        eta_i = eta_i - S4_t0 @ (S4_t0.T @ eta_i)
        eta_norm_before = float(np.sqrt(eta_r @ eta_r + eta_i @ eta_i))
        if eta_norm_before == 0.0:
            raise SystemExit("EXECUTION_FAILED: existing transverse direction has zero norm")
        eta = (eta_r + 1j * eta_i) / eta_norm_before
        real_columns = np.column_stack([eta.real, eta.imag])
        if np.linalg.matrix_rank(real_columns) != 2:
            raise SystemExit("EXECUTION_FAILED: Tperp real span is not 2D")
        Tbasis, _ = np.linalg.qr(real_columns)
        Tbasis = Tbasis[:, :2]
        Tprojector = Tbasis @ Tbasis.T
        eta_all[seed] = eta
        basis_all[seed] = Tbasis
        projector_all[seed] = Tprojector
        eta_s4_residual = float(np.sqrt(
            np.linalg.norm(S4_t0.T @ eta.real) ** 2 +
            np.linalg.norm(S4_t0.T @ eta.imag) ** 2
        ))
        basis_s4_residual = float(np.linalg.norm(S4_t0.T @ Tbasis))
        basis_orth_error = float(np.linalg.norm(Tbasis.T @ Tbasis - np.eye(2)))
        projector_idem_error = float(np.linalg.norm(Tprojector @ Tprojector - Tprojector))

        for epsilon in CFG["transverse"]["epsilons"]:
            Zb = Z0.copy()
            wpb = wp0.copy()
            Zt = Z0 + float(epsilon) * eta
            Zt = Zt / np.linalg.norm(Zt)
            wpt = wp0.copy()
            for _ in range(record_step - t0):
                Zb, wpb = transverse.evolve(sys_lr, Zb, wpb)
                Zt, wpt = transverse.evolve(sys_lr, Zt, wpt)
            baseline_match = np.array_equal(Zb, states[record_step])
            all_baseline_states_match = all_baseline_states_match and baseline_match
            S4_record = transverse.s4_basis(sys_lr, B0, Zb)
            diff = Zt - Zb
            total_difference = float(np.linalg.norm(diff))
            transverse_difference = float(np.linalg.norm(transverse.perp(S4_record, diff)))
            amplification = transverse_difference / float(epsilon)
            computed = {
                "baseline_norm": f"{np.linalg.norm(Zb):.10e}",
                "perturbed_norm": f"{np.linalg.norm(Zt):.10e}",
                "total_difference": f"{total_difference:.10e}",
                "transverse_difference": f"{transverse_difference:.10e}",
                "normalized_transverse_amplification": f"{amplification:.10e}",
                "norm_error": f"{abs(np.linalg.norm(Zt) - 1):.10e}",
                "conservation_error": f"{abs(np.linalg.norm(Zb) - 1):.10e}",
            }
            expected = expected_4200[(seed, float(epsilon))]
            matches = all(computed[k] == expected[k] for k in computed)
            all_record_strings_match = all_record_strings_match and matches
            rows.append({
                "seed": seed,
                "epsilon": f17(epsilon),
                "t0": t0,
                "verification_step": record_step,
                "direction_prng_seed": CFG["transverse"]["direction_prng_seed"],
                "eta_complex_norm": f17(np.linalg.norm(eta)),
                "eta_S4_perp_residual": f17(eta_s4_residual),
                "Tperp_basis_rank": 2,
                "Tperp_basis_orthogonality_error": f17(basis_orth_error),
                "Tperp_projector_idempotence_error": f17(projector_idem_error),
                "Tperp_vs_S4_projection_error": f17(basis_s4_residual),
                "baseline_state_bitwise_matches_replay_at_verification_step": str(baseline_match),
                "stage_a0_verification_record_string_match": str(matches),
                "recomputed_total_difference": computed["total_difference"],
                "stage_a0_total_difference": expected["total_difference"],
                "recomputed_transverse_difference": computed["transverse_difference"],
                "stage_a0_transverse_difference": expected["transverse_difference"],
                "recomputed_normalized_amplification": computed["normalized_transverse_amplification"],
                "stage_a0_normalized_amplification": expected["normalized_transverse_amplification"],
            })

    if not all_record_strings_match or not all_baseline_states_match:
        raise SystemExit(
            "TRANSVERSE_RECONSTRUCTION_MISMATCH: "
            f"records={all_record_strings_match}, baseline={all_baseline_states_match}"
        )
    np.save(RAW / "Tperp_eta_complex.npy", eta_all)
    np.save(RAW / "Tperp_real_basis.npy", basis_all)
    np.save(RAW / "Tperp_projectors.npy", projector_all)
    np.save(RAW / "S4_t0_basis.npy", S4_t0)
    write_csv(RAW / "transverse_direction_reconstruction.csv", list(rows[0]), rows)
    summary = {
        "stage": "A2c",
        "status": "TRANSVERSE_RECONSTRUCTED",
        "t0": t0,
        "crossing": replay["crossing"],
        "direction_prng_seed": CFG["transverse"]["direction_prng_seed"],
        "unique_seed_directions": CFG["transverse"]["seeds"],
        "epsilons": CFG["transverse"]["epsilons"],
        "eta_rule": "original eta_r/eta_i projected separately into S4(t0)^perp and jointly normalized",
        "comparison_subspace_rule": "QR basis of span(eta.real, eta.imag)",
        "new_arbitrary_direction_added": False,
        "s4_saved_projector_error": s4_saved_error,
        "stage_a0_first_record_strings_all_match": all_record_strings_match,
        "baseline_states_bitwise_match_replay": all_baseline_states_match,
    }
    (LOGS / "transverse_reconstruction_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("TRANSVERSE_RECONSTRUCTED: 3 existing directions × 4 existing eps; Stage A0 records exact")


if __name__ == "__main__":
    main()
