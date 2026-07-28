"""C0 gate and Stage G-B three-mode System A comparison."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from system_A_stage_G_copy import KAPPA_VALUES, Params, run_series


STAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[5]
F_ROWS = (
    REPO_ROOT
    / "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "codex_systemA_scattering_audit_v1/"
    "stage_F_original_systemA_integration/data/"
    "stage_F_R_sweep_results.csv"
)
OUTPUT = STAGE_ROOT / "data" / "stage_G_collision_results.csv"
FINAL_STATES = STAGE_ROOT / "data" / "stage_G_final_states.npz"
GATE_LOG = STAGE_ROOT / "logs" / "stage_G_C0_reproduction.json"
GATE_TOLERANCE = 5.0e-12
COMPARE_FIELDS = (
    "c_A",
    "c_B",
    "c_mean",
    "theta_eff",
    "R_eff",
    "T_eff",
    "L_A",
    "L_B",
    "N_eff_A",
    "N_eff_B",
    "spectral_similarity_A_to_initial_A",
    "spectral_similarity_A_to_initial_B",
    "spectral_similarity_B_to_initial_A",
    "spectral_similarity_B_to_initial_B",
    "path_A_to_A_norm",
    "path_B_to_A_norm",
    "path_B_to_B_norm",
    "path_A_to_B_norm",
    "interference_A",
    "interference_B",
    "raw_norm_A",
    "raw_norm_B",
    "next_state_norm_A",
    "next_state_norm_B",
)


def read_stage_f_rows() -> list[dict]:
    with F_ROWS.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def reference_run(
    rows: list[dict],
    *,
    reflection: float,
    scattering_mode: str,
    normalization_mode: str,
    kappa: float,
) -> list[dict]:
    selected = [
        row
        for row in rows
        if abs(float(row["R0"]) - reflection) <= 1.0e-15
        and row["scattering_mode"] == scattering_mode
        and row["normalization_mode"] == normalization_mode
        and abs(float(row["kappa"]) - kappa) <= 1.0e-15
    ]
    return sorted(selected, key=lambda row: int(row["collision_index"]))


def compare_runs(
    candidate: list[dict], reference: list[dict]
) -> dict:
    if len(candidate) != len(reference):
        return {
            "passed": False,
            "reason": "row_count_mismatch",
            "candidate_count": len(candidate),
            "reference_count": len(reference),
        }
    maximums = {field: 0.0 for field in COMPARE_FIELDS}
    for left, right in zip(candidate, reference):
        if int(left["collision_index"]) != int(right["collision_index"]):
            raise ValueError("collision index mismatch")
        for field in COMPARE_FIELDS:
            maximums[field] = max(
                maximums[field],
                abs(float(left[field]) - float(right[field])),
            )
    return {
        "passed": max(maximums.values()) <= GATE_TOLERANCE,
        "tolerance": GATE_TOLERANCE,
        "maximum_absolute_errors": maximums,
        "overall_maximum_absolute_error": max(maximums.values()),
    }


def main() -> None:
    params = Params()
    stage_f_rows = read_stage_f_rows()

    normalized_gate = run_series(
        params,
        reflection_baseline=0.55,
        scattering_mode="C0",
        normalization_mode="existing_normalization",
        kappa=0.01,
    )
    gate_reference = reference_run(
        stage_f_rows,
        reflection=0.55,
        scattering_mode="C0",
        normalization_mode="existing_normalization",
        kappa=0.01,
    )
    gate = compare_runs(normalized_gate.rows, gate_reference)
    GATE_LOG.write_text(
        json.dumps(gate, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if not gate["passed"]:
        raise SystemExit("Stage G C0 reproduction failed; stop before G-B")

    rows: list[dict] = []
    final_states: dict[str, np.ndarray] = {}
    for row in normalized_gate.rows:
        rows.append(
            {
                "series_role": (
                    "C0_existing_normalization_reproduction"
                ),
                **row,
            }
        )
    final_states["C0__existing_normalization__kappa_0p01__R0_0p55__A"] = (
        normalized_gate.final_a
    )
    final_states["C0__existing_normalization__kappa_0p01__R0_0p55__B"] = (
        normalized_gate.final_b
    )

    stage_f_regressions = []
    for reflection in (0.55, 0.70):
        for scattering_mode in ("C0", "reversed_C1", "relational_C1"):
            for kappa in KAPPA_VALUES:
                result = run_series(
                    params,
                    reflection_baseline=reflection,
                    scattering_mode=scattering_mode,
                    normalization_mode="raw_update",
                    kappa=kappa,
                )
                for row in result.rows:
                    rows.append({"series_role": "primary", **row})
                key = (
                    f"{scattering_mode}__raw_update__"
                    f"kappa_{kappa:g}__R0_{reflection:g}"
                ).replace(".", "p")
                final_states[f"{key}__A"] = result.final_a
                final_states[f"{key}__B"] = result.final_b
                if scattering_mode in ("C0", "reversed_C1"):
                    reference = reference_run(
                        stage_f_rows,
                        reflection=reflection,
                        scattering_mode=scattering_mode,
                        normalization_mode="raw_update",
                        kappa=kappa,
                    )
                    comparison = compare_runs(result.rows, reference)
                    stage_f_regressions.append(
                        {
                            "R0": reflection,
                            "scattering_mode": scattering_mode,
                            "kappa": kappa,
                            **comparison,
                        }
                    )
                    if not comparison["passed"]:
                        raise SystemExit(
                            "Stage F mode regression failed; "
                            "stop before completing G-B"
                        )

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(FINAL_STATES, **final_states)
    result = {
        "status": "complete",
        "C0_existing_normalization_gate": gate,
        "primary_run_count": 18,
        "primary_row_count": 18 * params.recursive_collision_count,
        "total_row_count": len(rows),
        "stage_F_regression_count": len(stage_f_regressions),
        "stage_F_regression_maximum_error": max(
            row["overall_maximum_absolute_error"]
            for row in stage_f_regressions
        ),
        "stage_F_regressions": stage_f_regressions,
    }
    (STAGE_ROOT / "logs" / "stage_G_repeated_run.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "C0_gate_maximum_error": gate[
                    "overall_maximum_absolute_error"
                ],
                "primary_run_count": result["primary_run_count"],
                "total_row_count": result["total_row_count"],
                "stage_F_regression_maximum_error": result[
                    "stage_F_regression_maximum_error"
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
