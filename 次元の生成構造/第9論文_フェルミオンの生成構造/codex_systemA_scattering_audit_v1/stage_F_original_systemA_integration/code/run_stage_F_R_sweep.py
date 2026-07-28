"""Run Stage F-B on the seven reflection points in the 20260713 System A."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from system_A_stage_F_copy import (
    KAPPA_VALUES,
    ORIGINAL_R_VALUES,
    Params,
    run_series,
)


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
OUTPUT = DATA_DIR / "stage_F_R_sweep_results.csv"
SUMMARY_OUTPUT = DATA_DIR / "stage_F_run_summary.csv"
FINAL_STATES_OUTPUT = DATA_DIR / "stage_F_final_states.npz"


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to create empty output: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict]) -> dict:
    min_l = min(rows, key=lambda row: abs(row["L_A"] - row["L_B"]))
    min_n = min(
        rows, key=lambda row: abs(row["N_eff_A"] - row["N_eff_B"])
    )
    return {
        "run_id": rows[0]["run_id"],
        "scattering_mode": rows[0]["scattering_mode"],
        "normalization_mode": rows[0]["normalization_mode"],
        "kappa": rows[0]["kappa"],
        "R0": rows[0]["R0"],
        "collision_count": len(rows),
        "min_L_difference": abs(min_l["L_A"] - min_l["L_B"]),
        "min_L_difference_collision": min_l["collision_index"],
        "min_N_eff_difference": abs(
            min_n["N_eff_A"] - min_n["N_eff_B"]
        ),
        "min_N_eff_difference_collision": min_n["collision_index"],
        "R_eff_min": min(row["R_eff"] for row in rows),
        "R_eff_max": max(row["R_eff"] for row in rows),
        "c_A_min": min(row["c_A"] for row in rows),
        "c_A_max": max(row["c_A"] for row in rows),
        "c_B_min": min(row["c_B"] for row in rows),
        "c_B_max": max(row["c_B"] for row in rows),
        "c_mean_min": min(row["c_mean"] for row in rows),
        "c_mean_max": max(row["c_mean"] for row in rows),
        "raw_norm_min": min(
            min(row["raw_norm_A"], row["raw_norm_B"]) for row in rows
        ),
        "raw_norm_max": max(
            max(row["raw_norm_A"], row["raw_norm_B"]) for row in rows
        ),
        "normalization_scale_min": min(
            min(
                row["normalization_scale_A"],
                row["normalization_scale_B"],
            )
            for row in rows
        ),
        "normalization_scale_max": max(
            max(
                row["normalization_scale_A"],
                row["normalization_scale_B"],
            )
            for row in rows
        ),
        "maximum_numerical_residual": max(
            max(
                row["unitarity_residual"],
                row["coefficient_orthogonality_residual"],
                row["path_sum_residual_A"],
                row["path_sum_residual_B"],
                row["total_norm_conservation_residual"],
                row["demodulation_reconstruction_residual"],
                row["parity_projection_sum_residual"],
            )
            for row in rows
        ),
        "nan_inf_count": sum(row["nan_inf_count"] for row in rows),
        "theta_range_violation_count": sum(
            int(bool(row["theta_range_violation"])) for row in rows
        ),
    }


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    params = Params()
    all_rows: list[dict] = []
    summaries: list[dict] = []
    final_states: dict[str, np.ndarray] = {}
    for reflection in ORIGINAL_R_VALUES:
        for scattering_mode in ("C0", "reversed_C1"):
            for normalization_mode in (
                "existing_normalization",
                "raw_update",
            ):
                for kappa in KAPPA_VALUES:
                    result = run_series(
                        params,
                        reflection_baseline=reflection,
                        scattering_mode=scattering_mode,
                        normalization_mode=normalization_mode,
                        kappa=kappa,
                    )
                    all_rows.extend(result.rows)
                    summaries.append(summarize(result.rows))
                    key = (
                        f"{scattering_mode}__{normalization_mode}__"
                        f"kappa_{kappa:g}__R0_{reflection:g}"
                    ).replace(".", "p")
                    final_states[f"{key}__A"] = result.final_a
                    final_states[f"{key}__B"] = result.final_b
    write_csv(OUTPUT, all_rows)
    write_csv(SUMMARY_OUTPUT, summaries)
    np.savez_compressed(FINAL_STATES_OUTPUT, **final_states)
    result = {
        "status": "complete",
        "experiment": "F-B",
        "R_values": list(ORIGINAL_R_VALUES),
        "run_count": len(summaries),
        "row_count": len(all_rows),
        "collision_count_per_run": params.recursive_collision_count,
        "theta_range_violation_count": sum(
            int(row["theta_range_violation_count"]) for row in summaries
        ),
        "nan_inf_count": sum(
            int(row["nan_inf_count"]) for row in summaries
        ),
        "maximum_numerical_residual": max(
            float(row["maximum_numerical_residual"]) for row in summaries
        ),
    }
    (STAGE_ROOT / "logs" / "stage_F_FB_run.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
