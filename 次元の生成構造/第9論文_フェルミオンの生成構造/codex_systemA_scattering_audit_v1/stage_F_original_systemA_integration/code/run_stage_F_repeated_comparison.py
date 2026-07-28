"""Run Stage F-A at the existing N_A=1, N_B=63, R0=0.55 condition."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from system_A_stage_F_copy import KAPPA_VALUES, Params, run_series


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
LOG_DIR = STAGE_ROOT / "logs"
OUTPUT = DATA_DIR / "stage_F_repeated_collision_results.csv"


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to create empty output: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    params = Params()
    rows: list[dict] = []
    final_states: dict[str, np.ndarray] = {}
    for scattering_mode in ("C0", "reversed_C1"):
        for normalization_mode in (
            "existing_normalization",
            "raw_update",
        ):
            for kappa in KAPPA_VALUES:
                result = run_series(
                    params,
                    reflection_baseline=0.55,
                    scattering_mode=scattering_mode,
                    normalization_mode=normalization_mode,
                    kappa=kappa,
                )
                rows.extend(result.rows)
                key = (
                    f"{scattering_mode}__{normalization_mode}__"
                    f"kappa_{kappa:g}__R0_0p55"
                ).replace(".", "p")
                final_states[f"{key}__A"] = result.final_a
                final_states[f"{key}__B"] = result.final_b
    write_csv(OUTPUT, rows)
    np.savez_compressed(
        LOG_DIR / "stage_F_FA_final_states.npz", **final_states
    )
    maximum_residual = max(
        max(
            float(row["unitarity_residual"]),
            float(row["coefficient_orthogonality_residual"]),
            float(row["path_sum_residual_A"]),
            float(row["path_sum_residual_B"]),
            float(row["total_norm_conservation_residual"]),
            float(row["demodulation_reconstruction_residual"]),
            float(row["parity_projection_sum_residual"]),
        )
        for row in rows
    )
    result = {
        "status": "complete",
        "experiment": "F-A",
        "run_count": 12,
        "row_count": len(rows),
        "collision_count_per_run": params.recursive_collision_count,
        "maximum_numerical_residual": maximum_residual,
        "theta_range_violation_count": sum(
            int(bool(row["theta_range_violation"])) for row in rows
        ),
        "nan_inf_count": sum(int(row["nan_inf_count"]) for row in rows),
    }
    (LOG_DIR / "stage_F_FA_run.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
