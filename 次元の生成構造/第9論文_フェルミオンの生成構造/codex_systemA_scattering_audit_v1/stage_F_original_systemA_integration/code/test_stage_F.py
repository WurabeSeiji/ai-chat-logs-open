"""Deterministic unit tests for the Stage F integration copy."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from parity_demodulation import channel_parity, norm2
from state_dependent_scattering import effective_angle, scatter_once
from system_A_stage_F_copy import (
    Params,
    initial_state_pair,
    make_grids,
    run_series,
)


STAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[5]
OUTPUT = STAGE_ROOT / "logs" / "stage_F_test_results.json"
TOLERANCE = 5.0e-12


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _assert_below(name: str, value: float, limit: float = TOLERANCE) -> None:
    if value > limit:
        raise AssertionError(f"{name}: {value} > {limit}")


def run_tests() -> dict:
    original_13 = (
        REPO_ROOT
        / "次元の生成構造/第9論文_フェルミオンの生成構造/"
        "対照実験_波束収縮_実行環境_v1/20260713/"
        "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"
    )
    copy_13 = (
        STAGE_ROOT
        / "source_copy/"
        "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1_ORIGINAL_SNAPSHOT.py"
    )
    original_15 = (
        REPO_ROOT
        / "次元の生成構造/第9論文_フェルミオンの生成構造/"
        "対照実験_波束収縮_実行環境_v1/20260715/"
        "run_system_A_localization_exchange_R_sweep_preliminary_v1.py"
    )
    copy_15 = (
        STAGE_ROOT
        / "source_copy/"
        "run_system_A_localization_exchange_R_sweep_preliminary_v1_ORIGINAL_SNAPSHOT.py"
    )
    if _sha(original_13) != _sha(copy_13) or _sha(original_15) != _sha(copy_15):
        raise AssertionError("source snapshot hash mismatch")

    params = Params()
    u, eta = make_grids(params)
    initial_a, initial_b = initial_state_pair(params)
    parity_a = channel_parity(initial_a, u, eta, params.p0)
    parity_b = channel_parity(initial_b, u, eta, params.p0)
    _assert_below("initial A odd parity", abs(parity_a.c_pi + 1.0))
    _assert_below("initial B odd parity", abs(parity_b.c_pi + 1.0))
    _assert_below(
        "initial demodulation",
        max(
            parity_a.reconstruction_residual,
            parity_b.reconstruction_residual,
            parity_a.projection_sum_residual,
            parity_b.projection_sum_residual,
        ),
    )

    c0 = effective_angle("C0", 0.55, 1.0, -1.0, -1.0)
    reversed_c1 = effective_angle(
        "reversed_C1", 0.55, 1.0, -1.0, -1.0
    )
    if not reversed_c1.R_eff > 0.55:
        raise AssertionError("reversed Candidate 1 direction is incorrect")
    _assert_below("C0 baseline", abs(c0.R_eff - 0.55))
    for endpoint in (0.0, 1.0):
        endpoint_angle = effective_angle(
            "reversed_C1", endpoint, 1.0, -1.0, -1.0
        )
        _assert_below(
            f"endpoint {endpoint}", abs(endpoint_angle.R_eff - endpoint)
        )

    collision = scatter_once(
        initial_a,
        initial_b,
        scattering_mode="C0",
        reflection_baseline=0.55,
        kappa=1.0,
        u=u,
        eta=eta,
        p0=params.p0,
    )
    _assert_below("unitarity", collision.angle.unitarity_residual)
    _assert_below(
        "coefficient orthogonality",
        collision.angle.coefficient_orthogonality_residual,
    )
    _assert_below(
        "path sum",
        max(
            collision.path_sum_residual_a,
            collision.path_sum_residual_b,
        ),
    )
    _assert_below("total norm", collision.total_norm_residual)

    normalized = run_series(
        params,
        reflection_baseline=0.55,
        scattering_mode="C0",
        normalization_mode="existing_normalization",
        kappa=0.01,
        collision_count=4,
    )
    raw = run_series(
        params,
        reflection_baseline=0.55,
        scattering_mode="C0",
        normalization_mode="raw_update",
        kappa=0.01,
        collision_count=4,
    )
    difference = math.sqrt(
        norm2(normalized.final_a - raw.final_a)
        + norm2(normalized.final_b - raw.final_b)
    )
    _assert_below("normalization/raw baseline difference", difference)
    if {
        row["normalization_mode"] for row in normalized.rows + raw.rows
    } != {"existing_normalization", "raw_update"}:
        raise AssertionError("normalization series were not kept distinct")

    results = {
        "status": "pass",
        "tolerance": TOLERANCE,
        "source_snapshot_hashes_match": True,
        "initial_c_A": parity_a.c_pi,
        "initial_c_B": parity_b.c_pi,
        "C0_R_eff": c0.R_eff,
        "reversed_C1_R_eff_kappa1": reversed_c1.R_eff,
        "normalization_raw_state_difference_4_collisions": difference,
        "maximum_checked_residual": max(
            parity_a.reconstruction_residual,
            parity_b.reconstruction_residual,
            collision.angle.unitarity_residual,
            collision.angle.coefficient_orthogonality_residual,
            collision.path_sum_residual_a,
            collision.path_sum_residual_b,
            collision.total_norm_residual,
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(results, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return results


if __name__ == "__main__":
    print(json.dumps(run_tests(), ensure_ascii=False, indent=2))
