"""Reproduce the original N_A=1, N_B=63, R=0.55 C0 series."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

from system_A_stage_F_copy import (
    Params,
    distribution_similarity,
    harmonic_distribution,
    initial_state_pair,
    run_series,
    state_metrics,
)


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
REPORT_DIR = STAGE_ROOT / "reports"
SNAPSHOT = (
    STAGE_ROOT
    / "source_copy/"
    "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1_ORIGINAL_SNAPSHOT.py"
)
OUTPUT = DATA_DIR / "stage_F_reproduction_gate.csv"
SUMMARY_OUTPUT = DATA_DIR / "stage_F_reproduction_gate_summary.json"
ABS_TOLERANCE = 5.0e-12
REL_TOLERANCE = 5.0e-10
FIELDS = (
    "L_A",
    "L_B",
    "N_eff_A",
    "N_eff_B",
    "spectral_similarity_A_to_initial_A",
    "spectral_similarity_A_to_initial_B",
    "spectral_similarity_B_to_initial_A",
    "spectral_similarity_B_to_initial_B",
    "channel_norm_A",
    "channel_norm_B",
)


def load_snapshot() -> Any:
    spec = importlib.util.spec_from_file_location(
        "stage_f_original_system_a_snapshot", SNAPSHOT
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load original System A snapshot")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def canonical_rows(src: Any) -> list[dict[str, float | int]]:
    params = src.Params()
    reflection = 0.55
    delta = src.delta_from_reflection_rate(reflection)
    t, r, _, _ = src.scattering_coefficients(delta)
    a = src.make_state(
        params, 1, params.q_A, params.m_A, True, params.A_A
    )
    b = src.make_state(
        params,
        params.high_n,
        params.q_B,
        params.m_B,
        True,
        params.A_B,
    )
    h_a0 = src.harmonic_distribution(params, a)
    h_b0 = src.harmonic_distribution(params, b)
    rows: list[dict[str, float | int]] = []
    for collision_index in range(params.recursive_collision_count + 1):
        h_a = src.harmonic_distribution(params, a)
        h_b = src.harmonic_distribution(params, b)
        rows.append(
            {
                "collision_index": collision_index,
                "L_A": src.localization(a),
                "L_B": src.localization(b),
                "N_eff_A": src.effective_n(h_a)[0],
                "N_eff_B": src.effective_n(h_b)[0],
                "spectral_similarity_A_to_initial_A": (
                    distribution_similarity(h_a, h_a0)
                ),
                "spectral_similarity_A_to_initial_B": (
                    distribution_similarity(h_a, h_b0)
                ),
                "spectral_similarity_B_to_initial_A": (
                    distribution_similarity(h_b, h_a0)
                ),
                "spectral_similarity_B_to_initial_B": (
                    distribution_similarity(h_b, h_b0)
                ),
                "channel_norm_A": src.norm2(a),
                "channel_norm_B": src.norm2(b),
            }
        )
        if collision_index < params.recursive_collision_count:
            a_next = src.normalize(r * a + t * b)
            b_next = src.normalize(t * a + r * b)
            a, b = a_next, b_next
    return rows


def integrated_rows() -> list[dict[str, float | int]]:
    params = Params()
    initial_a, initial_b = initial_state_pair(params)
    initial_h_a = harmonic_distribution(params, initial_a)
    initial_h_b = harmonic_distribution(params, initial_b)
    initial_metrics_a = state_metrics(
        params, initial_a, initial_h_a, initial_h_b
    )
    initial_metrics_b = state_metrics(
        params, initial_b, initial_h_a, initial_h_b
    )
    rows: list[dict[str, float | int]] = [
        {
            "collision_index": 0,
            "L_A": initial_metrics_a["L"],
            "L_B": initial_metrics_b["L"],
            "N_eff_A": initial_metrics_a["N_eff"],
            "N_eff_B": initial_metrics_b["N_eff"],
            "spectral_similarity_A_to_initial_A": initial_metrics_a[
                "similarity_to_initial_A"
            ],
            "spectral_similarity_A_to_initial_B": initial_metrics_a[
                "similarity_to_initial_B"
            ],
            "spectral_similarity_B_to_initial_A": initial_metrics_b[
                "similarity_to_initial_A"
            ],
            "spectral_similarity_B_to_initial_B": initial_metrics_b[
                "similarity_to_initial_B"
            ],
            "channel_norm_A": initial_metrics_a["norm2"],
            "channel_norm_B": initial_metrics_b["norm2"],
        }
    ]
    result = run_series(
        params,
        reflection_baseline=0.55,
        scattering_mode="C0",
        normalization_mode="existing_normalization",
        kappa=0.01,
    )
    for row in result.rows:
        rows.append(
            {
                "collision_index": row["collision_index"],
                "L_A": row["L_A"],
                "L_B": row["L_B"],
                "N_eff_A": row["N_eff_A"],
                "N_eff_B": row["N_eff_B"],
                "spectral_similarity_A_to_initial_A": row[
                    "spectral_similarity_A_to_initial_A"
                ],
                "spectral_similarity_A_to_initial_B": row[
                    "spectral_similarity_A_to_initial_B"
                ],
                "spectral_similarity_B_to_initial_A": row[
                    "spectral_similarity_B_to_initial_A"
                ],
                "spectral_similarity_B_to_initial_B": row[
                    "spectral_similarity_B_to_initial_B"
                ],
                "channel_norm_A": row["next_state_norm_A"],
                "channel_norm_B": row["next_state_norm_B"],
            }
        )
    return rows


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    canonical = canonical_rows(load_snapshot())
    integrated = integrated_rows()
    if len(canonical) != len(integrated):
        raise SystemExit("C0 gate failed: row count mismatch")

    comparison_rows = []
    maximum_absolute = {field: 0.0 for field in FIELDS}
    maximum_relative = {field: 0.0 for field in FIELDS}
    for reference, candidate in zip(canonical, integrated):
        output_row: dict[str, Any] = {
            "collision_index": reference["collision_index"]
        }
        for field in FIELDS:
            ref_value = float(reference[field])
            got_value = float(candidate[field])
            absolute = abs(got_value - ref_value)
            relative = absolute / max(abs(ref_value), 1.0e-300)
            maximum_absolute[field] = max(
                maximum_absolute[field], absolute
            )
            maximum_relative[field] = max(
                maximum_relative[field], relative
            )
            output_row[f"canonical_{field}"] = ref_value
            output_row[f"integrated_{field}"] = got_value
            output_row[f"absolute_error_{field}"] = absolute
            output_row[f"relative_error_{field}"] = relative
        comparison_rows.append(output_row)

    passed = all(
        maximum_absolute[field] <= ABS_TOLERANCE
        or maximum_relative[field] <= REL_TOLERANCE
        for field in FIELDS
    )
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0]))
        writer.writeheader()
        writer.writerows(comparison_rows)
    summary = {
        "status": "pass" if passed else "fail",
        "canonical_source": str(SNAPSHOT.relative_to(STAGE_ROOT)),
        "canonical_note": (
            "No pre-existing exact R=0.55 series artifact was found. "
            "The canonical series was executed from the byte-identical "
            "20260713 source snapshot inside Stage F."
        ),
        "condition": {
            "N_A": 1,
            "N_B": 63,
            "R": 0.55,
            "collision_count": 128,
            "normalization": "existing channel-wise normalization",
        },
        "absolute_tolerance": ABS_TOLERANCE,
        "relative_tolerance": REL_TOLERANCE,
        "maximum_absolute_error": maximum_absolute,
        "maximum_relative_error": maximum_relative,
    }
    SUMMARY_OUTPUT.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    max_abs = max(maximum_absolute.values())
    max_rel = max(maximum_relative.values())
    report = f"""# 01 C0再現ゲート

## 結果

`{summary["status"]}`。`N_A=1, N_B=63, R=0.55`、128回衝突、既存チャネル別正規化を再現した。

## 正本

既存ディレクトリにはR=0.55の各衝突を保存した正確なCSVが見つからなかった。このため、丸めた記憶値は判定に使わず、20260713正本のバイト同一スナップショットをStage F内で実行して正本系列を生成した。

## 誤差

- 最大絶対誤差: `{max_abs:.17g}`
- 最大相対誤差: `{max_rel:.17g}`
- 絶対許容値: `{ABS_TOLERANCE:.1e}`
- 相対許容値: `{REL_TOLERANCE:.1e}`

比較対象は `L_A`, `L_B`, `N_eff_A`, `N_eff_B`、初期A/Bへの4つのスペクトル余弦類似度、A/Bチャネルノルムである。

`B_to_A_transfer`を記載する場合、その意味は `spectral_similarity_to_initial_B; not path flux` であり、経路フラックスではない。
"""
    (REPORT_DIR / "01_C0_reproduction_gate.md").write_text(
        report, encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit("C0 reproduction gate failed; stop before F-A")


if __name__ == "__main__":
    main()
