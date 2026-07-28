"""Run C0 versus reversed Candidate 1 for 32 repeated System A collisions."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

import numpy as np

from system_A_experimental_copy import (
    ExperimentConfig,
    norm2,
    run_repeated,
)


STAGE_E_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = STAGE_E_ROOT / "data"
FIGURE_ROOT = STAGE_E_ROOT / "figures"
os.environ.setdefault(
    "MPLCONFIGDIR", str(DATA_ROOT / ".matplotlib-cache")
)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CONFIG = ExperimentConfig()
KERNELS = ("C0", "C1_reversed")
KAPPA_VALUES = (0.01, 0.1, 1.0)
CASES = ("F_x_F", "B_x_B", "F_x_B")


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _max_abs(rows: list[dict], *fields: str) -> float:
    return max(abs(float(row[field])) for row in rows for field in fields)


def _run_summary(rows: list[dict]) -> list[dict]:
    grouped = {}
    for row in rows:
        grouped.setdefault(row["run_id"], []).append(row)
    summaries = []
    for run_id, records in sorted(grouped.items()):
        records.sort(key=lambda row: int(row["collision"]))
        first = records[0]
        last = records[-1]
        summaries.append(
            {
                "run_id": run_id,
                "kernel": first["kernel"],
                "case_id": first["case_id"],
                "kappa_requested": first["kappa_requested"],
                "collision_count": len(records),
                "R_0": first["R_0"],
                "R_eff_first": first["R_eff"],
                "R_eff_last": last["R_eff"],
                "R_eff_min": min(float(row["R_eff"]) for row in records),
                "R_eff_max": max(float(row["R_eff"]) for row in records),
                "R_eff_spread": max(
                    float(row["R_eff"]) for row in records
                )
                - min(float(row["R_eff"]) for row in records),
                "c_A_first": first["c_A"],
                "c_B_first": first["c_B"],
                "c_A_last": last["c_A"],
                "c_B_last": last["c_B"],
                "c_mean_max_abs": max(
                    abs(float(row["c_mean"])) for row in records
                ),
                "c_A_min": min(float(row["c_A"]) for row in records),
                "c_A_max": max(float(row["c_A"]) for row in records),
                "c_B_min": min(float(row["c_B"]) for row in records),
                "c_B_max": max(float(row["c_B"]) for row in records),
                "a_raw_norm2_min": min(
                    float(row["a_raw_norm2"]) for row in records
                ),
                "a_raw_norm2_max": max(
                    float(row["a_raw_norm2"]) for row in records
                ),
                "b_raw_norm2_min": min(
                    float(row["b_raw_norm2"]) for row in records
                ),
                "b_raw_norm2_max": max(
                    float(row["b_raw_norm2"]) for row in records
                ),
                "max_abs_I_A": max(
                    abs(float(row["I_A"])) for row in records
                ),
                "max_abs_I_B": max(
                    abs(float(row["I_B"])) for row in records
                ),
                "max_unitarity_residual": max(
                    abs(float(row["unitarity_residual"]))
                    for row in records
                ),
                "max_path_sum_residual": max(
                    max(
                        abs(float(row["path_sum_residual_A"])),
                        abs(float(row["path_sum_residual_B"])),
                    )
                    for row in records
                ),
                "max_pair_norm_residual": max(
                    abs(float(row["pair_norm_conservation_residual"]))
                    for row in records
                ),
                "final_p_B_A": last["p_B_A_raw"],
                "final_p_F_A": last["p_F_A_raw"],
                "final_p_B_B": last["p_B_B_raw"],
                "final_p_F_B": last["p_F_B_raw"],
                "final_origin_A_weight_in_A": last[
                    "origin_A_weight_in_A_raw"
                ],
                "final_origin_A_weight_in_B": last[
                    "origin_A_weight_in_B_raw"
                ],
            }
        )
    return summaries


def _comparison_rows(
    rows: list[dict],
    final_states: dict[str, np.ndarray],
) -> list[dict]:
    lookup = {
        (
            row["kernel"],
            float(row["kappa_requested"]),
            row["case_id"],
            int(row["collision"]),
        ): row
        for row in rows
    }
    comparisons = []
    for kappa in KAPPA_VALUES:
        for case_id in CASES:
            per_collision = []
            for collision in range(1, CONFIG.collision_count + 1):
                c0 = lookup[("C0", kappa, case_id, collision)]
                new = lookup[
                    ("C1_reversed", kappa, case_id, collision)
                ]
                per_collision.append(
                    {
                        "R_eff_difference": (
                            float(new["R_eff"]) - float(c0["R_eff"])
                        ),
                        "path_a_to_a_difference": (
                            float(new["path_a_to_a_norm2"])
                            - float(c0["path_a_to_a_norm2"])
                        ),
                        "path_b_to_a_difference": (
                            float(new["path_b_to_a_norm2"])
                            - float(c0["path_b_to_a_norm2"])
                        ),
                        "p_B_A_difference": (
                            float(new["p_B_A_raw"])
                            - float(c0["p_B_A_raw"])
                        ),
                        "p_F_A_difference": (
                            float(new["p_F_A_raw"])
                            - float(c0["p_F_A_raw"])
                        ),
                    }
                )
            c0_id = f"C0__{case_id}__kappa_{kappa:g}"
            new_id = f"C1_reversed__{case_id}__kappa_{kappa:g}"
            final_difference = math.sqrt(
                norm2(
                    final_states[f"{new_id}__A"]
                    - final_states[f"{c0_id}__A"]
                )
                + norm2(
                    final_states[f"{new_id}__B"]
                    - final_states[f"{c0_id}__B"]
                )
            )
            first_c0 = lookup[("C0", kappa, case_id, 1)]
            first_new = lookup[
                ("C1_reversed", kappa, case_id, 1)
            ]
            comparisons.append(
                {
                    "case_id": case_id,
                    "kappa": kappa,
                    "R_0": CONFIG.reflection_baseline,
                    "C0_R_eff": first_c0["R_eff"],
                    "C1_reversed_R_eff": first_new["R_eff"],
                    "R_eff_shift": (
                        float(first_new["R_eff"])
                        - float(first_c0["R_eff"])
                    ),
                    "R_eff_relation": (
                        "enhanced"
                        if float(first_new["R_eff"])
                        > CONFIG.reflection_baseline + 1.0e-14
                        else "suppressed"
                        if float(first_new["R_eff"])
                        < CONFIG.reflection_baseline - 1.0e-14
                        else "baseline"
                    ),
                    "max_abs_R_eff_difference_over_32": max(
                        abs(item["R_eff_difference"])
                        for item in per_collision
                    ),
                    "max_abs_path_a_to_a_difference_over_32": max(
                        abs(item["path_a_to_a_difference"])
                        for item in per_collision
                    ),
                    "max_abs_path_b_to_a_difference_over_32": max(
                        abs(item["path_b_to_a_difference"])
                        for item in per_collision
                    ),
                    "max_abs_p_B_A_difference_over_32": max(
                        abs(item["p_B_A_difference"])
                        for item in per_collision
                    ),
                    "max_abs_p_F_A_difference_over_32": max(
                        abs(item["p_F_A_difference"])
                        for item in per_collision
                    ),
                    "final_state_pair_L2_difference": final_difference,
                }
            )
    return comparisons


def _build_summary(
    rows: list[dict],
    run_summaries: list[dict],
    comparisons: list[dict],
) -> dict:
    expected_rows = (
        len(KERNELS)
        * len(KAPPA_VALUES)
        * len(CASES)
        * CONFIG.collision_count
    )
    if len(rows) != expected_rows:
        raise ValueError(f"row count {len(rows)} != {expected_rows}")
    new_comparisons = {
        (row["case_id"], float(row["kappa"])): row
        for row in comparisons
    }
    sign_failures = []
    for kappa in KAPPA_VALUES:
        expected = {
            "F_x_F": "enhanced",
            "B_x_B": "suppressed",
            "F_x_B": "baseline",
        }
        for case_id, relation in expected.items():
            if new_comparisons[(case_id, kappa)]["R_eff_relation"] != relation:
                sign_failures.append(
                    {
                        "case_id": case_id,
                        "kappa": kappa,
                        "expected": relation,
                        "actual": new_comparisons[
                            (case_id, kappa)
                        ]["R_eff_relation"],
                    }
                )
    c0_summaries = [
        row for row in run_summaries if row["kernel"] == "C0"
    ]
    new_summaries = [
        row
        for row in run_summaries
        if row["kernel"] == "C1_reversed"
    ]
    c0_by_case = {}
    for row in c0_summaries:
        c0_by_case.setdefault(row["case_id"], []).append(row)
    c0_kappa_spread = max(
        max(float(row["R_eff_first"]) for row in records)
        - min(float(row["R_eff_first"]) for row in records)
        for records in c0_by_case.values()
    )
    return {
        "schema": "stage_E_reversed_candidate1_repeated_scattering_v1",
        "status": "PASS" if not sign_failures else "FAIL",
        "purpose": (
            "Implement the reversed Candidate 1 main hypothesis in an "
            "independent System A experimental copy and compare it with C0."
        ),
        "implementation_formula": (
            "theta_eff = theta_0 - kappa*rho(theta_0)*(c_A+c_B)/2"
        ),
        "configuration": {
            "physical_grid": [CONFIG.u_grid_n, CONFIG.eta_grid_n],
            "K": CONFIG.k_components,
            "R_0": CONFIG.reflection_baseline,
            "collision_count": CONFIG.collision_count,
            "kappa_values": list(KAPPA_VALUES),
            "cases": list(CASES),
            "kernels": list(KERNELS),
            "carrier_eta_lineages": {
                "origin_A": {"q": 1.0, "eta_mode": 1},
                "origin_B": {"q": -1.0, "eta_mode": 2},
            },
            "raw_output_is_next_physical_state": True,
            "channel_normalization_applied": False,
        },
        "row_count": len(rows),
        "expected_row_count": expected_rows,
        "run_count": len(run_summaries),
        "hypothesis_sign_check": {
            "status": "PASS" if not sign_failures else "FAIL",
            "failure_count": len(sign_failures),
            "failures": sign_failures,
        },
        "C1_reversed_effective_R": [
            {
                "case_id": row["case_id"],
                "kappa": row["kappa"],
                "R_eff": row["C1_reversed_R_eff"],
                "R_eff_shift": row["R_eff_shift"],
                "relation_to_R_0": row["R_eff_relation"],
            }
            for row in comparisons
        ],
        "maximum_residuals": {
            "unitarity": _max_abs(rows, "unitarity_residual"),
            "coefficient_orthogonality": _max_abs(
                rows, "coefficient_orthogonality_residual"
            ),
            "path_sum": _max_abs(
                rows, "path_sum_residual_A", "path_sum_residual_B"
            ),
            "pair_norm_conservation": _max_abs(
                rows, "pair_norm_conservation_residual"
            ),
            "lineage_reconstruction": _max_abs(
                rows,
                "input_reconstruction_residual_A",
                "input_reconstruction_residual_B",
                "raw_reconstruction_residual_A",
                "raw_reconstruction_residual_B",
            ),
            "interference": _max_abs(rows, "I_A", "I_B"),
            "raw_channel_norm_difference_from_one": max(
                abs(float(row[field]) - 1.0)
                for row in rows
                for field in ("a_raw_norm2", "b_raw_norm2")
            ),
        },
        "repeated_readout": {
            "max_pure_R_eff_drift": max(
                float(row["R_eff_spread"])
                for row in new_summaries
                if row["case_id"] in ("F_x_F", "B_x_B")
            ),
            "max_FxB_abs_c_mean": max(
                float(row["c_mean_max_abs"])
                for row in new_summaries
                if row["case_id"] == "F_x_B"
            ),
            "max_FxB_R_eff_departure_from_R_0": max(
                max(
                    abs(float(row["R_eff_min"]) - CONFIG.reflection_baseline),
                    abs(float(row["R_eff_max"]) - CONFIG.reflection_baseline),
                )
                for row in new_summaries
                if row["case_id"] == "F_x_B"
            ),
        },
        "C0_control": {
            "max_R_eff_spread_across_requested_kappa": c0_kappa_spread,
            "duplicated_per_kappa_for_direct_pairing": True,
        },
        "comparison_rows": comparisons,
        "original_System_A_modified": False,
        "N_system_integrated": False,
        "formal_integration_performed": False,
    }


def _select(
    rows: list[dict],
    *,
    kernel: str,
    case_id: str,
    kappa: float,
) -> list[dict]:
    selected = [
        row
        for row in rows
        if row["kernel"] == kernel
        and row["case_id"] == case_id
        and float(row["kappa_requested"]) == kappa
    ]
    return sorted(selected, key=lambda row: int(row["collision"]))


def _save_figures(rows: list[dict]) -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    colors = {0.01: "tab:green", 0.1: "tab:orange", 1.0: "tab:red"}

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)
    for axis, case_id in zip(axes, CASES):
        baseline = _select(
            rows, kernel="C0", case_id=case_id, kappa=0.01
        )
        axis.plot(
            [row["collision"] for row in baseline],
            [row["R_eff"] for row in baseline],
            color="black",
            linestyle="--",
            label="C0",
        )
        for kappa in KAPPA_VALUES:
            selected = _select(
                rows,
                kernel="C1_reversed",
                case_id=case_id,
                kappa=kappa,
            )
            axis.plot(
                [row["collision"] for row in selected],
                [row["R_eff"] for row in selected],
                color=colors[kappa],
                label=f"reversed C1, k={kappa:g}",
            )
        axis.set_title(case_id.replace("_x_", " x "))
        axis.set_xlabel("collision")
        axis.set_ylim(0.3, 1.0)
    axes[0].set_ylabel("effective reflection probability")
    axes[-1].legend(fontsize=8, loc="best")
    fig.suptitle("Repeated state-dependent scattering versus C0")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "R_eff_by_collision.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14, 4.4),
        sharey=True,
        constrained_layout=True,
    )
    for axis, case_id in zip(axes, CASES):
        selected = _select(
            rows,
            kernel="C1_reversed",
            case_id=case_id,
            kappa=1.0,
        )
        collision = [row["collision"] for row in selected]
        axis.plot(collision, [row["c_A"] for row in selected], label="c_A")
        axis.plot(collision, [row["c_B"] for row in selected], label="c_B")
        axis.plot(
            collision,
            [row["c_mean"] for row in selected],
            linestyle="--",
            label="mean",
        )
        axis.set_title(case_id.replace("_x_", " x "))
        axis.set_xlabel("collision")
        axis.set_ylim(-1.05, 1.05)
    axes[0].set_ylabel("input parity readout")
    axes[-1].legend(fontsize=8)
    fig.suptitle("Parity recomputed from source-resolved kernels (kappa=1)")
    fig.savefig(FIGURE_ROOT / "parity_readout_by_collision.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)
    for axis, case_id in zip(axes, CASES):
        for kappa in KAPPA_VALUES:
            selected = _select(
                rows,
                kernel="C1_reversed",
                case_id=case_id,
                kappa=kappa,
            )
            axis.plot(
                [row["collision"] for row in selected],
                [row["path_a_to_a_norm2"] for row in selected],
                color=colors[kappa],
                label=f"A->A, k={kappa:g}",
            )
            axis.plot(
                [row["collision"] for row in selected],
                [row["path_b_to_a_norm2"] for row in selected],
                color=colors[kappa],
                linestyle=":",
                label=f"B->A, k={kappa:g}",
            )
        axis.set_title(case_id.replace("_x_", " x "))
        axis.set_xlabel("collision")
        axis.set_ylim(0.0, 1.0)
    axes[0].set_ylabel("squared path norm")
    axes[-1].legend(fontsize=7)
    fig.suptitle("Path norms under the reversed Candidate 1")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "path_norms_by_collision.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9, 5))
    for kernel, linestyle in (("C0", "--"), ("C1_reversed", "-")):
        selected = _select(
            rows, kernel=kernel, case_id="F_x_B", kappa=1.0
        )
        axis.plot(
            [row["collision"] for row in selected],
            [row["origin_A_weight_in_A_raw"] for row in selected],
            linestyle=linestyle,
            label=f"{kernel}: origin A weight in channel A",
        )
        axis.plot(
            [row["collision"] for row in selected],
            [row["origin_A_weight_in_B_raw"] for row in selected],
            linestyle=linestyle,
            alpha=0.7,
            label=f"{kernel}: origin A weight in channel B",
        )
    axis.set_xlabel("collision")
    axis.set_ylabel("lineage weight")
    axis.set_ylim(0.0, 1.0)
    axis.set_title("F x B lineage transport; C0 and reversed C1 overlap")
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "FxB_lineage_transport.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    selected = _select(
        rows,
        kernel="C1_reversed",
        case_id="F_x_F",
        kappa=1.0,
    )
    collision = [row["collision"] for row in selected]
    axes[0].plot(
        collision, [row["a_raw_norm2"] for row in selected], label="A raw"
    )
    axes[0].plot(
        collision, [row["b_raw_norm2"] for row in selected], label="B raw"
    )
    axes[0].set_title("raw channel norms")
    axes[0].set_xlabel("collision")
    axes[0].legend()
    axes[1].plot(
        collision,
        [abs(float(row["I_A"])) for row in selected],
        label="|I_A|",
    )
    axes[1].plot(
        collision,
        [abs(float(row["I_B"])) for row in selected],
        label="|I_B|",
    )
    axes[1].set_yscale("log")
    axes[1].set_title("path interference magnitude")
    axes[1].set_xlabel("collision")
    axes[1].legend()
    fig.suptitle("Raw propagation diagnostics (F x F, kappa=1)")
    fig.tight_layout()
    fig.savefig(
        FIGURE_ROOT / "raw_norm_and_interference_diagnostics.png",
        dpi=180,
    )
    plt.close(fig)


def run() -> dict:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    final_states = {}
    for kernel_name in KERNELS:
        for kappa in KAPPA_VALUES:
            for case_id in CASES:
                rows, final_a, final_b = run_repeated(
                    case_id, kernel_name, kappa, CONFIG
                )
                all_rows.extend(rows)
                run_id = rows[0]["run_id"]
                final_states[f"{run_id}__A"] = final_a
                final_states[f"{run_id}__B"] = final_b

    run_summaries = _run_summary(all_rows)
    comparisons = _comparison_rows(all_rows, final_states)
    summary = _build_summary(all_rows, run_summaries, comparisons)
    if summary["status"] != "PASS":
        raise RuntimeError("main hypothesis sign check failed")

    _write_csv(DATA_ROOT / "stage_E_collision_results.csv", all_rows)
    _write_csv(DATA_ROOT / "stage_E_run_summary.csv", run_summaries)
    _write_csv(DATA_ROOT / "stage_E_C0_vs_reversed_C1.csv", comparisons)
    (DATA_ROOT / "stage_E_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    np.savez_compressed(
        DATA_ROOT / "stage_E_final_raw_states.npz",
        **final_states,
    )
    _save_figures(all_rows)
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, ensure_ascii=False))
