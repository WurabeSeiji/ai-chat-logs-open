from __future__ import annotations

import csv
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_abc_multigauge_interference_readout_v2 import (
    Gauge,
    Params as BaseParams,
    compute_verdicts as compute_single_collision_verdicts,
    default_gauges,
    readout_all,
    simulate_single_collision,
    summarize_stage_readouts,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_interference_readout_robustness_sweep_result_v2"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class Params(BaseParams):
    min_gauge_count: int = 8
    robustness_readout_tol: float = 1.0e-9
    robustness_r_gauge_tol: float = 1.0e-9
    robustness_tr_tol: float = 1.0e-10


def prefixed_gauges(case: str, gauges: Iterable[Gauge]) -> List[Gauge]:
    return [Gauge(**{**asdict(gauge), "name": f"{case}:{gauge.name}"}) for gauge in gauges]


def make_phase_center_gauges(params: Params) -> List[Gauge]:
    gauges: List[Gauge] = []
    h = 4.0e-4
    for delta_phi in [0.0, 0.13, 0.37, 0.91, 1.41]:
        for delta_chi in [-8.0e-4, 0.0, 8.0e-4]:
            for delta_tau in [-8.0e-4, 0.0, 8.0e-4]:
                gauges.append(
                    Gauge(
                        f"phase_center_phi_{delta_phi:.2f}_chi_{delta_chi:+.4f}_tau_{delta_tau:+.4f}",
                        delta_chi=delta_chi,
                        delta_tau=delta_tau,
                        delta_phi=delta_phi,
                        h_chi=h,
                        h_tau=h,
                        nh_chi_c=params.Nh_chi_C,
                        nh_tau_c=params.Nh_tau_C,
                        c_gain=params.A_C,
                    )
                )
    return gauges


def make_width_gain_gauges(params: Params) -> List[Gauge]:
    gauges: List[Gauge] = []
    for nh in [599, 799, 999, 1199]:
        for h in [2.5e-4, 5.0e-4, 7.5e-4]:
            for gain in [600.0, 1000.0, 1600.0]:
                gauges.append(
                    Gauge(
                        f"width_gain_nh_{nh}_h_{h:.4f}_gain_{gain:.0f}",
                        delta_phi=0.23,
                        h_chi=h,
                        h_tau=h,
                        nh_chi_c=nh,
                        nh_tau_c=nh,
                        c_gain=gain,
                    )
                )
    return gauges


def make_near_lobe_gauges(params: Params) -> List[Gauge]:
    gauges: List[Gauge] = []
    h = 3.5e-4
    offsets = [-1.45e-3, -7.5e-4, 0.0, 7.5e-4, 1.45e-3]
    for delta_chi in offsets:
        for delta_tau in offsets:
            gauges.append(
                Gauge(
                    f"near_lobe_chi_{delta_chi:+.5f}_tau_{delta_tau:+.5f}",
                    delta_chi=delta_chi,
                    delta_tau=delta_tau,
                    delta_phi=0.57,
                    h_chi=h,
                    h_tau=h,
                    nh_chi_c=params.Nh_chi_C,
                    nh_tau_c=params.Nh_tau_C,
                    c_gain=params.A_C,
                )
            )
    return gauges


def make_mixed_grid_gauges(params: Params) -> List[Gauge]:
    gauges: List[Gauge] = []
    for index, delta_phi in enumerate([0.0, 0.19, 0.43, 0.83]):
        for delta_chi in [-6.0e-4, 6.0e-4]:
            for delta_tau in [-6.0e-4, 6.0e-4]:
                nh = [699, 899, 999, 1099][index]
                gain = [750.0, 950.0, 1250.0, 1500.0][index]
                h = [3.0e-4, 4.5e-4, 5.5e-4, 6.5e-4][index]
                gauges.append(
                    Gauge(
                        f"mixed_{index}_chi_{delta_chi:+.4f}_tau_{delta_tau:+.4f}",
                        delta_chi=delta_chi,
                        delta_tau=delta_tau,
                        delta_phi=delta_phi,
                        h_chi=h,
                        h_tau=h,
                        nh_chi_c=nh,
                        nh_tau_c=nh,
                        c_gain=gain,
                    )
                )
    return gauges


def gauge_cases(params: Params) -> Dict[str, List[Gauge]]:
    return {
        "baseline_default": prefixed_gauges("baseline_default", default_gauges(params)),
        "phase_center_grid": prefixed_gauges("phase_center_grid", make_phase_center_gauges(params)),
        "width_gain_grid": prefixed_gauges("width_gain_grid", make_width_gain_gauges(params)),
        "near_lobe_offset": prefixed_gauges("near_lobe_offset", make_near_lobe_gauges(params)),
        "mixed_readout_grid": prefixed_gauges("mixed_readout_grid", make_mixed_grid_gauges(params)),
    }


def bool_all(verdicts: Dict[str, Any], keys: Iterable[str]) -> bool:
    return all(bool(verdicts[key]) for key in keys)


def compute_case_result(
    params: Params,
    case_name: str,
    gauges: List[Gauge],
    stages: List[Dict[str, Any]],
    events: List[Dict[str, Any]],
) -> Dict[str, Any]:
    gauge_rows = readout_all(stages, gauges, params)
    summaries = summarize_stage_readouts(gauge_rows)
    verdicts = compute_single_collision_verdicts(params, stages, events, gauge_rows, summaries)
    multigauge_used = len(gauges) >= params.min_gauge_count
    strict_case_valid = bool_all(
        verdicts,
        [
            "baseline_collision_valid",
            "label_modes_preserved",
            "closure_preserved",
            "p_reconstructed_all_gauges",
            "E_reconstructed_all_gauges",
            "R_reconstructed_all_gauges",
            "p_reflection_valid",
            "E_preserved",
            "R_preserved",
            "R_gauge_stable",
            "t_R_separation_valid",
        ],
    )
    case_valid = bool(strict_case_valid and multigauge_used)
    case_summary = {
        "case": case_name,
        "gauge_count": len(gauges),
        "multigauge_used": multigauge_used,
        "p_max_abs_error": float(verdicts["p_max_abs_error"]),
        "E_max_abs_error": float(verdicts["E_max_abs_error"]),
        "R_max_abs_error": float(verdicts["R_max_abs_error"]),
        "R_max_gauge_std": float(verdicts["R_max_gauge_std"]),
        "separation_ratio_time": float(verdicts["separation_ratio_time"]),
        "closure_residual_abs": float(verdicts["closure_residual_abs"]),
        "p_reflection_error_A": float(verdicts["p_reflection_error_A"]),
        "p_reflection_error_B": float(verdicts["p_reflection_error_B"]),
        "E_preservation_error_A": float(verdicts["E_preservation_error_A"]),
        "E_preservation_error_B": float(verdicts["E_preservation_error_B"]),
        "R_preservation_error_A": float(verdicts["R_preservation_error_A"]),
        "R_preservation_error_B": float(verdicts["R_preservation_error_B"]),
        "case_valid": case_valid,
    }
    for row in gauge_rows:
        row["case"] = case_name
    for row in summaries:
        row["case"] = case_name
    return {
        "case": case_name,
        "gauges": [asdict(gauge) for gauge in gauges],
        "case_summary": case_summary,
        "stage_summaries": summaries,
        "gauge_rows": gauge_rows,
        "verdicts": verdicts,
    }


def aggregate_verdict(case_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    case_summaries = [entry["case_summary"] for entry in case_results]
    return {
        "case_count": len(case_summaries),
        "total_gauge_count": int(sum(int(row["gauge_count"]) for row in case_summaries)),
        "all_cases_valid": all(bool(row["case_valid"]) for row in case_summaries),
        "max_p_abs_error_all_cases": float(max(float(row["p_max_abs_error"]) for row in case_summaries)),
        "max_E_abs_error_all_cases": float(max(float(row["E_max_abs_error"]) for row in case_summaries)),
        "max_R_abs_error_all_cases": float(max(float(row["R_max_abs_error"]) for row in case_summaries)),
        "max_R_gauge_std_all_cases": float(max(float(row["R_max_gauge_std"]) for row in case_summaries)),
        "max_separation_ratio_time_all_cases": float(
            max(float(row["separation_ratio_time"]) for row in case_summaries)
        ),
        "single_gauge_only_used": False,
        "robustness_sweep_valid": all(bool(row["case_valid"]) for row in case_summaries),
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(case_summaries: List[Dict[str, Any]]) -> None:
    cases = [str(row["case"]) for row in case_summaries]
    xs = np.arange(len(cases))
    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.22
    ax.bar(xs - width, [float(row["p_max_abs_error"]) for row in case_summaries], width, label="p error")
    ax.bar(xs, [float(row["E_max_abs_error"]) for row in case_summaries], width, label="E error")
    ax.bar(xs + width, [float(row["R_max_abs_error"]) for row in case_summaries], width, label="R error")
    ax.axhline(1.0e-9, color="black", linewidth=0.8, linestyle="--", label="1e-9")
    ax.set_yscale("symlog", linthresh=1e-20)
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("max absolute error")
    ax.set_title("multigauge readout robustness: p/E/R reconstruction")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best", ncol=4)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_errors_v2.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(xs - width / 2.0, [float(row["R_max_gauge_std"]) for row in case_summaries], width, label="R gauge std")
    ax.bar(
        xs + width / 2.0,
        [float(row["separation_ratio_time"]) for row in case_summaries],
        width,
        label="Var(R)/Var(t)",
    )
    ax.axhline(1.0e-10, color="black", linewidth=0.8, linestyle="--", label="1e-10")
    ax.set_yscale("symlog", linthresh=1e-32)
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("stability metric")
    ax.set_title("multigauge readout robustness: R stability and t/R separation")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best", ncol=3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_stability_v2.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Interference Readout Robustness Sweep v1",
        "",
        "## Purpose",
        "",
        "This sweep changes only the readout-gauge family while keeping the same one-collision ABC trajectory.",
        "It checks whether p-like, E-like, and R-like readouts are stable under reference phase, readout center, width, and gain changes.",
        "",
        "## Aggregate Verdict",
        "",
    ]
    for key, value in result["aggregate_verdict"].items():
        lines.append(f"- {key}: `{value}`")
    lines.extend(
        [
            "",
            "## Case Summary",
            "",
            "| case | gauges | p max err | E max err | R max err | R std | Var(R)/Var(t) | valid |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | {row['gauge_count']} | {row['p_max_abs_error']:.16e} | "
            f"{row['E_max_abs_error']:.16e} | {row['R_max_abs_error']:.16e} | "
            f"{row['R_max_gauge_std']:.16e} | {row['separation_ratio_time']:.16e} | "
            f"`{row['case_valid']}` |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_interference_readout_robustness_sweep_result_v2.json` |",
            "| case CSV | `abc_multigauge_interference_readout_robustness_sweep_cases_v2.csv` |",
            "| gauge CSV | `abc_multigauge_interference_readout_robustness_sweep_gauge_rows_v2.csv` |",
            "| stage summary CSV | `abc_multigauge_interference_readout_robustness_sweep_stage_summary_v2.csv` |",
            "| error plot | `abc_multigauge_interference_readout_robustness_sweep_errors_v2.png` |",
            "| stability plot | `abc_multigauge_interference_readout_robustness_sweep_stability_v2.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_report_v2.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    stages, events = simulate_single_collision(params)
    cases = gauge_cases(params)
    case_results = [
        compute_case_result(params, case_name, gauges, stages, events)
        for case_name, gauges in cases.items()
    ]
    case_summaries = [entry["case_summary"] for entry in case_results]
    aggregate = aggregate_verdict(case_results)
    return {
        "experiment": "abc_multigauge_interference_readout_robustness_sweep_v2",
        "parameters": asdict(params),
        "stages": stages,
        "events": events,
        "case_summaries": case_summaries,
        "aggregate_verdict": aggregate,
        "case_results": [
            {
                "case": entry["case"],
                "gauges": entry["gauges"],
                "case_summary": entry["case_summary"],
                "verdicts": entry["verdicts"],
            }
            for entry in case_results
        ],
        "note": (
            "The trajectory is fixed. Only multigauge readout settings are swept. "
            "A case is valid only when the full gauge family reconstructs p/E/R and preserves t/R separation."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    case_results_full = []
    params = Params(**result["parameters"])
    stages = result["stages"]
    events = result["events"]
    for case_name, gauges in gauge_cases(params).items():
        case_results_full.append(compute_case_result(params, case_name, gauges, stages, events))

    gauge_rows: List[Dict[str, Any]] = []
    stage_summaries: List[Dict[str, Any]] = []
    for entry in case_results_full:
        gauge_rows.extend(entry["gauge_rows"])
        stage_summaries.extend(entry["stage_summaries"])

    (OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_result_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_cases_v2.csv", result["case_summaries"])
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_gauge_rows_v2.csv", gauge_rows)
    write_csv(
        OUT_DIR / "abc_multigauge_interference_readout_robustness_sweep_stage_summary_v2.csv",
        stage_summaries,
    )
    make_plots(result["case_summaries"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
