from __future__ import annotations

import csv
import json
import os
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_abc_multigauge_interference_readout_v2 import default_gauges, summarize_stage_readouts
from run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v2 import (
    Params,
    case_summary,
    readout_all_float,
    simulate,
    stage_quantities,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXTREME_R_CASES: List[Dict[str, float]] = [
    {"A_A": 1.0, "A_B": 0.125, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 1.0, "A_B": 0.25, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 1.0, "A_B": 0.5, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 1.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 1.0, "A_B": 4.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 1.0, "A_B": 8.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 0.25, "A_B": 1.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 0.5, "A_B": 1.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 2.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 4.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -0.5},
    {"A_A": 8.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -0.5},
]


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    return float(max(abs(float(row[key])) for row in rows))


def run_case(index: int, settings: Dict[str, float]) -> Dict[str, Any]:
    params = Params(**settings)
    gauges = default_gauges(params)
    stages, events = simulate(params)
    gauge_rows = readout_all_float(stages, gauges, params)
    summaries = summarize_stage_readouts(gauge_rows)
    case = (
        f"c{index:02d}_A{params.A_A:.3f}_B{params.A_B:.3f}_"
        f"u{params.q_A0:.2f}_v{params.q_B0:.2f}"
    )
    for row in gauge_rows:
        row["case"] = case
        row["A_A"] = params.A_A
        row["A_B"] = params.A_B
        row["q_A0"] = params.q_A0
        row["q_B0"] = params.q_B0
    quantities = stage_quantities(case, summaries)
    summary = case_summary(params, case, gauge_rows, summaries, quantities, events)
    dynamic_range = max(params.A_A**2, params.A_B**2) / min(params.A_A**2, params.A_B**2)
    summary["R_dynamic_range"] = dynamic_range
    summary["boundary_case_valid"] = bool(
        summary["collision_cell_reached"]
        and summary["individual_readout_valid"]
        and summary["P_R_preserved"]
        and summary["K_R_phase_preserved"]
        and summary["relative_gradient_flipped"]
        and summary["E_tau_R_preserved"]
        and summary["R_total_preserved"]
    )
    return {
        "case": case,
        "parameters": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "stages": stages,
        "events": events,
        "gauge_rows": gauge_rows,
        "stage_summaries": summaries,
        "stage_quantities": quantities,
        "case_summary": summary,
    }


def aggregate(case_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    rows = [entry["case_summary"] for entry in case_results]
    return {
        "case_count": len(rows),
        "all_cases_valid": all(bool(row["boundary_case_valid"]) for row in rows),
        "collision_reached_all_cases": all(bool(row["collision_cell_reached"]) for row in rows),
        "individual_readout_valid_all_cases": all(bool(row["individual_readout_valid"]) for row in rows),
        "P_R_preserved_all_cases": all(bool(row["P_R_preserved"]) for row in rows),
        "K_R_phase_preserved_all_cases": all(bool(row["K_R_phase_preserved"]) for row in rows),
        "relative_gradient_flipped_all_cases": all(bool(row["relative_gradient_flipped"]) for row in rows),
        "E_tau_R_preserved_all_cases": all(bool(row["E_tau_R_preserved"]) for row in rows),
        "R_total_preserved_all_cases": all(bool(row["R_total_preserved"]) for row in rows),
        "max_R_dynamic_range": float(max(float(row["R_dynamic_range"]) for row in rows)),
        "min_R_ratio_B_over_A": float(min(float(row["R_ratio_B_over_A"]) for row in rows)),
        "max_R_ratio_B_over_A": float(max(float(row["R_ratio_B_over_A"]) for row in rows)),
        "max_p_abs_error": float(max(float(row["p_max_abs_error"]) for row in rows)),
        "max_E_abs_error": float(max(float(row["E_max_abs_error"]) for row in rows)),
        "max_R_abs_error": float(max(float(row["R_max_abs_error"]) for row in rows)),
        "max_R_gauge_std": float(max(float(row["R_max_gauge_std"]) for row in rows)),
        "max_within_particle_separation_ratio_time": float(
            max(float(row["within_particle_separation_ratio_time"]) for row in rows)
        ),
        "max_P_R_conservation_error": float(max(float(row["P_R_conservation_error"]) for row in rows)),
        "max_K_R_phase_conservation_error": float(
            max(float(row["K_R_phase_conservation_error"]) for row in rows)
        ),
        "max_relative_flip_error": float(max(float(row["relative_flip_error"]) for row in rows)),
        "max_E_tau_R_conservation_error": float(max(float(row["E_tau_R_conservation_error"]) for row in rows)),
        "max_R_total_conservation_error": float(max(float(row["R_total_conservation_error"]) for row in rows)),
        "single_gauge_only_used": False,
        "extreme_R_sweep_valid": all(bool(row["boundary_case_valid"]) for row in rows),
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
    xs = np.array([float(row["R_ratio_B_over_A"]) for row in case_summaries])
    p_errors = np.array([float(row["P_R_conservation_error"]) for row in case_summaries])
    k_errors = np.array([float(row["K_R_phase_conservation_error"]) for row in case_summaries])
    r_errors = np.array([float(row["R_max_abs_error"]) for row in case_summaries])
    floor = 1.0e-18
    order = np.argsort(xs)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(xs[order], np.maximum(p_errors[order], floor), marker="o", label="R*p error")
    ax.plot(xs[order], np.maximum(k_errors[order], floor), marker="o", label="R*p^2 error")
    ax.plot(xs[order], np.maximum(r_errors[order], floor), marker="o", label="R readout error")
    ax.axhline(1.0e-9, color="black", linestyle="--", linewidth=0.8, label="1e-9")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("R_B / R_A")
    ax.set_ylabel("absolute error")
    ax.set_title("generalized elastic collision extreme R sweep")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_errors_v2.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(xs[order], [float(case_summaries[i]["q_A_after"]) for i in order], marker="o", label="p_A after")
    ax.plot(xs[order], [float(case_summaries[i]["q_B_after"]) for i in order], marker="o", label="p_B after")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("R_B / R_A")
    ax.set_ylabel("post-collision phase gradient")
    ax.set_title("post-collision readout gradients across R ratios")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_outputs_v2.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Generalized Elastic Collision Extreme R Sweep v1",
        "",
        "## Purpose",
        "",
        "This experiment sweeps the readout R ratio over extreme asymmetric amplitude conditions.",
        "It checks whether the generalized R-weighted elastic map and multigauge readout survive large R contrast.",
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
            "| case | R_B/R_A | dynamic range | q after A/B | R*p err | R*p^2 err | valid |",
            "|---|---:|---:|---|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | {row['R_ratio_B_over_A']:.16e} | {row['R_dynamic_range']:.16e} | "
            f"{row['q_A_after']:.8g} / {row['q_B_after']:.8g} | "
            f"{row['P_R_conservation_error']:.16e} | "
            f"{row['K_R_phase_conservation_error']:.16e} | "
            f"`{row['boundary_case_valid']}` |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2.json` |",
            "| case CSV | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_cases_v2.csv` |",
            "| stage quantity CSV | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_stage_quantities_v2.csv` |",
            "| gauge CSV | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_gauge_rows_v2.csv` |",
            "| error plot | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_errors_v2.png` |",
            "| output plot | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_outputs_v2.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_report_v2.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_results = [run_case(index + 1, settings) for index, settings in enumerate(EXTREME_R_CASES)]
    case_summaries = [entry["case_summary"] for entry in case_results]
    return {
        "experiment": "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v2",
        "extreme_R_cases": EXTREME_R_CASES,
        "case_summaries": case_summaries,
        "aggregate_verdict": aggregate(case_results),
        "case_results": [
            {
                "case": entry["case"],
                "parameters": entry["parameters"],
                "events": entry["events"],
                "case_summary": entry["case_summary"],
            }
            for entry in case_results
        ],
        "note": "R is read as amplitude squared. This sweep tests large readout-R contrast, not standard mass units.",
    }


def write_outputs(result: Dict[str, Any]) -> None:
    full_case_results = [run_case(index + 1, settings) for index, settings in enumerate(EXTREME_R_CASES)]
    gauge_rows: List[Dict[str, Any]] = []
    stage_quantities_rows: List[Dict[str, Any]] = []
    for entry in full_case_results:
        gauge_rows.extend(entry["gauge_rows"])
        stage_quantities_rows.extend(entry["stage_quantities"])
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_cases_v2.csv", result["case_summaries"])
    write_csv(
        OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_stage_quantities_v2.csv",
        stage_quantities_rows,
    )
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_gauge_rows_v2.csv", gauge_rows)
    make_plots(result["case_summaries"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
