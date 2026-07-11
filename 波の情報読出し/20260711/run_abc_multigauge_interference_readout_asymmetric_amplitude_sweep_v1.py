from __future__ import annotations

import csv
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_multigauge_interference_readout_v1 import (
    Gauge,
    Params as BaseParams,
    compute_verdicts as compute_single_collision_verdicts,
    default_gauges,
    readout_all,
    simulate_single_collision,
    summarize_stage_readouts,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1"
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
    pass


AMPLITUDE_CASES: List[Tuple[float, float]] = [
    (1.0, 1.0),
    (1.0, 1.10),
    (1.0, 1.25),
    (1.0, 1.50),
    (1.0, 2.00),
    (1.0, 3.00),
    (1.50, 1.0),
    (2.00, 1.0),
]


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def summary_lookup(summaries: List[Dict[str, Any]]) -> Dict[tuple[str, str], Dict[str, Any]]:
    return {(str(row["stage"]), str(row["particle"])): row for row in summaries}


def stage_quantities(case_name: str, summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    lookup = summary_lookup(summaries)
    stages = ["initial", "collision_map", "final"]
    rows: List[Dict[str, Any]] = []
    for stage in stages:
        a = lookup[(stage, "A")]
        b = lookup[(stage, "B")]
        r_a = float(a["R_mean"])
        r_b = float(b["R_mean"])
        p_a = float(a["p_mean"])
        p_b = float(b["p_mean"])
        e_a = float(a["E_mean"])
        e_b = float(b["E_mean"])
        rows.append(
            {
                "case": case_name,
                "stage": stage,
                "R_A": r_a,
                "R_B": r_b,
                "p_A": p_a,
                "p_B": p_b,
                "E_A": e_a,
                "E_B": e_b,
                "P_phase_total": p_a + p_b,
                "P_R_weighted_total": r_a * p_a + r_b * p_b,
                "E_phase_total": e_a + e_b,
                "E_R_weighted_total": r_a * e_a + r_b * e_b,
                "R_total": r_a + r_b,
            }
        )
    return rows


def within_particle_tr_separation(summaries: List[Dict[str, Any]]) -> float:
    ratios: List[float] = []
    for particle in ["A", "B"]:
        selected = [row for row in summaries if row["particle"] == particle]
        r_values = np.array([float(row["R_mean"]) for row in selected])
        t_values = np.array([float(row["t_mean"]) for row in selected])
        var_r = float(np.var(r_values))
        var_t = float(np.var(t_values))
        ratios.append(float(var_r / var_t) if var_t > 0.0 else float("inf"))
    return float(max(ratios))


def case_from_quantities(
    params: Params,
    case_name: str,
    quantities: List[Dict[str, Any]],
    verdicts: Dict[str, Any],
    summaries: List[Dict[str, Any]],
) -> Dict[str, Any]:
    by_stage = {str(row["stage"]): row for row in quantities}
    initial = by_stage["initial"]
    after = by_stage["collision_map"]
    final = by_stage["final"]
    weighted_p_initial = float(initial["P_R_weighted_total"])
    weighted_p_after = float(after["P_R_weighted_total"])
    weighted_p_final = float(final["P_R_weighted_total"])
    phase_p_initial = float(initial["P_phase_total"])
    phase_p_after = float(after["P_phase_total"])
    weighted_e_initial = float(initial["E_R_weighted_total"])
    weighted_e_final = float(final["E_R_weighted_total"])
    r_total_initial = float(initial["R_total"])
    r_total_final = float(final["R_total"])
    weighted_p_collision_error = abs(weighted_p_after - weighted_p_initial)
    weighted_p_final_error = abs(weighted_p_final - weighted_p_initial)
    phase_p_collision_error = abs(phase_p_after - phase_p_initial)
    weighted_e_error = abs(weighted_e_final - weighted_e_initial)
    r_total_error = abs(r_total_final - r_total_initial)
    asymmetric = abs(params.A_A - params.A_B) > 1.0e-12
    detects_asymmetry = bool(asymmetric and weighted_p_collision_error > params.conservation_tol)
    equal_case_preserved = bool((not asymmetric) and weighted_p_collision_error <= params.conservation_tol)
    tr_ratio_within_particle = within_particle_tr_separation(summaries)
    individual_multigauge_valid = all(
        bool(verdicts[key])
        for key in [
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
        ]
    ) and tr_ratio_within_particle <= params.tr_separation_threshold
    return {
        "case": case_name,
        "A_A": params.A_A,
        "A_B": params.A_B,
        "R_A_expected": params.A_A**2,
        "R_B_expected": params.A_B**2,
        "R_ratio_B_over_A": (params.A_B**2) / (params.A_A**2),
        "individual_multigauge_valid": bool(individual_multigauge_valid),
        "p_max_abs_error": float(verdicts["p_max_abs_error"]),
        "E_max_abs_error": float(verdicts["E_max_abs_error"]),
        "R_max_abs_error": float(verdicts["R_max_abs_error"]),
        "R_max_gauge_std": float(verdicts["R_max_gauge_std"]),
        "global_R_contrast_ratio_time": float(verdicts["separation_ratio_time"]),
        "within_particle_separation_ratio_time": tr_ratio_within_particle,
        "phase_p_total_initial": phase_p_initial,
        "phase_p_total_after": phase_p_after,
        "phase_p_collision_error": phase_p_collision_error,
        "weighted_p_total_initial": weighted_p_initial,
        "weighted_p_total_after": weighted_p_after,
        "weighted_p_total_final": weighted_p_final,
        "weighted_p_collision_error": weighted_p_collision_error,
        "weighted_p_final_error": weighted_p_final_error,
        "weighted_e_total_initial": weighted_e_initial,
        "weighted_e_total_final": weighted_e_final,
        "weighted_e_error": weighted_e_error,
        "R_total_initial": r_total_initial,
        "R_total_final": r_total_final,
        "R_total_error": r_total_error,
        "simple_q_flip_weighted_momentum_preserved": bool(
            weighted_p_collision_error <= params.conservation_tol
            and weighted_p_final_error <= params.conservation_tol
        ),
        "weighted_energy_preserved": bool(weighted_e_error <= params.conservation_tol),
        "R_total_preserved": bool(r_total_error <= params.conservation_tol),
        "asymmetric_case": asymmetric,
        "detects_asymmetry_boundary": detects_asymmetry or equal_case_preserved,
    }


def run_case(a_amp: float, b_amp: float) -> Dict[str, Any]:
    params = Params(A_A=a_amp, A_B=b_amp)
    gauges = default_gauges(params)
    stages, events = simulate_single_collision(params)
    gauge_rows = readout_all(stages, gauges, params)
    summaries = summarize_stage_readouts(gauge_rows)
    verdicts = compute_single_collision_verdicts(params, stages, events, gauge_rows, summaries)
    case_name = f"A_{a_amp:.2f}_B_{b_amp:.2f}"
    quantities = stage_quantities(case_name, summaries)
    case_summary = case_from_quantities(params, case_name, quantities, verdicts, summaries)
    for row in gauge_rows:
        row["case"] = case_name
        row["A_A"] = a_amp
        row["A_B"] = b_amp
    return {
        "case": case_name,
        "parameters": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "events": events,
        "stages": stages,
        "stage_summaries": summaries,
        "stage_quantities": quantities,
        "gauge_rows": gauge_rows,
        "verdicts": verdicts,
        "case_summary": case_summary,
    }


def aggregate(case_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    summaries = [entry["case_summary"] for entry in case_results]
    asymmetric_rows = [row for row in summaries if row["asymmetric_case"]]
    equal_rows = [row for row in summaries if not row["asymmetric_case"]]
    return {
        "case_count": len(summaries),
        "asymmetric_case_count": len(asymmetric_rows),
        "individual_multigauge_valid_all_cases": all(bool(row["individual_multigauge_valid"]) for row in summaries),
        "weighted_energy_preserved_all_cases": all(bool(row["weighted_energy_preserved"]) for row in summaries),
        "R_total_preserved_all_cases": all(bool(row["R_total_preserved"]) for row in summaries),
        "equal_case_weighted_momentum_preserved": all(
            bool(row["simple_q_flip_weighted_momentum_preserved"]) for row in equal_rows
        ),
        "asymmetric_cases_detect_weighted_momentum_failure": all(
            not bool(row["simple_q_flip_weighted_momentum_preserved"]) for row in asymmetric_rows
        ),
        "max_p_abs_error": float(max(float(row["p_max_abs_error"]) for row in summaries)),
        "max_E_abs_error": float(max(float(row["E_max_abs_error"]) for row in summaries)),
        "max_R_abs_error": float(max(float(row["R_max_abs_error"]) for row in summaries)),
        "max_R_gauge_std": float(max(float(row["R_max_gauge_std"]) for row in summaries)),
        "max_global_R_contrast_ratio_time": float(
            max(float(row["global_R_contrast_ratio_time"]) for row in summaries)
        ),
        "max_within_particle_separation_ratio_time": float(
            max(float(row["within_particle_separation_ratio_time"]) for row in summaries)
        ),
        "max_weighted_p_collision_error": float(max(float(row["weighted_p_collision_error"]) for row in summaries)),
        "asymmetric_amplitude_diagnostic_valid": all(
            [
                all(bool(row["individual_multigauge_valid"]) for row in summaries),
                all(bool(row["weighted_energy_preserved"]) for row in summaries),
                all(bool(row["R_total_preserved"]) for row in summaries),
                all(bool(row["simple_q_flip_weighted_momentum_preserved"]) for row in equal_rows),
                all(not bool(row["simple_q_flip_weighted_momentum_preserved"]) for row in asymmetric_rows),
            ]
        ),
    }


def make_plots(case_summaries: List[Dict[str, Any]]) -> None:
    ratios = np.array([float(row["R_ratio_B_over_A"]) for row in case_summaries])
    floor = 1.0e-18
    p_errors = np.maximum(
        np.array([float(row["weighted_p_collision_error"]) for row in case_summaries]), floor
    )
    r_errors = np.maximum(np.array([float(row["R_total_error"]) for row in case_summaries]), floor)
    e_errors = np.maximum(np.array([float(row["weighted_e_error"]) for row in case_summaries]), floor)

    order = np.argsort(ratios)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ratios[order], p_errors[order], marker="o", label="R-weighted P conservation error")
    ax.plot(ratios[order], e_errors[order], marker="o", label="R-weighted E conservation error")
    ax.plot(ratios[order], r_errors[order], marker="o", label="R total conservation error")
    ax.axhline(1.0e-9, color="black", linestyle="--", linewidth=0.8, label="1e-9")
    ax.set_yscale("log")
    ax.set_xlabel("R_B / R_A")
    ax.set_ylabel("absolute conservation error")
    ax.set_title("asymmetric amplitude diagnostic under simple q-flip")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_conservation_v1.png", dpi=160)
    plt.close(fig)

    cases = [str(row["case"]) for row in case_summaries]
    xs = np.arange(len(cases))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(xs, [float(row["R_ratio_B_over_A"]) for row in case_summaries])
    ax.axhline(1.0, color="black", linewidth=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("R_B / R_A")
    ax.set_title("mass-like R asymmetry cases")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Interference Readout Asymmetric Amplitude Sweep v1",
        "",
        "## Purpose",
        "",
        "This diagnostic keeps the simple q-flip collision map and changes only A/B representative amplitudes.",
        "It checks whether individual p/E/R readouts remain reconstructable, and whether R-weighted total momentum remains compatible with the simple equal-mass reflection map.",
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
            "| case | R_B/R_A | individual valid | within-particle Var(R)/Var(t) | weighted P error | weighted E error | R total error | simple q flip P preserved |",
            "|---|---:|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | {row['R_ratio_B_over_A']:.16e} | `{row['individual_multigauge_valid']}` | "
            f"{row['within_particle_separation_ratio_time']:.16e} | {row['weighted_p_collision_error']:.16e} | "
            f"{row['weighted_e_error']:.16e} | "
            f"{row['R_total_error']:.16e} | `{row['simple_q_flip_weighted_momentum_preserved']}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The individual multigauge p/E/R readout remains valid in every amplitude case.",
            "However, the R-weighted total momentum is conserved by the simple q-flip map only in the equal-amplitude case.",
            "The global R variance is not used as a failure condition here, because unequal A/B amplitudes are the diagnostic signal itself.",
            "Thus the experiment separates readout validity from collision-law compatibility.",
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1.json` |",
            "| case CSV | `abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.csv` |",
            "| stage quantity CSV | `abc_multigauge_interference_readout_asymmetric_amplitude_stage_quantities_v1.csv` |",
            "| gauge CSV | `abc_multigauge_interference_readout_asymmetric_amplitude_gauge_rows_v1.csv` |",
            "| conservation plot | `abc_multigauge_interference_readout_asymmetric_amplitude_conservation_v1.png` |",
            "| cases plot | `abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_results = [run_case(a_amp, b_amp) for a_amp, b_amp in AMPLITUDE_CASES]
    case_summaries = [entry["case_summary"] for entry in case_results]
    return {
        "experiment": "abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1",
        "amplitude_cases": [{"A_A": a, "A_B": b} for a, b in AMPLITUDE_CASES],
        "case_summaries": case_summaries,
        "aggregate_verdict": aggregate(case_results),
        "case_results": [
            {
                "case": entry["case"],
                "parameters": entry["parameters"],
                "events": entry["events"],
                "verdicts": entry["verdicts"],
                "case_summary": entry["case_summary"],
            }
            for entry in case_results
        ],
        "note": (
            "This is a diagnostic sweep. It does not require the simple q-flip map to be valid for unequal R. "
            "Instead it checks whether the multigauge readout detects the boundary of that equal-amplitude map."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    full_case_results = [run_case(a_amp, b_amp) for a_amp, b_amp in AMPLITUDE_CASES]
    gauge_rows: List[Dict[str, Any]] = []
    stage_quantities_rows: List[Dict[str, Any]] = []
    for entry in full_case_results:
        gauge_rows.extend(entry["gauge_rows"])
        stage_quantities_rows.extend(entry["stage_quantities"])

    (OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_cases_v1.csv", result["case_summaries"])
    write_csv(
        OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_stage_quantities_v1.csv",
        stage_quantities_rows,
    )
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_asymmetric_amplitude_gauge_rows_v1.csv", gauge_rows)
    make_plots(result["case_summaries"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
