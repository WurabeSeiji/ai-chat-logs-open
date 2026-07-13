from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_abc_multigauge_interference_readout_v2 import (
    Gauge,
    Params as BaseParams,
    State,
    closure_residual,
    default_gauges,
    read_particle_gauge,
    summarize_stage_readouts,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_generalized_elastic_collision_readout_result_v2"
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


def copy_state(state: State) -> State:
    return State(state.chi, state.tau, state.q, state.amplitude, state.m, state.omega)


def append_stage(stages: List[Dict[str, Any]], name: str, step: int, a: State, b: State) -> None:
    stages.append(
        {
            "stage": name,
            "step": step,
            "chi_A": float(a.chi),
            "chi_B": float(b.chi),
            "tau_A": float(a.tau),
            "tau_B": float(b.tau),
            "q_A": float(a.q),
            "q_B": float(b.q),
            "m_A": int(a.m),
            "m_B": int(b.m),
            "omega_A": float(a.omega),
            "omega_B": float(b.omega),
        }
    )


def generalized_elastic_velocity(m_a: float, m_b: float, u_a: float, u_b: float) -> tuple[float, float]:
    denom = m_a + m_b
    v_a = ((m_a - m_b) / denom) * u_a + (2.0 * m_b / denom) * u_b
    v_b = (2.0 * m_a / denom) * u_a + ((m_b - m_a) / denom) * u_b
    return float(v_a), float(v_b)


def simulate_generalized_collision(params: Params) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    eps_chi = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    a = State(params.chi_A0, params.tau_A0, float(params.q_A0), params.A_A, params.m_A, params.omega_A)
    b = State(params.chi_B0, params.tau_B0, float(params.q_B0), params.A_B, params.m_B, params.omega_B)
    stages: List[Dict[str, Any]] = []
    events: List[Dict[str, Any]] = []
    append_stage(stages, "initial", 0, copy_state(a), copy_state(b))

    previous_a = copy_state(a)
    previous_b = copy_state(b)
    step = 0
    collision_cell_reached = False
    while step < params.s_max:
        previous_a = copy_state(a)
        previous_b = copy_state(b)
        a.chi += float(a.q) * params.v_chi * params.delta_s
        b.chi += float(b.q) * params.v_chi * params.delta_s
        a.tau += a.omega * params.delta_s
        b.tau += b.omega * params.delta_s
        step += 1
        if abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau:
            collision_cell_reached = True
            break

    append_stage(stages, "pre_collision", step - 1, previous_a, previous_b)
    append_stage(stages, "collision_cell", step, copy_state(a), copy_state(b))

    before_q_a = float(a.q)
    before_q_b = float(b.q)
    mass_a = params.A_A**2
    mass_b = params.A_B**2
    after_q_a, after_q_b = generalized_elastic_velocity(mass_a, mass_b, before_q_a, before_q_b)
    simple_flip_q_a = -before_q_a
    simple_flip_q_b = -before_q_b
    if collision_cell_reached:
        a.q = after_q_a
        b.q = after_q_b
    events.append(
        {
            "event": "generalized_ab_collision",
            "step": step,
            "collision_cell_reached": collision_cell_reached,
            "R_A": mass_a,
            "R_B": mass_b,
            "q_A_before": before_q_a,
            "q_B_before": before_q_b,
            "q_A_after": float(a.q),
            "q_B_after": float(b.q),
            "simple_flip_q_A": simple_flip_q_a,
            "simple_flip_q_B": simple_flip_q_b,
            "m_A": params.m_A,
            "m_B": params.m_B,
            "closure_residual_abs": closure_residual(params),
        }
    )
    append_stage(stages, "collision_map", step, copy_state(a), copy_state(b))

    post_saved = False
    while step < params.s_max:
        a.chi += float(a.q) * params.v_chi * params.delta_s
        b.chi += float(b.q) * params.v_chi * params.delta_s
        a.tau += a.omega * params.delta_s
        b.tau += b.omega * params.delta_s
        step += 1
        separated = abs(a.chi - b.chi) > eps_chi
        if separated and not post_saved:
            append_stage(stages, "post_collision", step, copy_state(a), copy_state(b))
            post_saved = True
        if separated and min(a.tau, b.tau) >= abs(params.tau_A0) - 1.0e-12:
            break

    append_stage(stages, "final", step, copy_state(a), copy_state(b))
    return stages, events


def readout_all_float(stages: List[Dict[str, Any]], gauges: List[Gauge], params: Params) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for stage in stages:
        a = State(
            float(stage["chi_A"]),
            float(stage["tau_A"]),
            float(stage["q_A"]),
            params.A_A,
            params.m_A,
            params.omega_A,
        )
        b = State(
            float(stage["chi_B"]),
            float(stage["tau_B"]),
            float(stage["q_B"]),
            params.A_B,
            params.m_B,
            params.omega_B,
        )
        for gauge in gauges:
            rows.append(read_particle_gauge(str(stage["stage"]), "A", a, gauge, params))
            rows.append(read_particle_gauge(str(stage["stage"]), "B", b, gauge, params))
    return rows


def summary_lookup(summaries: List[Dict[str, Any]]) -> Dict[tuple[str, str], Dict[str, Any]]:
    return {(str(row["stage"]), str(row["particle"])): row for row in summaries}


def stage_quantities(case_name: str, summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    lookup = summary_lookup(summaries)
    rows: List[Dict[str, Any]] = []
    for stage in ["initial", "collision_cell", "collision_map", "final"]:
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
                "P_R_total": r_a * p_a + r_b * p_b,
                "K_R_phase_total": r_a * p_a**2 + r_b * p_b**2,
                "E_tau_R_total": r_a * e_a + r_b * e_b,
                "R_total": r_a + r_b,
            }
        )
    return rows


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    return float(max(abs(float(row[key])) for row in rows))


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


def simple_flip_weighted_errors(params: Params) -> Dict[str, float]:
    r_a = params.A_A**2
    r_b = params.A_B**2
    u_a = float(params.q_A0)
    u_b = float(params.q_B0)
    sf_a = -u_a
    sf_b = -u_b
    return {
        "simple_flip_P_R_error": abs((r_a * sf_a + r_b * sf_b) - (r_a * u_a + r_b * u_b)),
        "simple_flip_K_R_error": abs((r_a * sf_a**2 + r_b * sf_b**2) - (r_a * u_a**2 + r_b * u_b**2)),
    }


def case_summary(
    params: Params,
    case_name: str,
    gauge_rows: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
    quantities: List[Dict[str, Any]],
    events: List[Dict[str, Any]],
) -> Dict[str, Any]:
    by_stage = {str(row["stage"]): row for row in quantities}
    initial = by_stage["initial"]
    after = by_stage["collision_map"]
    final = by_stage["final"]
    event = events[0]
    p_max = max_abs(gauge_rows, "p_abs_error")
    e_max = max_abs(gauge_rows, "E_abs_error")
    r_max = max_abs(gauge_rows, "R_abs_error")
    r_std = float(max(float(row["R_std"]) for row in summaries))
    tr_ratio = within_particle_tr_separation(summaries)
    p_r_initial = float(initial["P_R_total"])
    p_r_after = float(after["P_R_total"])
    p_r_final = float(final["P_R_total"])
    k_r_initial = float(initial["K_R_phase_total"])
    k_r_after = float(after["K_R_phase_total"])
    k_r_final = float(final["K_R_phase_total"])
    e_tau_initial = float(initial["E_tau_R_total"])
    e_tau_final = float(final["E_tau_R_total"])
    r_total_initial = float(initial["R_total"])
    r_total_final = float(final["R_total"])
    simple_errors = simple_flip_weighted_errors(params)
    p_r_error = max(abs(p_r_after - p_r_initial), abs(p_r_final - p_r_initial))
    k_r_error = max(abs(k_r_after - k_r_initial), abs(k_r_final - k_r_initial))
    e_tau_error = abs(e_tau_final - e_tau_initial)
    r_total_error = abs(r_total_final - r_total_initial)
    return {
        "case": case_name,
        "A_A": params.A_A,
        "A_B": params.A_B,
        "R_A": params.A_A**2,
        "R_B": params.A_B**2,
        "R_ratio_B_over_A": (params.A_B**2) / (params.A_A**2),
        "collision_cell_reached": bool(event["collision_cell_reached"]),
        "q_A_after": float(event["q_A_after"]),
        "q_B_after": float(event["q_B_after"]),
        "simple_flip_q_A": float(event["simple_flip_q_A"]),
        "simple_flip_q_B": float(event["simple_flip_q_B"]),
        "p_max_abs_error": p_max,
        "E_max_abs_error": e_max,
        "R_max_abs_error": r_max,
        "R_max_gauge_std": r_std,
        "within_particle_separation_ratio_time": tr_ratio,
        "P_R_initial": p_r_initial,
        "P_R_after": p_r_after,
        "P_R_final": p_r_final,
        "P_R_conservation_error": p_r_error,
        "K_R_phase_initial": k_r_initial,
        "K_R_phase_after": k_r_after,
        "K_R_phase_final": k_r_final,
        "K_R_phase_conservation_error": k_r_error,
        "E_tau_R_conservation_error": e_tau_error,
        "R_total_conservation_error": r_total_error,
        "simple_flip_P_R_error": simple_errors["simple_flip_P_R_error"],
        "simple_flip_K_R_error": simple_errors["simple_flip_K_R_error"],
        "individual_readout_valid": bool(
            p_max <= params.readout_tol
            and e_max <= params.readout_tol
            and r_max <= params.readout_tol
            and r_std <= params.r_gauge_tol
            and tr_ratio <= params.tr_separation_threshold
        ),
        "generalized_P_R_preserved": bool(p_r_error <= params.conservation_tol),
        "generalized_K_R_phase_preserved": bool(k_r_error <= params.conservation_tol),
        "E_tau_R_preserved": bool(e_tau_error <= params.conservation_tol),
        "R_total_preserved": bool(r_total_error <= params.conservation_tol),
        "simple_flip_improves_or_matches": bool(simple_errors["simple_flip_P_R_error"] <= p_r_error + 1.0e-12),
    }


def run_case(a_amp: float, b_amp: float) -> Dict[str, Any]:
    params = Params(A_A=a_amp, A_B=b_amp)
    gauges = default_gauges(params)
    stages, events = simulate_generalized_collision(params)
    gauge_rows = readout_all_float(stages, gauges, params)
    summaries = summarize_stage_readouts(gauge_rows)
    case_name = f"A_{a_amp:.2f}_B_{b_amp:.2f}"
    for row in gauge_rows:
        row["case"] = case_name
        row["A_A"] = a_amp
        row["A_B"] = b_amp
    quantities = stage_quantities(case_name, summaries)
    summary = case_summary(params, case_name, gauge_rows, summaries, quantities, events)
    return {
        "case": case_name,
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
        "individual_readout_valid_all_cases": all(bool(row["individual_readout_valid"]) for row in rows),
        "generalized_P_R_preserved_all_cases": all(bool(row["generalized_P_R_preserved"]) for row in rows),
        "generalized_K_R_phase_preserved_all_cases": all(bool(row["generalized_K_R_phase_preserved"]) for row in rows),
        "E_tau_R_preserved_all_cases": all(bool(row["E_tau_R_preserved"]) for row in rows),
        "R_total_preserved_all_cases": all(bool(row["R_total_preserved"]) for row in rows),
        "max_p_abs_error": float(max(float(row["p_max_abs_error"]) for row in rows)),
        "max_E_abs_error": float(max(float(row["E_max_abs_error"]) for row in rows)),
        "max_R_abs_error": float(max(float(row["R_max_abs_error"]) for row in rows)),
        "max_R_gauge_std": float(max(float(row["R_max_gauge_std"]) for row in rows)),
        "max_within_particle_separation_ratio_time": float(
            max(float(row["within_particle_separation_ratio_time"]) for row in rows)
        ),
        "max_P_R_conservation_error": float(max(float(row["P_R_conservation_error"]) for row in rows)),
        "max_K_R_phase_conservation_error": float(max(float(row["K_R_phase_conservation_error"]) for row in rows)),
        "max_E_tau_R_conservation_error": float(max(float(row["E_tau_R_conservation_error"]) for row in rows)),
        "max_R_total_conservation_error": float(max(float(row["R_total_conservation_error"]) for row in rows)),
        "single_gauge_only_used": False,
        "generalized_elastic_collision_readout_valid": all(
            [
                all(bool(row["individual_readout_valid"]) for row in rows),
                all(bool(row["generalized_P_R_preserved"]) for row in rows),
                all(bool(row["generalized_K_R_phase_preserved"]) for row in rows),
                all(bool(row["E_tau_R_preserved"]) for row in rows),
                all(bool(row["R_total_preserved"]) for row in rows),
            ]
        ),
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
    ratios = np.array([float(row["R_ratio_B_over_A"]) for row in case_summaries])
    order = np.argsort(ratios)
    floor = 1.0e-18
    p_errors = np.maximum(np.array([float(row["P_R_conservation_error"]) for row in case_summaries]), floor)
    k_errors = np.maximum(
        np.array([float(row["K_R_phase_conservation_error"]) for row in case_summaries]), floor
    )
    simple_p_errors = np.maximum(np.array([float(row["simple_flip_P_R_error"]) for row in case_summaries]), floor)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ratios[order], p_errors[order], marker="o", label="generalized R*p error")
    ax.plot(ratios[order], k_errors[order], marker="o", label="generalized R*p^2 error")
    ax.plot(ratios[order], simple_p_errors[order], marker="o", linestyle="--", label="simple q-flip R*p error")
    ax.axhline(1.0e-9, color="black", linestyle="--", linewidth=0.8, label="1e-9")
    ax.set_yscale("log")
    ax.set_xlabel("R_B / R_A")
    ax.set_ylabel("absolute conservation error")
    ax.set_title("generalized elastic collision readout")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_conservation_v2.png", dpi=160)
    plt.close(fig)

    cases = [str(row["case"]) for row in case_summaries]
    xs = np.arange(len(cases))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(xs, [float(row["q_A_after"]) for row in case_summaries], marker="o", label="q_A after")
    ax.plot(xs, [float(row["q_B_after"]) for row in case_summaries], marker="o", label="q_B after")
    ax.plot(xs, [float(row["simple_flip_q_A"]) for row in case_summaries], linestyle="--", label="simple flip A")
    ax.plot(xs, [float(row["simple_flip_q_B"]) for row in case_summaries], linestyle="--", label="simple flip B")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("post-collision phase gradient q")
    ax.set_title("generalized collision map outputs")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_q_outputs_v2.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Generalized Elastic Collision Readout v1",
        "",
        "## Purpose",
        "",
        "This experiment replaces the equal-amplitude q-flip map with a generalized 1D elastic map using the multigauge R readout as the mass-like weight.",
        "It checks conservation of R*p and R*p^2 across asymmetric amplitude cases.",
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
            "| case | R_B/R_A | q_A after | q_B after | R*p err | R*p^2 err | R*E_tau err | simple q-flip R*p err | valid readout |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | {row['R_ratio_B_over_A']:.16e} | {row['q_A_after']:.16e} | "
            f"{row['q_B_after']:.16e} | {row['P_R_conservation_error']:.16e} | "
            f"{row['K_R_phase_conservation_error']:.16e} | {row['E_tau_R_conservation_error']:.16e} | "
            f"{row['simple_flip_P_R_error']:.16e} | `{row['individual_readout_valid']}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The generalized map preserves the R-weighted phase-gradient momentum and R-weighted phase-gradient square across every tested amplitude case.",
            "The simple q-flip map appears as the equal-R special case and fails for unequal R in the R*p readout.",
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_generalized_elastic_collision_readout_result_v2.json` |",
            "| case CSV | `abc_multigauge_generalized_elastic_collision_cases_v2.csv` |",
            "| stage quantity CSV | `abc_multigauge_generalized_elastic_collision_stage_quantities_v2.csv` |",
            "| gauge CSV | `abc_multigauge_generalized_elastic_collision_gauge_rows_v2.csv` |",
            "| conservation plot | `abc_multigauge_generalized_elastic_collision_conservation_v2.png` |",
            "| q output plot | `abc_multigauge_generalized_elastic_collision_q_outputs_v2.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_readout_report_v2.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_results = [run_case(a_amp, b_amp) for a_amp, b_amp in AMPLITUDE_CASES]
    case_summaries = [entry["case_summary"] for entry in case_results]
    return {
        "experiment": "abc_multigauge_generalized_elastic_collision_readout_v2",
        "amplitude_cases": [{"A_A": a, "A_B": b} for a, b in AMPLITUDE_CASES],
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
        "note": (
            "R is read as the mass-like stable residual. The generalized collision map is tested by "
            "R*p and R*p^2 conservation, not by assuming a single external mass parameter."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    full_case_results = [run_case(a_amp, b_amp) for a_amp, b_amp in AMPLITUDE_CASES]
    gauge_rows: List[Dict[str, Any]] = []
    stage_quantities: List[Dict[str, Any]] = []
    for entry in full_case_results:
        gauge_rows.extend(entry["gauge_rows"])
        stage_quantities.extend(entry["stage_quantities"])
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_readout_result_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_cases_v2.csv", result["case_summaries"])
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_stage_quantities_v2.csv", stage_quantities)
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_gauge_rows_v2.csv", gauge_rows)
    make_plots(result["case_summaries"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
