from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

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
OUT_DIR = BASE_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v2"
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
    q_A0: float = 1.0
    q_B0: float = -1.0


VELOCITY_CASES: List[Dict[str, float]] = [
    {"A_A": 1.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -1.0},
    {"A_A": 1.0, "A_B": 1.0, "q_A0": 1.4, "q_B0": -0.6},
    {"A_A": 1.0, "A_B": 1.0, "q_A0": 0.5, "q_B0": -1.7},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.0, "q_B0": -1.0},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.4, "q_B0": -0.6},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.2, "q_B0": 0.2},
    {"A_A": 2.0, "A_B": 1.0, "q_A0": 0.8, "q_B0": -1.5},
    {"A_A": 1.5, "A_B": 1.0, "q_A0": 1.8, "q_B0": -0.2},
    {"A_A": 1.0, "A_B": 3.0, "q_A0": 0.8, "q_B0": -0.4},
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


def generalized_elastic_velocity(r_a: float, r_b: float, u_a: float, u_b: float) -> tuple[float, float]:
    denom = r_a + r_b
    v_a = ((r_a - r_b) / denom) * u_a + (2.0 * r_b / denom) * u_b
    v_b = (2.0 * r_a / denom) * u_a + ((r_b - r_a) / denom) * u_b
    return float(v_a), float(v_b)


def simulate(params: Params) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
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

    r_a = params.A_A**2
    r_b = params.A_B**2
    q_a_before = float(a.q)
    q_b_before = float(b.q)
    q_a_after, q_b_after = generalized_elastic_velocity(r_a, r_b, q_a_before, q_b_before)
    if collision_cell_reached:
        a.q = q_a_after
        b.q = q_b_after

    events.append(
        {
            "event": "generalized_velocity_sweep_collision",
            "step": step,
            "collision_cell_reached": collision_cell_reached,
            "R_A": r_a,
            "R_B": r_b,
            "q_A_before": q_a_before,
            "q_B_before": q_b_before,
            "q_A_after": float(a.q),
            "q_B_after": float(b.q),
            "relative_q_before": q_a_before - q_b_before,
            "relative_q_after": float(a.q) - float(b.q),
            "closure_residual_abs": closure_residual(params),
            "m_A": params.m_A,
            "m_B": params.m_B,
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


def stage_quantities(case: str, summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
                "case": case,
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
                "relative_p": p_a - p_b,
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


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    return float(max(abs(float(row[key])) for row in rows))


def case_summary(
    params: Params,
    case: str,
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
    k_initial = float(initial["K_R_phase_total"])
    k_after = float(after["K_R_phase_total"])
    k_final = float(final["K_R_phase_total"])
    e_initial = float(initial["E_tau_R_total"])
    e_final = float(final["E_tau_R_total"])
    r_initial = float(initial["R_total"])
    r_final = float(final["R_total"])
    relative_initial = float(initial["relative_p"])
    relative_after = float(after["relative_p"])
    p_r_error = max(abs(p_r_after - p_r_initial), abs(p_r_final - p_r_initial))
    k_error = max(abs(k_after - k_initial), abs(k_final - k_initial))
    e_error = abs(e_final - e_initial)
    r_error = abs(r_final - r_initial)
    relative_flip_error = abs(relative_after + relative_initial)
    return {
        "case": case,
        "A_A": params.A_A,
        "A_B": params.A_B,
        "q_A0": params.q_A0,
        "q_B0": params.q_B0,
        "R_A": params.A_A**2,
        "R_B": params.A_B**2,
        "R_ratio_B_over_A": (params.A_B**2) / (params.A_A**2),
        "collision_cell_reached": bool(event["collision_cell_reached"]),
        "q_A_after": float(event["q_A_after"]),
        "q_B_after": float(event["q_B_after"]),
        "relative_q_before": float(event["relative_q_before"]),
        "relative_q_after": float(event["relative_q_after"]),
        "relative_flip_error": relative_flip_error,
        "p_max_abs_error": p_max,
        "E_max_abs_error": e_max,
        "R_max_abs_error": r_max,
        "R_max_gauge_std": r_std,
        "within_particle_separation_ratio_time": tr_ratio,
        "P_R_conservation_error": p_r_error,
        "K_R_phase_conservation_error": k_error,
        "E_tau_R_conservation_error": e_error,
        "R_total_conservation_error": r_error,
        "individual_readout_valid": bool(
            p_max <= params.readout_tol
            and e_max <= params.readout_tol
            and r_max <= params.readout_tol
            and r_std <= params.r_gauge_tol
            and tr_ratio <= params.tr_separation_threshold
        ),
        "P_R_preserved": bool(p_r_error <= params.conservation_tol),
        "K_R_phase_preserved": bool(k_error <= params.conservation_tol),
        "relative_gradient_flipped": bool(relative_flip_error <= params.conservation_tol),
        "E_tau_R_preserved": bool(e_error <= params.conservation_tol),
        "R_total_preserved": bool(r_error <= params.conservation_tol),
    }


def run_case(index: int, settings: Dict[str, float]) -> Dict[str, Any]:
    params = Params(**settings)
    gauges = default_gauges(params)
    stages, events = simulate(params)
    gauge_rows = readout_all_float(stages, gauges, params)
    summaries = summarize_stage_readouts(gauge_rows)
    case = (
        f"c{index:02d}_A{params.A_A:.2f}_B{params.A_B:.2f}_"
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
        "collision_reached_all_cases": all(bool(row["collision_cell_reached"]) for row in rows),
        "individual_readout_valid_all_cases": all(bool(row["individual_readout_valid"]) for row in rows),
        "P_R_preserved_all_cases": all(bool(row["P_R_preserved"]) for row in rows),
        "K_R_phase_preserved_all_cases": all(bool(row["K_R_phase_preserved"]) for row in rows),
        "relative_gradient_flipped_all_cases": all(bool(row["relative_gradient_flipped"]) for row in rows),
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
        "max_K_R_phase_conservation_error": float(
            max(float(row["K_R_phase_conservation_error"]) for row in rows)
        ),
        "max_relative_flip_error": float(max(float(row["relative_flip_error"]) for row in rows)),
        "max_E_tau_R_conservation_error": float(max(float(row["E_tau_R_conservation_error"]) for row in rows)),
        "max_R_total_conservation_error": float(max(float(row["R_total_conservation_error"]) for row in rows)),
        "single_gauge_only_used": False,
        "velocity_sweep_generalized_collision_valid": all(
            [
                all(bool(row["collision_cell_reached"]) for row in rows),
                all(bool(row["individual_readout_valid"]) for row in rows),
                all(bool(row["P_R_preserved"]) for row in rows),
                all(bool(row["K_R_phase_preserved"]) for row in rows),
                all(bool(row["relative_gradient_flipped"]) for row in rows),
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
    cases = [str(row["case"]) for row in case_summaries]
    xs = np.arange(len(cases))
    floor = 1.0e-18
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(
        xs,
        np.maximum([float(row["P_R_conservation_error"]) for row in case_summaries], floor),
        marker="o",
        label="R*p error",
    )
    ax.plot(
        xs,
        np.maximum([float(row["K_R_phase_conservation_error"]) for row in case_summaries], floor),
        marker="o",
        label="R*p^2 error",
    )
    ax.plot(
        xs,
        np.maximum([float(row["relative_flip_error"]) for row in case_summaries], floor),
        marker="o",
        label="relative p flip error",
    )
    ax.axhline(1.0e-9, color="black", linestyle="--", linewidth=0.8, label="1e-9")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("absolute error")
    ax.set_title("generalized elastic collision velocity sweep")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_errors_v2.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(xs, [float(row["q_A0"]) for row in case_summaries], marker="o", linestyle="--", label="q_A before")
    ax.plot(xs, [float(row["q_B0"]) for row in case_summaries], marker="o", linestyle="--", label="q_B before")
    ax.plot(xs, [float(row["q_A_after"]) for row in case_summaries], marker="o", label="q_A after")
    ax.plot(xs, [float(row["q_B_after"]) for row in case_summaries], marker="o", label="q_B after")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("phase gradient")
    ax.set_title("pre/post collision phase gradients")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_q_outputs_v2.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Generalized Elastic Collision Velocity Sweep v1",
        "",
        "## Purpose",
        "",
        "This experiment tests the generalized R-weighted elastic collision map under asymmetric initial phase gradients.",
        "It checks that the construction is not limited to the initial +1/-1 counter-propagating condition.",
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
            "| case | R_B/R_A | q before A/B | q after A/B | R*p err | R*p^2 err | relative flip err | valid |",
            "|---|---:|---|---|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | {row['R_ratio_B_over_A']:.16e} | "
            f"{row['q_A0']:.6g} / {row['q_B0']:.6g} | "
            f"{row['q_A_after']:.6g} / {row['q_B_after']:.6g} | "
            f"{row['P_R_conservation_error']:.16e} | {row['K_R_phase_conservation_error']:.16e} | "
            f"{row['relative_flip_error']:.16e} | `{row['individual_readout_valid']}` |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v2.json` |",
            "| case CSV | `abc_multigauge_generalized_elastic_collision_velocity_sweep_cases_v2.csv` |",
            "| stage quantity CSV | `abc_multigauge_generalized_elastic_collision_velocity_sweep_stage_quantities_v2.csv` |",
            "| gauge CSV | `abc_multigauge_generalized_elastic_collision_velocity_sweep_gauge_rows_v2.csv` |",
            "| error plot | `abc_multigauge_generalized_elastic_collision_velocity_sweep_errors_v2.png` |",
            "| q output plot | `abc_multigauge_generalized_elastic_collision_velocity_sweep_q_outputs_v2.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_report_v2.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_results = [run_case(index + 1, settings) for index, settings in enumerate(VELOCITY_CASES)]
    case_summaries = [entry["case_summary"] for entry in case_results]
    return {
        "experiment": "abc_multigauge_generalized_elastic_collision_velocity_sweep_v2",
        "velocity_cases": VELOCITY_CASES,
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
            "The sweep tests non-unit and same-direction catch-up collisions. "
            "Only cases with q_A0 > q_B0 are used so that the interaction cell is reached."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    full_case_results = [run_case(index + 1, settings) for index, settings in enumerate(VELOCITY_CASES)]
    gauge_rows: List[Dict[str, Any]] = []
    stage_quantities_rows: List[Dict[str, Any]] = []
    for entry in full_case_results:
        gauge_rows.extend(entry["gauge_rows"])
        stage_quantities_rows.extend(entry["stage_quantities"])
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_cases_v2.csv", result["case_summaries"])
    write_csv(
        OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_stage_quantities_v2.csv",
        stage_quantities_rows,
    )
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_velocity_sweep_gauge_rows_v2.csv", gauge_rows)
    make_plots(result["case_summaries"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
