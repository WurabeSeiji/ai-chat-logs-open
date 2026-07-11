from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_abc_multigauge_interference_readout_v1 import (
    Gauge,
    Params as BaseParams,
    State,
    closure_residual,
    default_gauges,
    read_particle_gauge,
    summarize_stage_readouts,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_result_v1"
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
    target_ab_collisions: int = 6
    chi_min: float = -0.25
    chi_max: float = 0.25
    s_max: int = 12000


MULTI_CASES: List[Dict[str, float]] = [
    {"A_A": 1.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -1.0},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.4, "q_B0": -0.6},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.2, "q_B0": 0.2},
    {"A_A": 1.5, "A_B": 1.0, "q_A0": 1.8, "q_B0": -0.2},
]


def copy_state(state: State) -> State:
    return State(state.chi, state.tau, state.q, state.amplitude, state.m, state.omega)


def append_stage(
    rows: List[Dict[str, Any]],
    stage: str,
    step: int,
    a: State,
    b: State,
    ab_collision_count: int,
    wall_reflection_count: int,
) -> None:
    rows.append(
        {
            "stage": stage,
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
            "ab_collision_count": ab_collision_count,
            "wall_reflection_count": wall_reflection_count,
        }
    )


def generalized_elastic_velocity(r_a: float, r_b: float, u_a: float, u_b: float) -> tuple[float, float]:
    denom = r_a + r_b
    v_a = ((r_a - r_b) / denom) * u_a + (2.0 * r_b / denom) * u_b
    v_b = (2.0 * r_a / denom) * u_a + ((r_b - r_a) / denom) * u_b
    return float(v_a), float(v_b)


def reflect_at_walls(state: State, params: Params) -> bool:
    if state.chi <= params.chi_min:
        overshoot = params.chi_min - state.chi
        state.chi = params.chi_min + overshoot
        state.q = abs(float(state.q))
        return True
    if state.chi >= params.chi_max:
        overshoot = state.chi - params.chi_max
        state.chi = params.chi_max - overshoot
        state.q = -abs(float(state.q))
        return True
    return False


def in_ab_cell(a: State, b: State, eps_chi: float, eps_tau: float) -> bool:
    return abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau


def simulate(params: Params) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    eps_chi = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    a = State(params.chi_A0, params.tau_A0, float(params.q_A0), params.A_A, params.m_A, params.omega_A)
    b = State(params.chi_B0, params.tau_B0, float(params.q_B0), params.A_B, params.m_B, params.omega_B)
    stages: List[Dict[str, Any]] = []
    events: List[Dict[str, Any]] = []
    ab_collision_count = 0
    wall_reflection_count = 0
    inside_ab_cell = False
    append_stage(stages, "initial", 0, a, b, ab_collision_count, wall_reflection_count)

    for step in range(1, params.s_max + 1):
        a.chi += float(a.q) * params.v_chi * params.delta_s
        b.chi += float(b.q) * params.v_chi * params.delta_s
        a.tau += a.omega * params.delta_s
        b.tau += b.omega * params.delta_s

        wall_a = reflect_at_walls(a, params)
        wall_b = reflect_at_walls(b, params)
        if wall_a or wall_b:
            wall_reflection_count += int(wall_a) + int(wall_b)
            events.append(
                {
                    "event": "wall_reflection",
                    "step": step,
                    "wall_A": wall_a,
                    "wall_B": wall_b,
                    "q_A": float(a.q),
                    "q_B": float(b.q),
                    "wall_reflection_count": wall_reflection_count,
                }
            )
            append_stage(
                stages,
                f"wall_reflection_{wall_reflection_count}",
                step,
                copy_state(a),
                copy_state(b),
                ab_collision_count,
                wall_reflection_count,
            )

        currently_inside = in_ab_cell(a, b, eps_chi, eps_tau)
        if currently_inside and not inside_ab_cell:
            before_a = copy_state(a)
            before_b = copy_state(b)
            before_q_a = float(a.q)
            before_q_b = float(b.q)
            append_stage(
                stages,
                f"ab_collision_{ab_collision_count + 1}_before",
                step,
                before_a,
                before_b,
                ab_collision_count,
                wall_reflection_count,
            )
            r_a = params.A_A**2
            r_b = params.A_B**2
            after_q_a, after_q_b = generalized_elastic_velocity(r_a, r_b, before_q_a, before_q_b)
            a.q = after_q_a
            b.q = after_q_b
            ab_collision_count += 1
            events.append(
                {
                    "event": "generalized_ab_collision",
                    "collision_index": ab_collision_count,
                    "step": step,
                    "chi_A": float(a.chi),
                    "chi_B": float(b.chi),
                    "tau_A": float(a.tau),
                    "tau_B": float(b.tau),
                    "R_A": r_a,
                    "R_B": r_b,
                    "q_A_before": before_q_a,
                    "q_A_after": float(a.q),
                    "q_B_before": before_q_b,
                    "q_B_after": float(b.q),
                    "relative_q_before": before_q_a - before_q_b,
                    "relative_q_after": float(a.q) - float(b.q),
                    "m_A": params.m_A,
                    "m_B": params.m_B,
                    "closure_residual_abs": closure_residual(params),
                }
            )
            append_stage(
                stages,
                f"ab_collision_{ab_collision_count}_after",
                step,
                copy_state(a),
                copy_state(b),
                ab_collision_count,
                wall_reflection_count,
            )
            inside_ab_cell = True
        elif not currently_inside:
            inside_ab_cell = False

        if ab_collision_count >= params.target_ab_collisions and not inside_ab_cell:
            append_stage(stages, "final", step, a, b, ab_collision_count, wall_reflection_count)
            break
    else:
        append_stage(stages, "final", params.s_max, a, b, ab_collision_count, wall_reflection_count)
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


def collision_rows(events: List[Dict[str, Any]], summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    lookup = summary_lookup(summaries)
    rows: List[Dict[str, Any]] = []
    for event in [entry for entry in events if entry["event"] == "generalized_ab_collision"]:
        index = int(event["collision_index"])
        before_a = lookup[(f"ab_collision_{index}_before", "A")]
        before_b = lookup[(f"ab_collision_{index}_before", "B")]
        after_a = lookup[(f"ab_collision_{index}_after", "A")]
        after_b = lookup[(f"ab_collision_{index}_after", "B")]
        r_a_before = float(before_a["R_mean"])
        r_b_before = float(before_b["R_mean"])
        p_a_before = float(before_a["p_mean"])
        p_b_before = float(before_b["p_mean"])
        p_a_after = float(after_a["p_mean"])
        p_b_after = float(after_b["p_mean"])
        e_tau_before = r_a_before * float(before_a["E_mean"]) + r_b_before * float(before_b["E_mean"])
        e_tau_after = r_a_before * float(after_a["E_mean"]) + r_b_before * float(after_b["E_mean"])
        p_r_before = r_a_before * p_a_before + r_b_before * p_b_before
        p_r_after = r_a_before * p_a_after + r_b_before * p_b_after
        k_r_before = r_a_before * p_a_before**2 + r_b_before * p_b_before**2
        k_r_after = r_a_before * p_a_after**2 + r_b_before * p_b_after**2
        rows.append(
            {
                "collision_index": index,
                "R_A": r_a_before,
                "R_B": r_b_before,
                "p_A_before": p_a_before,
                "p_B_before": p_b_before,
                "p_A_after": p_a_after,
                "p_B_after": p_b_after,
                "P_R_before": p_r_before,
                "P_R_after": p_r_after,
                "P_R_error": abs(p_r_after - p_r_before),
                "K_R_before": k_r_before,
                "K_R_after": k_r_after,
                "K_R_error": abs(k_r_after - k_r_before),
                "relative_before": p_a_before - p_b_before,
                "relative_after": p_a_after - p_b_after,
                "relative_flip_error": abs((p_a_after - p_b_after) + (p_a_before - p_b_before)),
                "E_tau_R_error": abs(e_tau_after - e_tau_before),
                "R_A_error": abs(float(after_a["R_mean"]) - r_a_before),
                "R_B_error": abs(float(after_b["R_mean"]) - r_b_before),
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


def case_summary(
    params: Params,
    case: str,
    events: List[Dict[str, Any]],
    gauge_rows: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
    collisions: List[Dict[str, Any]],
) -> Dict[str, Any]:
    ab_events = [entry for entry in events if entry["event"] == "generalized_ab_collision"]
    closure_values = [float(entry["closure_residual_abs"]) for entry in ab_events]
    return {
        "case": case,
        "A_A": params.A_A,
        "A_B": params.A_B,
        "q_A0": params.q_A0,
        "q_B0": params.q_B0,
        "target_ab_collisions": params.target_ab_collisions,
        "ab_collision_count": len(ab_events),
        "wall_reflection_count": len([entry for entry in events if entry["event"] == "wall_reflection"]),
        "completed_target_collisions": bool(len(ab_events) >= params.target_ab_collisions),
        "p_max_abs_error": max_abs(gauge_rows, "p_abs_error"),
        "E_max_abs_error": max_abs(gauge_rows, "E_abs_error"),
        "R_max_abs_error": max_abs(gauge_rows, "R_abs_error"),
        "R_max_gauge_std": float(max(float(row["R_std"]) for row in summaries)),
        "within_particle_separation_ratio_time": within_particle_tr_separation(summaries),
        "closure_max_residual_abs": float(max(closure_values)) if closure_values else float("inf"),
        "max_P_R_error": max_abs(collisions, "P_R_error"),
        "max_K_R_error": max_abs(collisions, "K_R_error"),
        "max_relative_flip_error": max_abs(collisions, "relative_flip_error"),
        "max_E_tau_R_error": max_abs(collisions, "E_tau_R_error"),
        "max_R_A_error": max_abs(collisions, "R_A_error"),
        "max_R_B_error": max_abs(collisions, "R_B_error"),
        "label_preserved_all": all(int(entry["m_A"]) == params.m_A and int(entry["m_B"]) == params.m_B for entry in ab_events),
    }


def with_case_verdict(params: Params, summary: Dict[str, Any]) -> Dict[str, Any]:
    result = dict(summary)
    result.update(
        {
            "individual_readout_valid": bool(
                summary["p_max_abs_error"] <= params.readout_tol
                and summary["E_max_abs_error"] <= params.readout_tol
                and summary["R_max_abs_error"] <= params.readout_tol
                and summary["R_max_gauge_std"] <= params.r_gauge_tol
                and summary["within_particle_separation_ratio_time"] <= params.tr_separation_threshold
            ),
            "closure_preserved": bool(summary["closure_max_residual_abs"] <= params.closure_tol),
            "P_R_preserved_each_collision": bool(summary["max_P_R_error"] <= params.conservation_tol),
            "K_R_preserved_each_collision": bool(summary["max_K_R_error"] <= 2.0e-9),
            "relative_flip_each_collision": bool(summary["max_relative_flip_error"] <= params.conservation_tol),
            "E_tau_R_preserved_each_collision": bool(summary["max_E_tau_R_error"] <= params.conservation_tol),
            "R_preserved_each_collision": bool(
                summary["max_R_A_error"] <= params.conservation_tol
                and summary["max_R_B_error"] <= params.conservation_tol
            ),
        }
    )
    result["case_valid"] = all(
        bool(result[key])
        for key in [
            "completed_target_collisions",
            "label_preserved_all",
            "individual_readout_valid",
            "closure_preserved",
            "P_R_preserved_each_collision",
            "K_R_preserved_each_collision",
            "relative_flip_each_collision",
            "E_tau_R_preserved_each_collision",
            "R_preserved_each_collision",
        ]
    )
    return result


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
    collisions = collision_rows(events, summaries)
    for row in collisions:
        row["case"] = case
        row["A_A"] = params.A_A
        row["A_B"] = params.A_B
        row["q_A0"] = params.q_A0
        row["q_B0"] = params.q_B0
    summary = with_case_verdict(params, case_summary(params, case, events, gauge_rows, summaries, collisions))
    return {
        "case": case,
        "parameters": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "stages": stages,
        "events": events,
        "gauge_rows": gauge_rows,
        "stage_summaries": summaries,
        "collision_rows": collisions,
        "case_summary": summary,
    }


def aggregate(case_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    rows = [entry["case_summary"] for entry in case_results]
    return {
        "case_count": len(rows),
        "all_cases_valid": all(bool(row["case_valid"]) for row in rows),
        "completed_target_collisions_all_cases": all(bool(row["completed_target_collisions"]) for row in rows),
        "individual_readout_valid_all_cases": all(bool(row["individual_readout_valid"]) for row in rows),
        "closure_preserved_all_cases": all(bool(row["closure_preserved"]) for row in rows),
        "P_R_preserved_each_collision_all_cases": all(bool(row["P_R_preserved_each_collision"]) for row in rows),
        "K_R_preserved_each_collision_all_cases": all(bool(row["K_R_preserved_each_collision"]) for row in rows),
        "relative_flip_each_collision_all_cases": all(bool(row["relative_flip_each_collision"]) for row in rows),
        "E_tau_R_preserved_each_collision_all_cases": all(bool(row["E_tau_R_preserved_each_collision"]) for row in rows),
        "R_preserved_each_collision_all_cases": all(bool(row["R_preserved_each_collision"]) for row in rows),
        "max_ab_collision_count": int(max(int(row["ab_collision_count"]) for row in rows)),
        "max_wall_reflection_count": int(max(int(row["wall_reflection_count"]) for row in rows)),
        "max_p_abs_error": float(max(float(row["p_max_abs_error"]) for row in rows)),
        "max_E_abs_error": float(max(float(row["E_max_abs_error"]) for row in rows)),
        "max_R_abs_error": float(max(float(row["R_max_abs_error"]) for row in rows)),
        "max_R_gauge_std": float(max(float(row["R_max_gauge_std"]) for row in rows)),
        "max_within_particle_separation_ratio_time": float(
            max(float(row["within_particle_separation_ratio_time"]) for row in rows)
        ),
        "max_P_R_error": float(max(float(row["max_P_R_error"]) for row in rows)),
        "max_K_R_error": float(max(float(row["max_K_R_error"]) for row in rows)),
        "max_relative_flip_error": float(max(float(row["max_relative_flip_error"]) for row in rows)),
        "max_E_tau_R_error": float(max(float(row["max_E_tau_R_error"]) for row in rows)),
        "max_R_A_error": float(max(float(row["max_R_A_error"]) for row in rows)),
        "max_R_B_error": float(max(float(row["max_R_B_error"]) for row in rows)),
        "single_gauge_only_used": False,
        "generalized_multi_collision_valid": all(bool(row["case_valid"]) for row in rows),
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(case_summaries: List[Dict[str, Any]], collision_rows_all: List[Dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    for case in sorted({str(row["case"]) for row in collision_rows_all}):
        selected = [row for row in collision_rows_all if row["case"] == case]
        xs = [int(row["collision_index"]) for row in selected]
        ys = np.maximum([float(row["P_R_error"]) for row in selected], 1.0e-18)
        ax.plot(xs, ys, marker="o", label=case)
    ax.axhline(1.0e-9, color="black", linestyle="--", linewidth=0.8, label="1e-9")
    ax.set_yscale("log")
    ax.set_xlabel("AB collision index")
    ax.set_ylabel("R*p conservation error")
    ax.set_title("generalized multi-collision R*p conservation")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_p_errors_v1.png", dpi=160)
    plt.close(fig)

    cases = [str(row["case"]) for row in case_summaries]
    xs = np.arange(len(cases))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(xs - 0.2, [float(row["max_P_R_error"]) for row in case_summaries], 0.2, label="R*p")
    ax.bar(xs, [float(row["max_K_R_error"]) for row in case_summaries], 0.2, label="R*p^2")
    ax.bar(xs + 0.2, [float(row["max_relative_flip_error"]) for row in case_summaries], 0.2, label="relative flip")
    ax.axhline(1.0e-9, color="black", linestyle="--", linewidth=0.8, label="1e-9")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, rotation=25, ha="right")
    ax.set_ylabel("max error")
    ax.set_title("generalized multi-collision max errors")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_summary_v1.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Generalized Elastic Collision Multi Collision v1",
        "",
        "## Purpose",
        "",
        "This experiment repeats the generalized R-weighted elastic collision map across multiple AB collisions with wall returns.",
        "It checks whether R*p, R*p^2, relative phase-gradient flip, R*E_tau, and R stability survive repeated collision cycles.",
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
            "| case | AB collisions | wall reflections | max R*p err | max R*p^2 err | max relative err | valid |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case']} | {row['ab_collision_count']} | {row['wall_reflection_count']} | "
            f"{row['max_P_R_error']:.16e} | {row['max_K_R_error']:.16e} | "
            f"{row['max_relative_flip_error']:.16e} | `{row['case_valid']}` |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_generalized_elastic_collision_multi_collision_result_v1.json` |",
            "| case CSV | `abc_multigauge_generalized_elastic_collision_multi_collision_cases_v1.csv` |",
            "| collision CSV | `abc_multigauge_generalized_elastic_collision_multi_collision_readouts_v1.csv` |",
            "| gauge CSV | `abc_multigauge_generalized_elastic_collision_multi_collision_gauge_rows_v1.csv` |",
            "| P error plot | `abc_multigauge_generalized_elastic_collision_multi_collision_p_errors_v1.png` |",
            "| summary plot | `abc_multigauge_generalized_elastic_collision_multi_collision_summary_v1.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_results = [run_case(index + 1, settings) for index, settings in enumerate(MULTI_CASES)]
    case_summaries = [entry["case_summary"] for entry in case_results]
    return {
        "experiment": "abc_multigauge_generalized_elastic_collision_multi_collision_v1",
        "multi_cases": MULTI_CASES,
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
            "AB collision conservation is evaluated only across AB collision before/after stages. "
            "Wall reflections are used to return particles for repeated encounters and are counted separately."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    full_case_results = [run_case(index + 1, settings) for index, settings in enumerate(MULTI_CASES)]
    gauge_rows: List[Dict[str, Any]] = []
    collision_rows_all: List[Dict[str, Any]] = []
    for entry in full_case_results:
        gauge_rows.extend(entry["gauge_rows"])
        collision_rows_all.extend(entry["collision_rows"])
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_cases_v1.csv", result["case_summaries"])
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_readouts_v1.csv", collision_rows_all)
    write_csv(OUT_DIR / "abc_multigauge_generalized_elastic_collision_multi_collision_gauge_rows_v1.csv", gauge_rows)
    make_plots(result["case_summaries"], collision_rows_all)
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
