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
    readout_all,
    summarize_stage_readouts,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_interference_readout_multi_collision_result_v1"
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
    target_ab_collisions: int = 8
    chi_min: float = -0.25
    chi_max: float = 0.25
    s_max: int = 5000


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
            "chi_A": a.chi,
            "chi_B": b.chi,
            "tau_A": a.tau,
            "tau_B": b.tau,
            "q_A": a.q,
            "q_B": b.q,
            "m_A": a.m,
            "m_B": b.m,
            "omega_A": a.omega,
            "omega_B": b.omega,
            "ab_collision_count": ab_collision_count,
            "wall_reflection_count": wall_reflection_count,
        }
    )


def in_ab_cell(a: State, b: State, eps_chi: float, eps_tau: float) -> bool:
    return abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau


def reflect_at_walls(state: State, params: Params) -> bool:
    if state.chi <= params.chi_min:
        overshoot = params.chi_min - state.chi
        state.chi = params.chi_min + overshoot
        state.q = abs(state.q)
        return True
    if state.chi >= params.chi_max:
        overshoot = state.chi - params.chi_max
        state.chi = params.chi_max - overshoot
        state.q = -abs(state.q)
        return True
    return False


def simulate(params: Params) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    eps_chi = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    a = State(params.chi_A0, params.tau_A0, params.q_A0, params.A_A, params.m_A, params.omega_A)
    b = State(params.chi_B0, params.tau_B0, params.q_B0, params.A_B, params.m_B, params.omega_B)
    stages: List[Dict[str, Any]] = []
    events: List[Dict[str, Any]] = []
    ab_collision_count = 0
    wall_reflection_count = 0
    inside_ab_cell = False
    append_stage(stages, "initial", 0, a, b, ab_collision_count, wall_reflection_count)

    for step in range(1, params.s_max + 1):
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
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
                    "q_A": a.q,
                    "q_B": b.q,
                    "m_A": a.m,
                    "m_B": b.m,
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
            before_q_a = a.q
            before_q_b = b.q
            append_stage(
                stages,
                f"ab_collision_{ab_collision_count + 1}_before",
                step,
                before_a,
                before_b,
                ab_collision_count,
                wall_reflection_count,
            )
            a.q = -a.q
            b.q = -b.q
            ab_collision_count += 1
            residual = closure_residual(params)
            events.append(
                {
                    "event": "ab_collision",
                    "collision_index": ab_collision_count,
                    "step": step,
                    "chi_A": a.chi,
                    "chi_B": b.chi,
                    "tau_A": a.tau,
                    "tau_B": b.tau,
                    "q_A_before": before_q_a,
                    "q_A_after": a.q,
                    "q_B_before": before_q_b,
                    "q_B_after": b.q,
                    "m_A": a.m,
                    "m_B": b.m,
                    "closure_residual_abs": residual,
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


def group_by(rows: Iterable[Dict[str, Any]], *keys: str) -> Dict[tuple[Any, ...], List[Dict[str, Any]]]:
    grouped: Dict[tuple[Any, ...], List[Dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        grouped.setdefault(key, []).append(row)
    return grouped


def summary_lookup(summaries: List[Dict[str, Any]]) -> Dict[tuple[str, str], Dict[str, Any]]:
    return {(str(row["stage"]), str(row["particle"])): row for row in summaries}


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    return float(max(abs(float(row[key])) for row in rows))


def collision_readout_rows(
    params: Params,
    events: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    lookup = summary_lookup(summaries)
    rows: List[Dict[str, Any]] = []
    for event in [entry for entry in events if entry["event"] == "ab_collision"]:
        index = int(event["collision_index"])
        for particle in ["A", "B"]:
            before = lookup[(f"ab_collision_{index}_before", particle)]
            after = lookup[(f"ab_collision_{index}_after", particle)]
            rows.append(
                {
                    "collision_index": index,
                    "particle": particle,
                    "p_before": float(before["p_mean"]),
                    "p_after": float(after["p_mean"]),
                    "p_reflection_error": float(abs(float(after["p_mean"]) + float(before["p_mean"]))),
                    "E_before": float(before["E_mean"]),
                    "E_after": float(after["E_mean"]),
                    "E_preservation_error": float(abs(float(after["E_mean"]) - float(before["E_mean"]))),
                    "R_before": float(before["R_mean"]),
                    "R_after": float(after["R_mean"]),
                    "R_preservation_error": float(abs(float(after["R_mean"]) - float(before["R_mean"]))),
                    "R_std_before": float(before["R_std"]),
                    "R_std_after": float(after["R_std"]),
                    "q_before": int(event[f"q_{particle}_before"]),
                    "q_after": int(event[f"q_{particle}_after"]),
                }
            )
    return rows


def compute_verdicts(
    params: Params,
    stages: List[Dict[str, Any]],
    events: List[Dict[str, Any]],
    gauge_rows: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
    collision_rows: List[Dict[str, Any]],
) -> Dict[str, Any]:
    ab_events = [entry for entry in events if entry["event"] == "ab_collision"]
    closure_values = [float(entry["closure_residual_abs"]) for entry in ab_events]
    completed_target_collisions = len(ab_events) >= params.target_ab_collisions
    q_reversed_each_collision = all(
        int(event["q_A_after"]) == -int(event["q_A_before"])
        and int(event["q_B_after"]) == -int(event["q_B_before"])
        for event in ab_events
    )
    label_preserved_all = all(int(event["m_A"]) == params.m_A and int(event["m_B"]) == params.m_B for event in ab_events)
    p_max_error = max_abs(gauge_rows, "p_abs_error")
    e_max_error = max_abs(gauge_rows, "E_abs_error")
    r_max_error = max_abs(gauge_rows, "R_abs_error")
    r_max_gauge_std = float(max(float(row["R_std"]) for row in summaries))
    max_p_reflection_error = float(max(row["p_reflection_error"] for row in collision_rows))
    max_e_preservation_error = float(max(row["E_preservation_error"] for row in collision_rows))
    max_r_preservation_error = float(max(row["R_preservation_error"] for row in collision_rows))
    r_all = np.array([float(row["R_read"]) for row in gauge_rows])
    t_all = np.array([float(row["t_read"]) for row in gauge_rows])
    var_r = float(np.var(r_all))
    var_t = float(np.var(t_all))
    separation_ratio_time = float(var_r / var_t) if var_t > 0.0 else float("inf")
    verdicts = {
        "target_ab_collisions": params.target_ab_collisions,
        "ab_collision_count": len(ab_events),
        "wall_reflection_count": len([entry for entry in events if entry["event"] == "wall_reflection"]),
        "completed_target_collisions": completed_target_collisions,
        "q_reversed_each_collision": q_reversed_each_collision,
        "label_preserved_all": label_preserved_all,
        "closure_max_residual_abs": float(max(closure_values)) if closure_values else float("inf"),
        "closure_preserved_all": bool(closure_values and max(closure_values) <= params.closure_tol),
        "p_reconstructed_all_gauges": bool(p_max_error <= params.readout_tol),
        "E_reconstructed_all_gauges": bool(e_max_error <= params.readout_tol),
        "R_reconstructed_all_gauges": bool(r_max_error <= params.readout_tol),
        "p_max_abs_error": p_max_error,
        "E_max_abs_error": e_max_error,
        "R_max_abs_error": r_max_error,
        "max_p_reflection_error": max_p_reflection_error,
        "p_reflection_each_collision": bool(max_p_reflection_error <= params.conservation_tol),
        "max_E_preservation_error": max_e_preservation_error,
        "E_preserved_each_collision": bool(max_e_preservation_error <= params.conservation_tol),
        "max_R_preservation_error": max_r_preservation_error,
        "R_preserved_each_collision": bool(max_r_preservation_error <= params.conservation_tol),
        "R_max_gauge_std": r_max_gauge_std,
        "R_gauge_stable": bool(r_max_gauge_std <= params.r_gauge_tol),
        "var_R_all": var_r,
        "var_t_all": var_t,
        "separation_ratio_time": separation_ratio_time,
        "t_R_separation_valid": bool(separation_ratio_time <= params.tr_separation_threshold),
        "single_gauge_only_used": False,
    }
    verdicts["multi_collision_multigauge_valid"] = all(
        bool(verdicts[key])
        for key in [
            "completed_target_collisions",
            "q_reversed_each_collision",
            "label_preserved_all",
            "closure_preserved_all",
            "p_reconstructed_all_gauges",
            "E_reconstructed_all_gauges",
            "R_reconstructed_all_gauges",
            "p_reflection_each_collision",
            "E_preserved_each_collision",
            "R_preserved_each_collision",
            "R_gauge_stable",
            "t_R_separation_valid",
        ]
    )
    return verdicts


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(collision_rows: List[Dict[str, Any]], verdicts: Dict[str, Any]) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for particle, color in [("A", "tab:blue"), ("B", "tab:orange")]:
        selected = [row for row in collision_rows if row["particle"] == particle]
        xs = [int(row["collision_index"]) for row in selected]
        axes[0].plot(xs, [row["p_before"] for row in selected], marker="o", linestyle="--", color=color, label=f"p before {particle}")
        axes[0].plot(xs, [row["p_after"] for row in selected], marker="o", color=color, label=f"p after {particle}")
        axes[1].plot(xs, [row["E_before"] for row in selected], marker="o", linestyle="--", color=color, label=f"E before {particle}")
        axes[1].plot(xs, [row["E_after"] for row in selected], marker="o", color=color, label=f"E after {particle}")
        axes[2].plot(xs, [row["R_before"] for row in selected], marker="o", linestyle="--", color=color, label=f"R before {particle}")
        axes[2].plot(xs, [row["R_after"] for row in selected], marker="o", color=color, label=f"R after {particle}")
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set_ylabel("p_read")
    axes[1].set_ylabel("E_read")
    axes[2].set_ylabel("R_read")
    axes[2].set_xlabel("AB collision index")
    for axis in axes:
        axis.ticklabel_format(useOffset=False)
        axis.grid(True, alpha=0.25)
        axis.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_multi_collision_invariants_v1.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(["Var(R)", "Var(t)", "Var(R)/Var(t)"], [
        float(verdicts["var_R_all"]),
        float(verdicts["var_t_all"]),
        float(verdicts["separation_ratio_time"]),
    ])
    ax.set_yscale("symlog", linthresh=1e-20)
    ax.set_title("multi-collision t/R separation readout")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_interference_readout_multi_collision_tr_separation_v1.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    verdicts = result["verdicts"]
    lines = [
        "# ABC Multigauge Interference Readout Multi Collision Result v1",
        "",
        "## Verdict",
        "",
    ]
    for key in [
        "target_ab_collisions",
        "ab_collision_count",
        "wall_reflection_count",
        "completed_target_collisions",
        "q_reversed_each_collision",
        "label_preserved_all",
        "closure_preserved_all",
        "p_reconstructed_all_gauges",
        "E_reconstructed_all_gauges",
        "R_reconstructed_all_gauges",
        "p_reflection_each_collision",
        "E_preserved_each_collision",
        "R_preserved_each_collision",
        "R_gauge_stable",
        "t_R_separation_valid",
        "single_gauge_only_used",
        "multi_collision_multigauge_valid",
    ]:
        lines.append(f"- {key}: `{verdicts[key]}`")
    lines.extend(
        [
            "",
            "## Key Numerical Values",
            "",
            f"- p_max_abs_error: `{verdicts['p_max_abs_error']:.16e}`",
            f"- E_max_abs_error: `{verdicts['E_max_abs_error']:.16e}`",
            f"- R_max_abs_error: `{verdicts['R_max_abs_error']:.16e}`",
            f"- max_p_reflection_error: `{verdicts['max_p_reflection_error']:.16e}`",
            f"- max_E_preservation_error: `{verdicts['max_E_preservation_error']:.16e}`",
            f"- max_R_preservation_error: `{verdicts['max_R_preservation_error']:.16e}`",
            f"- R_max_gauge_std: `{verdicts['R_max_gauge_std']:.16e}`",
            f"- closure_max_residual_abs: `{verdicts['closure_max_residual_abs']:.16e}`",
            f"- separation_ratio_time: `{verdicts['separation_ratio_time']:.16e}`",
            "",
            "## AB Collision Readout Checks",
            "",
            "| index | particle | p before | p after | p flip error | E error | R error | R std before/after |",
            "|---:|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["collision_readouts"]:
        lines.append(
            f"| {row['collision_index']} | {row['particle']} | {row['p_before']:.16e} | {row['p_after']:.16e} | "
            f"{row['p_reflection_error']:.16e} | {row['E_preservation_error']:.16e} | "
            f"{row['R_preservation_error']:.16e} | {row['R_std_before']:.3e} / {row['R_std_after']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_interference_readout_multi_collision_result_v1.json` |",
            "| timeline CSV | `abc_multigauge_interference_readout_multi_collision_timeline_v1.csv` |",
            "| events CSV | `abc_multigauge_interference_readout_multi_collision_events_v1.csv` |",
            "| gauge CSV | `abc_multigauge_interference_readout_multi_collision_gauge_sweep_v1.csv` |",
            "| stage summary CSV | `abc_multigauge_interference_readout_multi_collision_stage_summary_v1.csv` |",
            "| collision readout CSV | `abc_multigauge_interference_readout_multi_collision_readouts_v1.csv` |",
            "| invariant plot | `abc_multigauge_interference_readout_multi_collision_invariants_v1.png` |",
            "| t/R plot | `abc_multigauge_interference_readout_multi_collision_tr_separation_v1.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_interference_readout_multi_collision_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    gauges = default_gauges(params)
    stages, events = simulate(params)
    gauge_rows = readout_all(stages, gauges, params)
    summaries = summarize_stage_readouts(gauge_rows)
    collision_rows = collision_readout_rows(params, events, summaries)
    verdicts = compute_verdicts(params, stages, events, gauge_rows, summaries, collision_rows)
    return {
        "experiment": "abc_multigauge_interference_readout_multi_collision_v1",
        "parameters": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "stages": stages,
        "events": events,
        "stage_summaries": summaries,
        "collision_readouts": collision_rows,
        "verdicts": verdicts,
        "note": "Repeated AB collisions test whether the v1 multigauge p/E/R readout remains stable under collision repetition and wall reflections.",
    }


def write_outputs(result: Dict[str, Any]) -> None:
    (OUT_DIR / "abc_multigauge_interference_readout_multi_collision_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_multi_collision_timeline_v1.csv", result["stages"])
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_multi_collision_events_v1.csv", result["events"])
    params = Params(**result["parameters"])
    gauges = [Gauge(**data) for data in result["gauges"]]
    gauge_rows = readout_all(result["stages"], gauges, params)
    write_csv(OUT_DIR / "abc_multigauge_interference_readout_multi_collision_gauge_sweep_v1.csv", gauge_rows)
    write_csv(
        OUT_DIR / "abc_multigauge_interference_readout_multi_collision_stage_summary_v1.csv",
        result["stage_summaries"],
    )
    write_csv(
        OUT_DIR / "abc_multigauge_interference_readout_multi_collision_readouts_v1.csv",
        result["collision_readouts"],
    )
    make_plots(result["collision_readouts"], result["verdicts"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["verdicts"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
