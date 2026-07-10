from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_multi_collision_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class Params:
    A_A: float = 1.0
    A_B: float = 1.0
    A_C: float = 1000.0
    Nh_chi_A: int = 99
    Nh_chi_B: int = 99
    Nh_tau_A: int = 99
    Nh_tau_B: int = 99
    Nh_chi_C: int = 999
    Nh_tau_C: int = 999
    chi_A0: float = -0.2
    chi_B0: float = 0.2
    tau_A0: float = -0.2
    tau_B0: float = -0.2
    q_A0: int = 1
    q_B0: int = -1
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    delta_s: float = 0.001
    s_max: int = 5000
    target_ab_collisions: int = 8
    chi_min: float = -0.25
    chi_max: float = 0.25
    m_A: int = 1
    m_B: int = 2
    closure_tol: float = 1e-12


@dataclass
class State:
    chi: float
    tau: float
    q: int
    amplitude: float
    m: int


def harmonic_component_count(nh: int) -> int:
    return (nh + 1) // 2


def closure_coefficients(amplitude: float, nh_chi: int, nh_tau: int) -> np.ndarray:
    count = harmonic_component_count(nh_chi) * harmonic_component_count(nh_tau)
    return np.full(count, amplitude / count, dtype=complex)


def closure_residual(params: Params) -> float:
    coeffs = np.concatenate(
        [
            closure_coefficients(params.A_A, params.Nh_chi_A, params.Nh_tau_A),
            closure_coefficients(params.A_B, params.Nh_chi_B, params.Nh_tau_B),
            closure_coefficients(params.A_C, params.Nh_chi_C, params.Nh_tau_C),
        ]
    )
    return float(abs(np.sum(coeffs**2) + np.sum((1j * coeffs) ** 2)))


def append_timeline(rows: List[Dict[str, object]], stage: str, step: int, a: State, b: State, ab_count: int, wall_count: int) -> None:
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
            "A_A": a.amplitude,
            "A_B": b.amplitude,
            "ab_collision_count": ab_count,
            "wall_reflection_count": wall_count,
        }
    )


def in_ab_cell(a: State, b: State, eps_chi: float, eps_tau: float) -> bool:
    return abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau


def reflect_at_walls(state: State, params: Params) -> bool:
    reflected = False
    if state.chi <= params.chi_min:
        overshoot = params.chi_min - state.chi
        state.chi = params.chi_min + overshoot
        state.q = abs(state.q)
        reflected = True
    elif state.chi >= params.chi_max:
        overshoot = state.chi - params.chi_max
        state.chi = params.chi_max - overshoot
        state.q = -abs(state.q)
        reflected = True
    return reflected


def run() -> Dict[str, object]:
    params = Params()
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)

    a = State(params.chi_A0, params.tau_A0, params.q_A0, params.A_A, params.m_A)
    b = State(params.chi_B0, params.tau_B0, params.q_B0, params.A_B, params.m_B)
    timeline: List[Dict[str, object]] = []
    events: List[Dict[str, object]] = []
    closure_events: List[float] = []

    ab_collision_count = 0
    wall_reflection_count = 0
    inside_ab_cell = False
    append_timeline(timeline, "initial", 0, a, b, ab_collision_count, wall_reflection_count)

    for step in range(1, params.s_max + 1):
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s

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
                }
            )
            append_timeline(timeline, "wall_reflection", step, a, b, ab_collision_count, wall_reflection_count)

        currently_inside = in_ab_cell(a, b, eps_chi_ab, eps_tau_ab)
        if currently_inside and not inside_ab_cell:
            before_q_a = a.q
            before_q_b = b.q
            a.q = -a.q
            b.q = -b.q
            ab_collision_count += 1
            residual = closure_residual(params)
            closure_events.append(residual)
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
                    "q_B_before": before_q_b,
                    "q_A_after": a.q,
                    "q_B_after": b.q,
                    "m_A": a.m,
                    "m_B": b.m,
                    "closure_residual_abs": residual,
                }
            )
            append_timeline(timeline, "ab_collision", step, a, b, ab_collision_count, wall_reflection_count)
            inside_ab_cell = True
        elif not currently_inside:
            inside_ab_cell = False

        if step % 10 == 0:
            append_timeline(timeline, "step", step, a, b, ab_collision_count, wall_reflection_count)

        if ab_collision_count >= params.target_ab_collisions and not inside_ab_cell:
            append_timeline(timeline, "final", step, a, b, ab_collision_count, wall_reflection_count)
            break
    else:
        append_timeline(timeline, "final", params.s_max, a, b, ab_collision_count, wall_reflection_count)

    ab_events = [event for event in events if event["event"] == "ab_collision"]
    label_preserved_all = all(event["m_A"] == params.m_A and event["m_B"] == params.m_B for event in ab_events)
    q_reversed_each_collision = all(
        event["q_A_after"] == -event["q_A_before"] and event["q_B_after"] == -event["q_B_before"]
        for event in ab_events
    )
    amplitude_preserved = math.isclose(a.amplitude, params.A_A) and math.isclose(b.amplitude, params.A_B)
    closure_preserved_all = all(value <= params.closure_tol for value in closure_events)
    completed_target_collisions = ab_collision_count >= params.target_ab_collisions
    final_labels_preserved = a.m == params.m_A and b.m == params.m_B
    final_inside_bounds = params.chi_min <= a.chi <= params.chi_max and params.chi_min <= b.chi <= params.chi_max
    multi_collision_valid = all(
        [
            completed_target_collisions,
            label_preserved_all,
            final_labels_preserved,
            q_reversed_each_collision,
            amplitude_preserved,
            closure_preserved_all,
            final_inside_bounds,
        ]
    )

    result = {
        "parameters": asdict(params),
        "cell_widths": {
            "epsilon_chi_AB": eps_chi_ab,
            "epsilon_tau_AB": eps_tau_ab,
        },
        "summary": {
            "target_ab_collisions": params.target_ab_collisions,
            "ab_collision_count": ab_collision_count,
            "wall_reflection_count": wall_reflection_count,
            "completed_target_collisions": completed_target_collisions,
            "label_preserved_all": label_preserved_all,
            "q_reversed_each_collision": q_reversed_each_collision,
            "amplitude_preserved": amplitude_preserved,
            "closure_preserved_all": closure_preserved_all,
            "final_inside_bounds": final_inside_bounds,
            "multi_collision_valid": multi_collision_valid,
        },
        "events": events,
        "timeline": timeline,
    }
    return result


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "multi_collision_result_v1.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    event_path = OUT_DIR / "multi_collision_events_v1.csv"
    events = result["events"]
    event_keys = sorted({key for event in events for key in event.keys()})
    with event_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=event_keys)
        writer.writeheader()
        writer.writerows(events)

    timeline_path = OUT_DIR / "multi_collision_timeline_v1.csv"
    timeline = result["timeline"]
    with timeline_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timeline[0].keys()))
        writer.writeheader()
        writer.writerows(timeline)

    summary = result["summary"]
    lines = [
        "# Multi Collision Result v1",
        "",
        "## Verdict",
        "",
        f"- target_ab_collisions: `{summary['target_ab_collisions']}`",
        f"- ab_collision_count: `{summary['ab_collision_count']}`",
        f"- wall_reflection_count: `{summary['wall_reflection_count']}`",
        f"- completed_target_collisions: `{summary['completed_target_collisions']}`",
        f"- label_preserved_all: `{summary['label_preserved_all']}`",
        f"- q_reversed_each_collision: `{summary['q_reversed_each_collision']}`",
        f"- amplitude_preserved: `{summary['amplitude_preserved']}`",
        f"- closure_preserved_all: `{summary['closure_preserved_all']}`",
        f"- multi_collision_valid: `{summary['multi_collision_valid']}`",
        "",
        "## AB Collision Events",
        "",
        "| index | step | chi_A | chi_B | q_A before/after | q_B before/after | m_A | m_B | closure residual |",
        "|---:|---:|---:|---:|---|---|---:|---:|---:|",
    ]
    for event in [entry for entry in events if entry["event"] == "ab_collision"]:
        lines.append(
            f"| {event['collision_index']} | {event['step']} | {event['chi_A']:.12g} | {event['chi_B']:.12g} | "
            f"{event['q_A_before']} -> {event['q_A_after']} | {event['q_B_before']} -> {event['q_B_after']} | "
            f"{event['m_A']} | {event['m_B']} | {event['closure_residual_abs']:.12g} |"
        )
    (OUT_DIR / "multi_collision_report_v1.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    steps = [row["step"] for row in timeline]
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(steps, [row["chi_A"] for row in timeline], label="chi_A")
    axes[0].plot(steps, [row["chi_B"] for row in timeline], label="chi_B")
    axes[0].axhline(result["parameters"]["chi_min"], color="black", linestyle="--", linewidth=1)
    axes[0].axhline(result["parameters"]["chi_max"], color="black", linestyle="--", linewidth=1)
    axes[0].set_ylabel("position phase")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].step(steps, [row["q_A"] for row in timeline], where="post", label="q_A")
    axes[1].step(steps, [row["q_B"] for row in timeline], where="post", label="q_B")
    axes[1].set_xlabel("step")
    axes[1].set_ylabel("direction readout")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "multi_collision_trajectory_v1.png", dpi=180)
    plt.close(fig)

    collision_events = [entry for entry in events if entry["event"] == "ab_collision"]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(
        [entry["collision_index"] for entry in collision_events],
        [entry["closure_residual_abs"] for entry in collision_events],
        marker="o",
    )
    ax.set_xlabel("AB collision index")
    ax.set_ylabel("closure residual abs")
    ax.set_title("Closure residual across repeated collisions")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "multi_collision_closure_v1.png", dpi=180)
    plt.close(fig)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print(f"result_dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
