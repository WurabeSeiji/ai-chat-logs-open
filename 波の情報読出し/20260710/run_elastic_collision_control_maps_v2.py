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
OUT_DIR = BASE_DIR / "elastic_collision_control_maps_result_v2"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fermionic_reflection_q(q_initial: float, delta_f: float = math.pi) -> float:
    reflection_rate = math.sin(delta_f / 2.0) ** 2
    transmission_rate = math.cos(delta_f / 2.0) ** 2
    return q_initial * (transmission_rate - reflection_rate)


@dataclass
class Params:
    A_A: float = 1.0
    A_B: float = 1.0
    A_C: float = 1000.0
    Nh_chi_A: int = 99
    Nh_chi_B: int = 99
    Nh_tau_A: int = 99
    Nh_tau_B: int = 99
    chi_A0: float = -0.2
    chi_B0: float = 0.2
    tau0: float = 0.2
    q_A0: int = 1
    q_B0: int = -1
    delta_s: float = 0.01
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    s_max: int = 10000
    m_A: int = 1
    m_B: int = 2
    grid_n: int = 4096


@dataclass
class State:
    name: str
    chi: float
    tau: float
    q: int
    m: int


def phase_grid(n: int) -> np.ndarray:
    return np.linspace(-math.pi, math.pi, n, endpoint=False)


def eta_overlap(m_particle: int, m_read: int, grid: np.ndarray) -> complex:
    return complex(np.mean(np.exp(1j * m_particle * grid) * np.exp(-1j * m_read * grid)))


def detect_mode(state: State, modes: List[int], grid: np.ndarray) -> int:
    observations = {mode: eta_overlap(state.m, mode, grid) for mode in modes}
    return max(observations, key=lambda mode: abs(observations[mode]))


def append_timeline(rows: List[Dict[str, object]], map_name: str, stage: str, step: int, a: State, b: State) -> None:
    rows.append(
        {
            "map": map_name,
            "stage": stage,
            "step": step,
            "chi_A": a.chi,
            "chi_B": b.chi,
            "tau_A": a.tau,
            "tau_B": b.tau,
            "q_A": a.q,
            "q_B": b.q,
            "m_A_state": a.m,
            "m_B_state": b.m,
        }
    )


def apply_collision_map(map_name: str, a: State, b: State) -> None:
    if map_name == "reflection":
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)
    elif map_name == "transmission":
        pass
    elif map_name == "label_exchange_reflection":
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)
        a.m, b.m = b.m, a.m
    elif map_name == "transmission_with_label_exchange":
        a.m, b.m = b.m, a.m
    else:
        raise ValueError(f"unknown map: {map_name}")


def run_case(params: Params, map_name: str, grid: np.ndarray) -> Dict[str, object]:
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    modes = [params.m_A, params.m_B]
    a = State("A", params.chi_A0, -params.tau0, params.q_A0, params.m_A)
    b = State("B", params.chi_B0, -params.tau0, params.q_B0, params.m_B)
    timeline: List[Dict[str, object]] = []
    append_timeline(timeline, map_name, "initial", 0, a, b)

    step = 0
    collision_cell_reached = False
    while abs(a.chi - b.chi) >= eps_chi_ab or abs(a.tau - b.tau) >= eps_tau_ab:
        if step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        step += 1
        append_timeline(timeline, map_name, "approach", step, a, b)
    else:
        collision_cell_reached = True

    append_timeline(timeline, map_name, "collision_cell", step, a, b)
    if collision_cell_reached:
        apply_collision_map(map_name, a, b)
    append_timeline(timeline, map_name, "collision_map", step, a, b)

    post_step = 0
    post_completed = False
    while min(a.tau, b.tau) < params.tau0 - 1e-12:
        if post_step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        post_step += 1
        append_timeline(timeline, map_name, "post_collision", step + post_step, a, b)
    else:
        post_completed = collision_cell_reached

    append_timeline(timeline, map_name, "final", step + post_step, a, b)

    final_left, final_right = sorted([a, b], key=lambda state: state.chi)
    final_left_mode = detect_mode(final_left, modes, grid)
    final_right_mode = detect_mode(final_right, modes, grid)

    q_reversed = a.q == -params.q_A0 and b.q == -params.q_B0
    identity_labels_preserved = a.m == params.m_A and b.m == params.m_B
    identity_labels_swapped = a.m == params.m_B and b.m == params.m_A
    reflection_spatial_pattern = final_left_mode == params.m_A and final_right_mode == params.m_B
    transmission_spatial_pattern = final_left_mode == params.m_B and final_right_mode == params.m_A
    reflection_position_order = a.chi < b.chi
    transmission_position_order = a.chi > b.chi
    reflection_valid = (
        collision_cell_reached
        and post_completed
        and q_reversed
        and identity_labels_preserved
        and reflection_spatial_pattern
        and reflection_position_order
    )

    return {
        "map": map_name,
        "collision_cell_reached": collision_cell_reached,
        "post_collision_propagation_completed": post_completed,
        "q_A_final": a.q,
        "q_B_final": b.q,
        "q_reversed": q_reversed,
        "A_identity_mode_final": a.m,
        "B_identity_mode_final": b.m,
        "identity_labels_preserved": identity_labels_preserved,
        "identity_labels_swapped": identity_labels_swapped,
        "final_chi_A": a.chi,
        "final_chi_B": b.chi,
        "reflection_position_order": reflection_position_order,
        "transmission_position_order": transmission_position_order,
        "left_slot_identity": final_left.name,
        "right_slot_identity": final_right.name,
        "left_slot_mode": final_left_mode,
        "right_slot_mode": final_right_mode,
        "reflection_spatial_pattern": reflection_spatial_pattern,
        "transmission_spatial_pattern": transmission_spatial_pattern,
        "reflection_valid": reflection_valid,
        "collision_step": step,
        "post_collision_steps": post_step,
        "timeline": timeline,
    }


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "control_maps_result_v2.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    case_rows = [{key: value for key, value in row.items() if key != "timeline"} for row in result["cases"]]
    csv_path = OUT_DIR / "control_maps_cases_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(case_rows[0].keys()))
        writer.writeheader()
        writer.writerows(case_rows)

    timeline_rows: List[Dict[str, object]] = []
    for row in result["cases"]:
        timeline_rows.extend(row["timeline"])
    timeline_path = OUT_DIR / "control_maps_timeline_v2.csv"
    with timeline_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timeline_rows[0].keys()))
        writer.writeheader()
        writer.writerows(timeline_rows)

    lines = [
        "# Control Maps Result v1",
        "",
        "## Verdict",
        "",
        f"- tested_maps: `{result['summary']['tested_maps']}`",
        f"- reflection_valid_maps: `{result['summary']['reflection_valid_maps']}`",
        "",
        "## Cases",
        "",
        "| map | q_reversed | identity_labels_preserved | identity_labels_swapped | left_slot_mode | right_slot_mode | reflection_spatial_pattern | transmission_spatial_pattern | reflection_valid |",
        "|---|---|---|---|---:|---:|---|---|---|",
    ]
    for row in result["cases"]:
        lines.append(
            f"| {row['map']} | {row['q_reversed']} | {row['identity_labels_preserved']} | "
            f"{row['identity_labels_swapped']} | {row['left_slot_mode']} | {row['right_slot_mode']} | "
            f"{row['reflection_spatial_pattern']} | {row['transmission_spatial_pattern']} | {row['reflection_valid']} |"
        )
    (OUT_DIR / "control_maps_report_v2.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 5))
    for row in result["cases"]:
        timeline = row["timeline"]
        steps = [entry["step"] for entry in timeline]
        ax.plot(steps, [entry["chi_A"] for entry in timeline], label=f"{row['map']} A")
        ax.plot(steps, [entry["chi_B"] for entry in timeline], linestyle="--", label=f"{row['map']} B")
    ax.set_xlabel("step")
    ax.set_ylabel("position phase")
    ax.set_title("Control map trajectories")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "control_maps_trajectories_v2.png", dpi=180)
    plt.close(fig)

    labels = [row["map"] for row in result["cases"]]
    x = range(len(labels))
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar([i - 0.2 for i in x], [1 if row["q_reversed"] else 0 for row in result["cases"]], width=0.2, label="q reversed")
    ax.bar(x, [1 if row["reflection_spatial_pattern"] else 0 for row in result["cases"]], width=0.2, label="reflection spatial pattern")
    ax.bar([i + 0.2 for i in x], [1 if row["reflection_valid"] else 0 for row in result["cases"]], width=0.2, label="reflection valid")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_yticks([0, 1])
    ax.set_ylabel("pass/fail")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "control_maps_verdict_v2.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    params = Params()
    grid = np.linspace(-math.pi, math.pi, params.grid_n, endpoint=False)
    map_names = [
        "reflection",
        "transmission",
        "label_exchange_reflection",
        "transmission_with_label_exchange",
    ]
    cases = [run_case(params, map_name, grid) for map_name in map_names]
    result = {
        "parameters": asdict(params),
        "summary": {
            "tested_maps": len(cases),
            "reflection_valid_maps": [row["map"] for row in cases if row["reflection_valid"]],
        },
        "cases": cases,
    }
    return result


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print(f"result_dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
