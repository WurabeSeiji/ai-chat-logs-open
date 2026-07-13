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
OUT_DIR = BASE_DIR / "elastic_collision_observation_perturbation_result_v2"
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
    Nh_chi_C: int = 999
    Nh_tau_A: int = 99
    Nh_tau_B: int = 99
    Nh_tau_C: int = 999
    chi_A0: float = -0.2
    chi_B0: float = 0.2
    tau_A0: float = -0.2
    tau_B0: float = -0.2
    tau_target: float = 0.2
    q_A0: int = 1
    q_B0: int = -1
    delta_s: float = 0.001
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    s_max: int = 20000
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


@dataclass
class PerturbationCase:
    name: str
    delta_chi_A_factor: float
    delta_chi_B_factor: float
    delta_tau_A_factor: float
    delta_tau_B_factor: float


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


def append_timeline(rows: List[Dict[str, object]], case_name: str, stage: str, step: int, a: State, b: State) -> None:
    rows.append(
        {
            "case": case_name,
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
        }
    )


def in_collision_cell(a: State, b: State, eps_chi: float, eps_tau: float) -> bool:
    return abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau


def run_case(params: Params, case: PerturbationCase) -> Dict[str, object]:
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    eps_chi_c = math.pi / (params.Nh_chi_C + 1)
    eps_tau_c = math.pi / (params.Nh_tau_C + 1)

    delta_chi_A = case.delta_chi_A_factor * eps_chi_c
    delta_chi_B = case.delta_chi_B_factor * eps_chi_c
    delta_tau_A = case.delta_tau_A_factor * eps_tau_c
    delta_tau_B = case.delta_tau_B_factor * eps_tau_c

    observation_bound_respected = (
        abs(delta_chi_A) <= eps_chi_c
        and abs(delta_chi_B) <= eps_chi_c
        and abs(delta_tau_A) <= eps_tau_c
        and abs(delta_tau_B) <= eps_tau_c
    )

    a = State(params.chi_A0, params.tau_A0, params.q_A0, params.A_A, params.m_A)
    b = State(params.chi_B0, params.tau_B0, params.q_B0, params.A_B, params.m_B)
    timeline: List[Dict[str, object]] = []
    append_timeline(timeline, case.name, "initial", 0, a, b)

    a.chi += delta_chi_A
    b.chi += delta_chi_B
    a.tau += delta_tau_A
    b.tau += delta_tau_B
    append_timeline(timeline, case.name, "after_observation_perturbation", 0, a, b)

    step = 0
    collision_cell_reached = False
    spatial_crossed_without_time_cell = False
    previous_gap = b.chi - a.chi
    min_abs_chi_gap = abs(previous_gap)
    min_abs_tau_gap = abs(a.tau - b.tau)

    while step < params.s_max:
        if in_collision_cell(a, b, eps_chi_ab, eps_tau_ab):
            collision_cell_reached = True
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        step += 1
        append_timeline(timeline, case.name, "approach", step, a, b)

        current_gap = b.chi - a.chi
        min_abs_chi_gap = min(min_abs_chi_gap, abs(current_gap))
        min_abs_tau_gap = min(min_abs_tau_gap, abs(a.tau - b.tau))
        if in_collision_cell(a, b, eps_chi_ab, eps_tau_ab):
            collision_cell_reached = True
            break
        if previous_gap * current_gap < 0:
            spatial_crossed_without_time_cell = True
            break
        previous_gap = current_gap

    append_timeline(timeline, case.name, "collision_cell", step, a, b)

    if collision_cell_reached:
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)
    append_timeline(timeline, case.name, "collision_map", step, a, b)

    post_step = 0
    post_completed = False
    if collision_cell_reached:
        while post_step < params.s_max:
            separated = a.chi < b.chi and abs(a.chi - b.chi) > eps_chi_ab
            if separated and min(a.tau, b.tau) >= params.tau_target - 1e-12:
                post_completed = True
                break
            a.chi += a.q * params.v_chi * params.delta_s
            b.chi += b.q * params.v_chi * params.delta_s
            a.tau += params.omega_A * params.delta_s
            b.tau += params.omega_B * params.delta_s
            post_step += 1
            append_timeline(timeline, case.name, "post_collision", step + post_step, a, b)
    append_timeline(timeline, case.name, "final", step + post_step, a, b)

    q_reversed = a.q == -params.q_A0 and b.q == -params.q_B0
    label_preserved = a.m == params.m_A and b.m == params.m_B
    amplitude_preserved = math.isclose(a.amplitude, params.A_A) and math.isclose(b.amplitude, params.A_B)
    separated_after_collision = collision_cell_reached and a.chi < b.chi and abs(a.chi - b.chi) > eps_chi_ab
    residual = closure_residual(params)
    closure_preserved = residual <= params.closure_tol
    collision_map_valid = all(
        [
            collision_cell_reached,
            post_completed,
            q_reversed,
            label_preserved,
            amplitude_preserved,
            separated_after_collision,
            closure_preserved,
        ]
    )
    model_valid = observation_bound_respected and collision_map_valid

    return {
        "case": case.name,
        "delta_chi_A": delta_chi_A,
        "delta_chi_B": delta_chi_B,
        "delta_tau_A": delta_tau_A,
        "delta_tau_B": delta_tau_B,
        "max_chi_perturbation_over_epsilon_C": max(abs(delta_chi_A), abs(delta_chi_B)) / eps_chi_c,
        "max_tau_perturbation_over_epsilon_C": max(abs(delta_tau_A), abs(delta_tau_B)) / eps_tau_c,
        "tau_difference_after_perturbation": abs((params.tau_A0 + delta_tau_A) - (params.tau_B0 + delta_tau_B)),
        "tau_difference_over_epsilon_AB": abs((params.tau_A0 + delta_tau_A) - (params.tau_B0 + delta_tau_B)) / eps_tau_ab,
        "chi_difference_after_perturbation": abs((params.chi_A0 + delta_chi_A) - (params.chi_B0 + delta_chi_B)),
        "observation_bound_respected": observation_bound_respected,
        "epsilon_chi_C": eps_chi_c,
        "epsilon_tau_C": eps_tau_c,
        "epsilon_chi_AB": eps_chi_ab,
        "epsilon_tau_AB": eps_tau_ab,
        "min_abs_chi_gap": min_abs_chi_gap,
        "min_abs_tau_gap": min_abs_tau_gap,
        "collision_cell_reached": collision_cell_reached,
        "spatial_crossed_without_time_cell": spatial_crossed_without_time_cell,
        "collision_step": step,
        "post_collision_steps": post_step,
        "post_collision_propagation_completed": post_completed,
        "q_reversed": q_reversed,
        "label_preserved": label_preserved,
        "amplitude_preserved": amplitude_preserved,
        "separated_after_collision": separated_after_collision,
        "closure_preserved": closure_preserved,
        "closure_residual_abs": residual,
        "collision_map_valid": collision_map_valid,
        "model_valid": model_valid,
        "timeline": timeline,
    }


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "observation_perturbation_result_v2.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    case_rows = [{key: value for key, value in row.items() if key != "timeline"} for row in result["cases"]]
    csv_path = OUT_DIR / "observation_perturbation_cases_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(case_rows[0].keys()))
        writer.writeheader()
        writer.writerows(case_rows)

    timeline_rows: List[Dict[str, object]] = []
    for row in result["cases"]:
        timeline_rows.extend(row["timeline"])
    timeline_path = OUT_DIR / "observation_perturbation_timeline_v2.csv"
    with timeline_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timeline_rows[0].keys()))
        writer.writeheader()
        writer.writerows(timeline_rows)

    lines = [
        "# Observation Perturbation Result v1",
        "",
        "## Verdict",
        "",
        f"- total_cases: `{result['summary']['total_cases']}`",
        f"- model_valid_cases: `{result['summary']['model_valid_cases']}`",
        f"- collision_map_valid_cases: `{result['summary']['collision_map_valid_cases']}`",
        f"- invalid_cases: `{result['summary']['invalid_cases']}`",
        "",
        "## Cases",
        "",
        "| case | max_chi/eps_C | max_tau/eps_C | tau_diff/eps_AB | observation_bound_respected | collision_cell_reached | collision_map_valid | model_valid |",
        "|---|---:|---:|---:|---|---|---|---|",
    ]
    for row in result["cases"]:
        lines.append(
            f"| {row['case']} | {row['max_chi_perturbation_over_epsilon_C']:.12g} | "
            f"{row['max_tau_perturbation_over_epsilon_C']:.12g} | "
            f"{row['tau_difference_over_epsilon_AB']:.12g} | "
            f"{row['observation_bound_respected']} | {row['collision_cell_reached']} | "
            f"{row['collision_map_valid']} | {row['model_valid']} |"
        )
    (OUT_DIR / "observation_perturbation_report_v2.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    labels = [row["case"] for row in result["cases"]]
    x = list(range(len(labels)))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar([i - 0.25 for i in x], [1 if row["observation_bound_respected"] else 0 for row in result["cases"]], width=0.25, label="within C bound")
    ax.bar(x, [1 if row["collision_map_valid"] else 0 for row in result["cases"]], width=0.25, label="collision map valid")
    ax.bar([i + 0.25 for i in x], [1 if row["model_valid"] else 0 for row in result["cases"]], width=0.25, label="model valid")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticks([0, 1])
    ax.set_ylabel("pass/fail")
    ax.set_title("Observation perturbation verdict")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "observation_perturbation_verdict_v2.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for row in result["cases"]:
        ax.scatter(
            row["max_tau_perturbation_over_epsilon_C"],
            row["tau_difference_over_epsilon_AB"],
            c="green" if row["model_valid"] else ("orange" if row["collision_map_valid"] else "red"),
        )
        ax.annotate(row["case"], (row["max_tau_perturbation_over_epsilon_C"], row["tau_difference_over_epsilon_AB"]), fontsize=8)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1, label="C perturbation bound")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="AB time-cell bound")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("max tau perturbation / epsilon_tau_C")
    ax.set_ylabel("tau difference after perturbation / epsilon_tau_AB")
    ax.set_title("Observation perturbation thresholds")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "observation_perturbation_thresholds_v2.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    params = Params()
    cases = [
        PerturbationCase("none", 0.0, 0.0, 0.0, 0.0),
        PerturbationCase("common_spatial_within_C", 0.8, 0.8, 0.0, 0.0),
        PerturbationCase("differential_spatial_within_C", 1.0, -1.0, 0.0, 0.0),
        PerturbationCase("differential_time_within_C", 0.0, 0.0, 1.0, -1.0),
        PerturbationCase("mixed_corner_within_C", 1.0, -1.0, 1.0, -1.0),
        PerturbationCase("over_C_time_still_in_AB", 0.0, 0.0, 3.0, -3.0),
        PerturbationCase("over_C_time_breaks_AB", 0.0, 0.0, 8.0, -8.0),
        PerturbationCase("over_C_spatial_only", -10.0, 10.0, 0.0, 0.0),
    ]
    rows = [run_case(params, case) for case in cases]
    result = {
        "parameters": asdict(params),
        "summary": {
            "total_cases": len(rows),
            "model_valid_cases": sum(1 for row in rows if row["model_valid"]),
            "collision_map_valid_cases": sum(1 for row in rows if row["collision_map_valid"]),
            "invalid_cases": sum(1 for row in rows if not row["model_valid"]),
        },
        "cases": rows,
    }
    return result


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print(f"result_dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
