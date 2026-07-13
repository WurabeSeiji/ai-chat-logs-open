from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_asymmetry_sweep_result_v2"
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
    v_chi_A: float = 1.0
    v_chi_B: float = 1.0
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
            "A_A": a.amplitude,
            "A_B": b.amplitude,
            "m_A": a.m,
            "m_B": b.m,
        }
    )


def in_collision_cell(a: State, b: State, eps_chi: float, eps_tau: float) -> bool:
    return abs(a.chi - b.chi) < eps_chi and abs(a.tau - b.tau) < eps_tau


def run_case(case_name: str, params: Params) -> Dict[str, object]:
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    a = State(params.chi_A0, params.tau_A0, params.q_A0, params.A_A, params.m_A)
    b = State(params.chi_B0, params.tau_B0, params.q_B0, params.A_B, params.m_B)
    initial_amplitude_a = a.amplitude
    initial_amplitude_b = b.amplitude
    timeline: List[Dict[str, object]] = []
    append_timeline(timeline, case_name, "initial", 0, a, b)

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
        a.chi += a.q * params.v_chi_A * params.delta_s
        b.chi += b.q * params.v_chi_B * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        step += 1
        append_timeline(timeline, case_name, "approach", step, a, b)

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

    append_timeline(timeline, case_name, "collision_cell", step, a, b)

    if collision_cell_reached:
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)
    append_timeline(timeline, case_name, "collision_map", step, a, b)

    post_step = 0
    post_completed = False
    if collision_cell_reached:
        while post_step < params.s_max:
            separated = a.chi < b.chi and abs(a.chi - b.chi) > eps_chi_ab
            if separated and min(a.tau, b.tau) >= params.tau_target - 1e-12:
                post_completed = True
                break
            a.chi += a.q * params.v_chi_A * params.delta_s
            b.chi += b.q * params.v_chi_B * params.delta_s
            a.tau += params.omega_A * params.delta_s
            b.tau += params.omega_B * params.delta_s
            post_step += 1
            append_timeline(timeline, case_name, "post_collision", step + post_step, a, b)
    append_timeline(timeline, case_name, "final", step + post_step, a, b)

    q_reversed = a.q == -params.q_A0 and b.q == -params.q_B0
    label_preserved = a.m == params.m_A and b.m == params.m_B
    amplitude_preserved = math.isclose(a.amplitude, initial_amplitude_a) and math.isclose(b.amplitude, initial_amplitude_b)
    separated_after_collision = collision_cell_reached and a.chi < b.chi and abs(a.chi - b.chi) > eps_chi_ab
    residual = closure_residual(params)
    closure_preserved = residual <= params.closure_tol
    case_valid = all(
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

    return {
        "case": case_name,
        "A_A": params.A_A,
        "A_B": params.A_B,
        "Nh_chi_A": params.Nh_chi_A,
        "Nh_chi_B": params.Nh_chi_B,
        "Nh_tau_A": params.Nh_tau_A,
        "Nh_tau_B": params.Nh_tau_B,
        "tau_A0": params.tau_A0,
        "tau_B0": params.tau_B0,
        "omega_A": params.omega_A,
        "omega_B": params.omega_B,
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
        "case_valid": case_valid,
        "timeline": timeline,
    }


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "asymmetry_sweep_result_v2.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    case_rows = [{key: value for key, value in row.items() if key != "timeline"} for row in result["cases"]]
    csv_path = OUT_DIR / "asymmetry_sweep_cases_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(case_rows[0].keys()))
        writer.writeheader()
        writer.writerows(case_rows)

    timeline_rows: List[Dict[str, object]] = []
    for row in result["cases"]:
        timeline_rows.extend(row["timeline"])
    timeline_path = OUT_DIR / "asymmetry_sweep_timeline_v2.csv"
    with timeline_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timeline_rows[0].keys()))
        writer.writeheader()
        writer.writerows(timeline_rows)

    lines = [
        "# Asymmetry Sweep Result v1",
        "",
        "## Verdict",
        "",
        f"- total_cases: `{result['summary']['total_cases']}`",
        f"- valid_cases: `{result['summary']['valid_cases']}`",
        f"- invalid_cases: `{result['summary']['invalid_cases']}`",
        "",
        "## Cases",
        "",
        "| case | A_A | A_B | Nh_chi_A | Nh_chi_B | tau_gap_initial | omega_A | omega_B | min_tau_gap | collision_cell_reached | spatial_crossed_without_time_cell | case_valid |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for row in result["cases"]:
        lines.append(
            f"| {row['case']} | {row['A_A']:.12g} | {row['A_B']:.12g} | {row['Nh_chi_A']} | {row['Nh_chi_B']} | "
            f"{abs(row['tau_A0'] - row['tau_B0']):.12g} | {row['omega_A']:.12g} | {row['omega_B']:.12g} | "
            f"{row['min_abs_tau_gap']:.12g} | {row['collision_cell_reached']} | "
            f"{row['spatial_crossed_without_time_cell']} | {row['case_valid']} |"
        )
    (OUT_DIR / "asymmetry_sweep_report_v2.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    labels = [row["case"] for row in result["cases"]]
    x = list(range(len(labels)))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar([i - 0.2 for i in x], [1 if row["collision_cell_reached"] else 0 for row in result["cases"]], width=0.2, label="cell reached")
    ax.bar(x, [1 if row["q_reversed"] else 0 for row in result["cases"]], width=0.2, label="q reversed")
    ax.bar([i + 0.2 for i in x], [1 if row["case_valid"] else 0 for row in result["cases"]], width=0.2, label="case valid")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticks([0, 1])
    ax.set_ylabel("pass/fail")
    ax.set_title("Asymmetry sweep verdict")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "asymmetry_sweep_verdict_v2.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for row in result["cases"]:
        ax.scatter(
            row["min_abs_tau_gap"] / row["epsilon_tau_AB"],
            row["min_abs_chi_gap"] / row["epsilon_chi_AB"],
            c="green" if row["case_valid"] else "red",
            label=row["case"],
        )
        ax.annotate(row["case"], (row["min_abs_tau_gap"] / row["epsilon_tau_AB"], row["min_abs_chi_gap"] / row["epsilon_chi_AB"]), fontsize=8)
    ax.axvline(1, color="black", linestyle="--", linewidth=1)
    ax.axhline(1, color="black", linestyle=":", linewidth=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("min tau gap / epsilon_tau_AB")
    ax.set_ylabel("min chi gap / epsilon_chi_AB")
    ax.set_title("Simultaneous cell-entry condition")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "asymmetry_sweep_cell_gaps_v2.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    base = Params()
    cases = [
        ("symmetric_baseline", base),
        ("amplitude_asymmetry", replace(base, A_A=1.0, A_B=2.0)),
        ("spatial_harmonic_asymmetry", replace(base, Nh_chi_A=199, Nh_chi_B=99)),
        ("temporal_harmonic_asymmetry", replace(base, Nh_tau_A=199, Nh_tau_B=99)),
        ("small_time_phase_offset", replace(base, tau_A0=-0.2, tau_B0=-0.18)),
        ("large_time_phase_offset", replace(base, tau_A0=-0.2, tau_B0=-0.15)),
        ("small_omega_mismatch", replace(base, omega_A=1.0, omega_B=1.05)),
        ("large_omega_mismatch", replace(base, omega_A=1.0, omega_B=1.5)),
        ("combined_asymmetry_pass", replace(base, A_A=1.5, A_B=0.8, Nh_chi_A=199, Nh_chi_B=99, Nh_tau_A=199, Nh_tau_B=99, tau_A0=-0.2, tau_B0=-0.22, omega_A=1.0, omega_B=1.1)),
        ("combined_asymmetry_fail", replace(base, A_A=1.5, A_B=0.8, Nh_chi_A=199, Nh_chi_B=99, Nh_tau_A=199, Nh_tau_B=99, tau_A0=-0.2, tau_B0=-0.15, omega_A=1.0, omega_B=1.5)),
    ]
    rows = [run_case(name, params) for name, params in cases]
    result = {
        "parameters": asdict(base),
        "summary": {
            "total_cases": len(rows),
            "valid_cases": sum(1 for row in rows if row["case_valid"]),
            "invalid_cases": sum(1 for row in rows if not row["case_valid"]),
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
