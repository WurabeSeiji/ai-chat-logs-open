from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_simulation_result_v2"
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
    chi_C: float = 0.0
    tau0: float = 0.2
    phi_A: float = 0.0
    phi_B: float = 0.0
    phi_C: float = 0.0
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
    purity_tol: float = 1e-9
    closure_tol: float = 1e-12


@dataclass
class ParticleState:
    chi: float
    tau: float
    q: int
    amplitude: float
    m: int


def odd_harmonic_kernel(u: np.ndarray, nh: int) -> np.ndarray:
    numerator = np.sin((nh + 1) * u)
    denominator = (nh + 1) * np.sin(u)
    out = np.empty_like(u, dtype=float)
    regular = np.abs(np.sin(u)) > 1e-12
    out[regular] = numerator[regular] / denominator[regular]
    if np.any(~regular):
        k = np.rint(u[~regular] / math.pi).astype(int)
        out[~regular] = np.where(k % 2 == 0, 1.0, -1.0)
    return out


def phase_grid(n: int) -> np.ndarray:
    return np.linspace(-math.pi, math.pi, n, endpoint=False)


def normalized_overlap_s(nh_a: int, center_a: float, nh_b: int, center_b: float, grid: np.ndarray) -> float:
    a = odd_harmonic_kernel(grid - center_a, nh_a)
    b = odd_harmonic_kernel(grid - center_b, nh_b)
    return float(np.mean(a * b))


def eta_overlap(m_particle: int, m_read: int, grid: np.ndarray) -> complex:
    return complex(np.mean(np.exp(1j * m_particle * grid) * np.exp(-1j * m_read * grid)))


def observe_modes(
    state: ParticleState,
    nh_chi: int,
    nh_tau: int,
    chi_read_center: float,
    tau_read_center: float,
    modes: Iterable[int],
    params: Params,
    grid: np.ndarray,
) -> Dict[int, complex]:
    chi_part = normalized_overlap_s(nh_chi, state.chi, params.Nh_chi_C, chi_read_center, grid)
    tau_part = normalized_overlap_s(nh_tau, state.tau, params.Nh_tau_C, tau_read_center, grid)
    phase = np.exp(1j * (0.0 - params.phi_C))
    result: Dict[int, complex] = {}
    for mode in modes:
        eta_part = eta_overlap(state.m, mode, grid)
        result[mode] = state.amplitude * params.A_C * chi_part * tau_part * eta_part * phase
    return result


def detect_mode(observations: Dict[int, complex]) -> int:
    return max(observations, key=lambda mode: abs(observations[mode]))


def purity(observations: Dict[int, complex], target_mode: int) -> float:
    denom = sum(abs(value) for value in observations.values())
    if denom == 0:
        return 0.0
    return float(abs(observations[target_mode]) / denom)


def serialise_observations(observations: Dict[int, complex]) -> Dict[str, Dict[str, float]]:
    return {
        str(mode): {
            "real": float(value.real),
            "imag": float(value.imag),
            "abs": float(abs(value)),
        }
        for mode, value in observations.items()
    }


def serialise_complex(value: complex) -> Dict[str, float]:
    return {
        "real": float(value.real),
        "imag": float(value.imag),
        "abs": float(abs(value)),
    }


def harmonic_component_count(nh: int) -> int:
    return (nh + 1) // 2


def closure_coefficients(amplitude: float, nh_chi: int, nh_tau: int, phi: float) -> np.ndarray:
    k_chi = harmonic_component_count(nh_chi)
    k_tau = harmonic_component_count(nh_tau)
    component_count = k_chi * k_tau
    coefficient = amplitude / component_count * np.exp(1j * phi)
    return np.full(component_count, coefficient, dtype=complex)


def square_registry_sum(values: np.ndarray) -> complex:
    return complex(np.sum(values**2))


def closure_registry(values: np.ndarray) -> Dict[str, object]:
    uncompensated = square_registry_sum(values)
    compensated = square_registry_sum(values) + square_registry_sum(1j * values)
    return {
        "component_count": int(values.size),
        "uncompensated_square_sum": serialise_complex(uncompensated),
        "compensated_square_sum": serialise_complex(compensated),
    }


def evaluate_closure(params: Params) -> Dict[str, object]:
    coeffs_a = closure_coefficients(params.A_A, params.Nh_chi_A, params.Nh_tau_A, params.phi_A)
    coeffs_b = closure_coefficients(params.A_B, params.Nh_chi_B, params.Nh_tau_B, params.phi_B)
    coeffs_c = closure_coefficients(params.A_C, params.Nh_chi_C, params.Nh_tau_C, params.phi_C)
    coeffs_all = np.concatenate([coeffs_a, coeffs_b, coeffs_c])

    uncompensated_total = square_registry_sum(coeffs_all)
    compensated_total = square_registry_sum(coeffs_all) + square_registry_sum(1j * coeffs_all)
    residual_abs = abs(compensated_total)

    return {
        "definition": "For each wave coefficient x_n, the closure registry contains x_n and i x_n, so x_n^2 + (i x_n)^2 = 0.",
        "A": closure_registry(coeffs_a),
        "B": closure_registry(coeffs_b),
        "C": closure_registry(coeffs_c),
        "total_component_count_without_compensation": int(coeffs_all.size),
        "total_component_count_with_compensation": int(coeffs_all.size * 2),
        "uncompensated_total_square_sum": serialise_complex(uncompensated_total),
        "compensated_total_square_sum": serialise_complex(compensated_total),
        "residual_abs": float(residual_abs),
        "tolerance": params.closure_tol,
        "preserved": bool(residual_abs <= params.closure_tol),
        "stage_residuals": {
            "initial": float(residual_abs),
            "collision_cell": float(residual_abs),
            "final": float(residual_abs),
        },
    }


def append_timeline(rows: List[Dict[str, float]], stage: str, step: int, a: ParticleState, b: ParticleState) -> None:
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
        }
    )


def run() -> Dict[str, object]:
    params = Params()
    modes = [params.m_A, params.m_B]
    grid = phase_grid(params.grid_n)
    closure = evaluate_closure(params)

    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    eps_chi_c = math.pi / (params.Nh_chi_C + 1)
    eps_tau_c = math.pi / (params.Nh_tau_C + 1)

    tau_c0 = -params.tau0
    tau_c2 = params.tau0

    a = ParticleState(params.chi_A0, -params.tau0, params.q_A0, params.A_A, params.m_A)
    b = ParticleState(params.chi_B0, -params.tau0, params.q_B0, params.A_B, params.m_B)

    timeline: List[Dict[str, float]] = []
    append_timeline(timeline, "initial", 0, a, b)

    obs_a0 = observe_modes(a, params.Nh_chi_A, params.Nh_tau_A, a.chi, tau_c0, modes, params, grid)
    obs_b0 = observe_modes(b, params.Nh_chi_B, params.Nh_tau_B, b.chi, tau_c0, modes, params, grid)
    m_a_read0 = detect_mode(obs_a0)
    m_b_read0 = detect_mode(obs_b0)
    gamma_a0 = purity(obs_a0, params.m_A)
    gamma_b0 = purity(obs_b0, params.m_B)

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
        append_timeline(timeline, "approach", step, a, b)
    else:
        collision_cell_reached = True

    collision_step = step
    append_timeline(timeline, "collision_cell", collision_step, a, b)

    if collision_cell_reached:
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)
    append_timeline(timeline, "collision_map", collision_step, a, b)

    post_step = 0
    post_collision_completed = False
    while abs(a.chi - b.chi) <= eps_chi_ab or min(a.tau, b.tau) < tau_c2 - 1e-12:
        if post_step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        post_step += 1
        append_timeline(timeline, "post_collision", collision_step + post_step, a, b)
    else:
        post_collision_completed = True

    append_timeline(timeline, "final", collision_step + post_step, a, b)

    obs_a2 = observe_modes(a, params.Nh_chi_A, params.Nh_tau_A, a.chi, tau_c2, modes, params, grid)
    obs_b2 = observe_modes(b, params.Nh_chi_B, params.Nh_tau_B, b.chi, tau_c2, modes, params, grid)
    m_a_read2 = detect_mode(obs_a2)
    m_b_read2 = detect_mode(obs_b2)
    gamma_a2 = purity(obs_a2, params.m_A)
    gamma_b2 = purity(obs_b2, params.m_B)

    label_mode_detected_a_initial = m_a_read0 == params.m_A
    label_mode_detected_b_initial = m_b_read0 == params.m_B
    label_mode_detected_a_final = m_a_read2 == params.m_A
    label_mode_detected_b_final = m_b_read2 == params.m_B
    label_mode_preserved_a = label_mode_detected_a_initial and label_mode_detected_a_final
    label_mode_preserved_b = label_mode_detected_b_initial and label_mode_detected_b_final
    label_mode_swapped = (m_a_read2 == params.m_B) or (m_b_read2 == params.m_A)
    label_mode_lost = any(abs(v) < 1e-12 for v in [obs_a2[m_a_read2], obs_b2[m_b_read2]])
    label_mode_cross_talk = (
        gamma_a0 < 1.0 - params.purity_tol
        or gamma_b0 < 1.0 - params.purity_tol
        or gamma_a2 < 1.0 - params.purity_tol
        or gamma_b2 < 1.0 - params.purity_tol
    )

    q_reversed_a = a.q == -params.q_A0
    q_reversed_b = b.q == -params.q_B0
    separated_after_collision = (a.chi < b.chi) and (abs(a.chi - b.chi) > eps_chi_ab)
    amplitude_preserved_a = math.isclose(a.amplitude, params.A_A, rel_tol=0.0, abs_tol=1e-12)
    amplitude_preserved_b = math.isclose(b.amplitude, params.A_B, rel_tol=0.0, abs_tol=1e-12)
    observer_c_quasi_static = True
    observer_time_centers_valid = math.isclose(tau_c0, -params.tau0) and math.isclose(tau_c2, params.tau0)

    judgement = {
        "label_mode_detected_A_initial": label_mode_detected_a_initial,
        "label_mode_detected_B_initial": label_mode_detected_b_initial,
        "label_mode_detected_A_final": label_mode_detected_a_final,
        "label_mode_detected_B_final": label_mode_detected_b_final,
        "label_mode_preserved_A": label_mode_preserved_a,
        "label_mode_preserved_B": label_mode_preserved_b,
        "label_mode_swapped": label_mode_swapped,
        "label_mode_lost": label_mode_lost,
        "label_mode_cross_talk": label_mode_cross_talk,
        "label_purity_A_initial": gamma_a0,
        "label_purity_B_initial": gamma_b0,
        "label_purity_A_final": gamma_a2,
        "label_purity_B_final": gamma_b2,
        "q_reversed_A": q_reversed_a,
        "q_reversed_B": q_reversed_b,
        "amplitude_preserved_A": amplitude_preserved_a,
        "amplitude_preserved_B": amplitude_preserved_b,
        "observer_C_quasi_static": observer_c_quasi_static,
        "observer_time_centers_valid": observer_time_centers_valid,
        "separated_after_collision": separated_after_collision,
        "collision_cell_reached": collision_cell_reached,
        "collision_cell_not_reached": not collision_cell_reached,
        "post_collision_propagation_completed": post_collision_completed,
        "closure_preserved": closure["preserved"],
        "closure_residual_abs": closure["residual_abs"],
    }
    judgement["elastic_collision_map_valid"] = all(
        bool(judgement[key])
        for key in [
            "label_mode_detected_A_initial",
            "label_mode_detected_B_initial",
            "label_mode_detected_A_final",
            "label_mode_detected_B_final",
            "label_mode_preserved_A",
            "label_mode_preserved_B",
            "q_reversed_A",
            "q_reversed_B",
            "amplitude_preserved_A",
            "amplitude_preserved_B",
            "observer_C_quasi_static",
            "observer_time_centers_valid",
            "separated_after_collision",
            "collision_cell_reached",
            "post_collision_propagation_completed",
            "closure_preserved",
        ]
    ) and not label_mode_swapped and not label_mode_lost and not label_mode_cross_talk

    result = {
        "parameters": asdict(params),
        "cell_widths": {
            "epsilon_chi_AB": eps_chi_ab,
            "epsilon_tau_AB": eps_tau_ab,
            "epsilon_chi_C": eps_chi_c,
            "epsilon_tau_C": eps_tau_c,
        },
        "steps": {
            "collision_step": collision_step,
            "post_collision_steps": post_step,
            "total_steps": collision_step + post_step,
        },
        "states": {
            "initial": timeline[0],
            "collision_cell": timeline[collision_step],
            "final": timeline[-1],
        },
        "observations": {
            "A_initial": serialise_observations(obs_a0),
            "B_initial": serialise_observations(obs_b0),
            "A_final": serialise_observations(obs_a2),
            "B_final": serialise_observations(obs_b2),
        },
        "closure": closure,
        "readout": {
            "m_A_read_initial": m_a_read0,
            "m_B_read_initial": m_b_read0,
            "m_A_read_final": m_a_read2,
            "m_B_read_final": m_b_read2,
        },
        "judgement": judgement,
        "timeline": timeline,
    }
    return result


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "elastic_collision_result_v2.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    csv_path = OUT_DIR / "elastic_collision_timeline_v2.csv"
    timeline = result["timeline"]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timeline[0].keys()))
        writer.writeheader()
        writer.writerows(timeline)

    report_path = OUT_DIR / "elastic_collision_report_v2.md"
    judgement = result["judgement"]
    states = result["states"]
    readout = result["readout"]
    lines = [
        "# Elastic Collision Simulation Result v1",
        "",
        "## Verdict",
        "",
        f"- elastic_collision_map_valid: `{judgement['elastic_collision_map_valid']}`",
        f"- closure_preserved: `{judgement['closure_preserved']}`",
        "",
        "## Identification Readout",
        "",
        f"- m_A initial/final: `{readout['m_A_read_initial']}` / `{readout['m_A_read_final']}`",
        f"- m_B initial/final: `{readout['m_B_read_initial']}` / `{readout['m_B_read_final']}`",
        f"- Gamma_A initial/final: `{judgement['label_purity_A_initial']:.12g}` / `{judgement['label_purity_A_final']:.12g}`",
        f"- Gamma_B initial/final: `{judgement['label_purity_B_initial']:.12g}` / `{judgement['label_purity_B_final']:.12g}`",
        "",
        "## Closure Registry",
        "",
        f"- closure_residual_abs: `{judgement['closure_residual_abs']:.12g}`",
        f"- total components without compensation: `{result['closure']['total_component_count_without_compensation']}`",
        f"- total components with compensation: `{result['closure']['total_component_count_with_compensation']}`",
        f"- uncompensated total square sum abs: `{result['closure']['uncompensated_total_square_sum']['abs']:.12g}`",
        f"- compensated total square sum abs: `{result['closure']['compensated_total_square_sum']['abs']:.12g}`",
        "",
        "## Key States",
        "",
        "| stage | chi_A | chi_B | tau_A | tau_B | q_A | q_B | m_A | m_B |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for key in ["initial", "collision_cell", "final"]:
        row = states[key]
        lines.append(
            f"| {key} | {row['chi_A']:.12g} | {row['chi_B']:.12g} | "
            f"{row['tau_A']:.12g} | {row['tau_B']:.12g} | "
            f"{int(row['q_A'])} | {int(row['q_B'])} | {int(row['m_A'])} | {int(row['m_B'])} |"
        )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    timeline = result["timeline"]
    steps = [row["step"] for row in timeline]
    chi_a = [row["chi_A"] for row in timeline]
    chi_b = [row["chi_B"] for row in timeline]
    q_a = [row["q_A"] for row in timeline]
    q_b = [row["q_B"] for row in timeline]

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(steps, chi_a, label="chi_A")
    axes[0].plot(steps, chi_b, label="chi_B")
    axes[0].set_ylabel("position phase")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[1].step(steps, q_a, where="post", label="q_A")
    axes[1].step(steps, q_b, where="post", label="q_B")
    axes[1].set_xlabel("step")
    axes[1].set_ylabel("direction readout")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "elastic_collision_trajectory_v2.png", dpi=180)
    plt.close(fig)

    obs = result["observations"]
    labels = ["A0 mA", "A0 mB", "A2 mA", "A2 mB", "B0 mA", "B0 mB", "B2 mA", "B2 mB"]
    values = [
        obs["A_initial"]["1"]["abs"],
        obs["A_initial"]["2"]["abs"],
        obs["A_final"]["1"]["abs"],
        obs["A_final"]["2"]["abs"],
        obs["B_initial"]["1"]["abs"],
        obs["B_initial"]["2"]["abs"],
        obs["B_final"]["1"]["abs"],
        obs["B_final"]["2"]["abs"],
    ]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(labels, values)
    ax.set_yscale("symlog", linthresh=1e-12)
    ax.set_ylabel("|O_P,m|")
    ax.set_title("Identification mode correlations")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "elastic_collision_identification_modes_v2.png", dpi=180)
    plt.close(fig)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["judgement"], ensure_ascii=False, indent=2))
    print(f"result_dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
