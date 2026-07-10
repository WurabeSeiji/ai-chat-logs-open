from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_observer_sweep_result_v1"
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
    Nh_chi_C: int = 999
    Nh_tau_A: int = 99
    Nh_tau_B: int = 99
    Nh_tau_C: int = 999
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
    closure_tol: float = 1e-12


@dataclass
class State:
    chi: float
    tau: float
    q: int
    amplitude: float
    m: int


def phase_grid(n: int) -> np.ndarray:
    return np.linspace(-math.pi, math.pi, n, endpoint=False)


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


def normalized_overlap_s(nh_a: int, center_a: float, nh_b: int, center_b: float, grid: np.ndarray) -> float:
    return float(np.mean(odd_harmonic_kernel(grid - center_a, nh_a) * odd_harmonic_kernel(grid - center_b, nh_b)))


def eta_overlap(m_particle: int, m_read: int, grid: np.ndarray) -> complex:
    return complex(np.mean(np.exp(1j * m_particle * grid) * np.exp(-1j * m_read * grid)))


def observe_modes(
    state: State,
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
    result: Dict[int, complex] = {}
    for mode in modes:
        result[mode] = state.amplitude * params.A_C * chi_part * tau_part * eta_overlap(state.m, mode, grid)
    return result


def detect_mode(observations: Dict[int, complex]) -> int:
    return max(observations, key=lambda mode: abs(observations[mode]))


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


def run_motion(params: Params) -> Tuple[State, State, bool, bool, int, int]:
    a = State(params.chi_A0, -params.tau0, params.q_A0, params.A_A, params.m_A)
    b = State(params.chi_B0, -params.tau0, params.q_B0, params.A_B, params.m_B)
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)

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
    else:
        collision_cell_reached = True

    if collision_cell_reached:
        a.q = -a.q
        b.q = -b.q

    post_step = 0
    post_collision_completed = False
    while abs(a.chi - b.chi) <= eps_chi_ab or min(a.tau, b.tau) < params.tau0 - 1e-12:
        if post_step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        post_step += 1
    else:
        post_collision_completed = True

    return a, b, collision_cell_reached, post_collision_completed, step, post_step


def run_case(params: Params, grid: np.ndarray) -> Dict[str, object]:
    modes = [params.m_A, params.m_B]
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    eps_chi_c = math.pi / (params.Nh_chi_C + 1)
    eps_tau_c = math.pi / (params.Nh_tau_C + 1)

    loop_gain = (params.A_A + params.A_B) / params.A_C
    delta_chi_c_est = loop_gain * eps_chi_ab
    delta_tau_c_est = loop_gain * eps_tau_ab
    observer_c_quasi_static = (
        delta_chi_c_est <= eps_chi_c or math.isclose(delta_chi_c_est, eps_chi_c, rel_tol=1e-12, abs_tol=1e-15)
    ) and (
        delta_tau_c_est <= eps_tau_c or math.isclose(delta_tau_c_est, eps_tau_c, rel_tol=1e-12, abs_tol=1e-15)
    )

    initial_a = State(params.chi_A0, -params.tau0, params.q_A0, params.A_A, params.m_A)
    initial_b = State(params.chi_B0, -params.tau0, params.q_B0, params.A_B, params.m_B)
    obs_a0 = observe_modes(initial_a, params.Nh_chi_A, params.Nh_tau_A, initial_a.chi, -params.tau0, modes, params, grid)
    obs_b0 = observe_modes(initial_b, params.Nh_chi_B, params.Nh_tau_B, initial_b.chi, -params.tau0, modes, params, grid)

    final_a, final_b, collision_reached, post_completed, collision_step, post_steps = run_motion(params)
    obs_a2 = observe_modes(final_a, params.Nh_chi_A, params.Nh_tau_A, final_a.chi, params.tau0, modes, params, grid)
    obs_b2 = observe_modes(final_b, params.Nh_chi_B, params.Nh_tau_B, final_b.chi, params.tau0, modes, params, grid)

    label_preserved = detect_mode(obs_a0) == params.m_A and detect_mode(obs_b0) == params.m_B and detect_mode(obs_a2) == params.m_A and detect_mode(obs_b2) == params.m_B
    q_reversed = final_a.q == -params.q_A0 and final_b.q == -params.q_B0
    separated = final_a.chi < final_b.chi and abs(final_a.chi - final_b.chi) > eps_chi_ab
    residual = closure_residual(params)
    closure_preserved = residual <= params.closure_tol
    model_valid = all([label_preserved, q_reversed, separated, collision_reached, post_completed, observer_c_quasi_static, closure_preserved])

    return {
        "A_C": params.A_C,
        "R_C": params.A_C,
        "loop_gain": loop_gain,
        "epsilon_chi_AB": eps_chi_ab,
        "epsilon_chi_C": eps_chi_c,
        "delta_chi_C_est": delta_chi_c_est,
        "delta_tau_C_est": delta_tau_c_est,
        "observer_C_quasi_static": observer_c_quasi_static,
        "label_preserved": label_preserved,
        "q_reversed": q_reversed,
        "separated_after_collision": separated,
        "collision_cell_reached": collision_reached,
        "post_collision_propagation_completed": post_completed,
        "closure_preserved": closure_preserved,
        "closure_residual_abs": residual,
        "model_valid": model_valid,
        "collision_step": collision_step,
        "post_collision_steps": post_steps,
        "initial_signal_A_abs": float(abs(obs_a0[params.m_A])),
        "initial_signal_B_abs": float(abs(obs_b0[params.m_B])),
        "final_signal_A_abs": float(abs(obs_a2[params.m_A])),
        "final_signal_B_abs": float(abs(obs_b2[params.m_B])),
    }


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "observer_sweep_result_v1.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    rows = result["cases"]
    csv_path = OUT_DIR / "observer_sweep_cases_v1.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Observer C Condition Sweep Result v1",
        "",
        "## Verdict",
        "",
        f"- total_cases: `{result['summary']['total_cases']}`",
        f"- valid_cases: `{result['summary']['valid_cases']}`",
        f"- invalid_cases: `{result['summary']['invalid_cases']}`",
        f"- first_valid_A_C: `{result['summary']['first_valid_A_C']}`",
        "",
        "## Cases",
        "",
        "| A_C | loop_gain | delta_chi_C_est | epsilon_chi_C | observer_C_quasi_static | closure_preserved | model_valid |",
        "|---:|---:|---:|---:|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['A_C']:.12g} | {row['loop_gain']:.12g} | {row['delta_chi_C_est']:.12g} | "
            f"{row['epsilon_chi_C']:.12g} | {row['observer_C_quasi_static']} | "
            f"{row['closure_preserved']} | {row['model_valid']} |"
        )
    (OUT_DIR / "observer_sweep_report_v1.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 5))
    x = [row["A_C"] for row in rows]
    ax.plot(x, [row["loop_gain"] for row in rows], marker="o", label="G_loop")
    ax.plot(x, [row["delta_chi_C_est"] for row in rows], marker="s", label="estimated delta chi_C")
    ax.plot(x, [row["epsilon_chi_C"] for row in rows], linestyle="--", label="epsilon chi_C")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("observer amplitude A_C")
    ax.set_ylabel("condition value")
    ax.set_title("Observer C quasi-static condition")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "observer_sweep_conditions_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(x, [1 if row["observer_C_quasi_static"] else 0 for row in rows], marker="o", label="C quasi-static")
    ax.plot(x, [1 if row["model_valid"] else 0 for row in rows], marker="s", label="model valid")
    ax.set_xscale("log")
    ax.set_xlabel("observer amplitude A_C")
    ax.set_ylabel("pass/fail")
    ax.set_yticks([0, 1])
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "observer_sweep_validity_v1.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    base_params = Params()
    grid = phase_grid(base_params.grid_n)
    observer_amplitudes = [1, 2, 3, 5, 10, 15, 20, 30, 50, 100, 300, 1000, 3000]
    cases: List[Dict[str, object]] = []
    for amplitude in observer_amplitudes:
        cases.append(run_case(replace(base_params, A_C=float(amplitude)), grid))

    valid_cases = [row for row in cases if row["model_valid"]]
    result = {
        "parameters": asdict(base_params),
        "summary": {
            "total_cases": len(cases),
            "valid_cases": len(valid_cases),
            "invalid_cases": len(cases) - len(valid_cases),
            "first_valid_A_C": valid_cases[0]["A_C"] if valid_cases else None,
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
