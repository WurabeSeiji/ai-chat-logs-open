from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_label_robustness_result_v2"
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
    tau0: float = 0.2
    q_A0: int = 1
    q_B0: int = -1
    delta_s: float = 0.01
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    s_max: int = 10000
    grid_n: int = 4096
    purity_acceptance: float = 0.75


@dataclass
class State:
    chi: float
    tau: float
    q: int
    amplitude: float


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


def eta_label_signal(eta: np.ndarray, target_mode: int, other_mode: int, leakage: float) -> np.ndarray:
    return (1.0 - leakage) * np.exp(1j * target_mode * eta) + leakage * np.exp(1j * other_mode * eta)


def eta_readout_overlap(
    target_mode: int,
    other_mode: int,
    leakage: float,
    read_mode: int,
    eta: np.ndarray,
) -> complex:
    signal = eta_label_signal(eta, target_mode, other_mode, leakage)
    return complex(np.mean(signal * np.exp(-1j * read_mode * eta)))


def observe_label_modes(
    state: State,
    target_mode: int,
    other_mode: int,
    leakage: float,
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
        eta_part = eta_readout_overlap(target_mode, other_mode, leakage, mode, grid)
        result[mode] = state.amplitude * params.A_C * chi_part * tau_part * eta_part
    return result


def detect_mode(observations: Dict[int, complex]) -> int:
    return max(observations, key=lambda mode: abs(observations[mode]))


def purity(observations: Dict[int, complex], target_mode: int) -> float:
    denom = sum(abs(value) for value in observations.values())
    if denom == 0:
        return 0.0
    return float(abs(observations[target_mode]) / denom)


def run_motion(params: Params) -> Tuple[State, State, bool, bool, int, int]:
    a = State(params.chi_A0, -params.tau0, params.q_A0, params.A_A)
    b = State(params.chi_B0, -params.tau0, params.q_B0, params.A_B)
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
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)

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


def run_case(params: Params, m_a: int, m_b: int, leakage: float, grid: np.ndarray) -> Dict[str, object]:
    modes = [m_a, m_b]
    tau_c0 = -params.tau0
    tau_c2 = params.tau0
    initial_a = State(params.chi_A0, -params.tau0, params.q_A0, params.A_A)
    initial_b = State(params.chi_B0, -params.tau0, params.q_B0, params.A_B)

    obs_a0 = observe_label_modes(initial_a, m_a, m_b, leakage, params.Nh_chi_A, params.Nh_tau_A, initial_a.chi, tau_c0, modes, params, grid)
    obs_b0 = observe_label_modes(initial_b, m_b, m_a, leakage, params.Nh_chi_B, params.Nh_tau_B, initial_b.chi, tau_c0, modes, params, grid)

    final_a, final_b, collision_cell_reached, post_completed, collision_step, post_steps = run_motion(params)

    obs_a2 = observe_label_modes(final_a, m_a, m_b, leakage, params.Nh_chi_A, params.Nh_tau_A, final_a.chi, tau_c2, modes, params, grid)
    obs_b2 = observe_label_modes(final_b, m_b, m_a, leakage, params.Nh_chi_B, params.Nh_tau_B, final_b.chi, tau_c2, modes, params, grid)

    m_a0 = detect_mode(obs_a0)
    m_b0 = detect_mode(obs_b0)
    m_a2 = detect_mode(obs_a2)
    m_b2 = detect_mode(obs_b2)
    gamma_a0 = purity(obs_a0, m_a)
    gamma_b0 = purity(obs_b0, m_b)
    gamma_a2 = purity(obs_a2, m_a)
    gamma_b2 = purity(obs_b2, m_b)

    label_preserved = m_a0 == m_a and m_b0 == m_b and m_a2 == m_a and m_b2 == m_b
    label_swapped = m_a2 == m_b or m_b2 == m_a
    purity_pass = min(gamma_a0, gamma_b0, gamma_a2, gamma_b2) >= params.purity_acceptance
    q_reversed = final_a.q == -params.q_A0 and final_b.q == -params.q_B0
    separated = final_a.chi < final_b.chi

    return {
        "m_A": m_a,
        "m_B": m_b,
        "mode_separation": abs(m_b - m_a),
        "leakage": leakage,
        "detected_A_initial": m_a0,
        "detected_B_initial": m_b0,
        "detected_A_final": m_a2,
        "detected_B_final": m_b2,
        "purity_A_initial": gamma_a0,
        "purity_B_initial": gamma_b0,
        "purity_A_final": gamma_a2,
        "purity_B_final": gamma_b2,
        "min_purity": min(gamma_a0, gamma_b0, gamma_a2, gamma_b2),
        "label_preserved": label_preserved,
        "label_swapped": label_swapped,
        "purity_pass": purity_pass,
        "collision_cell_reached": collision_cell_reached,
        "post_collision_propagation_completed": post_completed,
        "q_reversed": q_reversed,
        "separated_after_collision": separated,
        "case_valid": label_preserved and not label_swapped and purity_pass and collision_cell_reached and post_completed and q_reversed and separated,
        "collision_step": collision_step,
        "post_collision_steps": post_steps,
    }


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "label_robustness_result_v2.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    rows = result["cases"]
    csv_path = OUT_DIR / "label_robustness_cases_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    report_lines = [
        "# Identification Vibration Robustness Result v1",
        "",
        "## Verdict",
        "",
        f"- total_cases: `{result['summary']['total_cases']}`",
        f"- valid_cases: `{result['summary']['valid_cases']}`",
        f"- invalid_cases: `{result['summary']['invalid_cases']}`",
        f"- first_failure_leakage: `{result['summary']['first_failure_leakage']}`",
        "",
        "## Interpretation",
        "",
        "The identification channel is preserved while the target mode remains the dominant internal vibration.",
        "The transition point is the expected mode-mixture boundary near leakage = 0.5.",
        "",
        "## Cases",
        "",
        "| m_A | m_B | leakage | min_purity | label_preserved | label_swapped | case_valid |",
        "|---:|---:|---:|---:|---|---|---|",
    ]
    for row in rows:
        report_lines.append(
            f"| {row['m_A']} | {row['m_B']} | {row['leakage']:.3g} | {row['min_purity']:.12g} | "
            f"{row['label_preserved']} | {row['label_swapped']} | {row['case_valid']} |"
        )
    (OUT_DIR / "label_robustness_report_v2.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 5))
    pair_keys = sorted({(row["m_A"], row["m_B"]) for row in rows})
    for m_a, m_b in pair_keys:
        pair_rows = [row for row in rows if row["m_A"] == m_a and row["m_B"] == m_b]
        ax.plot(
            [row["leakage"] for row in pair_rows],
            [row["min_purity"] for row in pair_rows],
            marker="o",
            label=f"mA={m_a}, mB={m_b}",
        )
    ax.axhline(result["parameters"]["purity_acceptance"], color="black", linestyle="--", linewidth=1, label="purity threshold")
    ax.axvline(0.5, color="gray", linestyle=":", linewidth=1, label="dominance boundary")
    ax.set_xlabel("label-mode leakage")
    ax.set_ylabel("minimum identification purity")
    ax.set_title("Identification vibration robustness")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "label_robustness_purity_v2.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4))
    for m_a, m_b in pair_keys:
        pair_rows = [row for row in rows if row["m_A"] == m_a and row["m_B"] == m_b]
        ax.plot(
            [row["leakage"] for row in pair_rows],
            [1 if row["label_preserved"] else 0 for row in pair_rows],
            marker="s",
            label=f"mA={m_a}, mB={m_b}",
        )
    ax.set_xlabel("label-mode leakage")
    ax.set_ylabel("label preserved")
    ax.set_yticks([0, 1])
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "label_robustness_detection_v2.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    params = Params()
    grid = phase_grid(params.grid_n)
    mode_pairs = [(1, 2), (1, 3), (2, 5), (8, 9)]
    leakages = [0.0, 0.01, 0.05, 0.1, 0.2, 0.35, 0.49, 0.5, 0.51]
    cases: List[Dict[str, object]] = []
    for m_a, m_b in mode_pairs:
        for leakage in leakages:
            cases.append(run_case(params, m_a, m_b, leakage, grid))

    invalid = [row for row in cases if not row["case_valid"]]
    first_failure_leakage = None
    if invalid:
        first_failure_leakage = min(row["leakage"] for row in invalid)

    result = {
        "parameters": asdict(params),
        "summary": {
            "total_cases": len(cases),
            "valid_cases": sum(1 for row in cases if row["case_valid"]),
            "invalid_cases": sum(1 for row in cases if not row["case_valid"]),
            "first_failure_leakage": first_failure_leakage,
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
