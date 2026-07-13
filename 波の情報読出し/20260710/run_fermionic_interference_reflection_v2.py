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
OUT_DIR = BASE_DIR / "fermionic_interference_reflection_result_v2"
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
    m_A: int = 1
    m_B: int = 2
    q_A0: float = 1.0
    q_B0: float = -1.0
    delta_s: float = 0.01
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    s_max: int = 10000
    chi_grid_n: int = 1024
    phase_sweep_n: int = 33
    rho_grid_n: int = 4096
    rho_domain_length: float = 160.0
    packet_center: float = -12.0
    packet_k0: float = 2.4
    packet_sigma: float = 2.0
    final_time: float = 10.0
    time_sample_n: int = 41
    scattering_grid_n: int = 8192
    scattering_domain_length: float = 200.0
    scattering_packet_center: float = -40.0
    scattering_packet_k0: float = 5.0
    scattering_packet_sigma: float = 2.0
    scattering_window_inner: float = 25.0
    scattering_window_outer: float = 35.0
    node_tol: float = 1e-24
    dynamic_node_tol: float = 1e-10
    phase_error_tol: float = 1e-12
    transmission_left_prob_tol: float = 1e-2
    reflection_left_prob_tol: float = 1e-6
    scattering_error_tol: float = 1e-8
    reversibility_tol: float = 1e-10
    closure_tol: float = 1e-12


@dataclass
class LocalParticleState:
    chi: float
    tau: float
    q_readout: float
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


def exchange_diagonal_metrics(params: Params, delta_f: float, include_exchange: bool = True) -> Dict[str, float]:
    grid = phase_grid(params.chi_grid_n)
    wave_a = params.A_A * odd_harmonic_kernel(grid - params.chi_A0, params.Nh_chi_A)
    wave_b = params.A_B * odd_harmonic_kernel(grid - params.chi_B0, params.Nh_chi_B)

    direct_diag = wave_a * wave_b
    exchange_diag = wave_b * wave_a

    if include_exchange:
        psi_diag = (direct_diag + np.exp(1j * delta_f) * exchange_diag) / math.sqrt(2.0)
    else:
        psi_diag = direct_diag / math.sqrt(2.0)

    ref_norm = float(np.sum(np.abs(direct_diag) ** 2))
    diag_norm = float(np.sum(np.abs(psi_diag) ** 2))
    relative_norm = diag_norm / ref_norm if ref_norm else 0.0
    expected = abs(1.0 + np.exp(1j * delta_f)) ** 2 / 2.0 if include_exchange else 0.5

    return {
        "delta_f": float(delta_f),
        "diagonal_relative_norm": float(relative_norm),
        "expected_relative_norm": float(expected),
        "abs_error": float(abs(relative_norm - expected)),
        "diagonal_max_abs": float(np.max(np.abs(psi_diag))),
    }


def eta_overlap(m_particle: int, m_read: int, n: int = 4096) -> complex:
    eta = phase_grid(n)
    return complex(np.mean(np.exp(1j * m_particle * eta) * np.exp(-1j * m_read * eta)))


def eta_readout(m_particle: int, modes: List[int]) -> Dict[int, float]:
    return {mode: float(abs(eta_overlap(m_particle, mode))) for mode in modes}


def detect_mode(readout: Dict[int, float]) -> int:
    return max(readout, key=lambda mode: readout[mode])


def detect_abs_mode(readout: Dict[int, complex]) -> int:
    return max(readout, key=lambda mode: abs(readout[mode]))


def serialise_complex(value: complex) -> Dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag), "abs": float(abs(value))}


def rho_grid(params: Params) -> tuple[np.ndarray, float]:
    dx = params.rho_domain_length / params.rho_grid_n
    rho = (np.arange(params.rho_grid_n) - params.rho_grid_n // 2) * dx
    return rho, dx


def normalize(psi: np.ndarray, dx: float) -> np.ndarray:
    norm = math.sqrt(float(np.sum(np.abs(psi) ** 2) * dx))
    return psi / norm


def gaussian_packet(rho: np.ndarray, center: float, k0: float, sigma: float) -> np.ndarray:
    envelope = np.exp(-((rho - center) ** 2) / (4.0 * sigma**2))
    carrier = np.exp(1j * k0 * rho)
    return envelope * carrier


def evolve_free(psi0: np.ndarray, t: float, dx: float) -> np.ndarray:
    n = psi0.size
    k = 2.0 * math.pi * np.fft.fftfreq(n, d=dx)
    psi_ordered = np.fft.ifftshift(psi0)
    evolved_ordered = np.fft.ifft(np.fft.fft(psi_ordered) * np.exp(-0.5j * (k**2) * t))
    return np.fft.fftshift(evolved_ordered)


def current_density(psi: np.ndarray, dx: float) -> np.ndarray:
    dpsi = np.gradient(psi, dx)
    return np.imag(np.conj(psi) * dpsi)


def side_metrics(rho: np.ndarray, dx: float, psi: np.ndarray) -> Dict[str, float]:
    density = np.abs(psi) ** 2
    current = current_density(psi, dx)
    left = rho < 0.0
    right = rho > 0.0
    left_prob = float(np.sum(density[left]) * dx)
    right_prob = float(np.sum(density[right]) * dx)
    left_mean = float(np.sum(rho[left] * density[left]) * dx / left_prob) if left_prob else 0.0
    right_mean = float(np.sum(rho[right] * density[right]) * dx / right_prob) if right_prob else 0.0
    left_current = float(np.sum(current[left]) * dx)
    right_current = float(np.sum(current[right]) * dx)
    node_index = int(np.argmin(np.abs(rho)))
    return {
        "left_prob": left_prob,
        "right_prob": right_prob,
        "left_mean_rho": left_mean,
        "right_mean_rho": right_mean,
        "left_current": left_current,
        "right_current": right_current,
        "node_abs": float(abs(psi[node_index])),
        "node_current": float(current[node_index]),
    }


def normalized_overlap_s(nh_a: int, center_a: float, nh_b: int, center_b: float, grid: np.ndarray) -> float:
    a = odd_harmonic_kernel(grid - center_a, nh_a)
    b = odd_harmonic_kernel(grid - center_b, nh_b)
    return float(np.mean(a * b))


def observe_label_modes(
    state: LocalParticleState,
    nh_chi: int,
    nh_tau: int,
    chi_read_center: float,
    tau_read_center: float,
    modes: List[int],
    params: Params,
    grid: np.ndarray,
) -> Dict[int, complex]:
    chi_part = normalized_overlap_s(nh_chi, state.chi, params.Nh_chi_C, chi_read_center, grid)
    tau_part = normalized_overlap_s(nh_tau, state.tau, params.Nh_tau_C, tau_read_center, grid)
    result: Dict[int, complex] = {}
    for mode in modes:
        eta_part = eta_overlap(state.m, mode, grid.size)
        result[mode] = state.amplitude * params.A_C * chi_part * tau_part * eta_part
    return result


def purity_abs(readout: Dict[int, complex], target_mode: int) -> float:
    denom = sum(abs(value) for value in readout.values())
    return float(abs(readout[target_mode]) / denom) if denom else 0.0


def append_ab_timeline(
    rows: List[Dict[str, float | str]],
    stage: str,
    step: int,
    a: LocalParticleState,
    b: LocalParticleState,
) -> None:
    rows.append(
        {
            "stage": stage,
            "step": step,
            "chi_A": a.chi,
            "chi_B": b.chi,
            "tau_A": a.tau,
            "tau_B": b.tau,
            "q_A_readout": a.q_readout,
            "q_B_readout": b.q_readout,
            "m_A": a.m,
            "m_B": b.m,
        }
    )


def run_relative_dynamics(params: Params) -> List[Dict[str, float | str]]:
    rho, dx = rho_grid(params)
    single0 = gaussian_packet(rho, params.packet_center, params.packet_k0, params.packet_sigma)
    single0 = normalize(single0, dx)

    mirror0 = gaussian_packet(rho, -params.packet_center, -params.packet_k0, params.packet_sigma)
    mirror0 = normalize(mirror0, dx)
    odd0 = normalize(single0 - mirror0, dx)

    rows: List[Dict[str, float | str]] = []
    for t in np.linspace(0.0, params.final_time, params.time_sample_n):
        single = evolve_free(single0, float(t), dx)
        odd = evolve_free(odd0, float(t), dx)

        single_metrics = side_metrics(rho, dx, single)
        odd_metrics = side_metrics(rho, dx, odd)

        rows.append({"model": "single_free_packet", "time": float(t), **single_metrics})
        rows.append({"model": "odd_node_halfline_read", "time": float(t), **odd_metrics})

    return rows


def parity_partner(rho: np.ndarray, psi: np.ndarray) -> np.ndarray:
    partner_rho = -rho
    real = np.interp(partner_rho, rho, psi.real, left=0.0, right=0.0)
    imag = np.interp(partner_rho, rho, psi.imag, left=0.0, right=0.0)
    return real + 1j * imag


def interaction_window(rho: np.ndarray, inner: float, outer: float) -> np.ndarray:
    distance = np.abs(rho)
    window = np.zeros_like(rho)
    window[distance <= inner] = 1.0
    edge = (distance > inner) & (distance < outer)
    window[edge] = 0.5 * (1.0 + np.cos(math.pi * (distance[edge] - inner) / (outer - inner)))
    return window


def apply_even_odd_phase(
    rho: np.ndarray,
    psi: np.ndarray,
    delta_f: float,
    inner: float,
    outer: float,
) -> np.ndarray:
    partner = parity_partner(rho, psi)
    window = interaction_window(rho, inner, outer)
    local_delta = delta_f * window
    t = np.exp(0.5j * local_delta) * np.cos(0.5 * local_delta)
    r = -1j * np.exp(0.5j * local_delta) * np.sin(0.5 * local_delta)
    return t * psi + r * partner


def run_single_sided_scattering(params: Params) -> List[Dict[str, float]]:
    dx = params.scattering_domain_length / params.scattering_grid_n
    rho = (np.arange(params.scattering_grid_n) - params.scattering_grid_n // 2) * dx
    hit_time = -params.scattering_packet_center / params.scattering_packet_k0
    initial = gaussian_packet(
        rho,
        params.scattering_packet_center,
        params.scattering_packet_k0,
        params.scattering_packet_sigma,
    )
    initial = normalize(initial, dx)
    collision_state = evolve_free(initial, hit_time, dx)

    rows: List[Dict[str, float]] = []
    for delta_f in np.linspace(0.0, math.pi, params.phase_sweep_n):
        after_interaction = apply_even_odd_phase(
            rho,
            collision_state,
            float(delta_f),
            params.scattering_window_inner,
            params.scattering_window_outer,
        )
        final = evolve_free(after_interaction, hit_time, dx)
        metrics = side_metrics(rho, dx, final)
        norm_total = metrics["left_prob"] + metrics["right_prob"]
        expected_reflection = math.sin(float(delta_f) / 2.0) ** 2
        expected_transmission = math.cos(float(delta_f) / 2.0) ** 2
        rows.append(
            {
                "delta_f": float(delta_f),
                "reflection_rate": metrics["left_prob"],
                "transmission_rate": metrics["right_prob"],
                "expected_reflection": expected_reflection,
                "expected_transmission": expected_transmission,
                "reflection_abs_error": abs(metrics["left_prob"] - expected_reflection),
                "transmission_abs_error": abs(metrics["right_prob"] - expected_transmission),
                "norm_total": norm_total,
                "norm_abs_error": abs(norm_total - 1.0),
                "hit_time": hit_time,
                "final_time": 2.0 * hit_time,
            }
        )
    return rows


def l2_relative_error(a: np.ndarray, b: np.ndarray, dx: float) -> float:
    numerator = math.sqrt(float(np.sum(np.abs(a - b) ** 2) * dx))
    denominator = math.sqrt(float(np.sum(np.abs(b) ** 2) * dx))
    return numerator / denominator if denominator else 0.0


def compensated_square_residual(values: np.ndarray) -> Dict[str, object]:
    uncompensated = complex(np.sum(values**2))
    compensated = uncompensated + complex(np.sum((1j * values) ** 2))
    return {
        "uncompensated_square_sum": serialise_complex(uncompensated),
        "compensated_square_sum": serialise_complex(compensated),
        "residual_abs": float(abs(compensated)),
    }


def run_reversibility_and_closure(params: Params) -> Dict[str, object]:
    dx = params.scattering_domain_length / params.scattering_grid_n
    rho = (np.arange(params.scattering_grid_n) - params.scattering_grid_n // 2) * dx
    hit_time = -params.scattering_packet_center / params.scattering_packet_k0
    initial = gaussian_packet(
        rho,
        params.scattering_packet_center,
        params.scattering_packet_k0,
        params.scattering_packet_sigma,
    )
    initial = normalize(initial, dx)
    collision_state = evolve_free(initial, hit_time, dx)
    after_pi = apply_even_odd_phase(
        rho,
        collision_state,
        math.pi,
        params.scattering_window_inner,
        params.scattering_window_outer,
    )
    after_pi_twice = apply_even_odd_phase(
        rho,
        after_pi,
        math.pi,
        params.scattering_window_inner,
        params.scattering_window_outer,
    )
    pi_twice_error = l2_relative_error(after_pi_twice, collision_state, dx)

    inverse_rows: List[Dict[str, float]] = []
    for delta_f in np.linspace(0.0, math.pi, params.phase_sweep_n):
        forward = apply_even_odd_phase(
            rho,
            collision_state,
            float(delta_f),
            params.scattering_window_inner,
            params.scattering_window_outer,
        )
        backward = apply_even_odd_phase(
            rho,
            forward,
            -float(delta_f),
            params.scattering_window_inner,
            params.scattering_window_outer,
        )
        inverse_rows.append(
            {
                "delta_f": float(delta_f),
                "inverse_relative_error": l2_relative_error(backward, collision_state, dx),
            }
        )

    closure_stages = {
        "initial": compensated_square_residual(initial * math.sqrt(dx)),
        "collision_before_map": compensated_square_residual(collision_state * math.sqrt(dx)),
        "after_pi_map": compensated_square_residual(after_pi * math.sqrt(dx)),
        "after_pi_twice": compensated_square_residual(after_pi_twice * math.sqrt(dx)),
    }
    closure_max_residual = max(float(stage["residual_abs"]) for stage in closure_stages.values())

    return {
        "pi_twice_relative_error": float(pi_twice_error),
        "inverse_max_relative_error": float(max(row["inverse_relative_error"] for row in inverse_rows)),
        "inverse_rows": inverse_rows,
        "closure_stages": closure_stages,
        "closure_max_residual": float(closure_max_residual),
    }


def local_direction_readout(q_initial: float, reflection_rate: float, transmission_rate: float) -> float:
    return q_initial * (transmission_rate - reflection_rate)


def run_ab_c_replacement_test(params: Params, reflection_rate: float, transmission_rate: float) -> Dict[str, object]:
    modes = [params.m_A, params.m_B]
    grid = phase_grid(params.chi_grid_n)
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    tau_c0 = -params.tau0
    tau_c2 = params.tau0

    a = LocalParticleState(params.chi_A0, -params.tau0, params.q_A0, params.A_A, params.m_A)
    b = LocalParticleState(params.chi_B0, -params.tau0, params.q_B0, params.A_B, params.m_B)

    timeline: List[Dict[str, float | str]] = []
    append_ab_timeline(timeline, "initial", 0, a, b)

    obs_a0 = observe_label_modes(a, params.Nh_chi_A, params.Nh_tau_A, a.chi, tau_c0, modes, params, grid)
    obs_b0 = observe_label_modes(b, params.Nh_chi_B, params.Nh_tau_B, b.chi, tau_c0, modes, params, grid)
    m_a_initial = detect_abs_mode(obs_a0)
    m_b_initial = detect_abs_mode(obs_b0)
    gamma_a_initial = purity_abs(obs_a0, params.m_A)
    gamma_b_initial = purity_abs(obs_b0, params.m_B)

    step = 0
    collision_cell_reached = False
    while abs(a.chi - b.chi) >= eps_chi_ab or abs(a.tau - b.tau) >= eps_tau_ab:
        if step >= params.s_max:
            break
        a.chi += a.q_readout * params.v_chi * params.delta_s
        b.chi += b.q_readout * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        step += 1
        append_ab_timeline(timeline, "approach", step, a, b)
    else:
        collision_cell_reached = True

    collision_step = step
    append_ab_timeline(timeline, "collision_cell", collision_step, a, b)

    q_a_before = a.q_readout
    q_b_before = b.q_readout
    a.q_readout = local_direction_readout(q_a_before, reflection_rate, transmission_rate)
    b.q_readout = local_direction_readout(q_b_before, reflection_rate, transmission_rate)
    append_ab_timeline(timeline, "local_exchange_interference_map", collision_step, a, b)

    post_step = 0
    post_collision_completed = False
    while abs(a.chi - b.chi) <= eps_chi_ab or min(a.tau, b.tau) < tau_c2 - 1e-12:
        if post_step >= params.s_max:
            break
        a.chi += a.q_readout * params.v_chi * params.delta_s
        b.chi += b.q_readout * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        post_step += 1
        append_ab_timeline(timeline, "post_collision", collision_step + post_step, a, b)
    else:
        post_collision_completed = True

    append_ab_timeline(timeline, "final", collision_step + post_step, a, b)

    obs_a2 = observe_label_modes(a, params.Nh_chi_A, params.Nh_tau_A, a.chi, tau_c2, modes, params, grid)
    obs_b2 = observe_label_modes(b, params.Nh_chi_B, params.Nh_tau_B, b.chi, tau_c2, modes, params, grid)
    m_a_final = detect_abs_mode(obs_a2)
    m_b_final = detect_abs_mode(obs_b2)
    gamma_a_final = purity_abs(obs_a2, params.m_A)
    gamma_b_final = purity_abs(obs_b2, params.m_B)

    direction_generated_a = abs(a.q_readout + params.q_A0) < params.scattering_error_tol
    direction_generated_b = abs(b.q_readout + params.q_B0) < params.scattering_error_tol
    separated_after_collision = a.chi < b.chi and abs(a.chi - b.chi) > eps_chi_ab
    label_preserved_a = m_a_initial == params.m_A and m_a_final == params.m_A
    label_preserved_b = m_b_initial == params.m_B and m_b_final == params.m_B
    cross_talk = (
        gamma_a_initial < 1.0 - params.scattering_error_tol
        or gamma_b_initial < 1.0 - params.scattering_error_tol
        or gamma_a_final < 1.0 - params.scattering_error_tol
        or gamma_b_final < 1.0 - params.scattering_error_tol
    )

    verdict = {
        "external_q_flip_used": False,
        "collision_cell_reached": collision_cell_reached,
        "post_collision_completed": post_collision_completed,
        "direction_generated_A": bool(direction_generated_a),
        "direction_generated_B": bool(direction_generated_b),
        "label_mode_preserved_A": bool(label_preserved_a),
        "label_mode_preserved_B": bool(label_preserved_b),
        "label_mode_cross_talk": bool(cross_talk),
        "separated_after_collision": bool(separated_after_collision),
        "integrated_ab_c_replacement_valid": bool(
            collision_cell_reached
            and post_collision_completed
            and direction_generated_a
            and direction_generated_b
            and label_preserved_a
            and label_preserved_b
            and not cross_talk
            and separated_after_collision
        ),
    }

    return {
        "method": "q readouts are computed from local exchange-interference rates: q_out = q_in * (T - R).",
        "reflection_rate_used": float(reflection_rate),
        "transmission_rate_used": float(transmission_rate),
        "cell_widths": {
            "epsilon_chi_AB": eps_chi_ab,
            "epsilon_tau_AB": eps_tau_ab,
        },
        "steps": {
            "collision_step": collision_step,
            "post_collision_steps": post_step,
            "total_steps": collision_step + post_step,
        },
        "q_readouts": {
            "A_before": float(q_a_before),
            "B_before": float(q_b_before),
            "A_after": float(a.q_readout),
            "B_after": float(b.q_readout),
        },
        "readout": {
            "m_A_initial": m_a_initial,
            "m_B_initial": m_b_initial,
            "m_A_final": m_a_final,
            "m_B_final": m_b_final,
            "gamma_A_initial": gamma_a_initial,
            "gamma_B_initial": gamma_b_initial,
            "gamma_A_final": gamma_a_final,
            "gamma_B_final": gamma_b_final,
        },
        "timeline": timeline,
        "verdict": verdict,
    }


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_phase_plot(rows: List[Dict[str, float]], path: Path) -> None:
    delta = [row["delta_f"] for row in rows]
    actual = [row["diagonal_relative_norm"] for row in rows]
    expected = [row["expected_relative_norm"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(delta, actual, marker="o", label="direct + exchange")
    ax.plot(delta, expected, linestyle="--", label="expected")
    ax.set_xlabel("fermionic core phase delta_F")
    ax.set_ylabel("diagonal relative norm")
    ax.set_title("Exchange-interference node at delta_F = pi")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_dynamics_plot(rows: List[Dict[str, float | str]], path: Path) -> None:
    by_model: Dict[str, List[Dict[str, float | str]]] = {}
    for row in rows:
        by_model.setdefault(str(row["model"]), []).append(row)

    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for model, model_rows in by_model.items():
        time = [float(row["time"]) for row in model_rows]
        left_prob = [float(row["left_prob"]) for row in model_rows]
        left_current = [float(row["left_current"]) for row in model_rows]
        axes[0].plot(time, left_prob, marker="o", markersize=3, label=model)
        axes[1].plot(time, left_current, marker="o", markersize=3, label=model)

    axes[0].set_ylabel("left probability")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("left current")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.suptitle("Transmission control vs odd-node reflection readout")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_scattering_plot(rows: List[Dict[str, float]], path: Path) -> None:
    delta = [row["delta_f"] for row in rows]
    reflection = [row["reflection_rate"] for row in rows]
    transmission = [row["transmission_rate"] for row in rows]
    expected_reflection = [row["expected_reflection"] for row in rows]
    expected_transmission = [row["expected_transmission"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(delta, reflection, marker="o", label="R dynamic")
    ax.plot(delta, transmission, marker="o", label="T dynamic")
    ax.plot(delta, expected_reflection, linestyle="--", label="sin^2(delta_F/2)")
    ax.plot(delta, expected_transmission, linestyle="--", label="cos^2(delta_F/2)")
    ax.set_xlabel("fermionic core phase delta_F")
    ax.set_ylabel("rate")
    ax.set_title("Single-sided packet scattering by exchange-interference phase")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_ab_c_replacement_plot(rows: List[Dict[str, float | str]], path: Path) -> None:
    steps = [float(row["step"]) for row in rows]
    chi_a = [float(row["chi_A"]) for row in rows]
    chi_b = [float(row["chi_B"]) for row in rows]
    q_a = [float(row["q_A_readout"]) for row in rows]
    q_b = [float(row["q_B_readout"]) for row in rows]

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(steps, chi_a, label="chi_A")
    axes[0].plot(steps, chi_b, label="chi_B")
    axes[0].set_ylabel("position phase")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].step(steps, q_a, where="post", label="q_A readout")
    axes[1].step(steps, q_b, where="post", label="q_B readout")
    axes[1].set_xlabel("step")
    axes[1].set_ylabel("direction readout")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.suptitle("AB/C cell replacement by local exchange-interference map")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_report(result: Dict[str, object]) -> str:
    verdict = result["verdict"]
    controls = result["controls"]
    paths = result["outputs"]
    phase = result["phase_sweep_summary"]
    dynamics = result["relative_dynamics_summary"]
    scattering = result["single_sided_scattering_summary"]
    labels = result["identification_readout"]
    reversibility = result["reversibility_closure_summary"]
    ab = result["ab_c_replacement_summary"]
    ab_verdict = ab["verdict"]

    return f"""# フェルミオン的逆相核による干渉反射 実行結果 v2

## 目的

前回の完全弾性反射実験で手続き的に置いていた `q_A -> -q_A`, `q_B -> -q_B` を使わず、片側から入射する単一波束が、交換干渉チャンネルの位相だけで反射へ変換されるかを確認した。

## 判定

| 項目 | 結果 |
|---|---:|
| 外部 `q` 反転代入 | `{str(verdict["external_q_flip_used"]).lower()}` |
| 片側入射初期条件 | `{str(verdict["single_sided_initial_packet"]).lower()}` |
| `Delta_F=0` で通過 | `{str(verdict["dynamic_transmission_at_zero"]).lower()}` |
| `Delta_F=pi/2` で半反射 | `{str(verdict["dynamic_half_phase_split"]).lower()}` |
| `Delta_F=pi` で完全反射 | `{str(verdict["dynamic_reflection_at_pi"]).lower()}` |
| 局所写像位相掃引が期待式と一致 | `{str(verdict["dynamic_phase_sweep_matches_expected"]).lower()}` |
| 局所写像ノルム保存 | `{str(verdict["dynamic_norm_preserved"]).lower()}` |
| `Delta_F=pi` の交換干渉節 | `{str(verdict["static_node_at_pi"]).lower()}` |
| 位相掃引が解析式と一致 | `{str(verdict["phase_sweep_matches_expected"]).lower()}` |
| 交換経路が節形成に必要 | `{str(verdict["exchange_path_required"]).lower()}` |
| 識別振動 A 保存 | `{str(verdict["label_mode_preserved_A"]).lower()}` |
| 識別振動 B 保存 | `{str(verdict["label_mode_preserved_B"]).lower()}` |
| `U(pi)^2` 可逆性 | `{str(verdict["pi_twice_reversible"]).lower()}` |
| `U(delta)U(-delta)` 可逆性 | `{str(verdict["inverse_sweep_reversible"]).lower()}` |
| 補償付き二乗閉鎖 | `{str(verdict["compensated_square_closure_preserved"]).lower()}` |
| AB/C セル差し替え | `{str(verdict["ab_c_replacement_valid"]).lower()}` |
| 最小機構判定 | `{str(verdict["mechanism_valid_minimal"]).lower()}` |

## 片側入射からの局所交換干渉反射

初期状態は左側だけに置いた単一局在波束であり、鏡像初期条件は使っていない。自由伝播で相互作用領域へ到達させ、局所窓内で偶・奇チャンネルに分解し、奇チャンネルへ内部核位相 `Delta_F` だけを与えて再合成した。

本実装は、連続ハミルトニアンを細かい時間刻みで積分する方式ではない。自由伝播で相互作用領域へ到達した波束に、局所交換干渉写像を一度作用させ、その後ふたたび自由伝播させる方式である。

| 量 | 値 |
|---|---:|
| `Delta_F=0` 反射率 | `{scattering["reflection_at_zero"]:.16e}` |
| `Delta_F=0` 通過率 | `{scattering["transmission_at_zero"]:.16e}` |
| `Delta_F=pi/2` 反射率 | `{scattering["reflection_at_half_pi"]:.16e}` |
| `Delta_F=pi/2` 通過率 | `{scattering["transmission_at_half_pi"]:.16e}` |
| `Delta_F=pi` 反射率 | `{scattering["reflection_at_pi"]:.16e}` |
| `Delta_F=pi` 通過率 | `{scattering["transmission_at_pi"]:.16e}` |
| 位相掃引最大誤差 | `{scattering["max_rate_abs_error"]:.16e}` |
| ノルム最大誤差 | `{scattering["max_norm_abs_error"]:.16e}` |

![single-sided scattering]({paths["scattering_plot"]})

## 交換干渉節

直接経路と交換経路を

```math
\\Psi_\\Delta(1,2)=\\frac{{1}}{{\\sqrt2}}\\left[P_A(1)P_B(2)+e^{{i\\Delta_F}}P_A(2)P_B(1)\\right]
```

として合成した。完全重なり対角 `1=2` では、`Delta_F=pi` で対角ノルムが消える。

| 量 | 値 |
|---|---:|
| `Delta_F=pi` 対角相対ノルム | `{phase["pi_diagonal_relative_norm"]:.16e}` |
| 位相掃引最大誤差 | `{phase["max_abs_error"]:.16e}` |
| 交換あり `Delta_F=pi` | `{controls["with_exchange_pi_diagonal_relative_norm"]:.16e}` |
| 交換なし `Delta_F=pi` | `{controls["without_exchange_pi_diagonal_relative_norm"]:.16e}` |

![phase sweep]({paths["phase_plot"]})

## 相対位相方向の反射読出し

対照として単一自由波束を同じ相対位相座標で進めると、左半線確率はほぼゼロへ移り、通過する。一方、逆相核に対応する奇関数節を置いた波は、原点節を保ったまま左半線読出しで戻り、左向き電流へ反転する。

| 量 | 値 |
|---|---:|
| 単一自由波束の最終左確率 | `{dynamics["single_final_left_prob"]:.16e}` |
| 奇関数節の最終左確率 | `{dynamics["odd_final_left_prob"]:.16e}` |
| 奇関数節の最終左電流 | `{dynamics["odd_final_left_current"]:.16e}` |
| 奇関数節の最大節振幅 | `{dynamics["odd_max_node_abs"]:.16e}` |
| 奇関数節の最大節電流 | `{dynamics["odd_max_node_current_abs"]:.16e}` |

![relative dynamics]({paths["dynamics_plot"]})

## 識別振動チャネル

識別振動 `eta` は反射生成チャネルではなく、保存される読出しチャネルとして扱った。これは、異なる `m_A,m_B` を交換消去チャネルへ直接混ぜると空間重なり節が壊れるためである。今回の実装では、反射はフェルミオン的逆相核の交換干渉で生成し、A/B の識別は `eta` 相関で読出す。

| 項目 | 読出し |
|---|---:|
| A の検出モード | `{labels["detected_A"]}` |
| B の検出モード | `{labels["detected_B"]}` |
| A ターゲット振幅 | `{labels["readout_A"][str(labels["target_A"])]:.16e}` |
| B ターゲット振幅 | `{labels["readout_B"][str(labels["target_B"])]:.16e}` |

## 可逆性と補償付き二乗閉鎖

局所交換干渉写像を二回適用し、波形が元へ戻るかを測定した。また、各サンプル係数 `x_n` に補償対 `i x_n` を付けた二乗閉鎖を、主要段階で評価した。

| 量 | 値 |
|---|---:|
| `U(pi)^2` 相対誤差 | `{reversibility["pi_twice_relative_error"]:.16e}` |
| `U(delta)U(-delta)` 最大相対誤差 | `{reversibility["inverse_max_relative_error"]:.16e}` |
| 補償付き二乗閉鎖 最大残差 | `{reversibility["closure_max_residual"]:.16e}` |

## AB/C 相互作用セル差し替え

前回の AB/C 完全弾性衝突シミュレーションにおけるセル内の直接 `q` 反転を使わず、局所交換干渉写像から得た `R,T` により、

```text
q_out = q_in * (T - R)
```

として進行方向読出し量を生成した。

| 項目 | 値 |
|---|---:|
| AB/C 差し替え成立 | `{str(ab_verdict["integrated_ab_c_replacement_valid"]).lower()}` |
| collision_cell_reached | `{str(ab_verdict["collision_cell_reached"]).lower()}` |
| post_collision_completed | `{str(ab_verdict["post_collision_completed"]).lower()}` |
| q_A before/after | `{ab["q_readouts"]["A_before"]:.16e}` / `{ab["q_readouts"]["A_after"]:.16e}` |
| q_B before/after | `{ab["q_readouts"]["B_before"]:.16e}` / `{ab["q_readouts"]["B_after"]:.16e}` |
| label A initial/final | `{ab["readout"]["m_A_initial"]}` / `{ab["readout"]["m_A_final"]}` |
| label B initial/final | `{ab["readout"]["m_B_initial"]}` / `{ab["readout"]["m_B_final"]}` |

![ab-c replacement]({paths["ab_c_replacement_plot"]})

## 出力ファイル

| 種類 | ファイル |
|---|---|
| JSON | [{Path(paths["json"]).name}]({paths["json"]}) |
| 位相掃引 CSV | [{Path(paths["phase_csv"]).name}]({paths["phase_csv"]}) |
| 相対時間発展 CSV | [{Path(paths["dynamics_csv"]).name}]({paths["dynamics_csv"]}) |
| 位相掃引図 | [{Path(paths["phase_plot"]).name}]({paths["phase_plot"]}) |
| 相対時間発展図 | [{Path(paths["dynamics_plot"]).name}]({paths["dynamics_plot"]}) |
| 片側入射散乱 CSV | [{Path(paths["scattering_csv"]).name}]({paths["scattering_csv"]}) |
| 片側入射散乱図 | [{Path(paths["scattering_plot"]).name}]({paths["scattering_plot"]}) |
| 可逆性 CSV | [{Path(paths["reversibility_csv"]).name}]({paths["reversibility_csv"]}) |
| AB/C 差し替え CSV | [{Path(paths["ab_c_replacement_csv"]).name}]({paths["ab_c_replacement_csv"]}) |
| AB/C 差し替え図 | [{Path(paths["ab_c_replacement_plot"]).name}]({paths["ab_c_replacement_plot"]}) |
"""


def run() -> Dict[str, object]:
    params = Params()

    deltas = np.linspace(0.0, math.pi, params.phase_sweep_n)
    phase_rows = [exchange_diagonal_metrics(params, float(delta), True) for delta in deltas]
    pi_with_exchange = exchange_diagonal_metrics(params, math.pi, True)
    pi_without_exchange = exchange_diagonal_metrics(params, math.pi, False)

    max_abs_error = max(row["abs_error"] for row in phase_rows)
    dynamics_rows = run_relative_dynamics(params)
    scattering_rows = run_single_sided_scattering(params)

    single_rows = [row for row in dynamics_rows if row["model"] == "single_free_packet"]
    odd_rows = [row for row in dynamics_rows if row["model"] == "odd_node_halfline_read"]
    single_final = single_rows[-1]
    odd_final = odd_rows[-1]
    odd_max_node_abs = max(abs(float(row["node_abs"])) for row in odd_rows)
    odd_max_node_current_abs = max(abs(float(row["node_current"])) for row in odd_rows)

    modes = [params.m_A, params.m_B]
    readout_a = eta_readout(params.m_A, modes)
    readout_b = eta_readout(params.m_B, modes)
    detected_a = detect_mode(readout_a)
    detected_b = detect_mode(readout_b)

    scatter_zero = min(scattering_rows, key=lambda row: abs(row["delta_f"] - 0.0))
    scatter_half_pi = min(scattering_rows, key=lambda row: abs(row["delta_f"] - math.pi / 2.0))
    scatter_pi = min(scattering_rows, key=lambda row: abs(row["delta_f"] - math.pi))
    reversibility_closure = run_reversibility_and_closure(params)
    ab_c_replacement = run_ab_c_replacement_test(
        params,
        scatter_pi["reflection_rate"],
        scatter_pi["transmission_rate"],
    )
    scattering_max_rate_abs_error = max(
        max(row["reflection_abs_error"], row["transmission_abs_error"]) for row in scattering_rows
    )
    scattering_max_norm_abs_error = max(row["norm_abs_error"] for row in scattering_rows)

    scatter_dx = params.scattering_domain_length / params.scattering_grid_n
    scatter_rho = (np.arange(params.scattering_grid_n) - params.scattering_grid_n // 2) * scatter_dx
    scatter_initial = gaussian_packet(
        scatter_rho,
        params.scattering_packet_center,
        params.scattering_packet_k0,
        params.scattering_packet_sigma,
    )
    scatter_initial = normalize(scatter_initial, scatter_dx)
    scatter_initial_metrics = side_metrics(scatter_rho, scatter_dx, scatter_initial)

    verdict = {
        "external_q_flip_used": False,
        "single_sided_initial_packet": bool(
            scatter_initial_metrics["left_prob"] > 1.0 - params.scattering_error_tol
            and scatter_initial_metrics["right_prob"] < params.scattering_error_tol
        ),
        "dynamic_transmission_at_zero": bool(
            scatter_zero["reflection_rate"] < params.scattering_error_tol
            and scatter_zero["transmission_rate"] > 1.0 - params.scattering_error_tol
        ),
        "dynamic_half_phase_split": bool(
            abs(scatter_half_pi["reflection_rate"] - 0.5) < params.scattering_error_tol
            and abs(scatter_half_pi["transmission_rate"] - 0.5) < params.scattering_error_tol
        ),
        "dynamic_reflection_at_pi": bool(
            scatter_pi["reflection_rate"] > 1.0 - params.scattering_error_tol
            and scatter_pi["transmission_rate"] < params.scattering_error_tol
        ),
        "dynamic_phase_sweep_matches_expected": bool(scattering_max_rate_abs_error < params.scattering_error_tol),
        "dynamic_norm_preserved": bool(scattering_max_norm_abs_error < params.scattering_error_tol),
        "static_node_at_pi": bool(pi_with_exchange["diagonal_relative_norm"] < params.node_tol),
        "phase_sweep_matches_expected": bool(max_abs_error < params.phase_error_tol),
        "exchange_path_required": bool(
            pi_with_exchange["diagonal_relative_norm"] < params.node_tol
            and pi_without_exchange["diagonal_relative_norm"] > 0.1
        ),
        "single_packet_transmits": bool(float(single_final["left_prob"]) < params.transmission_left_prob_tol),
        "odd_node_reflects": bool(
            abs(float(odd_final["left_prob"]) - 0.5) < params.reflection_left_prob_tol
            and float(odd_final["left_current"]) < 0.0
        ),
        "odd_node_current_zero": bool(odd_max_node_current_abs < params.dynamic_node_tol),
        "label_mode_preserved_A": bool(detected_a == params.m_A),
        "label_mode_preserved_B": bool(detected_b == params.m_B),
        "pi_twice_reversible": bool(
            float(reversibility_closure["pi_twice_relative_error"]) < params.reversibility_tol
        ),
        "inverse_sweep_reversible": bool(
            float(reversibility_closure["inverse_max_relative_error"]) < params.reversibility_tol
        ),
        "compensated_square_closure_preserved": bool(
            float(reversibility_closure["closure_max_residual"]) < params.closure_tol
        ),
        "ab_c_replacement_valid": bool(
            ab_c_replacement["verdict"]["integrated_ab_c_replacement_valid"]
        ),
    }
    core_verdict_keys = [
        "single_sided_initial_packet",
        "dynamic_transmission_at_zero",
        "dynamic_half_phase_split",
        "dynamic_reflection_at_pi",
        "dynamic_phase_sweep_matches_expected",
        "dynamic_norm_preserved",
        "static_node_at_pi",
        "phase_sweep_matches_expected",
        "exchange_path_required",
        "label_mode_preserved_A",
        "label_mode_preserved_B",
        "pi_twice_reversible",
        "inverse_sweep_reversible",
        "compensated_square_closure_preserved",
        "ab_c_replacement_valid",
    ]
    verdict["mechanism_valid_minimal"] = bool(
        all(verdict[key] for key in core_verdict_keys) and verdict["external_q_flip_used"] is False
    )

    outputs = {
        "json": "fermionic_interference_reflection_result_v2.json",
        "phase_csv": "fermionic_interference_phase_sweep_v2.csv",
        "phase_plot": "fermionic_interference_phase_sweep_v2.png",
        "dynamics_csv": "fermionic_interference_relative_dynamics_v2.csv",
        "dynamics_plot": "fermionic_interference_relative_dynamics_v2.png",
        "scattering_csv": "fermionic_interference_single_sided_scattering_v2.csv",
        "scattering_plot": "fermionic_interference_single_sided_scattering_v2.png",
        "reversibility_csv": "fermionic_interference_reversibility_sweep_v2.csv",
        "ab_c_replacement_csv": "fermionic_interference_ab_c_replacement_timeline_v2.csv",
        "ab_c_replacement_plot": "fermionic_interference_ab_c_replacement_v2.png",
        "report": "fermionic_interference_reflection_report_v2.md",
    }

    result: Dict[str, object] = {
        "experiment": "fermionic_interference_reflection_v2",
        "params": asdict(params),
        "phase_sweep_summary": {
            "pi_diagonal_relative_norm": pi_with_exchange["diagonal_relative_norm"],
            "max_abs_error": max_abs_error,
        },
        "controls": {
            "with_exchange_pi_diagonal_relative_norm": pi_with_exchange["diagonal_relative_norm"],
            "without_exchange_pi_diagonal_relative_norm": pi_without_exchange["diagonal_relative_norm"],
        },
        "relative_dynamics_summary": {
            "single_final_left_prob": float(single_final["left_prob"]),
            "odd_final_left_prob": float(odd_final["left_prob"]),
            "odd_final_left_current": float(odd_final["left_current"]),
            "odd_max_node_abs": odd_max_node_abs,
            "odd_max_node_current_abs": odd_max_node_current_abs,
        },
        "single_sided_scattering_summary": {
            "initial_left_prob": scatter_initial_metrics["left_prob"],
            "initial_right_prob": scatter_initial_metrics["right_prob"],
            "reflection_at_zero": scatter_zero["reflection_rate"],
            "transmission_at_zero": scatter_zero["transmission_rate"],
            "reflection_at_half_pi": scatter_half_pi["reflection_rate"],
            "transmission_at_half_pi": scatter_half_pi["transmission_rate"],
            "reflection_at_pi": scatter_pi["reflection_rate"],
            "transmission_at_pi": scatter_pi["transmission_rate"],
            "max_rate_abs_error": scattering_max_rate_abs_error,
            "max_norm_abs_error": scattering_max_norm_abs_error,
        },
        "reversibility_closure_summary": {
            "pi_twice_relative_error": reversibility_closure["pi_twice_relative_error"],
            "inverse_max_relative_error": reversibility_closure["inverse_max_relative_error"],
            "closure_max_residual": reversibility_closure["closure_max_residual"],
            "closure_stages": reversibility_closure["closure_stages"],
        },
        "ab_c_replacement_summary": {
            "method": ab_c_replacement["method"],
            "reflection_rate_used": ab_c_replacement["reflection_rate_used"],
            "transmission_rate_used": ab_c_replacement["transmission_rate_used"],
            "cell_widths": ab_c_replacement["cell_widths"],
            "steps": ab_c_replacement["steps"],
            "q_readouts": ab_c_replacement["q_readouts"],
            "readout": ab_c_replacement["readout"],
            "verdict": ab_c_replacement["verdict"],
        },
        "identification_readout": {
            "policy": "eta modes are preserved spectator/readout labels; exchange-interference reflection is generated in the fermionic core channel.",
            "target_A": params.m_A,
            "target_B": params.m_B,
            "detected_A": detected_a,
            "detected_B": detected_b,
            "readout_A": {str(k): v for k, v in readout_a.items()},
            "readout_B": {str(k): v for k, v in readout_b.items()},
        },
        "verdict": verdict,
        "output_dir": str(OUT_DIR.relative_to(BASE_DIR)),
        "outputs": outputs,
    }

    write_csv(OUT_DIR / outputs["phase_csv"], phase_rows)
    write_csv(OUT_DIR / outputs["dynamics_csv"], dynamics_rows)
    write_csv(OUT_DIR / outputs["scattering_csv"], scattering_rows)
    write_csv(OUT_DIR / outputs["reversibility_csv"], reversibility_closure["inverse_rows"])
    write_csv(OUT_DIR / outputs["ab_c_replacement_csv"], ab_c_replacement["timeline"])
    write_phase_plot(phase_rows, OUT_DIR / outputs["phase_plot"])
    write_dynamics_plot(dynamics_rows, OUT_DIR / outputs["dynamics_plot"])
    write_scattering_plot(scattering_rows, OUT_DIR / outputs["scattering_plot"])
    write_ab_c_replacement_plot(ab_c_replacement["timeline"], OUT_DIR / outputs["ab_c_replacement_plot"])

    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / outputs["report"]).write_text(build_report(result), encoding="utf-8")

    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))
