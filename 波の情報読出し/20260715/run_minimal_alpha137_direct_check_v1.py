#!/usr/bin/env python3
"""Direct one-point check for the minimal System A alpha test case."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


@dataclass
class Params:
    chi_grid_n: int = 512
    eta_grid_n: int = 16
    high_n: int = 63
    chi_center: float = 0.0
    p0: float = 1.0
    q_A: float = 1.0
    q_B: float = -1.0
    m_A: int = 1
    m_B: int = 2
    max_collision: int = 256


def normalize(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    if norm <= 0.0:
        raise ValueError("zero norm state")
    return v / norm


def norm2(v: np.ndarray) -> float:
    return float(np.vdot(v, v).real)


def make_grids(params: Params) -> Tuple[np.ndarray, np.ndarray]:
    chi = np.linspace(-math.pi, math.pi, params.chi_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, params.eta_grid_n, endpoint=False)
    return chi, eta


def make_state(
    params: Params,
    harmonics: Tuple[int, ...],
    weights: Tuple[float, ...],
    phases: Tuple[float, ...],
    wavelength_scales: Tuple[float, ...],
    q: float,
    m: int,
) -> np.ndarray:
    chi, eta = make_grids(params)
    u = chi - params.chi_center
    kernel = np.zeros_like(chi, dtype=float)
    for harmonic, weight, phase, wavelength_scale in zip(harmonics, weights, phases, wavelength_scales):
        kernel += float(weight) * np.cos((float(harmonic) / float(wavelength_scale)) * u + float(phase))
    weight_norm = math.sqrt(sum(float(weight) ** 2 for weight in weights))
    if weight_norm > 0.0:
        kernel = kernel / weight_norm
    phase_chi = np.exp(1j * q * params.p0 * u)
    eta_phase = np.exp(1j * m * eta)
    return normalize(((kernel * phase_chi)[:, None] * eta_phase[None, :]).reshape(-1))


def scattering_coefficients(reflection_rate: float) -> Tuple[complex, complex, float, float]:
    delta_f = 2.0 * math.asin(math.sqrt(reflection_rate))
    t = np.exp(0.5j * delta_f) * math.cos(0.5 * delta_f)
    r = -1j * np.exp(0.5j * delta_f) * math.sin(0.5 * delta_f)
    return complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


class Metrics:
    def __init__(self, params: Params):
        self.params = params
        self.freqs = np.fft.fftfreq(params.chi_grid_n, d=1.0 / params.chi_grid_n)
        rounded = np.rint(self.freqs).astype(int)
        self.indices_by_abs_n: Dict[int, List[int]] = {}
        for n_abs in range(min(params.chi_grid_n // 2, params.high_n + 2) + 1):
            self.indices_by_abs_n[n_abs] = [
                int(i) for i, freq in enumerate(rounded) if abs(int(freq)) == n_abs
            ]

    def harmonic_distribution(self, vector: np.ndarray) -> Dict[int, float]:
        denom = norm2(vector)
        arr = vector.reshape(self.params.chi_grid_n, self.params.eta_grid_n)
        transformed = np.fft.fft(arr, axis=0, norm="ortho")
        power = np.sum(np.abs(transformed) ** 2, axis=1)
        raw: Dict[int, float] = {}
        total = 0.0
        for n_abs, indices in self.indices_by_abs_n.items():
            amount = float(np.sum(power[indices]).real) / denom if indices else 0.0
            if amount > 1.0e-14:
                raw[n_abs] = max(amount, 0.0)
                total += raw[n_abs]
        return {k: v / total for k, v in raw.items()} if total > 0.0 else {0: 1.0}

    def localization(self, vector: np.ndarray) -> float:
        prob = np.abs(vector) ** 2 / norm2(vector)
        return float(np.sum(prob**2))


def effective_n(distribution: Dict[int, float]) -> float:
    return float(sum(float(n) * float(weight) for n, weight in distribution.items()))


def distribution_similarity(a: Dict[int, float], b: Dict[int, float]) -> float:
    keys = set(a) | set(b)
    va = np.asarray([float(a.get(k, 0.0)) for k in keys], dtype=float)
    vb = np.asarray([float(b.get(k, 0.0)) for k in keys], dtype=float)
    return float(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))


def default_states(params: Params) -> Tuple[np.ndarray, np.ndarray]:
    a = make_state(params, (1,), (1.0,), (0.0,), (1.0,), params.q_A, params.m_A)
    b = make_state(
        params,
        (1, 2),
        (1.0, 0.001),
        (0.0, 0.30 * 2.0 * math.pi),
        (1.0, 1.30),
        params.q_B,
        params.m_B,
    )
    return a, b


def alpha_inv_to_r(alpha_inv: float) -> float:
    return 1.0 - math.sqrt(4.0 * math.pi / alpha_inv)


def check_one(params: Params, reflection_rate: float) -> Dict[str, float]:
    metrics = Metrics(params)
    t, r, T, R_actual = scattering_coefficients(reflection_rate)
    a, b = default_states(params)
    h_b0 = metrics.harmonic_distribution(b)
    l_gap_best = (0, float("inf"))
    n_gap_best = (0, float("inf"))
    transfer_best = (0, -float("inf"))
    for collision in range(params.max_collision + 1):
        h_a = metrics.harmonic_distribution(a)
        h_b = metrics.harmonic_distribution(b)
        l_gap = abs(metrics.localization(a) - metrics.localization(b))
        n_gap = abs(effective_n(h_a) - effective_n(h_b))
        transfer = distribution_similarity(h_a, h_b0)
        if l_gap < l_gap_best[1]:
            l_gap_best = (collision, l_gap)
        if n_gap < n_gap_best[1]:
            n_gap_best = (collision, n_gap)
        if transfer > transfer_best[1]:
            transfer_best = (collision, transfer)
        if collision >= params.max_collision:
            break
        a_next = normalize(r * a + t * b)
        b_next = normalize(t * a + r * b)
        a, b = a_next, b_next
    return {
        "R_input": float(reflection_rate),
        "R_actual": float(R_actual),
        "T_actual": float(T),
        "g": float(1.0 - reflection_rate),
        "sqrt_g": math.sqrt(max(1.0 - reflection_rate, 0.0)),
        "sqrt_R": math.sqrt(max(reflection_rate, 0.0)),
        "theta_deg": math.degrees(math.asin(math.sqrt(max(1.0 - reflection_rate, 0.0)))),
        "alpha_inv_from_R": 4.0 * math.pi / ((1.0 - reflection_rate) ** 2),
        "L_gap_min": float(l_gap_best[1]),
        "L_gap_min_collision": int(l_gap_best[0]),
        "N_eff_gap_min": float(n_gap_best[1]),
        "N_eff_gap_min_collision": int(n_gap_best[0]),
        "max_B_to_A_transfer": float(transfer_best[1]),
        "B_to_A_transfer_collision": int(transfer_best[0]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha-inv", type=float, default=137.035999177)
    parser.add_argument("--R", type=float, help="override R directly")
    parser.add_argument("--max-collision", type=int, default=256)
    parser.add_argument("--chi-grid-n", type=int, default=512)
    parser.add_argument("--eta-grid-n", type=int, default=16)
    args = parser.parse_args()

    params = Params(
        chi_grid_n=args.chi_grid_n,
        eta_grid_n=args.eta_grid_n,
        max_collision=args.max_collision,
    )
    r_value = float(args.R) if args.R is not None else alpha_inv_to_r(args.alpha_inv)
    payload = {
        "case": "A=(1), B=(1 + 0.001*cos(2 chi / 1.3 + 0.3*2pi))",
        "alpha_inv_input": args.alpha_inv,
        "R_from_alpha_inv": alpha_inv_to_r(args.alpha_inv),
        "checked": check_one(params, r_value),
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
