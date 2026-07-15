#!/usr/bin/env python3
"""Minimal System A alpha-neighborhood R sweep.

This script intentionally keeps only the operations needed to check whether the
System A two-channel scattering score has a sharp bottom near alpha=137.
It does not generate plots, reports, or multi-case sweeps.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "minimal_alpha137_bugcheck_result_v1"


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


def inner(v: np.ndarray, w: np.ndarray) -> complex:
    return complex(np.vdot(v, w))


def make_grids(params: Params) -> Tuple[np.ndarray, np.ndarray]:
    chi = np.linspace(-math.pi, math.pi, params.chi_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, params.eta_grid_n, endpoint=False)
    return chi, eta


def packet_kernel(
    params: Params,
    harmonics: Tuple[int, ...],
    weights: Tuple[float, ...],
    phases: Tuple[float, ...],
    wavelength_scales: Tuple[float, ...],
) -> np.ndarray:
    chi, _eta = make_grids(params)
    u = chi - params.chi_center
    kernel = np.zeros_like(chi, dtype=float)
    for harmonic, weight, phase, wavelength_scale in zip(harmonics, weights, phases, wavelength_scales):
        effective_frequency = float(harmonic) / float(wavelength_scale)
        kernel += float(weight) * np.cos(effective_frequency * u + float(phase))
    weight_norm = math.sqrt(sum(float(weight) ** 2 for weight in weights))
    return kernel / weight_norm if weight_norm > 0.0 else kernel


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
    kernel = packet_kernel(params, harmonics, weights, phases, wavelength_scales)
    phase_chi = np.exp(1j * q * params.p0 * (chi - params.chi_center))
    eta_phase = np.exp(1j * m * eta)
    psi = (kernel * phase_chi)[:, None] * eta_phase[None, :]
    return normalize(psi.reshape(-1))


def scattering_coefficients(reflection_rate: float) -> Tuple[complex, complex, float, float]:
    if reflection_rate < 0.0 or reflection_rate > 1.0:
        raise ValueError(f"reflection rate must be in [0, 1]: {reflection_rate}")
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
        if denom <= 0.0:
            return {0: 1.0}
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
        if total <= 0.0:
            return {0: 1.0}
        return {k: v / total for k, v in raw.items()}

    def localization(self, vector: np.ndarray) -> float:
        denom = norm2(vector)
        if denom <= 0.0:
            return float("nan")
        prob = np.abs(vector) ** 2 / denom
        return float(np.sum(prob**2))


def effective_n(distribution: Dict[int, float]) -> float:
    return float(sum(float(n) * float(weight) for n, weight in distribution.items()))


def distribution_similarity(a: Dict[int, float], b: Dict[int, float]) -> float:
    keys = set(a) | set(b)
    va = np.asarray([float(a.get(k, 0.0)) for k in keys], dtype=float)
    vb = np.asarray([float(b.get(k, 0.0)) for k in keys], dtype=float)
    na = float(np.linalg.norm(va))
    nb = float(np.linalg.norm(vb))
    if na <= 0.0 or nb <= 0.0:
        return float("nan")
    return float(max(0.0, min(1.0, np.dot(va, vb) / (na * nb))))


def projection_weight(vector: np.ndarray, basis: np.ndarray) -> float:
    denom = norm2(vector)
    base_denom = norm2(basis)
    if denom <= 0.0 or base_denom <= 0.0:
        return float("nan")
    v = vector / math.sqrt(denom)
    b = basis / math.sqrt(base_denom)
    return float(abs(inner(b, v)) ** 2)


def default_states(params: Params) -> Tuple[np.ndarray, np.ndarray]:
    # Bug-check case:
    # A: base harmonic only.
    # B: base + second harmonic, second component has 0.1% amplitude,
    #    30% phase offset, and 30% wavelength scale offset.
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


def summarize_r(params: Params, metrics: Metrics, reflection_rate: float) -> Dict[str, float]:
    t, r, T, R_actual = scattering_coefficients(reflection_rate)
    a, b = default_states(params)
    initial_a = a.copy()
    initial_b = b.copy()
    h_a0 = metrics.harmonic_distribution(initial_a)
    h_b0 = metrics.harmonic_distribution(initial_b)

    l_records: List[Tuple[int, float]] = []
    n_records: List[Tuple[int, float]] = []
    transfer_records: List[Tuple[int, float]] = []
    for collision in range(params.max_collision + 1):
        h_a = metrics.harmonic_distribution(a)
        h_b = metrics.harmonic_distribution(b)
        l_gap = abs(metrics.localization(a) - metrics.localization(b))
        n_gap = abs(effective_n(h_a) - effective_n(h_b))
        b_to_a_transfer = distribution_similarity(h_a, h_b0)
        l_records.append((collision, l_gap))
        n_records.append((collision, n_gap))
        transfer_records.append((collision, b_to_a_transfer))
        if collision >= params.max_collision:
            break
        a_next = normalize(r * a + t * b)
        b_next = normalize(t * a + r * b)
        a, b = a_next, b_next

    l_collision, l_min = min(l_records, key=lambda item: item[1])
    n_collision, n_min = min(n_records, key=lambda item: item[1])
    transfer_collision, transfer_max = max(transfer_records, key=lambda item: item[1])
    return {
        "R_input": float(reflection_rate),
        "R_actual": float(R_actual),
        "T_actual": float(T),
        "g": float(1.0 - reflection_rate),
        "L_gap_min": float(l_min),
        "L_gap_min_collision": int(l_collision),
        "N_eff_gap_min": float(n_min),
        "N_eff_gap_min_collision": int(n_collision),
        "max_B_to_A_transfer": float(transfer_max),
        "B_to_A_transfer_collision": int(transfer_collision),
        "origin_A_at_initial_A": projection_weight(initial_a, initial_a),
        "origin_B_at_initial_B": projection_weight(initial_b, initial_b),
    }


def r_values(start: float, stop: float, step: float) -> List[float]:
    if step <= 0.0:
        raise ValueError("R step must be positive")
    if stop < start:
        raise ValueError("R stop must be greater than or equal to R start")
    count = int(math.floor((stop - start) / step + 1.0e-12))
    values = [round(start + i * step, 15) for i in range(count + 1)]
    if not values or abs(values[-1] - stop) > 0.5 * step:
        values.append(round(stop, 15))
    return sorted(set(values))


def add_joint_scores(rows: List[Dict[str, float]]) -> None:
    l_norm = max(max([float(row["L_gap_min"]) for row in rows] or [1.0]), 1.0e-300)
    n_norm = max(max([float(row["N_eff_gap_min"]) for row in rows] or [1.0]), 1.0e-300)
    for row in rows:
        row["joint_R_score"] = (
            float(row["L_gap_min"]) / l_norm
            + float(row["N_eff_gap_min"]) / n_norm
            + (1.0 - float(row["max_B_to_A_transfer"]))
        )


def alpha_to_r(alpha_inv: float) -> float:
    return 1.0 - math.sqrt(4.0 * math.pi / alpha_inv)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha-inv", type=float, default=137.035999177)
    parser.add_argument("--half-width", type=float, default=3.0e-7)
    parser.add_argument("--r-step", type=float, default=1.0e-8)
    parser.add_argument("--r-min", type=float)
    parser.add_argument("--r-max", type=float)
    parser.add_argument("--max-collision", type=int, default=256)
    parser.add_argument("--chi-grid-n", type=int, default=512)
    parser.add_argument("--eta-grid-n", type=int, default=16)
    parser.add_argument("--output-dir", default=str(OUT_DIR))
    args = parser.parse_args()

    r_alpha = alpha_to_r(args.alpha_inv)
    start = args.r_min if args.r_min is not None else r_alpha - args.half_width
    stop = args.r_max if args.r_max is not None else r_alpha + args.half_width
    params = Params(
        chi_grid_n=args.chi_grid_n,
        eta_grid_n=args.eta_grid_n,
        max_collision=args.max_collision,
    )
    metrics = Metrics(params)
    rows = [summarize_r(params, metrics, r_value) for r_value in r_values(start, stop, args.r_step)]
    add_joint_scores(rows)
    best = min(rows, key=lambda row: float(row["joint_R_score"]))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "minimal_alpha137_sweep_summary_v1.csv"
    json_path = out_dir / "minimal_alpha137_sweep_best_v1.json"
    fieldnames = list(rows[0].keys())
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "case": "A=(1), B=(1 + 0.001*cos(2 chi / 1.3 + 0.3*2pi))",
        "alpha_inv_input": args.alpha_inv,
        "R_from_alpha_inv": r_alpha,
        "R_min": start,
        "R_max": stop,
        "R_step": args.r_step,
        "row_count": len(rows),
        "best": best,
        "csv": str(csv_path),
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
