#!/usr/bin/env python3
"""Section 18.2-18.4: Harmonic Dependence and Internal State Comparison

Test whether:
  18.2 反復回数依存性 (iteration count dependency)
       - Structure present initially or emerges with iterations?
       - Is iteration an amplifier or creator of peaks?

  18.3 N=1 と N=2 の比較 (single vs double harmonic)
       - Does N=1 (no additional mode) prevent peak appearance?
       - Does N=2 cause immediate R_137 manifestation?

  18.4 内部状態の同一性検査 (internal state identity check)
       - Are the two peaks A_128 and A_137 same or different?
       - Same state, different readout windows?
       - A↔B symmetric?
"""

from __future__ import annotations

import csv
import json
import math
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np


DEFAULT_STEPS_MAX = 1024
METASTABLE_AMP_TARGET = 0.02

R_137_CANDIDATE = 0.697177902556148
R_128_CANDIDATE = 0.688363902556148


def normalize_pair(a: complex, b: complex) -> Tuple[complex, complex]:
    q = abs(a) ** 2 + abs(b) ** 2
    if q <= 0.0:
        raise ValueError("zero AB norm")
    scale = 1.0 / math.sqrt(q)
    return a * scale, b * scale


def state_from_s_phi(s: float, phi: float) -> Tuple[complex, complex]:
    if abs(s) >= 0.5:
        raise ValueError(f"|s| must be < 0.5: {s}")
    a = math.sqrt(0.5 + s)
    b = math.sqrt(0.5 - s) * complex(math.cos(phi), math.sin(phi))
    return normalize_pair(a, b)


def scattering_coefficients(reflection_rate: float) -> Tuple[complex, complex, float, float]:
    if reflection_rate < 0.0 or reflection_rate > 1.0:
        raise ValueError(f"reflection rate must be in [0, 1]: {reflection_rate}")
    delta_f = 2.0 * math.asin(math.sqrt(reflection_rate))
    half_delta = 0.5 * delta_f
    phase = complex(math.cos(half_delta), math.sin(half_delta))
    t = phase * math.cos(half_delta)
    r = -1j * phase * math.sin(half_delta)
    return complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


def prefix_error(prefix_sums: List[float], s_min: float, s_max: float) -> Tuple[float, float, float, float]:
    n = len(prefix_sums) - 1
    if n <= 0:
        raise ValueError("empty prefix")
    total = prefix_sums[n]
    s_mean = total / n
    s_amp = 0.5 * (s_max - s_min)
    mid = n // 2
    if mid == 0:
        s_drift = 0.0
    else:
        first_mean = prefix_sums[mid] / mid
        second_mean = (total - prefix_sums[mid]) / (n - mid)
        s_drift = abs(second_mean - first_mean)
    gray_error = abs(s_mean) + abs(s_amp - METASTABLE_AMP_TARGET) + s_drift
    return gray_error, s_mean, s_amp, s_drift


def run_time_series_harmonic_N(
    reflection_rate: float,
    phi: float,
    harmonic_orders: List[int],  # [1], [2], [1,2,3], etc.
    max_steps: int,
    s0: float = 1.0e-2,
) -> Dict[str, Any]:
    """
    Run AB scattering with harmonic mixing.

    harmonic_orders specifies which harmonics to include.
    For simplicity, we embed harmonics as phase shifts in the initial B.

    Returns time series with internal states at each step.
    """
    # Initialize with harmonic modulation
    a = math.sqrt(0.5 + s0)
    if len(harmonic_orders) == 1 and harmonic_orders[0] == 1:
        # Single harmonic baseline
        b = math.sqrt(0.5 - s0) * complex(math.cos(phi), math.sin(phi))
    else:
        # Mix harmonics by adding phase contributions
        b_real = math.sqrt(0.5 - s0) * math.cos(phi)
        b_imag = math.sqrt(0.5 - s0) * math.sin(phi)
        for n in harmonic_orders[1:]:
            phase_shift = (n - 1) * 0.1  # Arbitrary phase shift per harmonic
            b_real += 0.1 * math.cos(phi + phase_shift)
            b_imag += 0.1 * math.sin(phi + phase_shift)
        b = complex(b_real, b_imag)

    a, b = normalize_pair(a, b)
    t, r, _t_power, _r_power = scattering_coefficients(reflection_rate)

    prefix_sums = [0.0]
    s_min = float("inf")
    s_max = float("-inf")
    s_values = []
    states_a = [a]
    states_b = [b]

    for step in range(max_steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s_value = (p_a - p_b) / q
        s_values.append(s_value)
        s_min = min(s_min, s_value)
        s_max = max(s_max, s_value)
        prefix_sums.append(prefix_sums[-1] + s_value)

        if step < max_steps:
            a, b = normalize_pair(r * a + t * b, t * a + r * b)
            states_a.append(a)
            states_b.append(b)

    gray_error, s_mean, s_amp, s_drift = prefix_error(prefix_sums, s_min, s_max)

    return {
        "R": reflection_rate,
        "phi": phi,
        "s0": s0,
        "harmonic_orders": harmonic_orders,
        "max_steps": max_steps,
        "s_values": s_values,
        "states_a": states_a,
        "states_b": states_b,
        "prefix_sums": prefix_sums,
        "s_min": s_min,
        "s_max": s_max,
        "s_mean": s_mean,
        "s_amp": s_amp,
        "s_drift": s_drift,
        "gray_error": gray_error,
        "gray_depth": -math.log10(max(gray_error, 1.0e-300)),
    }


def iteration_count_sweep(
    R_values: List[float],
    phi: float = 0.0,
    step_counts: List[int] = None,
) -> List[Dict[str, Any]]:
    """
    18.2: Test whether peaks appear at initial iterations or emerge later.
    """
    if step_counts is None:
        step_counts = [1, 2, 3, 5, 10, 20, 50, 100, 250, 512, 1024]

    results = []
    for R in R_values:
        for max_step in step_counts:
            data = run_time_series_harmonic_N(R, phi, [1, 2], max_step)
            results.append({
                "R": R,
                "phi": phi,
                "max_steps": max_step,
                "gray_error": data["gray_error"],
                "gray_depth": data["gray_depth"],
                "s_mean": data["s_mean"],
                "s_amp": data["s_amp"],
            })

    return results


def harmonic_order_comparison(
    R_values: List[float],
    phi: float = 0.0,
    max_steps: int = 512,
) -> List[Dict[str, Any]]:
    """
    18.3: Compare N=1 (no additional mode) vs N=2, N=3, etc.
    """
    results = []
    harmonic_configs = [
        [1],           # N=1: single mode only
        [1, 2],        # N=2: base + first harmonic
        [1, 2, 3],     # N=3
        [1, 2, 4],     # N=4 (skip odd)
    ]

    for R in R_values:
        for harmonics in harmonic_configs:
            data = run_time_series_harmonic_N(R, phi, harmonics, max_steps)
            results.append({
                "R": R,
                "phi": phi,
                "harmonics": harmonics,
                "n_harmonics": len(harmonics),
                "gray_error": data["gray_error"],
                "gray_depth": data["gray_depth"],
                "s_mean": data["s_mean"],
                "s_amp": data["s_amp"],
            })

    return results


def internal_state_comparison(
    R_128: float = R_128_CANDIDATE,
    R_137: float = R_137_CANDIDATE,
    phi: float = 0.0,
    max_steps: int = 512,
) -> Dict[str, Any]:
    """
    18.4: Compare internal states at the two R values.

    Are the final attractor states A_128 and A_137 identical,
    symmetric (A↔B swap), or fundamentally different?
    """
    data_128 = run_time_series_harmonic_N(R_128, phi, [1, 2], max_steps)
    data_137 = run_time_series_harmonic_N(R_137, phi, [1, 2], max_steps)

    # Extract final state (or running average near end)
    tail_len = min(100, max_steps // 4)
    a_128_tail = data_128["states_a"][-tail_len:]
    b_128_tail = data_128["states_b"][-tail_len:]
    a_137_tail = data_137["states_a"][-tail_len:]
    b_137_tail = data_137["states_b"][-tail_len:]

    a_128_mean = np.mean([complex(c) for c in a_128_tail])
    b_128_mean = np.mean([complex(c) for c in b_128_tail])
    a_137_mean = np.mean([complex(c) for c in a_137_tail])
    b_137_mean = np.mean([complex(c) for c in b_137_tail])

    # Compute distances
    dist_a_to_a = abs(a_128_mean - a_137_mean)
    dist_b_to_b = abs(b_128_mean - b_137_mean)
    dist_a_to_b = abs(a_128_mean - b_137_mean)  # Check if swapped
    dist_b_to_a = abs(b_128_mean - a_137_mean)

    # Check closure residual
    def closure_residual(a_seq: List[complex], b_seq: List[complex]) -> List[float]:
        """Compute E_close = |A_n^2 + B_n^2| / (|A_n|^2 + |B_n|^2)"""
        residuals = []
        for a, b in zip(a_seq, b_seq):
            a, b = complex(a), complex(b)
            closure = abs(a**2 + b**2)
            norm = abs(a)**2 + abs(b)**2
            if norm > 0:
                residuals.append(closure / norm)
            else:
                residuals.append(0.0)
        return residuals

    closure_128 = closure_residual(a_128_tail, b_128_tail)
    closure_137 = closure_residual(a_137_tail, b_137_tail)

    return {
        "R_128": R_128,
        "R_137": R_137,
        "a_128_mean": complex(a_128_mean),
        "b_128_mean": complex(b_128_mean),
        "a_137_mean": complex(a_137_mean),
        "b_137_mean": complex(b_137_mean),
        "distance_a_to_a": dist_a_to_a,
        "distance_b_to_b": dist_b_to_b,
        "distance_a_to_b_swapped": dist_a_to_b,
        "distance_b_to_a_swapped": dist_b_to_a,
        "closure_mean_128": float(np.mean(closure_128)),
        "closure_mean_137": float(np.mean(closure_137)),
        "closure_std_128": float(np.std(closure_128)),
        "closure_std_137": float(np.std(closure_137)),
        "interpretation": _interpret_state_comparison(
            dist_a_to_a, dist_b_to_b, dist_a_to_b, dist_b_to_a
        ),
    }


def _interpret_state_comparison(
    dist_aa: float, dist_bb: float, dist_ab: float, dist_ba: float
) -> str:
    """Interpret the relationship between A_128 and A_137."""
    min_dist = min(dist_aa, dist_bb, dist_ab, dist_ba)

    if min_dist == dist_aa and dist_aa < 0.01:
        return "Identical_internal_states"
    elif min_dist == dist_ab and dist_ab < 0.01:
        return "Symmetric_A_B_swap"
    elif min_dist == dist_ba and dist_ba < 0.01:
        return "Symmetric_B_A_swap"
    else:
        return "Distinct_attractor_states"


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Section 18.2-18.4: Harmonic dependence and state comparison"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for results"
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Section 18.2: Iteration Count Dependency ===")
    R_sweep = np.linspace(0.686, 0.700, 29)
    iter_results = iteration_count_sweep(R_sweep.tolist())

    iter_csv = args.output_dir / "section18_iteration_dependency_v1.csv"
    with open(iter_csv, "w", newline="") as f:
        fieldnames = ["R", "phi", "max_steps", "gray_error", "gray_depth", "s_mean", "s_amp"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(iter_results)
    print(f"Saved to {iter_csv}")

    print("\n=== Section 18.3: Harmonic Order Comparison ===")
    harm_results = harmonic_order_comparison(R_sweep.tolist())

    harm_csv = args.output_dir / "section18_harmonic_comparison_v1.csv"
    with open(harm_csv, "w", newline="") as f:
        fieldnames = ["R", "phi", "harmonics", "n_harmonics", "gray_error", "gray_depth", "s_mean", "s_amp"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in harm_results:
            row["harmonics"] = str(row["harmonics"])
            writer.writerow(row)
    print(f"Saved to {harm_csv}")

    print("\n=== Section 18.4: Internal State Comparison ===")
    state_comp = internal_state_comparison()

    state_json = args.output_dir / "section18_state_comparison_v1.json"
    with open(state_json, "w") as f:
        json.dump({
            "R_128": state_comp["R_128"],
            "R_137": state_comp["R_137"],
            "a_128_mean": str(state_comp["a_128_mean"]),
            "b_128_mean": str(state_comp["b_128_mean"]),
            "a_137_mean": str(state_comp["a_137_mean"]),
            "b_137_mean": str(state_comp["b_137_mean"]),
            "distance_a_to_a": state_comp["distance_a_to_a"],
            "distance_b_to_b": state_comp["distance_b_to_b"],
            "distance_a_to_b_swapped": state_comp["distance_a_to_b_swapped"],
            "distance_b_to_a_swapped": state_comp["distance_b_to_a_swapped"],
            "closure_mean_128": state_comp["closure_mean_128"],
            "closure_mean_137": state_comp["closure_mean_137"],
            "closure_std_128": state_comp["closure_std_128"],
            "closure_std_137": state_comp["closure_std_137"],
            "interpretation": state_comp["interpretation"],
        }, f, indent=2)
    print(f"Saved to {state_json}")
    print(f"\nInterpretation: {state_comp['interpretation']}")


if __name__ == "__main__":
    main()
