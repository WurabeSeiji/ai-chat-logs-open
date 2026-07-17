#!/usr/bin/env python3
"""Section 18 Verification: Scattering Matrix Invariants

This script tests the hypothesis that R_128 and R_137 are common fixed points
by examining whether trace, determinant, or eigenvalue magnitudes coincide
at these R values across multiple harmonic conditions.

Based on Section 18.1 of 20260717 CHATGPT思考実験.md:

  18.1 一周期散乱行列の不変量比較

  複数条件 `c` に対して、二点近傍で以下を比較する。

  Tr M_c(R)
  det M_c(R)
  |λ_±,c(R)|
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
from scipy import linalg


DEFAULT_STEPS = 512
METASTABLE_AMP_TARGET = 0.02

# Target R values from 20260715 experiment
R_137_CANDIDATE = 0.697177902556148
R_128_CANDIDATE = 0.688363902556148

# Sweep ranges around these points
R_SWEEP_MIN = 0.686
R_SWEEP_MAX = 0.700
R_POINTS = 141  # ~1e-4 precision


def normalize_pair(a: complex, b: complex) -> Tuple[complex, complex]:
    q = abs(a) ** 2 + abs(b) ** 2
    if q <= 0.0:
        raise ValueError("zero AB norm")
    scale = 1.0 / math.sqrt(q)
    return a * scale, b * scale


def state_from_s_phi(s: float, phi: float) -> Tuple[complex, complex]:
    """Generate initial AB state from (s, phi) parameters."""
    if abs(s) >= 0.5:
        raise ValueError(f"|s| must be < 0.5: {s}")
    a = math.sqrt(0.5 + s)
    b = math.sqrt(0.5 - s) * complex(math.cos(phi), math.sin(phi))
    return normalize_pair(a, b)


def scattering_coefficients(reflection_rate: float) -> Tuple[complex, complex, float, float]:
    """Compute scattering matrix elements from reflection rate R."""
    if reflection_rate < 0.0 or reflection_rate > 1.0:
        raise ValueError(f"reflection rate must be in [0, 1]: {reflection_rate}")
    delta_f = 2.0 * math.asin(math.sqrt(reflection_rate))
    half_delta = 0.5 * delta_f
    phase = complex(math.cos(half_delta), math.sin(half_delta))
    t = phase * math.cos(half_delta)
    r = -1j * phase * math.sin(half_delta)
    return complex(t), complex(r), float(abs(t) ** 2), float(abs(r) ** 2)


def condition_id(phi: float, s0: float) -> str:
    return f"phi{phi / math.pi:.6g}_s{s0:.6g}"


def run_time_series(
    reflection_rate: float,
    phi: float,
    steps: int,
    s0: float = 1.0e-2,
) -> Dict[str, Any]:
    """Run AB scattering for fixed R, phi, return time series state and S values."""
    a, b = state_from_s_phi(s0, phi)
    t, r, _t_power, _r_power = scattering_coefficients(reflection_rate)

    states_a = [a]
    states_b = [b]
    s_values = []

    for step in range(steps):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s_value = (p_a - p_b) / q
        s_values.append(s_value)

        a, b = normalize_pair(r * a + t * b, t * a + r * b)
        states_a.append(a)
        states_b.append(b)

    return {
        "condition_id": condition_id(phi, s0),
        "R": reflection_rate,
        "phi": phi,
        "s0": s0,
        "states_a": states_a,
        "states_b": states_b,
        "s_values": s_values,
    }


def one_cycle_transfer_matrix(
    reflection_rate: float,
    phi_init: float,
    s0_init: float,
    num_steps: int,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Construct one-cycle transfer matrix M by tracking how the scattering
    affects the 4-component state (A_Re, A_Im, B_Re, B_Im).

    Returns:
      M: 4x4 transfer matrix (complex components as reals)
      metadata: dict with Tr, det, eigenvalues
    """
    a_init, b_init = state_from_s_phi(s0_init, phi_init)
    t, r, _, _ = scattering_coefficients(reflection_rate)

    # Extract real/imaginary components of initial state
    x0 = np.array([
        a_init.real, a_init.imag,
        b_init.real, b_init.imag
    ], dtype=np.complex128)

    # Run scattering and accumulate changes
    a, b = a_init, b_init
    states = [np.array([a.real, a.imag, b.real, b.imag])]

    for _ in range(num_steps):
        a, b = normalize_pair(r * a + t * b, t * a + r * b)
        states.append(np.array([a.real, a.imag, b.real, b.imag]))

    # Fit M such that x_{n+1} ≈ M @ x_n (least squares over all transitions)
    X_before = np.column_stack(states[:-1])  # 4 x num_steps
    X_after = np.column_stack(states[1:])    # 4 x num_steps

    M, residuals, rank, singular = np.linalg.lstsq(X_before.T, X_after.T, rcond=None)
    M = M.T  # 4x4 matrix

    # Compute invariants
    trace = np.trace(M)
    det = np.linalg.det(M)
    eigenvalues = np.linalg.eigvals(M)
    eigenvalue_mags = np.abs(eigenvalues)

    return M, {
        "R": reflection_rate,
        "condition_id": condition_id(phi_init, s0_init),
        "trace": float(np.real(trace)),
        "determinant": float(np.real(det)),
        "eigenvalues": [complex(ev) for ev in eigenvalues],
        "eigenvalue_mags": list(eigenvalue_mags),
        "max_eigenvalue_mag": float(np.max(eigenvalue_mags)),
        "min_eigenvalue_mag": float(np.min(eigenvalue_mags)),
    }


def floquet_discriminant(M: np.ndarray) -> float:
    """
    Floquet discriminant = Tr(M) for 2D periodic systems,
    generalized here for 4D.
    """
    return float(np.real(np.trace(M)))


def sweep_invariants(
    reflection_rates: List[float],
    harmonics: List[Tuple[float, float]],  # [(phi, s0), ...]
    num_steps: int,
) -> List[Dict[str, Any]]:
    """
    Sweep R values and compute invariants for each (R, harmonic) pair.

    Args:
      reflection_rates: R values to test
      harmonics: list of (phi, s0) initial conditions
      num_steps: iteration count per condition

    Returns:
      List of result dicts with Tr, det, |λ| for each condition
    """
    results = []

    for R in reflection_rates:
        for phi, s0 in harmonics:
            try:
                M, metadata = one_cycle_transfer_matrix(R, phi, s0, num_steps)
                results.append(metadata)
            except Exception as e:
                print(f"Error at R={R:.10f}, phi={phi:.6f}: {e}", file=sys.stderr)

    return results


def identify_common_points(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Find R values where invariants (trace, determinant, eigenvalue mags)
    are consistent across multiple harmonics.
    """
    # Group by R
    by_r = {}
    for res in results:
        R = res["R"]
        if R not in by_r:
            by_r[R] = []
        by_r[R].append(res)

    # For each R, check if invariants align
    common_candidates = []
    for R in sorted(by_r.keys()):
        conditions = by_r[R]
        if len(conditions) < 2:
            continue  # Need at least 2 harmonics to compare

        traces = [c["trace"] for c in conditions]
        dets = [c["determinant"] for c in conditions]
        max_mags = [c["max_eigenvalue_mag"] for c in conditions]

        trace_std = float(np.std(traces))
        det_std = float(np.std(dets))
        mag_std = float(np.std(max_mags))

        # Flag if variance is small (< 0.1%)
        trace_consensus = trace_std / (1.0 + abs(np.mean(traces))) < 0.001
        det_consensus = det_std / (1.0 + abs(np.mean(dets))) < 0.001
        mag_consensus = mag_std / (1.0 + np.mean(max_mags)) < 0.001

        common_candidates.append({
            "R": R,
            "num_conditions": len(conditions),
            "trace_mean": float(np.mean(traces)),
            "trace_std": trace_std,
            "trace_consensus": trace_consensus,
            "det_mean": float(np.mean(dets)),
            "det_std": det_std,
            "det_consensus": det_consensus,
            "max_mag_mean": float(np.mean(max_mags)),
            "max_mag_std": mag_std,
            "mag_consensus": mag_consensus,
        })

    return common_candidates


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Section 18.1: Verify common fixed points via transfer matrix invariants"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for CSV and JSON results"
    )
    parser.add_argument(
        "--num-steps",
        type=int,
        default=DEFAULT_STEPS,
        help="Number of iteration steps"
    )
    parser.add_argument(
        "--r-min",
        type=float,
        default=R_SWEEP_MIN,
        help="Minimum R value to sweep"
    )
    parser.add_argument(
        "--r-max",
        type=float,
        default=R_SWEEP_MAX,
        help="Maximum R value to sweep"
    )
    parser.add_argument(
        "--r-points",
        type=int,
        default=R_POINTS,
        help="Number of R points in sweep"
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # R sweep
    R_values = np.linspace(args.r_min, args.r_max, args.r_points)

    # Multiple harmonic conditions to test common-point hypothesis
    harmonics = [
        (0.0, 0.01),           # phi=0, s0=0.01
        (math.pi, 0.01),       # phi=π, s0=0.01
        (0.0, 0.005),          # phi=0, s0=0.005
        (math.pi / 2, 0.01),   # phi=π/2, s0=0.01
    ]

    print(f"Sweeping {len(R_values)} R values × {len(harmonics)} harmonics...")
    print(f"R range: [{args.r_min:.10f}, {args.r_max:.10f}]")
    print(f"Number of iteration steps: {args.num_steps}")

    # Run sweep
    results = sweep_invariants(R_values, harmonics, args.num_steps)

    print(f"Completed {len(results)} condition evaluations.")

    # Save detailed results
    csv_path = args.output_dir / "section18_invariants_detailed_v1.csv"
    with open(csv_path, "w", newline="") as f:
        fieldnames = [
            "R", "condition_id", "trace", "determinant",
            "max_eigenvalue_mag", "min_eigenvalue_mag",
            "eigenvalues_str"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for res in results:
            writer.writerow({
                "R": res["R"],
                "condition_id": res["condition_id"],
                "trace": res["trace"],
                "determinant": res["determinant"],
                "max_eigenvalue_mag": res["max_eigenvalue_mag"],
                "min_eigenvalue_mag": res["min_eigenvalue_mag"],
                "eigenvalues_str": str(res["eigenvalues"]),
            })
    print(f"Saved detailed results to {csv_path}")

    # Identify common points
    common = identify_common_points(results)

    # Save consensus points
    consensus_path = args.output_dir / "section18_common_points_v1.csv"
    with open(consensus_path, "w", newline="") as f:
        fieldnames = [
            "R", "num_conditions",
            "trace_mean", "trace_std", "trace_consensus",
            "det_mean", "det_std", "det_consensus",
            "max_mag_mean", "max_mag_std", "mag_consensus"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in common:
            writer.writerow(row)
    print(f"Saved consensus points to {consensus_path}")

    # JSON summary
    summary = {
        "model": "section18_transfer_matrix_invariants_v1",
        "parameters": {
            "r_min": args.r_min,
            "r_max": args.r_max,
            "r_points": args.r_points,
            "num_steps": args.num_steps,
            "harmonics": [{"phi": h[0], "s0": h[1]} for h in harmonics],
        },
        "results_count": len(results),
        "common_point_candidates": [
            {k: (bool(v) if isinstance(v, (bool, np.bool_)) else v) for k, v in row.items()}
            for row in common
        ],
    }

    json_path = args.output_dir / "section18_summary_v1.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"Saved summary to {json_path}")

    # Print highlights
    print("\n=== Common Point Candidates (trace_consensus or mag_consensus) ===")
    for row in common:
        if row["trace_consensus"] or row["mag_consensus"]:
            print(
                f"R={row['R']:.10f}  "
                f"trace_mean={row['trace_mean']:+.6f} (std={row['trace_std']:.2e})  "
                f"max_mag={row['max_mag_mean']:.6f} (std={row['max_mag_std']:.2e})"
            )


if __name__ == "__main__":
    main()
