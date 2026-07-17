#!/usr/bin/env python3
"""Phase 5: Eigenvalue and Fixed Point Analysis

Derive the stability conditions from transfer matrix eigenvalues.
Test whether R_137 and R_128 emerge as natural stable points from first principles.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np
from scipy import linalg


DEFAULT_STEPS = 256
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


def compute_transfer_matrix(
    reflection_rate: float,
    phi: float = 0.0,
    s0: float = 0.01,
    num_steps: int = DEFAULT_STEPS,
) -> np.ndarray:
    """Construct 4×4 transfer matrix from time evolution."""
    a_init, b_init = state_from_s_phi(s0, phi)
    t, r, _, _ = scattering_coefficients(reflection_rate)

    # Initial state in 4D
    x0 = np.array([a_init.real, a_init.imag, b_init.real, b_init.imag], dtype=np.float64)

    # Record states
    states = [x0.copy()]
    a, b = a_init, b_init

    for _ in range(num_steps):
        a, b = normalize_pair(r * a + t * b, t * a + r * b)
        x = np.array([a.real, a.imag, b.real, b.imag], dtype=np.float64)
        states.append(x)

    states = np.array(states)

    # Fit M: x_{n+1} ≈ M @ x_n (least squares)
    X_before = states[:-1].T  # 4 × num_steps
    X_after = states[1:].T    # 4 × num_steps

    M, residuals, rank, sv = np.linalg.lstsq(X_before.T, X_after.T, rcond=None)
    M = M.T  # 4×4

    return M


def analyze_eigenvalues(
    M: np.ndarray,
) -> Dict[str, Any]:
    """Analyze eigenvalues and stability."""
    eigenvalues = np.linalg.eigvals(M)
    eigenvalue_mags = np.abs(eigenvalues)
    eigenvalue_phases = np.angle(eigenvalues)

    trace = np.trace(M)
    determinant = np.linalg.det(M)

    # Floquet discriminant (stability criterion)
    floquet_discr = abs(trace) - 2

    # Lyapunov exponent (largest eigenvalue magnitude)
    lyapunov = np.log(np.max(eigenvalue_mags))

    # Resonance analysis: λⁿ ≈ 1
    resonances = []
    for n in range(1, 256):
        for i, lam in enumerate(eigenvalues):
            if abs(lam**n - 1) < 0.1:
                resonances.append({
                    "eigenvalue_index": i,
                    "n": n,
                    "lambda_n": complex(lam**n),
                    "distance_to_1": abs(lam**n - 1),
                })

    return {
        "trace": complex(trace),
        "determinant": complex(determinant),
        "eigenvalues": [complex(ev) for ev in eigenvalues],
        "eigenvalue_mags": list(eigenvalue_mags),
        "eigenvalue_phases": list(eigenvalue_phases),
        "max_mag": float(np.max(eigenvalue_mags)),
        "min_mag": float(np.min(eigenvalue_mags)),
        "floquet_discriminant": float(floquet_discr),
        "lyapunov_exponent": float(lyapunov),
        "is_stable": float(abs(trace)) < 2.0 or np.max(eigenvalue_mags) < 1.0,
        "resonances": resonances[:5],  # Top 5 resonances
    }


def sweep_R_eigenvalue(
    R_min: float = 0.686,
    R_max: float = 0.700,
    R_points: int = 141,
    num_steps: int = DEFAULT_STEPS,
) -> List[Dict[str, Any]]:
    """Sweep R and compute eigenvalue properties."""
    R_values = np.linspace(R_min, R_max, R_points)
    results = []

    print(f"Sweeping {len(R_values)} R values...")

    for i, R in enumerate(R_values):
        if i % 20 == 0:
            print(f"  {i}/{len(R_values)}")

        M = compute_transfer_matrix(R, num_steps=num_steps)
        analysis = analyze_eigenvalues(M)

        results.append({
            "R": R,
            "trace": str(analysis["trace"]),
            "determinant": str(analysis["determinant"]),
            "trace_real": float(analysis["trace"].real),
            "trace_imag": float(analysis["trace"].imag),
            "det_real": float(analysis["determinant"].real),
            "det_imag": float(analysis["determinant"].imag),
            "max_eigenvalue_mag": analysis["max_mag"],
            "min_eigenvalue_mag": analysis["min_mag"],
            "floquet_discriminant": analysis["floquet_discriminant"],
            "lyapunov_exponent": analysis["lyapunov_exponent"],
            "is_stable": analysis["is_stable"],
            "num_resonances": len(analysis["resonances"]),
        })

    return results


def identify_stable_points(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Find R values where stability conditions are satisfied."""
    # Floquet discriminant = 0 indicates stability boundary
    # |Tr(M)| = 2 is the stability criterion

    stable_candidates = []

    for r in results:
        floquet = abs(r["floquet_discriminant"])

        # Threshold for "near boundary"
        if floquet < 0.5:  # |Tr(M)| between 1.5 and 2.5
            stable_candidates.append({
                "R": r["R"],
                "floquet_discriminant": r["floquet_discriminant"],
                "trace_real": r["trace_real"],
                "max_mag": r["max_eigenvalue_mag"],
                "lyapunov": r["lyapunov_exponent"],
            })

    # Sort by floquet discriminant magnitude
    stable_candidates.sort(key=lambda x: abs(x["floquet_discriminant"]))

    return stable_candidates


def check_resonance_pattern(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze resonance patterns (λⁿ = 1)."""
    resonance_map = {}

    for r in results:
        R = r["R"]
        num_res = r["num_resonances"]

        if num_res > 0:
            if num_res not in resonance_map:
                resonance_map[num_res] = []
            resonance_map[num_res].append(R)

    # Find "cluster" of resonances
    resonance_summary = {
        f"{num_res}_resonances": R_list
        for num_res, R_list in sorted(resonance_map.items())
    }

    return resonance_summary


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Phase 5: Eigenvalue and Fixed Point Analysis")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Output directory")
    parser.add_argument("--r-min", type=float, default=0.686, help="Min R")
    parser.add_argument("--r-max", type=float, default=0.700, help="Max R")
    parser.add_argument("--r-points", type=int, default=141, help="Number of R points")
    parser.add_argument("--num-steps", type=int, default=DEFAULT_STEPS, help="Transfer matrix steps")

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Phase 5: Eigenvalue and Fixed Point Analysis")
    print(f"R range: [{args.r_min:.10f}, {args.r_max:.10f}]")
    print(f"R points: {args.r_points}")
    print()

    # Sweep eigenvalues
    results = sweep_R_eigenvalue(args.r_min, args.r_max, args.r_points, args.num_steps)

    # Save detailed results
    csv_path = args.output_dir / "phase5_eigenvalue_sweep_v1.csv"
    with open(csv_path, "w", newline="") as f:
        fieldnames = [
            "R", "trace_real", "trace_imag", "det_real", "det_imag",
            "max_eigenvalue_mag", "min_eigenvalue_mag",
            "floquet_discriminant", "lyapunov_exponent", "is_stable", "num_resonances"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        # Filter to only include the fieldnames
        for row in results:
            filtered = {k: v for k, v in row.items() if k in fieldnames}
            writer.writerow(filtered)
    print(f"Saved eigenvalue sweep to {csv_path}")

    # Identify stable points
    stable = identify_stable_points(results)

    stable_csv = args.output_dir / "phase5_stable_candidates_v1.csv"
    with open(stable_csv, "w", newline="") as f:
        fieldnames = ["R", "floquet_discriminant", "trace_real", "max_mag", "lyapunov"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(stable)
    print(f"Saved stable candidates to {stable_csv}")

    # Resonance analysis
    resonance_map = check_resonance_pattern(results)

    # Print summary
    print("\n=== Eigenvalue Analysis Summary ===")
    print(f"Total R points analyzed: {len(results)}")
    print(f"Stable candidates (|Floquet| < 0.5): {len(stable)}")

    if len(stable) > 0:
        print("\nTop stable candidates:")
        for i, s in enumerate(stable[:5]):
            print(f"  {i+1}. R={s['R']:.10f}  Floquet={s['floquet_discriminant']:+.6f}  max|λ|={s['max_mag']:.6f}")

    # Compare with R_137 and R_128
    print(f"\n=== Comparison with Targets ===")
    for target_name, target_R in [("R_137", R_137_CANDIDATE), ("R_128", R_128_CANDIDATE)]:
        # Find nearest result
        nearest = min(results, key=lambda r: abs(r["R"] - target_R))
        print(f"{target_name}={target_R:.10f}")
        print(f"  Nearest result R={nearest['R']:.10f}")
        print(f"    Floquet={nearest['floquet_discriminant']:+.6f}")
        print(f"    max|λ|={nearest['max_eigenvalue_mag']:.6f}")
        print(f"    Lyapunov={nearest['lyapunov_exponent']:.6f}")

    # Summary JSON
    summary = {
        "model": "phase5_eigenvalue_fixed_point_v1",
        "parameters": {
            "r_min": args.r_min,
            "r_max": args.r_max,
            "r_points": args.r_points,
            "num_steps": args.num_steps,
        },
        "results_count": len(results),
        "stable_candidates": len(stable),
        "top_stable_R": [s["R"] for s in stable[:3]],
    }

    summary_json = args.output_dir / "phase5_summary_v1.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved summary to {summary_json}")

    print("\nPhase 5 complete.")


if __name__ == "__main__":
    main()
