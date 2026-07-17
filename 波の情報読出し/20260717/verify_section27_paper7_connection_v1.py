#!/usr/bin/env python3
"""Section 27: Verify Paper7 R=3 Complete Inclusion Connection

This script tests the hypothesis that the two experimental peaks (R_128, R_137)
correspond to the dynamics realization of Paper7's 4D complete inclusion problem.

Paper7 defines:
  - State space: (a, b, c, d) ∈ ℝ⁴  (components of A, B ∈ ℂ)
  - Sphere radius: R = 3
  - Cell half-width: 1/2
  - Complete inclusion condition: Σ(|c_j| + 1/2)² ≤ 3²
  - Cell count: 137

Current experiment dynamics:
  - State components: (A_Re, A_Im, B_Re, B_Im) ∈ ℝ⁴
  - Closure condition: A² + B² = 0  (dynamical constraint)
  - Stability: reflected at certain R values

Verification strategy:
  1. Check if closure residual E_close(n, R) → 0 only at R_128, R_137
  2. Check if the effective state radius matches Paper7's R=3 interpretation
  3. Check if stable region boundaries align with complete inclusion
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np


DEFAULT_STEPS = 512
R_137_CANDIDATE = 0.697177902556148
R_128_CANDIDATE = 0.688363902556148

# Paper7 parameters
PAPER7_RADIUS = 3.0
PAPER7_CELL_HALFWIDTH = 0.5
PAPER7_TARGET_CELL_COUNT = 137


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


def run_time_series(
    reflection_rate: float,
    phi: float = 0.0,
    steps: int = DEFAULT_STEPS,
    s0: float = 1.0e-2,
) -> Dict[str, Any]:
    """Run AB scattering and record closure residuals at each step."""
    a, b = state_from_s_phi(s0, phi)
    t, r, _, _ = scattering_coefficients(reflection_rate)

    states_a = [a]
    states_b = [b]
    closure_residuals = []  # E_close = |A^2 + B^2| / (|A|^2 + |B|^2)
    state_radii = []

    for step in range(steps + 1):
        # Closure residual: measure departure from A^2 + B^2 = 0
        closure = abs(a**2 + b**2)
        norm = abs(a)**2 + abs(b)**2
        if norm > 0:
            e_close = closure / norm
        else:
            e_close = 0.0
        closure_residuals.append(e_close)

        # State radius in 4D: norm of (A_Re, A_Im, B_Re, B_Im)
        state_4d = np.array([a.real, a.imag, b.real, b.imag])
        radius = np.linalg.norm(state_4d)
        state_radii.append(radius)

        if step < steps:
            a, b = normalize_pair(r * a + t * b, t * a + r * b)
            states_a.append(a)
            states_b.append(b)

    return {
        "R": reflection_rate,
        "phi": phi,
        "s0": s0,
        "steps": steps,
        "states_a": states_a,
        "states_b": states_b,
        "closure_residuals": closure_residuals,
        "state_radii": state_radii,
        "closure_mean": float(np.mean(closure_residuals)),
        "closure_std": float(np.std(closure_residuals)),
        "closure_min": float(np.min(closure_residuals)),
        "radius_mean": float(np.mean(state_radii)),
        "radius_std": float(np.std(state_radii)),
        "radius_max": float(np.max(state_radii)),
    }


def sweep_closure_residuals(
    R_values: List[float],
    phi: float = 0.0,
    steps: int = DEFAULT_STEPS,
) -> List[Dict[str, Any]]:
    """Sweep R and compute closure residuals for each."""
    results = []
    for R in R_values:
        data = run_time_series(R, phi, steps)
        results.append({
            "R": R,
            "closure_mean": data["closure_mean"],
            "closure_std": data["closure_std"],
            "closure_min": data["closure_min"],
            "radius_mean": data["radius_mean"],
            "radius_std": data["radius_std"],
            "radius_max": data["radius_max"],
        })
    return results


def identify_closure_minima(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Find R values where closure residual E_close is minimized."""
    # Sort by closure_min
    sorted_res = sorted(results, key=lambda x: x["closure_min"])

    # Find local minima (if closure_min changes from decreasing to increasing)
    minima = []
    for i in range(1, len(sorted_res) - 1):
        prev_R = sorted_res[i-1]["R"]
        curr_R = sorted_res[i]["R"]
        next_R = sorted_res[i+1]["R"]
        prev_close = sorted_res[i-1]["closure_min"]
        curr_close = sorted_res[i]["closure_min"]
        next_close = sorted_res[i+1]["closure_min"]

        if curr_close < prev_close and curr_close < next_close:
            minima.append({
                "R": curr_R,
                "closure_min": curr_close,
                "closure_mean": sorted_res[i]["closure_mean"],
                "radius_mean": sorted_res[i]["radius_mean"],
            })

    return minima


def check_paper7_alignment(
    closure_minima: List[Dict[str, Any]],
    R_128_target: float = R_128_CANDIDATE,
    R_137_target: float = R_137_CANDIDATE,
) -> Dict[str, Any]:
    """
    Check if closure minima align with Paper7 predictions.

    Hypothesis: R_128 and R_137 are the two R values where
    the dynamical system most closely satisfies A^2 + B^2 = 0.
    """
    if len(closure_minima) == 0:
        return {
            "paper7_alignment": "NO_MINIMA_FOUND",
            "candidates": [],
        }

    # Sort by closure_min to find deepest minima
    sorted_minima = sorted(closure_minima, key=lambda x: x["closure_min"])

    # Check if top 2 minima align with R_128, R_137
    top_2 = sorted_minima[:2]

    matches_137 = []
    matches_128 = []

    for peak in top_2:
        if abs(peak["R"] - R_137_target) < 0.002:
            matches_137.append(peak)
        if abs(peak["R"] - R_128_target) < 0.002:
            matches_128.append(peak)

    if matches_137 and matches_128:
        alignment = "BOTH_PEAKS_ALIGNED"
    elif matches_137:
        alignment = "ONLY_137_ALIGNED"
    elif matches_128:
        alignment = "ONLY_128_ALIGNED"
    else:
        alignment = "NO_ALIGNMENT_WITH_TARGETS"

    return {
        "paper7_alignment": alignment,
        "target_R_137": R_137_target,
        "target_R_128": R_128_target,
        "top_2_minima": top_2,
        "candidates_near_137": matches_137,
        "candidates_near_128": matches_128,
    }


def radius_statistic_analysis(
    R_values: List[float],
    sweep_results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Analyze if state radius (norm of 4D state) is constrained
    to Paper7-like effective radius R~3.
    """
    radii = [r["radius_mean"] for r in sweep_results]
    radius_min = min(radii)
    radius_max = max(radii)
    radius_mean = float(np.mean(radii))
    radius_std = float(np.std(radii))

    # Check if any R shows especially small or large radius
    # (might indicate different attractor basins)
    r_and_rad = list(zip(R_values, radii))
    r_and_rad.sort(key=lambda x: x[1])

    return {
        "radius_range": [radius_min, radius_max],
        "radius_mean": radius_mean,
        "radius_std": radius_std,
        "paper7_radius": PAPER7_RADIUS,
        "radius_to_paper7_ratio": radius_mean / PAPER7_RADIUS,
        "extreme_radii": {
            "minimum_R": r_and_rad[0][0],
            "minimum_radius": r_and_rad[0][1],
            "maximum_R": r_and_rad[-1][0],
            "maximum_radius": r_and_rad[-1][1],
        },
    }


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Section 27: Verify Paper7 R=3 complete inclusion connection"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory"
    )
    parser.add_argument(
        "--r-min",
        type=float,
        default=0.685,
        help="Min R value"
    )
    parser.add_argument(
        "--r-max",
        type=float,
        default=0.702,
        help="Max R value"
    )
    parser.add_argument(
        "--r-points",
        type=int,
        default=171,
        help="Number of R points"
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_STEPS,
        help="Iterations per condition"
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    R_values = np.linspace(args.r_min, args.r_max, args.r_points)

    print(f"Section 27: Paper7 Connection Verification")
    print(f"R range: [{args.r_min:.10f}, {args.r_max:.10f}]")
    print(f"R points: {len(R_values)}")
    print(f"Iterations per point: {args.steps}")
    print()

    # Sweep closure residuals
    print("Computing closure residuals E_close = |A^2 + B^2| / (|A|^2 + |B|^2)...")
    sweep_results = sweep_closure_residuals(R_values.tolist(), steps=args.steps)

    # Save sweep results
    sweep_csv = args.output_dir / "section27_closure_residual_sweep_v1.csv"
    with open(sweep_csv, "w", newline="") as f:
        fieldnames = [
            "R", "closure_mean", "closure_std", "closure_min",
            "radius_mean", "radius_std", "radius_max"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sweep_results)
    print(f"Saved sweep results to {sweep_csv}")

    # Identify closure minima
    print("Identifying closure residual minima...")
    closure_minima = identify_closure_minima(sweep_results)
    print(f"Found {len(closure_minima)} local minima")

    # Check Paper7 alignment
    print("Checking Paper7 alignment...")
    alignment = check_paper7_alignment(closure_minima)

    alignment_json = args.output_dir / "section27_paper7_alignment_v1.json"
    with open(alignment_json, "w") as f:
        json.dump(alignment, f, indent=2, default=str)
    print(f"Saved alignment analysis to {alignment_json}")

    print(f"Result: {alignment['paper7_alignment']}")

    # Radius statistics
    print("Analyzing state radius statistics...")
    radius_stats = radius_statistic_analysis(R_values.tolist(), sweep_results)

    print(f"State radius range: [{radius_stats['radius_range'][0]:.6f}, {radius_stats['radius_range'][1]:.6f}]")
    print(f"Paper7 reference radius: {radius_stats['paper7_radius']:.1f}")
    print(f"Ratio (mean/Paper7): {radius_stats['radius_to_paper7_ratio']:.4f}")

    radius_json = args.output_dir / "section27_radius_analysis_v1.json"
    with open(radius_json, "w") as f:
        json.dump(radius_stats, f, indent=2)
    print(f"Saved radius analysis to {radius_json}")

    # Summary
    summary = {
        "model": "section27_paper7_connection_v1",
        "parameters": {
            "r_min": args.r_min,
            "r_max": args.r_max,
            "r_points": args.r_points,
            "steps": args.steps,
        },
        "paper7_alignment": alignment["paper7_alignment"],
        "top_2_minima": alignment.get("top_2_minima", []),
        "radius_analysis": radius_stats,
    }

    summary_json = args.output_dir / "section27_summary_v1.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved summary to {summary_json}")


if __name__ == "__main__":
    main()
