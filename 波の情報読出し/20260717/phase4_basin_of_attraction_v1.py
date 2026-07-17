#!/usr/bin/env python3
"""DEPRECATED Phase 4 v1 endpoint-comparison experiment.

This historical implementation keeps R fixed and compares two independently
evolved endpoints.  It does not measure attraction in R-space, and its
classifier always returns one of two labels.  It is retained only as provenance.
Use ``phase4_basin_of_attraction_v2.py`` for the H8-conditional basin analysis.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


DEFAULT_STEPS = 512
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


def evolve_trajectory(
    reflection_rate: float,
    phi: float,
    s0: float,
    steps: int,
) -> Dict[str, Any]:
    """Run AB scattering and return final state."""
    a, b = state_from_s_phi(s0, phi)
    t, r, _, _ = scattering_coefficients(reflection_rate)

    for _ in range(steps):
        a, b = normalize_pair(r * a + t * b, t * a + r * b)

    return {
        "R": reflection_rate,
        "phi": phi,
        "s0": s0,
        "final_a": a,
        "final_b": b,
        "final_a_norm": abs(a),
        "final_b_norm": abs(b),
        "final_a_phase": np.angle(a),
        "final_b_phase": np.angle(b),
    }


def classify_attractor(
    state_137: Dict[str, Any],
    state_128: Dict[str, Any],
    threshold: float = 0.1,
) -> Tuple[str, float]:
    """
    Classify whether the state converges to 137, 128, or neither.

    Returns:
      ("137", confidence), ("128", confidence), or ("ambiguous", None)
    """
    a_137 = complex(state_137["final_a"])
    b_137 = complex(state_137["final_b"])
    a_128 = complex(state_128["final_a"])
    b_128 = complex(state_128["final_b"])

    # Distance to 137 attractor
    dist_to_137 = abs(a_137 - a_128) + abs(b_137 - b_128)

    # Distance to 128 attractor (via symmetry: should be close to (b_137, a_137))
    dist_to_128_via_symmetry = abs(a_137 - b_128) + abs(b_137 - a_128)

    # A↔B swapped?
    is_swapped = dist_to_128_via_symmetry < dist_to_137

    if is_swapped:
        confidence = 1.0 - min(dist_to_128_via_symmetry, 1.0)
        return "128", confidence
    else:
        confidence = 1.0 - min(dist_to_137, 1.0)
        return "137", confidence


def sweep_basin(
    phi_points: int = 100,
    s0_points: int = 50,
    phi_range: Tuple[float, float] = (0.0, 2 * np.pi),
    s0_range: Tuple[float, float] = (0.001, 0.05),
    steps: int = DEFAULT_STEPS,
) -> List[Dict[str, Any]]:
    """Sweep (φ, s₀) parameter space and classify attractor."""
    phis = np.linspace(phi_range[0], phi_range[1], phi_points)
    s0s = np.linspace(s0_range[0], s0_range[1], s0_points)

    results = []
    total = phi_points * s0_points
    count = 0

    print(f"Sweeping {total} initial conditions...")

    for phi in phis:
        for s0 in s0s:
            count += 1
            if count % max(1, total // 20) == 0:
                print(f"  {count}/{total} ({100*count/total:.0f}%)")

            # Evolve at both R values
            state_137 = evolve_trajectory(R_137_CANDIDATE, phi, s0, steps)
            state_128 = evolve_trajectory(R_128_CANDIDATE, phi, s0, steps)

            # Classify
            attractor, confidence = classify_attractor(state_137, state_128)

            results.append({
                "phi": phi,
                "phi_over_pi": phi / np.pi,
                "s0": s0,
                "attractor": attractor,
                "confidence": confidence,
            })

    return results


def compute_basin_statistics(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute statistics about basin structure."""
    count_137 = sum(1 for r in results if r["attractor"] == "137")
    count_128 = sum(1 for r in results if r["attractor"] == "128")
    count_ambig = len(results) - count_137 - count_128

    # Check A↔B symmetry
    # Group by (φ mod π, s₀) and check if φ and φ+π have swapped attractors
    symmetry_violations = 0
    symmetry_checks = 0

    phi_groups = {}
    for r in results:
        key = (round(r["s0"], 6), round(r["phi_over_pi"] % 1, 3))
        if key not in phi_groups:
            phi_groups[key] = []
        phi_groups[key].append(r)

    for group in phi_groups.values():
        if len(group) >= 2:
            group.sort(key=lambda x: x["phi"])
            for i in range(len(group) - 1):
                # Check if φ and nearby φ have consistent attractor
                if abs(group[i]["phi"] - group[i+1]["phi"]) < 0.2:
                    if group[i]["attractor"] != group[i+1]["attractor"]:
                        symmetry_violations += 1
                    symmetry_checks += 1

    symmetry_ratio = symmetry_violations / max(1, symmetry_checks)

    return {
        "total_points": len(results),
        "count_137": count_137,
        "count_128": count_128,
        "count_ambiguous": count_ambig,
        "fraction_137": count_137 / len(results),
        "fraction_128": count_128 / len(results),
        "fraction_ambiguous": count_ambig / len(results),
        "symmetry_violations": symmetry_violations,
        "symmetry_checks": symmetry_checks,
        "symmetry_violation_ratio": symmetry_ratio,
    }


def plot_basin_heatmap(results: List[Dict[str, Any]], output_path: Path):
    """Plot basin as 2D heatmap."""
    # Convert to grid
    phis = sorted(set(r["phi"] for r in results))
    s0s = sorted(set(r["s0"] for r in results))

    phi_to_idx = {phi: i for i, phi in enumerate(phis)}
    s0_to_idx = {s0: i for i, s0 in enumerate(s0s)}

    # Color mapping: 0=128 (blue), 1=137 (red), 2=ambiguous (gray)
    grid = np.full((len(s0s), len(phis)), 2, dtype=int)
    confidence_grid = np.full((len(s0s), len(phis)), 0.5, dtype=float)

    for r in results:
        i_phi = phi_to_idx[r["phi"]]
        i_s0 = s0_to_idx[r["s0"]]
        if r["attractor"] == "137":
            grid[i_s0, i_phi] = 1
        elif r["attractor"] == "128":
            grid[i_s0, i_phi] = 0
        else:
            grid[i_s0, i_phi] = 2
        confidence_grid[i_s0, i_phi] = r["confidence"]

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Heatmap 1: Attractor classification
    cmap = ListedColormap(["blue", "red", "lightgray"])
    im1 = axes[0].imshow(
        grid,
        extent=[0, 2*np.pi, s0s[0], s0s[-1]],
        aspect="auto",
        origin="lower",
        cmap=cmap,
    )
    axes[0].set_xlabel(r"$\phi$ (radians)")
    axes[0].set_ylabel(r"$s_0$")
    axes[0].set_title("Basin of Attraction\nRed=R_137, Blue=R_128")
    axes[0].set_xticks([0, np.pi/2, np.pi, 1.5*np.pi, 2*np.pi])
    axes[0].set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])

    # Heatmap 2: Confidence
    im2 = axes[1].imshow(
        confidence_grid,
        extent=[0, 2*np.pi, s0s[0], s0s[-1]],
        aspect="auto",
        origin="lower",
        cmap="viridis",
        vmin=0, vmax=1,
    )
    axes[1].set_xlabel(r"$\phi$ (radians)")
    axes[1].set_ylabel(r"$s_0$")
    axes[1].set_title("Convergence Confidence")
    axes[1].set_xticks([0, np.pi/2, np.pi, 1.5*np.pi, 2*np.pi])
    axes[1].set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Saved basin plot to {output_path}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Phase 4: Basin of Attraction")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Output directory")
    parser.add_argument("--phi-points", type=int, default=100, help="Number of φ points")
    parser.add_argument("--s0-points", type=int, default=50, help="Number of s₀ points")
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS, help="Evolution steps")

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Phase 4: Basin of Attraction Analysis")
    print(f"Parameters: {args.phi_points} φ × {args.s0_points} s₀ = {args.phi_points * args.s0_points} initial conditions")
    print(f"R_137={R_137_CANDIDATE:.10f}")
    print(f"R_128={R_128_CANDIDATE:.10f}")
    print()

    # Run sweep
    results = sweep_basin(
        phi_points=args.phi_points,
        s0_points=args.s0_points,
        steps=args.steps,
    )

    # Save results
    csv_path = args.output_dir / "phase4_basin_sweep_v1.csv"
    with open(csv_path, "w", newline="") as f:
        fieldnames = ["phi", "phi_over_pi", "s0", "attractor", "confidence"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"Saved basin sweep to {csv_path}")

    # Compute statistics
    stats = compute_basin_statistics(results)

    stats_json = args.output_dir / "phase4_basin_statistics_v1.json"
    with open(stats_json, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Saved basin statistics to {stats_json}")

    # Print summary
    print("\n=== Basin Statistics ===")
    print(f"Total points: {stats['total_points']}")
    print(f"137 attractor: {stats['count_137']} ({stats['fraction_137']*100:.1f}%)")
    print(f"128 attractor: {stats['count_128']} ({stats['fraction_128']*100:.1f}%)")
    print(f"Ambiguous: {stats['count_ambiguous']} ({stats['fraction_ambiguous']*100:.1f}%)")
    print(f"Symmetry violation ratio: {stats['symmetry_violation_ratio']:.4f}")

    # Plot
    plot_path = args.output_dir / "phase4_basin_heatmap_v1.png"
    plot_basin_heatmap(results, plot_path)

    print("\nPhase 4 complete.")


if __name__ == "__main__":
    main()
