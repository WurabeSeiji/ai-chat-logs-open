#!/usr/bin/env python3
"""Phase 6: conditional N(R) = 4π/(1-R)² readout under H2 and D4.

Evaluate the independently selected empirical candidates R1 and R2, and keep
them distinct from inverse images obtained by inserting external diagnostics
into the same conditional formula.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_STEPS = 512
R1_CANDIDATE = 0.697177902556148
R2_CANDIDATE = 0.688363902556148

# External diagnostic values; they are not used to select R1 or R2.
ALPHA_INV_ZERO = 137.035999177
ALPHA_INV_MZ = 128.946


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


def compute_state_radius_distribution(
    reflection_rate: float,
    num_initial_conditions: int = 100,
    steps: int = DEFAULT_STEPS,
) -> Dict[str, Any]:
    """Compute effective state space radius for given R."""
    radii = []
    trajectories = []

    phi_vals = np.linspace(0, 2*np.pi, int(np.sqrt(num_initial_conditions)))
    s0_vals = np.linspace(0.001, 0.05, int(np.sqrt(num_initial_conditions)))

    for phi in phi_vals:
        for s0 in s0_vals:
            a, b = state_from_s_phi(s0, phi)
            t, r, _, _ = scattering_coefficients(reflection_rate)

            a_traj = [a]
            b_traj = [b]

            for _ in range(steps):
                a, b = normalize_pair(r * a + t * b, t * a + r * b)
                a_traj.append(a)
                b_traj.append(b)

            # Extract final 100 states for statistics
            final_states = list(zip(a_traj[-100:], b_traj[-100:]))
            trajectories.extend(final_states)

    # Compute radii in 4D space
    for a, b in trajectories:
        state_4d = np.array([a.real, a.imag, b.real, b.imag])
        radius = np.linalg.norm(state_4d)
        radii.append(radius)

    radii = np.array(radii)

    return {
        "R": reflection_rate,
        "num_states": len(radii),
        "radius_mean": float(np.mean(radii)),
        "radius_std": float(np.std(radii)),
        "radius_min": float(np.min(radii)),
        "radius_max": float(np.max(radii)),
        "radius_median": float(np.median(radii)),
    }


def compute_poincare_section_density(
    reflection_rate: float,
    num_samples: int = 10000,
    steps: int = DEFAULT_STEPS,
) -> Dict[str, Any]:
    """Count states on Poincaré section (e.g., Im(A) = 0)."""
    crossing_count = 0
    crossing_phases = []

    for i in range(num_samples):
        phi = 2 * np.pi * i / num_samples
        s0 = 0.01  # Fixed s0
        a, b = state_from_s_phi(s0, phi)
        t, r, _, _ = scattering_coefficients(reflection_rate)

        a_prev = a
        for step in range(steps):
            a, b = normalize_pair(r * a + t * b, t * a + r * b)

            # Check if crossed Im(A) = 0
            if (a_prev.imag < 0 and a.imag > 0) or (a_prev.imag > 0 and a.imag < 0):
                crossing_count += 1
                crossing_phases.append(np.angle(b))

            a_prev = a

    return {
        "R": reflection_rate,
        "num_crossings": crossing_count,
        "crossing_density": crossing_count / steps,
    }


def compute_N_theory(R: float) -> float:
    """Compute N(R) = 4π/(1-R)²."""
    return 4 * np.pi / (1 - R)**2


def compute_N_from_formula(R: float) -> float:
    """Evaluate the conditional inverse-coupling readout."""
    N = compute_N_theory(R)
    return N


def inverse_R_from_N(N: float) -> float:
    """Invert the same conditional readout; this is not an independent result."""
    if N <= 0.0:
        raise ValueError(f"N must be positive: {N}")
    return 1.0 - math.sqrt(4.0 * math.pi / N)


def sweep_N_R(
    R_min: float = 0.686,
    R_max: float = 0.700,
    R_points: int = 141,
) -> List[Dict[str, Any]]:
    """Sweep R and compute N(R)."""
    R_values = np.linspace(R_min, R_max, R_points)
    results = []

    print(f"Computing N(R) for {len(R_values)} R values...")

    for i, R in enumerate(R_values):
        if i % 20 == 0:
            print(f"  {i}/{len(R_values)}")

        N_theory = compute_N_theory(R)

        # Also compute radius distribution (lighter computation)
        radius_stats = compute_state_radius_distribution(R, num_initial_conditions=20)

        results.append({
            "R": R,
            "N_theory": N_theory,
            "inv_alpha": N_theory,
            "radius_mean": radius_stats["radius_mean"],
            "radius_std": radius_stats["radius_std"],
        })

    return results


def compare_candidates_and_diagnostics(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compare empirical R candidates with external diagnostic inverse images."""

    def nearest_sweep(target: float) -> Dict[str, float]:
        best = min(results, key=lambda row: abs(row["inv_alpha"] - target))
        return {
            "R": float(best["R"]),
            "N": float(best["N_theory"]),
            "difference": float(best["N_theory"] - target),
        }

    def candidate(candidate_R: float, diagnostic: float) -> Dict[str, float]:
        candidate_N = compute_N_theory(candidate_R)
        return {
            "R": candidate_R,
            "N": candidate_N,
            "diagnostic_N": diagnostic,
            "delta_N": candidate_N - diagnostic,
            "relative_difference": abs(candidate_N - diagnostic) / diagnostic,
        }

    return {
        "classification": "conditional_readout_under_H2_and_D4",
        "diagnostics": {
            "alpha_inverse_zero": {
                "N": ALPHA_INV_ZERO,
                "conditional_inverse_R": inverse_R_from_N(ALPHA_INV_ZERO),
                "nearest_sweep_point": nearest_sweep(ALPHA_INV_ZERO),
            },
            "alpha_inverse_MZ_squared": {
                "N": ALPHA_INV_MZ,
                "conditional_inverse_R": inverse_R_from_N(ALPHA_INV_MZ),
                "nearest_sweep_point": nearest_sweep(ALPHA_INV_MZ),
            },
        },
        "empirical_candidates": {
            "R1": candidate(R1_CANDIDATE, ALPHA_INV_ZERO),
            "R2": candidate(R2_CANDIDATE, ALPHA_INV_MZ),
        },
        "logical_status": (
            "Diagnostic inverse images are algebraic inversions of the same "
            "conditional formula, not independent derivations or selectors."
        ),
    }


def check_N_R_formula(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Check that the implementation reproduces its H2-D4 input formula."""
    # Check if N increases monotonically with R
    N_values = [r["N_theory"] for r in results]
    R_values = [r["R"] for r in results]

    is_monotonic = all(N_values[i] <= N_values[i+1] for i in range(len(N_values)-1))

    # Re-fit the generated data as N = C [1/(1-R)]^p.
    # Log-log regression
    log_inv_R = np.log(1 / (1 - np.array(R_values)))
    log_N = np.log(np.array(N_values))

    # Fit log(N) = log(C) + p log(1/(1-R)).
    coeffs = np.polyfit(log_inv_R, log_N, 1)
    fitted_exponent = coeffs[0]
    log_coeff = coeffs[1]
    coeff = np.exp(log_coeff)

    return {
        "is_monotonic_increasing": is_monotonic,
        "fit_exponent": float(fitted_exponent),
        "fit_coefficient": float(coeff),
        "expected_exponent": 2.0,
        "expected_coefficient": 4 * np.pi,
        "exponent_error": float(abs(fitted_exponent - 2.0)),
        "coefficient_error": float(abs(coeff - 4*np.pi) / (4*np.pi)),
    }


def plot_implementation_check(
    results: List[Dict[str, Any]],
    formula_check: Dict[str, Any],
    output: Path,
) -> None:
    """Visualize a same-formula re-fit as an implementation check only."""
    reflection_rates = np.array([row["R"] for row in results], dtype=float)
    generated_values = np.array([row["N_theory"] for row in results], dtype=float)
    inverse_gap = 1.0 / (1.0 - reflection_rates)
    fitted_exponent = float(formula_check["fit_exponent"])
    fitted_coefficient = float(formula_check["fit_coefficient"])
    fitted_values = fitted_coefficient * inverse_gap**fitted_exponent

    figure, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))

    direct_axis = axes[0]
    direct_axis.plot(
        reflection_rates,
        generated_values,
        color="#2554c7",
        linewidth=2.2,
        label=r"input: $N(R)=4\pi/(1-R)^2$",
    )
    direct_axis.set_xlabel(r"reflection coefficient $R$")
    direct_axis.set_ylabel(r"conditional readout $N(R)$")
    direct_axis.set_title("Direct evaluation of the H2-D4 readout")
    direct_axis.grid(alpha=0.25)
    direct_axis.legend(loc="upper left")

    fit_axis = axes[1]
    fit_axis.loglog(
        inverse_gap,
        generated_values,
        "o",
        markersize=4.2,
        color="#5b63f2",
        alpha=0.75,
        label="data generated from the input formula",
    )
    fit_axis.loglog(
        inverse_gap,
        fitted_values,
        color="#d62728",
        linewidth=2.2,
        label=rf"same-family re-fit: $C[1/(1-R)]^p$",
    )
    fit_axis.set_xlabel(r"$1/(1-R)$")
    fit_axis.set_ylabel(r"conditional readout $N(R)$")
    fit_axis.set_title("Re-fit of data generated from the same formula")
    fit_axis.grid(alpha=0.25, which="both")
    fit_axis.legend(loc="upper left")
    fit_axis.text(
        0.04,
        0.07,
        rf"$p_{{\rm fit}}={fitted_exponent:.15f}$" "\n"
        rf"$C_{{\rm fit}}={fitted_coefficient:.15f}$" "\n"
        rf"$|p_{{\rm fit}}-2|={formula_check['exponent_error']:.2e}$" "\n"
        rf"$|C_{{\rm fit}}-4\pi|/(4\pi)={formula_check['coefficient_error']:.2e}$",
        transform=fit_axis.transAxes,
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92},
    )

    figure.suptitle(
        "Numerical implementation consistency check of the conditional readout",
        fontsize=15,
    )
    figure.text(
        0.5,
        0.012,
        "The re-fit uses data generated from the same H2-D4 formula; it is not an independent test of H2 or D4.",
        ha="center",
        fontsize=10.5,
    )
    figure.tight_layout(rect=(0.0, 0.055, 1.0, 0.94))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Phase 6: conditional N(R) readout")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "phase6_results",
        help="Output directory",
    )
    parser.add_argument("--r-min", type=float, default=0.686, help="Min R")
    parser.add_argument("--r-max", type=float, default=0.700, help="Max R")
    parser.add_argument("--r-points", type=int, default=141, help="Number of R points")
    parser.add_argument("--figure-dir", type=Path, default=HERE / "figures", help="Figure directory")

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Phase 6: conditional N(R) = 4π/(1-R)² readout under H2 and D4")
    print(f"R range: [{args.r_min:.10f}, {args.r_max:.10f}]")
    print("Empirical candidates: R1 and R2")
    print("External diagnostics: alpha^-1(0)=137.035999177, alpha^-1(M_Z^2)=128.946")
    print()

    # Sweep N(R)
    results = sweep_N_R(args.r_min, args.r_max, args.r_points)

    # Save results
    csv_path = args.output_dir / "phase6_N_R_sweep_v1.csv"
    with open(csv_path, "w", newline="") as f:
        fieldnames = ["R", "N_theory", "inv_alpha", "radius_mean", "radius_std"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"Saved N(R) sweep to {csv_path}")

    # Compare independent candidates with conditional diagnostic inverse images.
    identification = compare_candidates_and_diagnostics(results)

    id_json = args.output_dir / "phase6_137_128_identification_v1.json"
    with open(id_json, "w") as f:
        json.dump(identification, f, indent=2)
    print(f"Saved identification to {id_json}")

    # Check formula
    formula_check = check_N_R_formula(results)

    formula_json = args.output_dir / "phase6_formula_check_v1.json"
    with open(formula_json, "w") as f:
        json.dump(formula_check, f, indent=2)
    print(f"Saved formula check to {formula_json}")

    figure_path = args.figure_dir / "fig6_formula_validation.png"
    plot_implementation_check(results, formula_check, figure_path)
    print(f"Saved implementation-check figure to {figure_path}")

    # Print summary
    print("\n=== Conditional N(R) Readout Summary ===")
    print(f"\nImplementation Check: conditional N(R) = 4π/(1-R)²")
    print(f"  Fit exponent: {formula_check['fit_exponent']:.6f} (expected: 2.0)")
    print(f"  Fit coefficient: {formula_check['fit_coefficient']:.6f} (expected: {4*np.pi:.6f})")
    print(f"  Exponent error: {formula_check['exponent_error']:.2e}")
    print(f"  Coefficient error: {formula_check['coefficient_error']:.4%}")

    for candidate_name, values in identification["empirical_candidates"].items():
        print(f"\nEmpirical candidate {candidate_name}:")
        print(f"  R: {values['R']:.15f}")
        print(f"  Conditional N(R): {values['N']:.12f}")
        print(f"  External diagnostic: {values['diagnostic_N']:.12f}")
        print(f"  Delta N: {values['delta_N']:+.12f}")
        print(f"  Relative difference: {values['relative_difference']:.8%}")

    for diagnostic_name, values in identification["diagnostics"].items():
        print(f"\nDiagnostic {diagnostic_name}:")
        print(f"  N: {values['N']:.12f}")
        print(f"  Conditional inverse image R: {values['conditional_inverse_R']:.15f}")

    # Summary JSON
    summary = {
        "model": "phase6_N_R_derivation_v1",
        "theory": "conditional N(R) = 4π/(1-R)² under H2 and D4",
        "parameters": {
            "r_min": args.r_min,
            "r_max": args.r_max,
            "r_points": args.r_points,
        },
        "formula_check": formula_check,
        "identification": identification,
        "conclusion": (
            "R1 maps close to the low-energy diagnostic, whereas R2 maps to "
            "129.394062925 rather than 128.946. Diagnostic inverse images are "
            "not independent derivations."
        ),
    }

    summary_json = args.output_dir / "phase6_summary_v1.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved summary to {summary_json}")

    print("\nPhase 6 complete.")


if __name__ == "__main__":
    main()
