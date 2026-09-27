#!/usr/bin/env python3
"""Analytic leading-order Einstein-Maxwell quasi-circular inspiral benchmark.

Purpose
-------
Create a reproducible *reference branch* for later comparison with the anonymous
state-interaction generator.  The generator must never call this code or use its
readouts as feedback.

Approximation level
-------------------
- Geometric/natural units: G = M = c = 1, epsilon_0 = 1/(4*pi).
- Point particles, nonspinning, adiabatic quasi-circular inspiral.
- Newtonian conservative gravity + Coulomb interaction through Z = 1-lambda_A lambda_B.
- Leading gravitational quadrupole radiation.
- Leading electromagnetic electric-dipole radiation.
- The resulting radial decay ODE is solved in closed form; no time-step orbit
  integration is used.
- The 1PN charged Kepler law from Zhang et al. (arXiv:2604.10654, Eq. 79/80)
  is evaluated only as a diagnostic truncation estimate and never fed back into
  the leading-order evolution.

Output
------
- charged_binary_analytic_reference_v1_raw.csv
- charged_binary_analytic_reference_v1_metadata.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


TWO_PI = 2.0 * math.pi


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def primitive_time(r: np.ndarray, A: float, B: float) -> np.ndarray:
    """Primitive of r^3/(A r+B) with respect to r."""
    r = np.asarray(r, dtype=float)
    if abs(A) < 1e-30:
        return r**4 / (4.0 * B)
    return (
        r**3 / (3.0 * A)
        - B * r**2 / (2.0 * A**2)
        + B**2 * r / A**3
        - B**3 * np.log(A * r + B) / A**4
    )


def primitive_phase(r: np.ndarray, A: float, B: float, Z: float) -> np.ndarray:
    """Primitive of sqrt(Z)*r^(3/2)/(A r+B) with respect to r."""
    r = np.asarray(r, dtype=float)
    if abs(A) < 1e-30:
        return math.sqrt(Z) * (2.0 / (5.0 * B)) * r ** 2.5
    sr = np.sqrt(r)
    term = (
        2.0 * r ** 1.5 / (3.0 * A)
        - 2.0 * B * sr / A**2
        + 2.0 * B ** 1.5 / A ** 2.5 * np.arctan(np.sqrt(A * r / B))
    )
    return math.sqrt(Z) * term


def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(y)
    dx = np.diff(x)
    out[1:] = np.cumsum(0.5 * (y[:-1] + y[1:]) * dx)
    return out


def compute_reference(
    n_rows: int = 40001,
    r0: float = 50.0,
    rf: float = 20.0,
    mass_ratio: float = 1.0,
    lambda_A: float = 0.30,
    lambda_B: float = -0.30,
):
    if not (r0 > rf > 0):
        raise ValueError("Require r0 > rf > 0")
    if mass_ratio <= 0:
        raise ValueError("mass_ratio=mA/mB must be positive")

    # M = mA + mB = 1.
    mB = 1.0 / (1.0 + mass_ratio)
    mA = mass_ratio * mB
    M = 1.0
    mu = mA * mB / M
    eta = mu / M
    chiA = mA / M
    chiB = mB / M

    # Geometric charge: q_i = lambda_i m_i in G=1 units.
    qA = lambda_A * mA
    qB = lambda_B * mB
    Q = lambda_A * lambda_B  # qA qB/(mA mB) because G=1.
    Z = 1.0 - Q
    if Z <= 0:
        raise ValueError("Z=1-lambda_A*lambda_B must be >0 for an attractive circular benchmark")
    delta_lambda = lambda_A - lambda_B

    # Zhang et al. Eq. (73), G=M=c=1.
    zeta0 = qA**2 / (mA * M) + qB**2 / (mB * M)

    # Leading adiabatic radial-decay coefficients:
    # dr/dt = -A/r^2 - B/r^3 = -(A r+B)/r^3
    # A: leading EM electric-dipole radiation
    # B: leading GW mass-quadrupole radiation
    A = (4.0 / 3.0) * mu * delta_lambda**2 * Z
    B = (64.0 / 5.0) * mu * Z**2

    # Sample monotonically from r0 down to rf. Closed-form t(r), phi(r).
    r = np.linspace(r0, rf, n_rows, dtype=float)
    Ft0 = float(primitive_time(np.array([r0]), A, B)[0])
    Fp0 = float(primitive_phase(np.array([r0]), A, B, Z)[0])
    t = Ft0 - primitive_time(r, A, B)
    phi = Fp0 - primitive_phase(r, A, B, Z)

    omega = np.sqrt(Z / r**3)
    f_orb = omega / TWO_PI
    f_em = f_orb
    f_gw = omega / math.pi  # dominant m=2 quadrupole frequency
    v_rel = r * omega
    pn_epsilon = Z / r  # v^2 at LO
    x_pn = omega ** (2.0 / 3.0) / Z ** (4.0 / 3.0)

    # Center-of-mass trajectories at leading order.
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    x_rel = r * cphi
    y_rel = r * sphi
    xA = chiB * x_rel
    yA = chiB * y_rel
    xB = -chiA * x_rel
    yB = -chiA * y_rel

    # Orbital binding energy and radiation powers.
    E_orb = -0.5 * mu * Z / r
    P_gw = (32.0 / 5.0) * mu**2 * Z**3 / r**5
    P_em = (2.0 / 3.0) * mu**2 * delta_lambda**2 * Z**2 / r**4
    P_total = P_gw + P_em
    ratio_em_gw = P_em / P_gw
    dr_dt = -(A / r**2 + B / r**3)

    E_gw_cum = cumulative_trapezoid(P_gw, t)
    E_em_cum = cumulative_trapezoid(P_em, t)
    E_total_cum = E_gw_cum + E_em_cum
    E_expected = E_orb[0] - E_orb
    energy_balance_residual = E_total_cum - E_expected

    # 1PN charged Kepler diagnostic from Zhang et al. Eq. (79)/(80).
    # This is diagnostic only; no 1PN quantity feeds the LO evolution.
    F1 = (
        eta
        - 3.0
        - 0.5 * (4.0 * eta + 1.0) * Q**2
        + (4.5 + eta) * Q
        - zeta0
    )
    omega_1pn_sq = (Z + F1 / r) / r**3
    omega_1pn_diag = np.sqrt(np.maximum(omega_1pn_sq, 0.0))
    r_1pn_from_omega = (Z / omega**2) ** (1.0 / 3.0) * (1.0 + x_pn * F1 / (3.0 * Z))
    rel_r_1pn_correction = r_1pn_from_omega / r - 1.0
    rel_omega_1pn_correction = omega_1pn_diag / omega - 1.0

    cycles = phi / TWO_PI

    # Power crossover P_EM=P_GW (if delta_lambda != 0).
    if abs(delta_lambda) > 0:
        r_cross = 48.0 * Z / (5.0 * delta_lambda**2)
    else:
        r_cross = math.inf

    data = {
        "index": np.arange(n_rows, dtype=int),
        "r": r,
        "t": t,
        "phi": phi,
        "cycles": cycles,
        "x_rel": x_rel,
        "y_rel": y_rel,
        "x_A": xA,
        "y_A": yA,
        "x_B": xB,
        "y_B": yB,
        "omega_lo": omega,
        "omega_1pn_diag": omega_1pn_diag,
        "f_orb": f_orb,
        "f_em_dipole": f_em,
        "f_gw_quadrupole": f_gw,
        "v_rel": v_rel,
        "pn_epsilon": pn_epsilon,
        "x_pn": x_pn,
        "E_orb": E_orb,
        "P_gw": P_gw,
        "P_em": P_em,
        "P_total": P_total,
        "P_em_over_P_gw": ratio_em_gw,
        "dr_dt": dr_dt,
        "E_gw_cumulative": E_gw_cum,
        "E_em_cumulative": E_em_cum,
        "E_total_cumulative": E_total_cum,
        "E_balance_expected": E_expected,
        "E_balance_residual": energy_balance_residual,
        "r_1pn_from_omega_diag": r_1pn_from_omega,
        "rel_r_1pn_correction": rel_r_1pn_correction,
        "rel_omega_1pn_correction": rel_omega_1pn_correction,
    }

    meta = {
        "model": "leading-order Einstein-Maxwell adiabatic quasi-circular inspiral reference",
        "units": "G=M=c=1, epsilon0=1/(4*pi)",
        "source_primary": {
            "title": "Post-Newtonian dynamics of charged compact binaries",
            "authors": "Zi-Han Zhang, Tan Liu, Shuai Zhang, Zong-Kuan Guo",
            "arxiv": "2604.10654",
            "date": "2026-04-12",
            "used_for": [
                "Q=qA qB/(G mA mB), Z=1-Q definitions (Eq. 63)",
                "zeta0 definition (Eq. 73)",
                "charged 1PN Kepler law diagnostic (Eq. 79/80)",
                "GW quadrupole and EM dipole radiation hierarchy / flux-balance framework",
            ],
        },
        "initial_conditions": {
            "M": M,
            "mA": mA,
            "mB": mB,
            "mass_ratio_mA_over_mB": mass_ratio,
            "eta": eta,
            "mu": mu,
            "lambda_A_q_over_m": lambda_A,
            "lambda_B_q_over_m": lambda_B,
            "qA": qA,
            "qB": qB,
            "Q": Q,
            "Z": Z,
            "delta_lambda": delta_lambda,
            "zeta0": zeta0,
            "r0": r0,
            "rf": rf,
            "phi0": 0.0,
            "n_rows": n_rows,
        },
        "analytic_coefficients": {
            "A_em_in_dr_dt": A,
            "B_gw_in_dr_dt": B,
            "F1_1pn_kepler_diagnostic": F1,
            "r_power_crossover_Pem_eq_Pgw": r_cross,
        },
        "closed_form_relations": {
            "omega": "sqrt(Z/r^3)",
            "E_orb": "-mu*Z/(2*r)",
            "P_gw": "(32/5)*mu^2*Z^3/r^5",
            "P_em": "(2/3)*mu^2*(lambda_A-lambda_B)^2*Z^2/r^4",
            "dr_dt": "-(4/3)*mu*DeltaLambda^2*Z/r^2 - (64/5)*mu*Z^2/r^3",
            "t_of_r": "F_t(r0)-F_t(r), with F_t'=r^3/(A*r+B)",
            "phi_of_r": "F_phi(r0)-F_phi(r), with F_phi'=sqrt(Z)*r^(3/2)/(A*r+B)",
        },
        "scope_warning": (
            "This is a controlled leading-order analytic benchmark, not the full general Einstein-Maxwell two-body solution. "
            "The 1PN charged Kepler expression is diagnostic only and is not fed back into the LO inspiral."
        ),
    }
    return data, meta


def write_csv(data: dict[str, np.ndarray], path: Path) -> None:
    names = list(data.keys())
    n = len(data[names[0]])
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(names)
        for i in range(n):
            row = []
            for name in names:
                v = data[name][i]
                if name == "index":
                    row.append(int(v))
                else:
                    row.append(f"{float(v):.17g}")
            w.writerow(row)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument("--rows", type=int, default=40001)
    ap.add_argument("--r0", type=float, default=50.0)
    ap.add_argument("--rf", type=float, default=20.0)
    ap.add_argument("--mass-ratio", type=float, default=1.0)
    ap.add_argument("--lambda-A", type=float, default=0.30)
    ap.add_argument("--lambda-B", type=float, default=-0.30)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    data, meta = compute_reference(
        n_rows=args.rows,
        r0=args.r0,
        rf=args.rf,
        mass_ratio=args.mass_ratio,
        lambda_A=args.lambda_A,
        lambda_B=args.lambda_B,
    )

    csv_path = args.outdir / "charged_binary_analytic_reference_v1_raw.csv"
    meta_path = args.outdir / "charged_binary_analytic_reference_v1_metadata.json"
    write_csv(data, csv_path)

    # Numerical self-checks, recorded in metadata.
    max_energy_res = float(np.max(np.abs(data["E_balance_residual"])))
    final_energy_res = float(data["E_balance_residual"][-1])
    meta["results_summary"] = {
        "t_final": float(data["t"][-1]),
        "phi_final": float(data["phi"][-1]),
        "orbital_cycles": float(data["cycles"][-1]),
        "omega_initial": float(data["omega_lo"][0]),
        "omega_final": float(data["omega_lo"][-1]),
        "P_em_over_P_gw_initial": float(data["P_em_over_P_gw"][0]),
        "P_em_over_P_gw_final": float(data["P_em_over_P_gw"][-1]),
        "E_gw_emitted": float(data["E_gw_cumulative"][-1]),
        "E_em_emitted": float(data["E_em_cumulative"][-1]),
        "E_total_emitted": float(data["E_total_cumulative"][-1]),
        "orbital_energy_loss_expected": float(data["E_balance_expected"][-1]),
        "max_abs_energy_balance_residual": max_energy_res,
        "final_energy_balance_residual": final_energy_res,
        "max_pn_epsilon": float(np.max(data["pn_epsilon"])),
        "max_abs_rel_r_1pn_diagnostic": float(np.max(np.abs(data["rel_r_1pn_correction"]))),
        "max_abs_rel_omega_1pn_diagnostic": float(np.max(np.abs(data["rel_omega_1pn_correction"]))),
    }
    meta["sha256"] = {"raw_csv": sha256_file(csv_path)}
    with meta_path.open("w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print(json.dumps(meta["results_summary"], indent=2))
    print(f"raw_csv={csv_path}")
    print(f"metadata={meta_path}")


if __name__ == "__main__":
    main()
