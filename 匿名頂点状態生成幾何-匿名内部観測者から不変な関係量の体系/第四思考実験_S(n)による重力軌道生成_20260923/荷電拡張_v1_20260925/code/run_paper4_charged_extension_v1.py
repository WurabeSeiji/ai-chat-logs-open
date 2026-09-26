#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 4 charged extension v1
============================

Purpose
-------
Minimal charged extension of the Paper-4 state-dependent orbit experiment.
Keep the same macro state/readout scaffold
    (chi, p, e, phi, t) -> X=(a,b) -> Y=(x,y)
and add only
    * discrete charge integers n_A, n_B,
    * one continuous coupling-scale Gamma = k e0^2/(G M^2),
    * Coulomb contribution to the Newtonian central coupling,
    * leading EM dipole radiation + leading mass-quadrupole GW radiation,
    * waveform readouts and simple anonymous invariant readouts.

This is a diagnostic scaffold, NOT a self-consistent Einstein-Maxwell solution.
The conservative angular gauge can be either:
    newtonian : dphi/dchi = 1
    paper4     : retain Paper-4 Schwarzschild/Darwin angular factor as a weak
                 gravity correction.  In a Coulomb-dominated orbit this is a
                 hybrid construction and is explicitly flagged as such.

Units
-----
G=c=M=1, M=m_A+m_B.
Mass fractions x_A=m_A/M, x_B=m_B/M, nu=x_A x_B.
Charges q_i=n_i e0, n_i integer.
Gamma = k e0^2/(G M^2).
The relative Newtonian central acceleration is
    rddot = -alpha r/r^3,
    alpha = 1 - Gamma*n_A*n_B/nu.
For opposite charges, alpha>1 and Coulomb is attractive.
The force ratio is
    |F_C/F_G| = Gamma*|n_A*n_B|/nu.

To make a Coulomb-dominated test remain nonrelativistic while preserving the
Paper-4 C1 orbital speed scale, use --scale-orbit-with-alpha.  Then p,a,b are
multiplied by alpha, so v ~ sqrt(alpha/p) remains comparable to the C1 value.

Radiation-reaction model
------------------------
Both loss channels are obtained from leading multipole fluxes evaluated on the
instantaneous Newtonian osculating orbit with central coupling alpha:

GW:
  P_GW = (nu^2/5) sum_ij (Q'''_ij)^2  [total, normalized]
  Jdot_GW,z = (2 nu^2/5) eps_zjk Q''_jl Q'''_kl
  Specific orbital E,h losses divide by nu.

EM dipole:
  beta_q = n_A*x_B - n_B*x_A
  P_EM = (2/3) Gamma beta_q^2 |rddot|^2 [total, normalized]
  Jdot_EM,z = (2/3) Gamma beta_q^2 (v x rddot)_z.
  Specific orbital E,h losses divide by nu.

Then, with p=h^2/alpha and Eorb=-alpha(1-e^2)/(2p),
  pdot = 2p hdot/h,
  edot = (p Edot + Eorb pdot)/(alpha e).
These are integrated versus chi by RK4, exactly as in Paper 4's osculating
p,e shell.  Radiation reaction can be toggled per channel.

Waveform readouts
-----------------
For a face-on observer at normalized distance D:
  GW basis from Q''; waveform coefficient is nu/D.
  EM basis from relative acceleration; waveform coefficient is
      sqrt(Gamma)*beta_q/D.
The common distance cancels in their coefficient ratio.

Outputs
-------
CSV time series + JSON summary in --out.
"""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

BASE_A = 30.0
BASE_B = 27.0


def coeffs(N: int) -> np.ndarray:
    c = [1.0]
    for m in range(N):
        c.append(c[-1] * (2 * m + 1) / (m + 1))
    return np.array(c, dtype=float)


def gN_scalar(chi: float, p: float, e: float, c: np.ndarray) -> float:
    u = (3.0 + e * math.cos(chi)) / p
    y = float(c[-1])
    for a in c[-2::-1]:
        y = float(a) + u * y
    return y


def rot(th: float) -> np.ndarray:
    c, s = math.cos(th), math.sin(th)
    return np.array([[c, -s], [s, c]], dtype=float)


def Phi(X: np.ndarray) -> np.ndarray:
    a, b = X[..., 0], X[..., 1]
    return np.column_stack((a * a - b * b, 2 * a * b))


def mass_fractions(q_mass: float):
    if q_mass <= 0:
        raise ValueError("mass ratio must be positive")
    xA = q_mass / (1.0 + q_mass)
    xB = 1.0 / (1.0 + q_mass)
    return xA, xB, xA * xB


def charge_couplings(q_mass: float, nA: int, nB: int, Gamma: float):
    xA, xB, nu = mass_fractions(q_mass)
    rho = Gamma * abs(nA * nB) / nu if nA * nB != 0 else 0.0
    alpha = 1.0 - Gamma * nA * nB / nu
    beta_q = nA * xB - nB * xA
    return {
        "xA": xA,
        "xB": xB,
        "nu": nu,
        "rho_coulomb_to_gravity": rho,
        "alpha": alpha,
        "beta_q": beta_q,
        "gamma_nprod": Gamma * nA * nB,
        "gamma_beta2": Gamma * beta_q * beta_q,
    }


def orbital_vectors(chi: float, p: float, e: float, phi: float, alpha: float, gphi: float):
    """Instantaneous osculating Newtonian kinematics, ignoring slow p,e derivatives.

    r = p/(1+e cos chi), h=sqrt(alpha p), chi_dot=h/r^2.
    radial speed is the generalized-Kepler value e sin(chi)*sqrt(alpha/p).
    tangential speed uses phi_dot=gphi*chi_dot, so v_t=r*phi_dot.
    """
    ce, se = math.cos(chi), math.sin(chi)
    den = 1.0 + e * ce
    r = p / den
    h = math.sqrt(alpha * p)
    chi_dot = h / (r * r)
    rdot = e * se * math.sqrt(alpha / p)
    phi_dot = gphi * chi_dot
    vt = r * phi_dot
    cp, sp = math.cos(phi), math.sin(phi)
    er = np.array([cp, sp], dtype=float)
    et = np.array([-sp, cp], dtype=float)
    rvec = r * er
    vvec = rdot * er + vt * et
    # Leading central acceleration.  The small Paper-4 angular gauge correction
    # is a readout/conservative-phase diagnostic, not included in the force here.
    avec = -alpha * rvec / (r ** 3)
    rv = float(np.dot(rvec, vvec))
    jvec = -alpha * (vvec / r ** 3 - 3.0 * rv * rvec / r ** 5)
    return {
        "r": r,
        "h": h,
        "chi_dot": chi_dot,
        "rdot": rdot,
        "phi_dot": phi_dot,
        "vt": vt,
        "rvec": rvec,
        "vvec": vvec,
        "avec": avec,
        "jvec": jvec,
        "speed2": float(np.dot(vvec, vvec)),
    }


def quadrupole_derivatives(rvec2, vvec2, avec2, jvec2):
    r = np.array([rvec2[0], rvec2[1], 0.0])
    v = np.array([vvec2[0], vvec2[1], 0.0])
    a = np.array([avec2[0], avec2[1], 0.0])
    j = np.array([jvec2[0], jvec2[1], 0.0])
    I = np.eye(3)
    qdd = (np.outer(a, r) + np.outer(r, a) + 2.0 * np.outer(v, v)
           - (2.0 / 3.0) * I * (np.dot(v, v) + np.dot(r, a)))
    qddd = (np.outer(j, r) + np.outer(r, j)
            + 3.0 * (np.outer(a, v) + np.outer(v, a))
            - (2.0 / 3.0) * I * (3.0 * np.dot(v, a) + np.dot(r, j)))
    return qdd, qddd


def radiation_fluxes(chi, p, e, phi, alpha, gphi, nu, Gamma, beta_q,
                     enable_gw=True, enable_em=True):
    """Return specific orbital Edot, hdot and detailed channel diagnostics."""
    ov = orbital_vectors(chi, p, e, phi, alpha, gphi)
    qdd, qddd = quadrupole_derivatives(ov["rvec"], ov["vvec"], ov["avec"], ov["jvec"])

    # Gravitational-wave channel: total normalized flux and specific losses.
    gw_power_total = (nu * nu / 5.0) * float(np.sum(qddd * qddd)) if enable_gw else 0.0
    # epsilon_zjk qdd_jl qddd_kl = qdd_xl*qddd_yl - qdd_yl*qddd_xl
    gw_jflux_total = ((2.0 * nu * nu / 5.0)
                      * float(np.sum(qdd[0, :] * qddd[1, :] - qdd[1, :] * qddd[0, :]))) if enable_gw else 0.0
    Edot_gw = -gw_power_total / nu if enable_gw else 0.0
    hdot_gw = -gw_jflux_total / nu if enable_gw else 0.0

    # Electric-dipole channel. beta_q is exchange-odd, beta_q^2 is anonymous.
    a2 = float(np.dot(ov["avec"], ov["avec"]))
    vxa = float(ov["vvec"][0] * ov["avec"][1] - ov["vvec"][1] * ov["avec"][0])
    em_power_total = ((2.0 / 3.0) * Gamma * beta_q * beta_q * a2) if enable_em else 0.0
    em_jflux_total = ((2.0 / 3.0) * Gamma * beta_q * beta_q * vxa) if enable_em else 0.0
    Edot_em = -em_power_total / nu if enable_em else 0.0
    hdot_em = -em_jflux_total / nu if enable_em else 0.0

    return {
        "Edot_gw": Edot_gw,
        "hdot_gw": hdot_gw,
        "Edot_em": Edot_em,
        "hdot_em": hdot_em,
        "gw_power_total": gw_power_total,
        "em_power_total": em_power_total,
        "gw_jflux_total": gw_jflux_total,
        "em_jflux_total": em_jflux_total,
        "qdd": qdd,
        "qddd": qddd,
        **ov,
    }


def rates_chi(chi, p, e, phi, alpha, nu, Gamma, beta_q, c, angular_mode,
              enable_gw_rr, enable_em_rr):
    if angular_mode == "newtonian":
        gphi = 1.0
    elif angular_mode == "paper4":
        gphi = gN_scalar(chi, p, e, c)
    else:
        raise ValueError(angular_mode)

    flux = radiation_fluxes(chi, p, e, phi, alpha, gphi, nu, Gamma, beta_q,
                            enable_gw=enable_gw_rr, enable_em=enable_em_rr)
    r, h = flux["r"], flux["h"]
    dt_dchi = r * r / h
    Edot = flux["Edot_gw"] + flux["Edot_em"]
    hdot = flux["hdot_gw"] + flux["hdot_em"]
    pdot = 2.0 * p * hdot / h
    Eorb = -alpha * (1.0 - e * e) / (2.0 * p)
    if e <= 1e-12:
        # Circular limit requires a separate regularized eccentricity variable.
        edot = 0.0
    else:
        edot = (p * Edot + Eorb * pdot) / (alpha * e)
    return pdot * dt_dchi, edot * dt_dchi, dt_dchi, gphi, flux


def rk4_macro(chi, p, e, phi, hchi, alpha, nu, Gamma, beta_q, c, angular_mode,
              enable_gw_rr, enable_em_rr):
    """RK4 on (p,e,t,phi) versus chi; phi derivative is gphi."""
    def f(x, y):
        pp, ee, tt, ph = y
        dp, de, dt, gp, _ = rates_chi(x, pp, ee, ph, alpha, nu, Gamma, beta_q, c,
                                      angular_mode, enable_gw_rr, enable_em_rr)
        return np.array([dp, de, dt, gp], dtype=float)

    y = np.array([p, e, 0.0, phi], dtype=float)
    k1 = f(chi, y)
    k2 = f(chi + 0.5 * hchi, y + 0.5 * hchi * k1)
    k3 = f(chi + 0.5 * hchi, y + 0.5 * hchi * k2)
    k4 = f(chi + hchi, y + hchi * k3)
    dy = (hchi / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return p + dy[0], e + dy[1], dy[2], phi + dy[3]


def fit_common_coefficient(data_channels, basis_channels):
    num = 0.0
    den = 0.0
    for d, b in zip(data_channels, basis_channels):
        d = np.asarray(d, dtype=float)
        b = np.asarray(b, dtype=float)
        num += float(np.dot(d, b))
        den += float(np.dot(b, b))
    return num / den if den > 0 else float("nan")


def charge_candidates(nu_fit, gamma_nprod_fit, gamma_beta2_fit, Gamma, nmax=4):
    if not (0 < nu_fit <= 0.25 + 1e-8) or Gamma <= 0:
        return []
    s = math.sqrt(max(0.0, 1.0 - 4.0 * nu_fit))
    xA, xB = 0.5 * (1.0 + s), 0.5 * (1.0 - s)
    target_prod = gamma_nprod_fit / Gamma
    target_beta2 = gamma_beta2_fit / Gamma
    rows = []
    for nA in range(-nmax, nmax + 1):
        for nB in range(-nmax, nmax + 1):
            if nA == 0 and nB == 0:
                continue
            prod = nA * nB
            beta = nA * xB - nB * xA
            err = abs(prod - target_prod) + abs(beta * beta - target_beta2)
            rows.append((err, nA, nB, prod, beta * beta))
    rows.sort(key=lambda z: z[0])
    return [
        {"nA": int(a), "nB": int(b), "nprod": int(prod), "beta2": float(be2), "score": float(err)}
        for err, a, b, prod, be2 in rows[:12]
    ]


def run(args):
    xA, xB, nu = mass_fractions(args.mass_ratio)
    if args.gamma is None:
        if args.coulomb_to_gravity < 0:
            raise ValueError("coulomb-to-gravity must be nonnegative")
        prod_abs = abs(args.nA * args.nB)
        if prod_abs == 0:
            Gamma = 0.0
        else:
            Gamma = args.coulomb_to_gravity * nu / prod_abs
    else:
        Gamma = args.gamma

    coupl = charge_couplings(args.mass_ratio, args.nA, args.nB, Gamma)
    alpha = coupl["alpha"]
    if alpha <= 0:
        raise RuntimeError(
            f"Net central interaction is not attractive: alpha={alpha}. "
            "A bound ellipse is not available in this v1 model."
        )

    e0 = math.sqrt(1.0 - (BASE_B / BASE_A) ** 2)
    p_base = BASE_A * (1.0 - e0 * e0)
    scale = (alpha if args.scale_orbit_with_alpha else 1.0) * args.orbit_dilation
    p0 = p_base * scale
    a0 = BASE_A * scale
    b0 = BASE_B * scale

    nsteps = int(round(args.steps_per_radial_orbit * args.radial_orbits))
    dchi = 2.0 * math.pi / args.steps_per_radial_orbit
    c = coeffs(args.order)

    chi = np.arange(nsteps + 1, dtype=float) * dchi
    p = np.empty(nsteps + 1); e = np.empty(nsteps + 1)
    phi = np.empty(nsteps + 1); t = np.empty(nsteps + 1)
    p[0], e[0], phi[0], t[0] = p0, e0, 0.0, 0.0

    for k in range(nsteps):
        pn, en, dt, phin = rk4_macro(
            chi[k], p[k], e[k], phi[k], dchi, alpha, nu, Gamma, coupl["beta_q"], c,
            args.angular_mode, not args.disable_gw_rr, not args.disable_em_rr
        )
        if not (pn > 0 and 0 < en < 1):
            raise RuntimeError(f"Invalid osculating elements at step {k}: p={pn}, e={en}")
        p[k + 1], e[k + 1], phi[k + 1], t[k + 1] = pn, en, phin, t[k] + dt

    r = p / (1.0 + e * np.cos(chi))
    theta = 0.5 * phi
    X = np.zeros((nsteps + 1, 2), dtype=float)
    X[0] = [math.sqrt(r[0]), 0.0]
    det_err = []
    map_err = []
    for k in range(nsteps):
        rho = math.sqrt(r[k + 1] / r[k])
        S = rot(theta[k + 1]) @ np.diag([rho, 1.0 / rho]) @ rot(-theta[k])
        X[k + 1] = S @ X[k]
        det_err.append(abs(np.linalg.det(S) - 1.0))
        target = math.sqrt(r[k + 1]) * np.array([math.cos(theta[k + 1]), math.sin(theta[k + 1])])
        map_err.append(float(np.linalg.norm(X[k + 1] - target)))
    Y = Phi(X)

    # Readout arrays, generated independently from the charge/mass coefficients.
    gw_plus = np.empty(nsteps + 1); gw_cross = np.empty(nsteps + 1)
    em_x = np.empty(nsteps + 1); em_y = np.empty(nsteps + 1)
    gw_bp = np.empty(nsteps + 1); gw_bc = np.empty(nsteps + 1)
    em_bx = np.empty(nsteps + 1); em_by = np.empty(nsteps + 1)
    speed = np.empty(nsteps + 1)
    pgw = np.empty(nsteps + 1); pem = np.empty(nsteps + 1)
    for k in range(nsteps + 1):
        gphi = 1.0 if args.angular_mode == "newtonian" else gN_scalar(chi[k], p[k], e[k], c)
        fl = radiation_fluxes(chi[k], p[k], e[k], phi[k], alpha, gphi, nu, Gamma, coupl["beta_q"],
                              enable_gw=True, enable_em=True)
        qdd = fl["qdd"]
        gw_bp[k] = qdd[0, 0] - qdd[1, 1]
        gw_bc[k] = 2.0 * qdd[0, 1]
        gw_plus[k] = (nu / args.distance) * gw_bp[k]
        gw_cross[k] = (nu / args.distance) * gw_bc[k]
        em_bx[k], em_by[k] = fl["avec"][0], fl["avec"][1]
        em_coeff = math.sqrt(max(Gamma, 0.0)) * coupl["beta_q"] / args.distance
        em_x[k], em_y[k] = em_coeff * em_bx[k], em_coeff * em_by[k]
        speed[k] = math.sqrt(fl["speed2"])
        pgw[k], pem[k] = fl["gw_power_total"], fl["em_power_total"]

    # First-orbit readout, avoiding late secular evolution bias.
    n1 = min(args.steps_per_radial_orbit + 1, nsteps + 1)
    gw_coeff_fit = fit_common_coefficient(
        [gw_plus[:n1], gw_cross[:n1]], [gw_bp[:n1], gw_bc[:n1]]
    )
    em_coeff_fit = fit_common_coefficient(
        [em_x[:n1], em_y[:n1]], [em_bx[:n1], em_by[:n1]]
    )
    nu_fit = abs(gw_coeff_fit) * args.distance
    gamma_beta2_fit = (abs(em_coeff_fit) * args.distance) ** 2
    if nsteps >= args.steps_per_radial_orbit:
        T1 = t[args.steps_per_radial_orbit] - t[0]
        alpha_fit = (2.0 * math.pi * a0 ** 1.5 / T1) ** 2
    else:
        T1 = float("nan")
        alpha_fit = float("nan")
    gamma_nprod_fit = nu_fit * (1.0 - alpha_fit) if math.isfinite(alpha_fit) else float("nan")

    # Local clock readout.  Unlike the period fit, this is not biased by secular
    # radiation-reaction evolution across a full orbit.  It uses the Paper-4
    # state/readout relation dt/dchi = r^2/sqrt(alpha p).
    dt_dchi_num = np.gradient(t[:n1], chi[:n1])
    alpha_local_series = r[:n1] ** 4 / (p[:n1] * dt_dchi_num ** 2)
    edge = max(2, min(20, n1 // 20))
    core = alpha_local_series[edge:-edge] if n1 > 2 * edge else alpha_local_series
    alpha_local_fit = float(np.median(core))
    gamma_nprod_local_fit = nu_fit * (1.0 - alpha_local_fit)
    candidates = charge_candidates(nu_fit, gamma_nprod_local_fit, gamma_beta2_fit, Gamma, args.candidate_nmax)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = f"charged_v1_qm{args.mass_ratio:g}_n{args.nA}_{args.nB}_rho{coupl['rho_coulomb_to_gravity']:.6g}"
    csv_path = out / f"{stem}.csv"
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "k","chi","t","p","e","phi","x","y","r","X0","X1","speed",
            "gw_plus","gw_cross","em_x","em_y","gw_power","em_power"
        ])
        for k in range(nsteps + 1):
            w.writerow([
                k, f"{chi[k]:.15g}", f"{t[k]:.15g}", f"{p[k]:.15g}", f"{e[k]:.15g}",
                f"{phi[k]:.15g}", f"{Y[k,0]:.15g}", f"{Y[k,1]:.15g}", f"{r[k]:.15g}",
                f"{X[k,0]:.15g}", f"{X[k,1]:.15g}", f"{speed[k]:.15g}",
                f"{gw_plus[k]:.15g}", f"{gw_cross[k]:.15g}", f"{em_x[k]:.15g}", f"{em_y[k]:.15g}",
                f"{pgw[k]:.15g}", f"{pem[k]:.15g}"
            ])

    # Basic diagnostic quantities.
    one_orbit = min(args.steps_per_radial_orbit, nsteps)
    avg_pgw = float(np.trapezoid(pgw[:one_orbit + 1], t[:one_orbit + 1]) / max(t[one_orbit] - t[0], 1e-300))
    avg_pem = float(np.trapezoid(pem[:one_orbit + 1], t[:one_orbit + 1]) / max(t[one_orbit] - t[0], 1e-300))
    power_ratio = avg_pem / avg_pgw if avg_pgw > 0 else float("inf")

    summary = {
        "model": "Paper-4 charged extension v1 diagnostic scaffold",
        "warnings": [
            "Not a self-consistent Einstein-Maxwell solution.",
            "Paper4 angular mode in a Coulomb-dominated orbit is a hybrid gauge, not a derived charged-GR orbit.",
            "Radiation reaction uses leading multipole fluxes on a Newtonian osculating orbit.",
            "Waveform amplitudes are normalized source readouts; absolute SI calibration requires a physical M and D.",
        ],
        "inputs": {
            "mass_ratio_mA_over_mB": args.mass_ratio,
            "nA": args.nA, "nB": args.nB,
            "Gamma": Gamma,
            "requested_coulomb_to_gravity": args.coulomb_to_gravity if args.gamma is None else None,
            "angular_mode": args.angular_mode,
            "scale_orbit_with_alpha": args.scale_orbit_with_alpha,
            "orbit_dilation": args.orbit_dilation,
            "gw_rr_enabled": not args.disable_gw_rr,
            "em_rr_enabled": not args.disable_em_rr,
            "distance_normalized": args.distance,
        },
        "derived_couplings": coupl,
        "initial_orbit": {
            "base_a": BASE_A, "base_b": BASE_B, "base_p": p_base, "e0": e0,
            "scale": scale, "a0": a0, "b0": b0, "p0": p0,
            "periapsis": p0 / (1.0 + e0),
            "apoapsis": p0 / (1.0 - e0),
        },
        "final_state": {
            "p": float(p[-1]), "e": float(e[-1]), "phi": float(phi[-1]), "t": float(t[-1]),
        },
        "diagnostics": {
            "max_speed_over_c": float(np.max(speed)),
            "nonrelativistic_warning": bool(np.max(speed) >= 0.3),
            "superluminal_invalid": bool(np.max(speed) >= 1.0),
            "first_orbit_period": float(T1),
            "first_orbit_avg_gw_power": avg_pgw,
            "first_orbit_avg_em_power": avg_pem,
            "first_orbit_em_to_gw_power_ratio": float(power_ratio),
            "max_det_S_minus_1": float(max(det_err, default=0.0)),
            "max_state_map_error": float(max(map_err, default=0.0)),
        },
        "readout": {
            "gw_coefficient_fit": float(gw_coeff_fit),
            "expected_gw_coefficient_nu_over_D": float(nu / args.distance),
            "nu_fit": float(nu_fit),
            "em_coefficient_fit": float(em_coeff_fit),
            "expected_em_coefficient_sqrtGamma_beta_over_D": float(math.sqrt(max(Gamma, 0.0)) * coupl["beta_q"] / args.distance),
            "gamma_beta2_fit": float(gamma_beta2_fit),
            "gamma_beta2_true": float(coupl["gamma_beta2"]),
            "alpha_period_fit": float(alpha_fit),
            "alpha_local_clock_fit": float(alpha_local_fit),
            "alpha_true": float(alpha),
            "gamma_nprod_period_fit": float(gamma_nprod_fit),
            "gamma_nprod_local_fit": float(gamma_nprod_local_fit),
            "gamma_nprod_true": float(coupl["gamma_nprod"]),
            "charge_candidates_using_known_Gamma": candidates,
        },
        "output_csv": csv_path.name,
    }
    json_path = out / f"{stem}_summary.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
        f.write("\n")
    return summary, csv_path, json_path


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass-ratio", type=float, default=1.0)
    ap.add_argument("--nA", type=int, default=1)
    ap.add_argument("--nB", type=int, default=-1)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--gamma", type=float, default=None,
                   help="Gamma = k e0^2/(G M^2)")
    g.add_argument("--coulomb-to-gravity", type=float, default=100.0,
                   help="target |F_C/F_G|; Gamma is derived from integer charges and mass ratio")
    ap.add_argument("--scale-orbit-with-alpha", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--angular-mode", choices=["newtonian", "paper4"], default="paper4")
    ap.add_argument("--orbit-dilation", type=float, default=1.0,
                    help="extra multiplicative dilation after optional alpha scaling; >1 makes RR more adiabatic and speeds smaller")
    ap.add_argument("--disable-gw-rr", action="store_true")
    ap.add_argument("--disable-em-rr", action="store_true")
    ap.add_argument("--order", type=int, default=12)
    ap.add_argument("--steps-per-radial-orbit", type=int, default=4000)
    ap.add_argument("--radial-orbits", type=float, default=3.0)
    ap.add_argument("--distance", type=float, default=1.0)
    ap.add_argument("--candidate-nmax", type=int, default=4)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "results"))
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    summary, csv_path, json_path = run(args)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(csv_path)
    print(json_path)
