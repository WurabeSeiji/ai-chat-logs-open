#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 4 Experiment 03: nine-case transferability test of one S(n) rule.

Purpose
-------
Keep the S(n) recurrence algorithm fixed and change only the initial orbital
parameters. Compare against the same Paper-3 PN reference equations for a
3 x 3 grid of semilatus rectum p and eccentricity e at fixed mass ratio q=1.

Grid
----
    p in {24.3, 60, 120}
    e in {0.20, 0.45, 0.55}
    q = 1, nu = 0.25

S(n) model
----------
Conservative angular gauge, order N:
    dphi/dchi = sum_m c_m ((3 + e cos chi)/p)^m
    c_0 = 1, c_{m+1} = ((2m+1)/(m+1)) c_m

Dissipative recurrence:
    local 2.5PN radiation reaction from the Paper-3 implementation,
    converted to osculating p,e drift with Newtonian relations and integrated
    by RK4 in chi.

Internal map:
    X_{k+1} = S_k X_k
    S_k = R(theta_{k+1}) diag(rho, rho^-1) R(-theta_k)
    theta = phi/2, rho = sqrt(r_{k+1}/r_k)
    det(S_k) = 1 analytically.

Outputs
-------
raw/pn_integrator_<case>.csv  : raw PN RK4 rows before plotting/interpolation
raw/pn_plot_<case>.csv        : uniform-time PN rows used in orbit plots
raw/sn_<case>.csv             : every S(n) step used in orbit plots, including
                                p,e and local RR diagnostic rows
summary/case_metrics.csv      : one row per case
summary/case_metrics.json     : same metrics + metadata
summary/cases.json            : exact case definitions
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import paper3_pn_core as p3

P_VALUES = (24.3, 60.0, 120.0)
E_VALUES = (0.20, 0.45, 0.55)
Q_MASS = 1.0
DEFAULT_ORDER = 12
DEFAULT_STEPS = 4000
DEFAULT_RADIAL_ORBITS = 3.0


def write_csv(path: Path, columns: dict[str, np.ndarray | list]):
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    arrays = [np.asarray(columns[k]) for k in names]
    n = len(arrays[0])
    if any(len(a) != n for a in arrays):
        raise ValueError("CSV columns have different lengths")
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(names)
        for i in range(n):
            row = []
            for a in arrays:
                v = a[i]
                if isinstance(v, (np.integer, int)):
                    row.append(str(int(v)))
                else:
                    row.append(f"{float(v):.15g}")
            w.writerow(row)


def coeffs(N: int) -> np.ndarray:
    c = [1.0]
    for m in range(N):
        c.append(c[-1] * (2 * m + 1) / (m + 1))
    return np.asarray(c, dtype=float)


def gN_scalar(chi: float, p: float, e: float, c: np.ndarray) -> float:
    u = (3.0 + e * math.cos(chi)) / p
    y = float(c[-1])
    for a in c[-2::-1]:
        y = float(a) + u * y
    return y


def rot(th: float) -> np.ndarray:
    cc, ss = math.cos(th), math.sin(th)
    return np.array([[cc, -ss], [ss, cc]], dtype=float)


def Phi(X: np.ndarray) -> np.ndarray:
    a, b = X[..., 0], X[..., 1]
    return np.column_stack((a * a - b * b, 2.0 * a * b))


def rr_rates_chi(chi: float, p: float, e: float, nu: float):
    ce, se = math.cos(chi), math.sin(chi)
    den = 1.0 + e * ce
    r = p / den
    h = math.sqrt(p)
    rdot = e * se / h
    vt = den / h
    v2 = rdot * rdot + vt * vt
    invr = 1.0 / r

    # Same local 2.5PN RR acceleration term as Paper 3.
    C = 1.6 * nu * invr ** 3
    an = rdot * (18.0 * v2 + (2.0 / 3.0) * invr - 25.0 * rdot * rdot)
    bv = -(6.0 * v2 - 2.0 * invr - 15.0 * rdot * rdot)

    Edot = C * (an * rdot + bv * v2)
    pdot = 2.0 * C * bv * p
    E = (e * e - 1.0) / (2.0 * p)
    edot = (p * Edot + E * pdot) / e
    dt_dchi = r * r / h

    return pdot * dt_dchi, edot * dt_dchi, dt_dchi, {
        "r_local": r,
        "rdot": rdot,
        "vt": vt,
        "v2": v2,
        "Edot": Edot,
        "pdot": pdot,
        "edot": edot,
        "C_rr": C,
        "an_rr": an,
        "bv_rr": bv,
    }


def rk4_elements(chi: float, p: float, e: float, hchi: float, nu: float):
    def f(x, pp, ee):
        dp, de, dt, _ = rr_rates_chi(x, pp, ee, nu)
        return np.array([dp, de, dt], dtype=float)

    y = np.array([p, e, 0.0], dtype=float)
    k1 = f(chi, y[0], y[1])
    y2 = y + 0.5 * hchi * k1
    k2 = f(chi + 0.5 * hchi, y2[0], y2[1])
    y3 = y + 0.5 * hchi * k2
    k3 = f(chi + 0.5 * hchi, y3[0], y3[1])
    y4 = y + hchi * k3
    k4 = f(chi + hchi, y4[0], y4[1])
    dy = (hchi / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return p + dy[0], e + dy[1], dy[2]


def periapsis_values(phi: np.ndarray, r: np.ndarray):
    ii = np.where((r[1:-1] < r[:-2]) & (r[1:-1] <= r[2:]))[0] + 1
    idx = np.r_[0, ii]
    return idx, r[idx], phi[idx]


def precession_from_rows(phi: np.ndarray, r: np.ndarray) -> float:
    _, _, pp = periapsis_values(phi, r)
    if len(pp) < 2:
        return float("nan")
    return float(np.mean(np.diff(pp) - 2.0 * math.pi))


def radial_rms_by_phi(pn_phi, pn_r, sn_phi, sn_r, a):
    hi = min(float(pn_phi[-1]), float(sn_phi[-1]))
    q = np.linspace(0.0, hi, 16000)
    rr = np.interp(q, pn_phi, pn_r)
    rs = np.interp(q, sn_phi, sn_r)
    rms = float(np.sqrt(np.mean((rr - rs) ** 2)))
    return rms, rms / a


def case_id(p: float, e: float) -> str:
    return f"P{str(p).replace('.', 'p')}_E{int(round(e * 100)):02d}"


def run_case(case, base: str, order: int, steps: int, radial_orbits: float):
    p0 = float(case["p0"])
    e0 = float(case["e0"])
    cid = case["case_id"]
    basep = Path(base)
    raw = basep / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    qmass = Q_MASS
    ma = qmass / (1.0 + qmass)
    mb = 1.0 / (1.0 + qmass)
    nu = ma * mb

    # Same ellipse definitions as Paper 3: p=a(1-e^2), b=a sqrt(1-e^2).
    a0 = p0 / (1.0 - e0 * e0)
    b0 = a0 * math.sqrt(1.0 - e0 * e0)
    T0 = 2.0 * math.pi * a0 ** 1.5

    # ---------------------------------------------------------------- PN reference
    init, _, init_info = p3.initial_state(a0, b0, nu, ma, mb)
    ts, st, dv, stop = p3.integrate(a0, b0, nu, ma, mb,
                                     radial_orbits * T0, rr=True, init=init)
    # Preserve raw RK4 rows before any interpolation/plotting.
    write_csv(raw / f"pn_integrator_{cid}.csv", {
        "step": np.arange(len(ts), dtype=int),
        "t": ts,
        "x": st[:, 0], "y": st[:, 1], "vx": st[:, 2], "vy": st[:, 3],
        "tau_A": st[:, 4], "tau_B": st[:, 5],
        "dxdt": dv[:, 0], "dydt": dv[:, 1], "ax": dv[:, 2], "ay": dv[:, 3],
    })

    tq = np.linspace(0.0, ts[-1], int(round(steps * radial_orbits)) + 1)
    yy = p3.interpolate(ts, st, dv, tq)
    xpn, ypn, vxpn, vypn = yy[:, 0], yy[:, 1], yy[:, 2], yy[:, 3]
    rpn = np.hypot(xpn, ypn)
    phipn = np.unwrap(np.arctan2(ypn, xpn))
    write_csv(raw / f"pn_plot_{cid}.csv", {
        "row": np.arange(len(tq), dtype=int),
        "t": tq, "x": xpn, "y": ypn, "vx": vxpn, "vy": vypn,
        "r": rpn, "phi": phipn,
    })

    # ---------------------------------------------------------------- S(n) + same 2.5PN local RR recurrence
    nsteps = int(round(steps * radial_orbits))
    dchi = 2.0 * math.pi / steps
    c = coeffs(order)
    chi = np.arange(nsteps + 1, dtype=float) * dchi
    pp = np.empty(nsteps + 1)
    ee = np.empty(nsteps + 1)
    phi = np.empty(nsteps + 1)
    tosc = np.empty(nsteps + 1)
    pp[0], ee[0], phi[0], tosc[0] = p0, e0, 0.0, 0.0

    # diagnostics are evaluated at row k before the k -> k+1 update.
    diag_names = ["r_local", "rdot", "vt", "v2", "Edot", "pdot", "edot", "C_rr", "an_rr", "bv_rr"]
    diag = {name: np.full(nsteps + 1, np.nan) for name in diag_names}

    for k in range(nsteps):
        _, _, _, d = rr_rates_chi(chi[k], pp[k], ee[k], nu)
        for name in diag_names:
            diag[name][k] = d[name]

        p1, e1, dt = rk4_elements(chi[k], pp[k], ee[k], dchi, nu)
        if not (p1 > 0.0 and 0.0 < e1 < 1.0):
            raise RuntimeError(f"{cid}: invalid elements at step {k}: p={p1}, e={e1}")
        pp[k + 1], ee[k + 1] = p1, e1
        tosc[k + 1] = tosc[k] + dt
        cm = chi[k] + 0.5 * dchi
        pm = 0.5 * (pp[k] + pp[k + 1])
        em = 0.5 * (ee[k] + ee[k + 1])
        phi[k + 1] = phi[k] + gN_scalar(cm, pm, em, c) * dchi

    _, _, _, dlast = rr_rates_chi(chi[-1], pp[-1], ee[-1], nu)
    for name in diag_names:
        diag[name][-1] = dlast[name]

    rsn_target = pp / (1.0 + ee * np.cos(chi))
    theta = 0.5 * phi
    X = np.zeros((nsteps + 1, 2), dtype=float)
    X[0] = [math.sqrt(rsn_target[0]), 0.0]
    det_err = np.full(nsteps + 1, np.nan)
    map_err = np.full(nsteps + 1, np.nan)

    for k in range(nsteps):
        rho = math.sqrt(rsn_target[k + 1] / rsn_target[k])
        S = rot(theta[k + 1]) @ np.diag([rho, 1.0 / rho]) @ rot(-theta[k])
        X[k + 1] = S @ X[k]
        det_err[k] = abs(np.linalg.det(S) - 1.0)
        target = math.sqrt(rsn_target[k + 1]) * np.array([
            math.cos(theta[k + 1]), math.sin(theta[k + 1])])
        map_err[k] = float(np.linalg.norm(X[k + 1] - target))

    Y = Phi(X)
    xsn, ysn = Y[:, 0], Y[:, 1]
    rsn = np.hypot(xsn, ysn)
    phisn = np.unwrap(np.arctan2(ysn, xsn))
    sn_cols = {
        "k": np.arange(nsteps + 1, dtype=int),
        "chi": chi, "t_osc": tosc, "p": pp, "e": ee,
        "x": xsn, "y": ysn, "r": rsn, "phi": phisn,
        "X0": X[:, 0], "X1": X[:, 1],
        "detS_error": det_err, "state_map_error": map_err,
    }
    sn_cols.update(diag)
    write_csv(raw / f"sn_{cid}.csv", sn_cols)

    # ---------------------------------------------------------------- metrics
    pn_prec = precession_from_rows(phipn, rpn)
    sn_prec = precession_from_rows(phisn, rsn)
    rms, rms_over_a = radial_rms_by_phi(phipn, rpn, phisn, rsn, a0)
    ipn, rp_pn, ph_pn = periapsis_values(phipn, rpn)
    isn, rp_sn, ph_sn = periapsis_values(phisn, rsn)

    nperi = int(min(len(rp_pn), len(rp_sn)))
    peri_initial_pn = float(rp_pn[0])
    peri_initial_sn = float(rp_sn[0])
    peri_last_pn = float(rp_pn[nperi - 1]) if nperi else float("nan")
    peri_last_sn = float(rp_sn[nperi - 1]) if nperi else float("nan")

    metrics = {
        "case_id": cid,
        "p0": p0, "e0": e0, "a0": a0, "b0": b0,
        "q_mass_ratio": qmass, "nu": nu,
        "order_N": order, "steps_per_radial_orbit": steps,
        "radial_orbits_requested": radial_orbits,
        "pn_stop_reason": stop,
        "pn_precession_rad_per_orbit": pn_prec,
        "sn_precession_rad_per_orbit": sn_prec,
        "precession_error_rad": float(sn_prec - pn_prec),
        "radial_rms_M": rms,
        "radial_rms_over_a": rms_over_a,
        "common_periapsis_count": nperi,
        "pn_peri_initial": peri_initial_pn,
        "pn_peri_last_common": peri_last_pn,
        "sn_peri_initial": peri_initial_sn,
        "sn_peri_last_common": peri_last_sn,
        "pn_peri_shrink_fraction": float(1.0 - peri_last_pn / peri_initial_pn) if nperi else float("nan"),
        "sn_peri_shrink_fraction": float(1.0 - peri_last_sn / peri_initial_sn) if nperi else float("nan"),
        "p_final_sn": float(pp[-1]),
        "e_final_sn": float(ee[-1]),
        "max_det_S_minus_1": float(np.nanmax(det_err)),
        "max_state_map_error": float(np.nanmax(map_err)),
        "pn_initial_v_peri": float(init_info["v_periapsis"]),
        "pn_initial_r_peri": float(init_info["r_periapsis"]),
        "pn_target_r_apo": float(init_info["r_apoapsis_target"]),
    }
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--order", type=int, default=DEFAULT_ORDER)
    ap.add_argument("--steps-per-radial-orbit", type=int, default=DEFAULT_STEPS)
    ap.add_argument("--radial-orbits", type=float, default=DEFAULT_RADIAL_ORBITS)
    ap.add_argument("--jobs", type=int, default=min(3, os.cpu_count() or 1))
    args = ap.parse_args()

    here = Path(__file__).resolve().parent.parent
    base = Path(args.out) if args.out else here
    (base / "raw").mkdir(parents=True, exist_ok=True)
    (base / "summary").mkdir(parents=True, exist_ok=True)

    cases = []
    for ip, p0 in enumerate(P_VALUES):
        for ie, e0 in enumerate(E_VALUES):
            cases.append({
                "case_id": case_id(p0, e0),
                "p_index": ip,
                "e_index": ie,
                "p0": p0,
                "e0": e0,
            })

    with open(base / "summary" / "cases.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "03 nine-case S(n) transferability",
            "fixed_rule": "same S(n) order and same local 2.5PN RR recurrence in all cases",
            "q_mass_ratio": Q_MASS,
            "order_N": args.order,
            "steps_per_radial_orbit": args.steps_per_radial_orbit,
            "radial_orbits": args.radial_orbits,
            "cases": cases,
        }, f, ensure_ascii=False, indent=2)

    metrics = []
    if args.jobs <= 1:
        for c in cases:
            print(f"running {c['case_id']} ...", flush=True)
            metrics.append(run_case(c, str(base), args.order, args.steps_per_radial_orbit, args.radial_orbits))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_case, c, str(base), args.order,
                              args.steps_per_radial_orbit, args.radial_orbits): c for c in cases}
            for fut in as_completed(futs):
                c = futs[fut]
                m = fut.result()
                print(f"done {c['case_id']}: RMS/a={m['radial_rms_over_a']:.5g}, "
                      f"dprec={m['precession_error_rad']:.5g}", flush=True)
                metrics.append(m)

    metrics.sort(key=lambda m: (P_VALUES.index(m["p0"]), E_VALUES.index(m["e0"])))

    fields = list(metrics[0].keys())
    with open(base / "summary" / "case_metrics.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for m in metrics:
            w.writerow(m)

    payload = {
        "experiment": "03 nine-case transferability",
        "interpretation_scope": "representation/transferability test; not a derivation of GR from S(n)",
        "same_Sn_rule_all_cases": True,
        "q_mass_ratio": Q_MASS,
        "order_N": args.order,
        "steps_per_radial_orbit": args.steps_per_radial_orbit,
        "radial_orbits": args.radial_orbits,
        "p_values": list(P_VALUES),
        "e_values": list(E_VALUES),
        "metrics": metrics,
    }
    with open(base / "summary" / "case_metrics.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
