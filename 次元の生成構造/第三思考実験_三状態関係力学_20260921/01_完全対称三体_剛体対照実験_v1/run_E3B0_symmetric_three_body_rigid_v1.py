#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3-B0 v1: fully symmetric three-state rigid-body control experiment.

Purpose
-------
Test the first genuine N=3 control point: if three anonymous states and all
pair couplings are exactly permutation-symmetric, does any body-specific
readout (mass-like directed ratio, charge-like pair difference, geometry)
appear spontaneously?

Candidate law (complete graph K3, common scalar coupling eps):
    X_i[k+1] + X_i[k-1]
      = tau X_i[k] + eps * sum_{j != i} X_j[k]
for i in {a,b,c}, X_i in R^2.

The initial configuration is an equilateral triangle centered at zero.  tau is
chosen for each eps so that the centroid-zero relative mode has the same closed
period n=31:
    lambda_rel = tau - eps = 2 cos(2 pi / n).
Thus tau = lambda_rel + eps.

This is a deliberately symmetric control, not a derived final three-body law.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SUMMARY_CSV = HERE / "E3B0_symmetric_three_body_summary_v1.csv"
TRAJ_CSV = HERE / "E3B0_symmetric_three_body_trajectories_v1.csv"
RESULTS_JSON = HERE / "E3B0_symmetric_three_body_results_v1.json"

N_PERIOD = 31
STEPS = N_PERIOD * 100
EPS_LIST = [0.0, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2]
BODY_NAMES = ["a", "b", "c"]
TOL = 1e-10

FAIL_CONDITIONS = {
    "E3B1_centroid_closure": "centroid/closure norm exceeds 1e-10 relative to state scale",
    "E3B2_rigid_pair_distances": "pair-distance equality or drift exceeds 1e-10",
    "E3B3_area_conservation": "oriented triangle-area relative drift exceeds 1e-10",
    "E3B4_readout_rigidity": "quadratic-readout triangle pair-distance equality or drift exceeds 1e-10",
    "E3B5_directed_couplings": "for eps>0, fitted six directed couplings differ from eps by more than 1e-10 relative",
    "E3B6_masslike_ratios": "for eps>0, any directed reverse-pair ratio differs from 1 by more than 1e-10",
    "E3B7_no_body_specific_split": "body-wise radius or one-step area readouts split by more than 1e-10 relative",
}


def rot(t: float) -> np.ndarray:
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]], dtype=float)


def cross(u: np.ndarray, v: np.ndarray) -> float:
    return float(u[0] * v[1] - u[1] * v[0])


def pair_distances(X: np.ndarray) -> np.ndarray:
    return np.array([
        np.linalg.norm(X[0] - X[1]),
        np.linalg.norm(X[1] - X[2]),
        np.linalg.norm(X[2] - X[0]),
    ], dtype=float)


def readout_points(X: np.ndarray) -> np.ndarray:
    w = X[:, 0] + 1j * X[:, 1]
    return w * w


def readout_pair_distances(z: np.ndarray) -> np.ndarray:
    return np.array([abs(z[0]-z[1]), abs(z[1]-z[2]), abs(z[2]-z[0])], dtype=float)


def triangle_area(X: np.ndarray) -> float:
    return 0.5 * cross(X[1]-X[0], X[2]-X[0])


def initial_triangle(theta: float):
    phases = np.array([0.0, 2*math.pi/3, 4*math.pi/3])
    X0 = np.c_[np.cos(phases), np.sin(phases)]
    Rm = rot(-theta)
    Xm1 = X0 @ Rm.T
    return Xm1, X0


def fit_directed_couplings(Xs: np.ndarray, tau: float):
    # For each receiver i, fit scalar coefficients c_{i<-j}, c_{i<-l} from
    # Y_i(k) = X_i(k+1)+X_i(k-1)-tau X_i(k)
    #        = c_ij X_j(k) + c_il X_l(k).
    fits = {}
    rel_resids = {}
    K = Xs.shape[0] - 2
    for i in range(3):
        js = [j for j in range(3) if j != i]
        A_rows, y_rows = [], []
        for k in range(1, K+1):
            y = Xs[k+1, i] + Xs[k-1, i] - tau * Xs[k, i]
            for comp in range(2):
                A_rows.append([Xs[k, js[0], comp], Xs[k, js[1], comp]])
                y_rows.append(y[comp])
        A = np.asarray(A_rows)
        y = np.asarray(y_rows)
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ coef
        rel = float(np.linalg.norm(pred-y) / max(np.linalg.norm(y), 1e-30))
        fits[(i, js[0])] = float(coef[0])
        fits[(i, js[1])] = float(coef[1])
        rel_resids[i] = rel
    return fits, rel_resids


def simulate(eps: float):
    theta = 2*math.pi/N_PERIOD
    lam_rel = 2*math.cos(theta)
    tau = lam_rel + eps
    Xm1, X0 = initial_triangle(theta)

    Xs = np.empty((STEPS+1, 3, 2), dtype=float)
    Xs[0] = X0
    # exact first step from recurrence using Xm1
    sums0 = np.sum(X0, axis=0)
    for i in range(3):
        others = sums0 - X0[i]
        Xs[1, i] = tau*X0[i] + eps*others - Xm1[i]
    for k in range(1, STEPS):
        total = np.sum(Xs[k], axis=0)
        for i in range(3):
            others = total - Xs[k, i]
            Xs[k+1, i] = tau*Xs[k, i] + eps*others - Xs[k-1, i]

    centroid = np.mean(Xs, axis=1)
    centroid_norm = np.linalg.norm(centroid, axis=1)
    state_scale = float(np.mean(np.linalg.norm(Xs, axis=2)))
    centroid_rel_max = float(np.max(centroid_norm) / max(state_scale, 1e-30))

    pd = np.array([pair_distances(X) for X in Xs])
    pd_eq = float(np.max(np.ptp(pd, axis=1)) / max(np.mean(pd), 1e-30))
    pd_drift = float(np.max(np.abs(pd - pd[0])) / max(np.mean(pd[0]), 1e-30))

    areas = np.array([triangle_area(X) for X in Xs])
    area_rel_drift = float(np.max(np.abs(areas-areas[0])) / max(abs(areas[0]), 1e-30))

    Zs = np.array([readout_points(X) for X in Xs])
    rpd = np.array([readout_pair_distances(z) for z in Zs])
    rpd_eq = float(np.max(np.ptp(rpd, axis=1)) / max(np.mean(rpd), 1e-30))
    rpd_drift = float(np.max(np.abs(rpd-rpd[0])) / max(np.mean(rpd[0]), 1e-30))

    radii = np.linalg.norm(Xs, axis=2)
    body_radius_split = float(np.max(np.ptp(radii, axis=1)) / max(np.mean(radii), 1e-30))
    qbody = np.empty((STEPS,3))
    for k in range(STEPS):
        for i in range(3):
            qbody[k,i] = cross(Xs[k,i], Xs[k+1,i])
    qbody_split = float(np.max(np.ptp(qbody, axis=1)) / max(np.mean(np.abs(qbody)), 1e-30))

    fits, fit_resids = fit_directed_couplings(Xs, tau)
    fit_vals = np.array(list(fits.values()))
    if eps > 0:
        coupling_rel_err = float(np.max(np.abs(fit_vals-eps))/eps)
        coupling_spread = float(np.ptp(fit_vals)/eps)
        ratios = np.array([
            fits[(0,1)]/fits[(1,0)],
            fits[(1,2)]/fits[(2,1)],
            fits[(2,0)]/fits[(0,2)],
        ])
        mass_ratio_max_dev = float(np.max(np.abs(ratios-1)))
        mass_ratio_cycle = float(np.prod(ratios))
    else:
        coupling_rel_err = coupling_spread = 0.0
        ratios = np.array([np.nan,np.nan,np.nan])
        mass_ratio_max_dev = 0.0
        mass_ratio_cycle = np.nan

    summary = {
        "eps": eps,
        "tau": tau,
        "lambda_relative": lam_rel,
        "lambda_centroid": tau + 2*eps,
        "centroid_rel_max": centroid_rel_max,
        "pair_distance_equality_error": pd_eq,
        "pair_distance_drift": pd_drift,
        "triangle_area_rel_drift": area_rel_drift,
        "readout_pair_distance_equality_error": rpd_eq,
        "readout_pair_distance_drift": rpd_drift,
        "body_radius_split": body_radius_split,
        "body_q_split": qbody_split,
        "directed_coupling_rel_error": coupling_rel_err,
        "directed_coupling_spread": coupling_spread,
        "fit_resid_max": float(max(fit_resids.values())),
        "mass_ratio_ab": None if eps == 0 else float(ratios[0]),
        "mass_ratio_bc": None if eps == 0 else float(ratios[1]),
        "mass_ratio_ca": None if eps == 0 else float(ratios[2]),
        "mass_ratio_max_dev": mass_ratio_max_dev,
        "mass_ratio_cycle_product": None if eps == 0 else mass_ratio_cycle,
    }
    directed = {f"{BODY_NAMES[i]}<-{BODY_NAMES[j]}": v for (i,j),v in fits.items()}
    return Xs, Zs, pd, rpd, areas, qbody, summary, directed


def main():
    all_summaries = []
    all_results = []
    traj_rows = []

    for eps in EPS_LIST:
        Xs, Zs, pd, rpd, areas, qbody, summary, directed = simulate(eps)
        all_summaries.append(summary)
        all_results.append({"summary": summary, "directed_couplings": directed})
        for k in range(STEPS+1):
            row = {"eps":eps,"step":k,"area":areas[k]}
            for i,name in enumerate(BODY_NAMES):
                row[f"{name}_x"] = Xs[k,i,0]
                row[f"{name}_y"] = Xs[k,i,1]
                row[f"{name}_readout_x"] = Zs[k,i].real
                row[f"{name}_readout_y"] = Zs[k,i].imag
            row["d_ab"], row["d_bc"], row["d_ca"] = pd[k]
            row["rd_ab"], row["rd_bc"], row["rd_ca"] = rpd[k]
            if k < STEPS:
                row["q_a"], row["q_b"], row["q_c"] = qbody[k]
            else:
                row["q_a"] = row["q_b"] = row["q_c"] = ""
            traj_rows.append(row)

    # verdicts across all eps
    nonzero = [s for s in all_summaries if s["eps"] > 0]
    verdict = {
        "E3B1_centroid_closure": max(s["centroid_rel_max"] for s in all_summaries) <= TOL,
        "E3B2_rigid_pair_distances": max(max(s["pair_distance_equality_error"],s["pair_distance_drift"]) for s in all_summaries) <= TOL,
        "E3B3_area_conservation": max(s["triangle_area_rel_drift"] for s in all_summaries) <= TOL,
        "E3B4_readout_rigidity": max(max(s["readout_pair_distance_equality_error"],s["readout_pair_distance_drift"]) for s in all_summaries) <= TOL,
        "E3B5_directed_couplings": max(s["directed_coupling_rel_error"] for s in nonzero) <= TOL,
        "E3B6_masslike_ratios": max(s["mass_ratio_max_dev"] for s in nonzero) <= TOL,
        "E3B7_no_body_specific_split": max(max(s["body_radius_split"],s["body_q_split"]) for s in all_summaries) <= TOL,
    }

    # summary csv
    cols = list(all_summaries[0].keys())
    with open(SUMMARY_CSV,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(all_summaries)

    # trajectory csv
    cols2 = list(traj_rows[0].keys())
    with open(TRAJ_CSV,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols2); w.writeheader(); w.writerows(traj_rows)

    out = {
        "experiment":"E3-B0 symmetric three-body rigid control v1",
        "status":"fully symmetric N=3 control; not final derived three-body law",
        "law":"X_i[k+1]+X_i[k-1]=tau X_i[k]+eps sum_{j!=i} X_j[k]",
        "period":N_PERIOD,
        "steps":STEPS,
        "eps_list":EPS_LIST,
        "fail_conditions":FAIL_CONDITIONS,
        "verdict":{k:("PASS" if v else "FAIL") for k,v in verdict.items()},
        "runs":all_results,
    }
    with open(RESULTS_JSON,"w",encoding="utf-8") as f:
        json.dump(out,f,ensure_ascii=False,indent=2)

    print(json.dumps({"verdict":out["verdict"],"summaries":all_summaries,"files":[str(SUMMARY_CSV),str(TRAJ_CSV),str(RESULTS_JSON)]},ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
