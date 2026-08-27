# -*- coding: utf-8 -*-
"""§1 of the memo check: does the N=5 metastable state have a finite recurrence in tau?
Input : data/N5_phase_by_edge_5000steps.csv, data/N5_phase_increments_5000steps.csv (from N5_complex_simplex_complete_analysis_20260826.zip),
        data/floquet_spectrum.csv (from N5_dynamics_followup_theorems_and_stability_20260826.zip)
Output: results/n5_recurrence.json
"""
import csv, json, math, os
import numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "data"); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
GAMMA = math.tan(math.pi / 144.0)

rows = list(csv.DictReader(open(os.path.join(D, "N5_phase_by_edge_5000steps.csv"), encoding="utf-8")))
T = 1 + max(int(r["step"]) for r in rows); E = 10
P = np.zeros((T, E)); A = np.zeros((T, E))
for r in rows:
    P[int(r["step"]), int(r["edge_index"])] = float(r["theta"]); A[int(r["step"]), int(r["edge_index"])] = float(r["amplitude"])
dP = np.angle(np.exp(1j * np.diff(P, axis=0)))
out = {"gamma": GAMMA, "parent_dphi_theory": 2 * math.atan(GAMMA * math.sqrt(14)), "windows": {}}
for w0, w1, name in [(0, 150, "latent"), (300, 450, "expansion"), (1000, 2000, "ordering"), (4000, 5000, "late_metastable")]:
    seg = dP[w0:w1]; m = seg.mean(axis=0); dphi = float(m.mean())
    out["windows"][name] = {"steps": [w0, w1], "edge_mean_increment_min": float(m.min()), "edge_mean_increment_max": float(m.max()),
                            "edge_spread": float(m.max() - m.min()), "time_std_max": float(seg.std(axis=0).max()),
                            "mean_dphi": dphi, "steps_per_turn": 2 * math.pi / dphi, "dphi_over_2pi": dphi / (2 * math.pi),
                            "sigma_eff": math.tan(dphi / 2) / GAMMA, "sigma_eff2": (math.tan(dphi / 2) / GAMMA) ** 2}
rel = np.angle(np.exp(1j * (P - P[:, [0]]))); late = rel[4000:]
out["late_relative_phase_over_pi"] = [float(x) for x in np.round(late.mean(axis=0) / math.pi, 6)]
out["late_relative_phase_time_std"] = [float(x) for x in late.std(axis=0)]
Z = A * np.exp(1j * P)
def dist(n, t0=4000, t1=4900):
    d = []; dg = []
    for t in range(t0, t1 - n):
        a, b = Z[t], Z[t + n]; d.append(np.linalg.norm(a - b) / np.linalg.norm(a))
        g = np.vdot(a, b); dg.append(np.linalg.norm(a - b * np.exp(-1j * np.angle(g))) / np.linalg.norm(a))
    return float(np.mean(d)), float(np.mean(dg))
exact = sorted(((dist(n)[0], n) for n in range(1, 400)), key=lambda x: x[0])[:5]
out["exact_recurrence_best"] = [{"n": n, "rel_dist": d} for d, n in exact]
out["gauge_removed_recurrence"] = {str(n): dist(n)[1] for n in (1, 2, 5, 36, 39, 50)}
m_late = float(dP[4000:].mean()); x = m_late / (2 * math.pi)
out["late_dphi_over_2pi_rational_approx"] = {f"q<={q}": [str(Fraction(x).limit_denominator(q)), abs(x - float(Fraction(x).limit_denominator(q)))] for q in (100, 1000, 10000)}
q = P[4000:] / (math.pi / 72); out["max_deviation_from_pi_over_72_grid"] = float(np.abs(q - np.round(q)).max())
inc = list(csv.DictReader(open(os.path.join(D, "N5_phase_increments_5000steps.csv"), encoding="utf-8")))
late_inc = [r for r in inc if 4000 <= int(r["step"]) < 5000]
out["chatgpt_columns_late"] = {"max_rat_err_q_le_256_mean": float(np.mean([float(r["max_rat_err_q_le_256"]) for r in late_inc])),
                               "best_denominators_last": late_inc[-1]["best_denominators"]}
fl = [r for r in csv.DictReader(open(os.path.join(D, "floquet_spectrum.csv"), encoding="utf-8")) if r["fd_eps"] == "3e-06"]
base = out["parent_dphi_theory"]
out["floquet_angles_over_parent_dphi"] = [{"rank": int(r["rank"]), "modulus": float(r["modulus"]), "angle_over_dphi": math.atan2(float(r["eig_im"]), float(r["eig_re"])) / base} for r in fl if float(r["eig_im"]) >= 0]
json.dump(out, open(os.path.join(R, "n5_recurrence.json"), "w"), indent=1)
w = out["windows"]["late_metastable"]
print(f"late: dphi={w['mean_dphi']:.7f} rad, edge spread={w['edge_spread']:.1e}, sigma_eff^2={w['sigma_eff2']:.6f}, steps/turn={w['steps_per_turn']:.4f}")
print("exact recurrence best:", out["exact_recurrence_best"][:3]); print("gauge-removed n=1:", out["gauge_removed_recurrence"]["1"])
print("neutral floquet angle/dphi:", [round(f["angle_over_dphi"], 4) for f in out["floquet_angles_over_parent_dphi"] if abs(f["modulus"] - 1) < 1e-6])
