#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deterministic mp-free hm parent generator, N=3..40.

Compatibility policy:
- N=3..16: reproduce the historical hm scale using frozen legacy norms.
  These constants are compatibility data only; no mp calculation is performed.
- N>=17: use the intrinsic hm normalization r^2 = 1/15 with no adjustment.
- Phase construction is analytic for every N:
  N=3 special Z3; even N one-factorization; odd N>=5 cyclic distance classes.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np

RBAR2 = 1.0 / 15.0
LEGACY_NORM = {
    3: 0.4472135954999693,
    4: 0.8164965809276707,
    5: 0.9271726499455791,
    6: 1.095445115010219,
    7: 1.1766901825131457,
    8: 1.2815849983171124,
    9: 1.352987915930318,
    10: 1.4857455929760082,
    11: 1.5611551506001249,
    12: 1.6270724736363549,
    13: 1.7103356194389499,
    14: 1.770758915601077,
    15: 1.8328264435092678,
    16: 1.9097313919229777,
}

def edges(N: int):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]

def classes(N: int):
    if N < 3:
        raise ValueError("N must be >= 3")
    if N == 3:
        return {(0, 2): 0, (1, 2): 1, (0, 1): 2}, 3, 60.0
    if N == 4:
        return {(0, 3): 0, (1, 2): 0, (0, 2): 1, (1, 3): 1,
                (0, 1): 2, (2, 3): 2}, 3, 60.0
    if N % 2 == 0:
        n = N - 1
        col = {}
        for rr in range(n):
            col[tuple(sorted((rr, N - 1)))] = rr
            for k in range(1, N // 2):
                col[tuple(sorted(((rr - k) % n, (rr + k) % n)))] = rr
        return col, n, 180.0 / (N - 1)
    q = (N - 1) // 2
    col = {}
    for d in range(1, q + 1):
        for i in range(N):
            col[tuple(sorted((i, (i + d) % N)))] = d - 1
    return col, q, 180.0 / q

def equimodular_base(N: int):
    E = edges(N)
    col, q, step = classes(N)
    theta = np.array([np.radians(step * col[e]) for e in E])
    v = np.sqrt(RBAR2) * np.exp(1j * theta)
    return v, E, col, q, step

def hm_parent(N: int):
    v, E, col, q, step = equimodular_base(N)
    if N <= 16:
        # Deliberately identical arithmetic to historical pass1:
        # v = sp.equimodular(N); v = v*NORM[N]/np.linalg.norm(v)
        v = v * LEGACY_NORM[N] / np.linalg.norm(v)
        scale_policy = "legacy_norm_compatibility"
    else:
        # No mp-derived adjustment. Intrinsic common rule.
        scale_policy = "r2_equals_1_over_15"
    return v, E, col, q, step, scale_policy

def diagnostics(N: int, v: np.ndarray, E):
    S = np.zeros(N, dtype=complex)
    for k, (i, j) in enumerate(E):
        S[i] += v[k] ** 2
        S[j] += v[k] ** 2
    global_sum = np.sum(v ** 2)
    r2 = float(np.vdot(v, v).real / len(v))
    M = len(E)
    A = np.zeros((M, M), dtype=float)
    for a in range(M):
        ea = set(E[a])
        for b in range(a + 1, M):
            if ea.intersection(E[b]):
                A[a,b] = A[b,a] = 1.0
    H = A * (np.conj(v)[:,None] * v[None,:])
    Hv = H @ v
    mu_H = float((np.vdot(v,Hv)/np.vdot(v,v)).real)
    residual_H = float(np.linalg.norm(Hv-mu_H*v)/np.linalg.norm(v))
    # For N>=4 local closure gives exact analytic eigenvalue mu_H=-2 r^2
    # in the interference-H representation. N=3 is the historical special Z3 case.
    return {
        "N": N,
        "M": len(E),
        "norm": float(np.linalg.norm(v)),
        "mean_amp2": r2,
        "amp_min": float(np.abs(v).min()),
        "amp_max": float(np.abs(v).max()),
        "global_sum_z2_abs": float(abs(global_sum)),
        "local_sum_z2_max_abs": float(np.max(np.abs(S))),
        "global_closure_normalized": float(abs(global_sum) / np.vdot(v, v).real),
        "local_closure_normalized": float(np.max(np.abs(S)) / r2),
        "mu_H": mu_H,
        "mu_H_over_r2": float(mu_H/r2),
        "selfconsistency_residual_H": residual_H,
    }

def save_one(root: Path, N: int):
    v, E, col, q, step, scale_policy = hm_parent(N)
    d = root / "data" / f"hm_N{N}"
    d.mkdir(parents=True, exist_ok=True)
    color = np.array([col[e] for e in E], dtype=int)
    design = "handmade_equimodular_" + ("Z3" if N == 3 else "1factor" if N % 2 == 0 else "distance_classes")
    diag = diagnostics(N, v, E)
    diag.update({"design": design, "q": q, "phase_step_deg": step, "scale_policy": scale_policy})
    np.savez_compressed(d / "parent_v.npz", v=v, edges=np.array(E, dtype=int), color=color,
                        theta=np.angle(v), design=design, r=float(np.sqrt(diag["mean_amp2"])),
                        mean_amp2=diag["mean_amp2"], scale_policy=scale_policy)
    with open(d / "parent_v.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["edge_index","i","j","class","theta_deg","a_Re","b_Im","abs_z"])
        for k, (i, j) in enumerate(E):
            w.writerow([k, i, j, int(color[k]), np.degrees(np.angle(v[k])) % 360,
                        v[k].real, v[k].imag, abs(v[k])])
    with open(d / "parent_checks.json", "w", encoding="utf-8") as f:
        json.dump(diag, f, indent=2, ensure_ascii=False)
    return diag

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-min", type=int, default=3)
    ap.add_argument("--n-max", type=int, default=40)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1])
    args = ap.parse_args()
    rows = [save_one(args.out, N) for N in range(args.n_min, args.n_max + 1)]
    with open(args.out / "generation_summary.csv", "w", newline="", encoding="utf-8") as f:
        keys = list(rows[0].keys()); w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"generated N={args.n_min}..{args.n_max} at {args.out}")

if __name__ == "__main__":
    main()
