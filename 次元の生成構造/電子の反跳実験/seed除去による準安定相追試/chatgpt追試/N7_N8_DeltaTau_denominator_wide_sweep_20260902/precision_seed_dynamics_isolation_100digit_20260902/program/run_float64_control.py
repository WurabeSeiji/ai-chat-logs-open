#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Condition A: IC64 + Dynamics64 (same-environment control rerun).

Dynamics identical to the parent sweep's run_sweep.py (instruction section 2):
H = A*(conj(z)[:,None]*z[None,:]); w,V = numpy.linalg.eigh(H);
z = V@(exp(-1j*(2*pi/D)*w)*(V.conj().T@z)). No seed/clip/renorm.

Also verifies reproduction against the existing sweep timeseries for the
same (N, D=N): stage A (500 steps) and stage B, comparing Hperp_frac.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, platform, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]
SWEEP = HERE.parent  # N7_N8_DeltaTau_denominator_wide_sweep_20260902
STEPS = 2000


def adjacency(N):
    ea, eb = np.triu_indices(N, k=1)
    M = len(ea)
    A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        share[e] = False
        A[e, share] = 1.0
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True, choices=[7, 8])
    args = ap.parse_args()
    N = args.n
    D = N
    M = N * (N - 1) // 2
    pf = SWEEP / "data" / f"N{N}" / "parent_v.npz"
    sha = hashlib.sha256(pf.read_bytes()).hexdigest()
    z0 = np.load(pf)["v"].astype(np.complex128)
    A = adjacency(N)
    p = z0.real.copy(); p /= np.linalg.norm(p)
    q = z0.imag.copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    H0 = float(np.vdot(z0, z0).real)
    dt = np.float64(2.0 * math.pi / D)

    out = HERE / "data" / f"N{N}_D{D}" / "A_IC64_DYN64"
    out.mkdir(parents=True, exist_ok=True)
    z = z0.copy()
    rows = []
    t0 = time.time()
    for t in range(STEPS + 1):
        h = float(np.vdot(z, z).real)
        zp = z - p * np.dot(p, z) - q * np.dot(q, z)
        hp = float(np.vdot(zp, zp).real)
        f = hp / h
        a2 = np.abs(z)
        s2 = a2 * a2
        pr = float((s2.sum() ** 2) / np.dot(s2, s2))
        rows.append([t, dt * t, f"{f:.17g}",
                     ("-inf" if f == 0.0 else f"{math.log10(f):.17g}"),
                     f"{1.0 - f:.17g}", f"{h:.17g}",
                     f"{abs(h - H0) / H0:.6g}",
                     f"{float(abs(z @ z) / h):.17g}",
                     f"{pr:.17g}", f"{pr / M:.17g}",
                     f"{a2.min():.17g}", f"{a2.max():.17g}",
                     f"{a2.std():.17g}", bool(np.all(np.isfinite(z)))])
        if t < STEPS:
            H = A * (np.conj(z)[:, None] * z[None, :])
            w, V = np.linalg.eigh(H)
            z = V @ (np.exp(-1j * dt * w) * (V.conj().T @ z))
    with open(out / "timeseries.csv", "w", newline="", encoding="utf-8") as f_:
        w_ = csv.writer(f_)
        w_.writerow(["step", "tau", "Hperp_frac", "log10_Hperp_frac",
                     "H_parallel_frac", "H_total", "H_total_rel_drift",
                     "global_closure", "PR", "PR_over_M",
                     "amp_min", "amp_max", "amp_std", "finite"])
        w_.writerows(rows)

    # reproduction check against existing sweep (stage A and stage B, D=N)
    repro = {}
    ours = np.array([float(r[2]) for r in rows])
    for stage in ("A", "B"):
        ref = SWEEP / "data" / f"N{N}" / f"D{D:04d}" / \
            f"timeseries_stage{stage}.csv"
        if not ref.exists():
            repro[f"stage{stage}"] = "reference_missing"
            continue
        with open(ref, newline="", encoding="utf-8") as f_:
            r = csv.reader(f_)
            hdr = next(r)
            i_f = hdr.index("Hperp_frac")
            refs = np.array([float(row[i_f]) for row in r])
        n = min(len(refs), len(ours))
        diff = np.abs(ours[:n] - refs[:n])
        repro[f"stage{stage}"] = {
            "compared_steps": int(n),
            "max_abs_diff_Hperp_frac": float(diff.max()),
            "bitwise_identical": bool(diff.max() == 0.0)}
    info = {"condition": "A_IC64_DYN64", "N": N, "D": D, "steps": STEPS,
            "parent_file": str(pf), "parent_sha256": sha,
            "dtype": "complex128/float64",
            "elapsed_sec": time.time() - t0,
            "environment": {"python": platform.python_version(),
                            "numpy": np.__version__,
                            "platform": platform.platform()},
            "reproduction_check": repro}
    with open(out / "run_info.json", "w", encoding="utf-8") as f_:
        json.dump(info, f_, indent=1)
    print(json.dumps(repro, indent=1))
    print(f"A N={N} done, {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
