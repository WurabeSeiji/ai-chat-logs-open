#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D-sweep runner: N=7,8 interference-preserving dynamics, Delta tau = 2*pi/D.

Dynamics fixed by CLAUDE_CODE_EXPERIMENT_INSTRUCTION_D_over_N_tau_sweep_20260902.md
section 2: H_ef = A_ef conj(z_e) z_f, one step z' = exp(-i (2pi/D) H) z via
numpy.linalg.eigh, H rebuilt from current z each step. float64/complex128.
No seed, no clipping, no renormalization, no noise, no state-dependent branching.

Stage A: D = 2..256 all integers plus 320, 384, 512; 500 steps (501 points).
Stage B: emphasis set; S_max = ceil(500*D/N) steps (same tau window as D=N/500).

Checkpoint: one run = one CSV under data/N{N}/D{D:04d}/; a run whose
run_info_stage{X}.json reports status "ok" (or a recorded failure) is skipped.
"""
from __future__ import annotations
import argparse, csv, json, math, os, sys, time, platform
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

STAGE_A_D = list(range(2, 257)) + [320, 384, 512]
STAGE_B_D = list(range(2, 33)) + [40, 48, 64, 80, 96, 112, 124, 128,
                                  160, 192, 224, 256, 320, 384, 512]
STAGE_A_STEPS = 500

WARN_DRIFT = 1e-10
FAIL_DRIFT = 1e-7


def adjacency(N: int) -> np.ndarray:
    ea, eb = np.triu_indices(N, k=1)
    M = len(ea)
    A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        share[e] = False
        A[e, share] = 1.0
    return A


def plane(z0: np.ndarray):
    p = z0.real.astype(np.float64).copy()
    p /= np.linalg.norm(p)
    q = z0.imag.astype(np.float64).copy()
    q -= np.dot(q, p) * p
    q /= np.linalg.norm(q)
    return p, q


def run_one(N: int, D: int, steps: int, z0: np.ndarray, A: np.ndarray,
            p: np.ndarray, q: np.ndarray, stage: str, outdir: Path) -> dict:
    M = len(z0)
    r2bar = float(np.vdot(z0, z0).real / M)
    H0 = float(np.vdot(z0, z0).real)
    dt = np.float64(2.0 * math.pi / D)
    chi_coef = 2.0 * r2bar * (N - 2)

    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"timeseries_stage{stage}.csv"
    info_path = outdir / f"run_info_stage{stage}.json"

    z = z0.copy()
    rows = []
    max_drift = 0.0
    max_herm = 0.0
    status = "ok"
    stopped_at = None
    t0 = time.time()
    for t in range(steps + 1):
        finite = bool(np.all(np.isfinite(z)))
        h = float(np.vdot(z, z).real)
        zp = z - p * np.dot(p, z) - q * np.dot(q, z)
        hp = float(np.vdot(zp, zp).real)
        f_perp = hp / h
        drift = abs(h - H0) / H0
        max_drift = max(max_drift, drift)
        a2 = np.abs(z)
        s2 = a2 * a2
        pr = float((s2.sum() ** 2) / np.dot(s2, s2))
        tau = float(dt * t)
        rows.append((N, D, D / N, t, tau, t / D, chi_coef * tau,
                     f_perp, 1.0 - f_perp, h, drift,
                     float(abs(z @ z) / h), pr, pr / M,
                     float(a2.min()), float(a2.max()), float(a2.std()),
                     finite, max_herm if t else 0.0))
        if not finite:
            status = "failure_nonfinite"
            stopped_at = t
            break
        if drift > FAIL_DRIFT:
            status = "failure_drift_gt_1e-7"
            stopped_at = t
            break
        if t < steps:
            H = A * (np.conj(z)[:, None] * z[None, :])
            herm = float(np.linalg.norm(H - H.conj().T) /
                         max(np.linalg.norm(H), 1e-300))
            max_herm = max(max_herm, herm)
            w, V = np.linalg.eigh(H)
            z = V @ (np.exp(-1j * dt * w) * (V.conj().T @ z))
    elapsed = time.time() - t0

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N", "D", "D_over_N", "step", "tau", "cycles", "chi",
                    "Hperp_frac", "H_parallel_frac", "H_total",
                    "H_total_rel_drift", "global_closure", "PR", "PR_over_M",
                    "amp_min", "amp_max", "amp_std", "finite", "herm_rel_max"])
        for r in rows:
            w.writerow([r[0], r[1], f"{r[2]:.10g}", r[3], f"{r[4]:.12g}",
                        f"{r[5]:.12g}", f"{r[6]:.12g}", f"{r[7]:.12g}",
                        f"{r[8]:.12g}", f"{r[9]:.15g}", f"{r[10]:.6g}",
                        f"{r[11]:.12g}", f"{r[12]:.12g}", f"{r[13]:.12g}",
                        f"{r[14]:.12g}", f"{r[15]:.12g}", f"{r[16]:.12g}",
                        r[17], f"{r[18]:.6g}"])
    info = {"N": N, "D": D, "stage": stage, "steps_requested": steps,
            "rows_written": len(rows), "status": status,
            "stopped_at": stopped_at, "max_rel_drift": max_drift,
            "drift_warning_gt_1e-10": bool(max_drift > WARN_DRIFT),
            "max_herm_rel": max_herm, "elapsed_sec": elapsed}
    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=1)
    return info


def already_done(outdir: Path, stage: str) -> bool:
    info_path = outdir / f"run_info_stage{stage}.json"
    if not info_path.exists():
        return False
    try:
        info = json.load(open(info_path))
    except Exception:
        return False
    return info.get("status", "").startswith(("ok", "failure"))


def audit_initial(N: int, z0: np.ndarray, p, q, sha: str) -> dict:
    M = len(z0)
    h = float(np.vdot(z0, z0).real)
    zp = z0 - p * np.dot(p, z0) - q * np.dot(q, z0)
    return {"N": N, "M_edge": M, "norm": float(np.linalg.norm(z0)),
            "mean_amp2": h / M, "H_total": h,
            "Hperp_frac_step0": float(np.vdot(zp, zp).real / h),
            "global_closure_step0": float(abs(z0 @ z0) / h),
            "parent_sha256": sha, "dtype": str(z0.dtype)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["A", "B"], required=True)
    ap.add_argument("--n-list", default="7,8")
    args = ap.parse_args()
    Ns = [int(x) for x in args.n_list.split(",")]

    audits = {}
    for N in Ns:
        pf = ROOT / "data" / f"N{N}" / "parent_v.npz"
        import hashlib
        sha = hashlib.sha256(pf.read_bytes()).hexdigest()
        z0 = np.load(pf)["v"].astype(np.complex128)
        A = adjacency(N)
        p, q = plane(z0)
        audits[N] = audit_initial(N, z0, p, q, sha)
        Ds = STAGE_A_D if args.stage == "A" else STAGE_B_D
        for D in Ds:
            steps = STAGE_A_STEPS if args.stage == "A" else \
                math.ceil(500 * D / N)
            outdir = ROOT / "data" / f"N{N}" / f"D{D:04d}"
            if already_done(outdir, args.stage):
                continue
            info = run_one(N, D, steps, z0, A, p, q, args.stage, outdir)
            print(f"stage{args.stage} N={N} D={D} steps={steps} "
                  f"status={info['status']} drift={info['max_rel_drift']:.3g} "
                  f"{info['elapsed_sec']:.1f}s", flush=True)
    ap_path = ROOT / "results" / f"initial_state_audit_stage{args.stage}.json"
    ap_path.parent.mkdir(exist_ok=True)
    with open(ap_path, "w", encoding="utf-8") as f:
        json.dump({"audits": audits,
                   "environment": {"python": platform.python_version(),
                                   "numpy": np.__version__,
                                   "platform": platform.platform(),
                                   "machine": platform.machine()}},
                  f, indent=1)
    print("STAGE", args.stage, "COMPLETE")


if __name__ == "__main__":
    main()
