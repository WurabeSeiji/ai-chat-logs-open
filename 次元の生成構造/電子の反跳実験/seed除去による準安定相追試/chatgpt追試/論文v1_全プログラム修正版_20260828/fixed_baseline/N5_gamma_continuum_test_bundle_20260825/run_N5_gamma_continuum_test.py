#!/usr/bin/env python3
import argparse
import importlib.util
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def load_engine(path):
    spec = importlib.util.spec_from_file_location("engine", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def trajectory(engine, n_den, tau_max=70.0, seed=0):
    gamma = math.tan(math.pi / n_den)  # 記録用
    engine.ANGLE = 2.0 * math.pi / n_den  # R1/R3(iii): 線形回転の刻み角を掃引
    dphi = engine.ANGLE

    N = 5
    M = N * (N - 1) // 2
    steps = int(math.ceil(tau_max / dphi))

    syslr = engine.LowRankSystem(N)
    rng = np.random.default_rng(40260722 + 1000 * N + seed)
    v, residual, sig = engine.make_parent(syslr, rng, iters=1200, tol=1e-12)

    Z = v.copy()
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    rows = []
    for t in range(steps + 1):
        h1 = abs(p @ Z) ** 2 + abs(q @ Z) ** 2
        norm2 = float(np.vdot(Z, Z).real)
        f = max(0.0, 1.0 - h1 / norm2)
        rows.append((t, t * dphi, f, abs(Z @ Z), norm2))
        if t < steps:
            syslr.set_state(Z)  # A4
            Z = syslr.linear_rotation_step(Z)  # R1

    df = pd.DataFrame(rows, columns=["step", "tau", "f", "abs_ztz", "norm2"])
    return df, residual, gamma, dphi

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", default="run_n_scaling_lowrank_v1_abprobe.py")
    ap.add_argument("--tau-max", type=float, default=70.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="N5_gamma_continuum_output")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    engine = load_engine(Path(args.engine))
    n_den_values = [144, 288, 576, 1152, 2304]

    trajs = {}
    meta = []
    for nd in n_den_values:
        tr, residual, gamma, dphi = trajectory(engine, nd, args.tau_max, args.seed)
        trajs[nd] = tr
        meta.append((nd, gamma, dphi, residual))

    # common tau comparison against finest step
    grid = np.linspace(0.0, args.tau_max, int(args.tau_max * 100) + 1)
    interp = {nd: np.interp(grid, tr["tau"], tr["f"]) for nd, tr in trajs.items()}
    ref = interp[2304]
    mask = (ref > 1e-10) & (ref < 0.8)

    stats = []
    for nd in n_den_values:
        x = interp[nd]
        stats.append((
            nd,
            float(np.sqrt(np.mean((x[mask] - ref[mask]) ** 2))),
            float(np.sqrt(np.mean(((x[mask] - ref[mask]) / ref[mask]) ** 2))),
            float(np.max(np.abs(x[mask] - ref[mask]))),
        ))
    stats_df = pd.DataFrame(
        stats,
        columns=["n_den", "abs_rmse_vs_2304", "rel_rmse_vs_2304", "max_abs_vs_2304"]
    )

    rates = []
    crossings = []
    for nd, tr in trajs.items():
        m = (tr["f"] > 1e-10) & (tr["f"] < 1e-3)
        rate = float(np.polyfit(tr.loc[m, "tau"], np.log(tr.loc[m, "f"]), 1)[0])
        rates.append((nd, rate))
        for thr in [1e-10, 1e-6, 1e-3, 0.05, 0.5]:
            ix = np.where(tr["f"].values > thr)[0]
            tau_cross = float(tr["tau"].iloc[ix[0]]) if len(ix) else np.nan
            crossings.append((nd, thr, tau_cross))

    rates_df = pd.DataFrame(rates, columns=["n_den", "rate_per_radian"])
    crossings_df = pd.DataFrame(crossings, columns=["n_den", "threshold", "tau_cross"])

    combined = []
    for nd, tr in trajs.items():
        tt = tr.copy()
        tt["n_den"] = nd
        tt["gamma"] = math.tan(math.pi / nd)
        tt["dphi"] = 2.0 * math.pi / nd
        combined.append(tt)
    combined_df = pd.concat(combined, ignore_index=True)

    combined_df.to_csv(out / "N5_gamma_continuum_all_timeseries.csv", index=False)
    stats_df.to_csv(out / "N5_gamma_continuum_stats.csv", index=False)
    rates_df.to_csv(out / "N5_gamma_continuum_rates.csv", index=False)
    crossings_df.to_csv(out / "N5_gamma_continuum_crossings.csv", index=False)

    plt.figure(figsize=(8, 5))
    for nd, tr in trajs.items():
        plt.semilogy(tr["tau"], np.maximum(tr["f"], 1e-18), label=f"nγ={nd}")
    plt.xlabel("cumulative phase tau [rad]")
    plt.ylabel("f outside parent plane")
    plt.title("N=5: small-step convergence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "N5_gamma_continuum_convergence.png", dpi=180)
    plt.close()

    print(stats_df)
    print(rates_df)
    print(crossings_df)

if __name__ == "__main__":
    main()
