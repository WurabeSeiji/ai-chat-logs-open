#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
大N自発的分裂走行（note図用・低ランク頂点分解エンジン使用）

    python3 run_spontaneous_splitting_largeN_v1.py N [delta] [--cap=N] [--seed=N]

親（自己無撞着円偏波固有モード）＋零閉鎖核種 δ から逐次再構成で走行し、
休眠フラクション f(τ) = 1 - h_plane1/h_total を毎ステップ記録する。
閾値 f > 0.05 の交差後も after ステップ走行して飽和域まで残す。

出力: largeN_splitting_result_v1/fcurve_N{n}_delta{delta}_seed{seed}.csv / .json
"""

import csv
import json
import math
import os
import sys
import time

import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
from run_n_scaling_lowrank_v1 import (
    GAMMA, LowRankSystem, make_parent, zero_closure_kernel_seed, progress,
)

RESULT_DIR = os.path.join(BASE_DIR, "largeN_splitting_result_v1")


def run(n, delta, seed, cap, after, tol=1e-8):
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + seed)
    t0 = time.time()
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=tol)
    t_parent = time.time() - t0
    progress(f"N={n} 親構成完了 残差={residual:.2e} 平面数={len(sig)} ({t_parent:.1f}s)")

    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    ztz0 = abs(complex(Z @ Z))

    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    f_hist = []
    max_ztz = 0.0
    crossed_at = None
    wp = rng.normal(size=sys_lr.m)
    t0 = time.time()
    for t in range(cap + 1):
        # 直交補への射影で休眠エネルギーを直接測る（1-h1/htot の桁落ち回避）
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        htot = float(np.real(np.conj(Z) @ Z))
        f = float(np.real(np.conj(Zp) @ Zp)) / htot
        f_hist.append(f)
        if crossed_at is None and f > 0.05:
            crossed_at = t
            progress(f"N={n} 閾値交差 τ={t}")
        if crossed_at is not None and t >= crossed_at + after:
            break
        if t % 500 == 0 and t > 0:
            progress(f"N={n} τ={t} f={f:.3e} ({(time.time()-t0)/t*1000:.1f} ms/step)")
        max_ztz = max(max_ztz, abs(complex(Z @ Z)))
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    t_run = time.time() - t0

    f_arr = np.array(f_hist)
    lo, hi = max(10.0 * f_arr[0], 1e-300), 1e-3
    mask = (f_arr > lo) & (f_arr < hi)
    idx = np.where(mask)[0]
    rate = float(np.polyfit(idx, np.log(f_arr[idx]), 1)[0]) if len(idx) >= 5 else None

    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = f"N{n:05d}_delta{delta:.0e}_seed{seed}"
    with open(os.path.join(RESULT_DIR, f"fcurve_{tag}.csv"), "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["tau", "f"])
        for t, f in enumerate(f_hist):
            wtr.writerow([t, f])
    summary = {
        "n": n, "m": sys_lr.m, "delta": delta, "seed": seed,
        "parent_residual": residual,
        "parent_rank_planes": int(len(sig)),
        "abs_ztz_initial": ztz0, "abs_ztz_max": max_ztz,
        "f_initial": float(f_arr[0]), "f_final": float(f_arr[-1]),
        "crossing_tau": crossed_at, "steps_run": len(f_hist) - 1,
        "onset_rate_per_step": rate,
        "steps_per_decade": (math.log(10.0) / rate) if rate else None,
        "t_parent_sec": t_parent, "t_run_sec": t_run,
    }
    with open(os.path.join(RESULT_DIR, f"summary_{tag}.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return summary


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[1:]
             if a.startswith("--") and "=" in a}
    n = int(args[0]) if args else 40
    delta = float(args[1]) if len(args) > 1 else 1e-8
    cap = int(flags.get("--cap", 8000))
    seed = int(flags.get("--seed", 0))
    after = int(flags.get("--after", 2000))
    tol = float(flags.get("--tol", 1e-8))
    run(n, delta, seed, cap, after, tol=tol)


if __name__ == "__main__":
    main()
