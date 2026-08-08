#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
大N自発的分裂走行 + σ_ratio(τ)=σ_2/σ_1 記録（別プログラム・原本は改変しない）

    python3 run_sigma_ratio_largeN_v1.py N [delta] [--cap=N] [--seed=N] [--after=N] [--tol=X] [--sub=K]

run_spontaneous_splitting_largeN_v1.py の run() ループを厳密に複製し、休眠フラクション
f(τ)=1-h_plane1/h_total を毎ステップ記録すると同時に、σ_ratio(τ)=σ_2/σ_1 を K ステップごとに
記録する。σ_ratio は読み取り専用 sigma_spectrum() で計算し Z・rng・wp を一切触らないため、
軌道（f も含む）は原本とビット単位で同一になる。

【内蔵テスト】出力の f 列が原本 fcurve と bit-exact 一致すれば、軌道が同一＝σ_ratio が
正しい軌道上の値であることの証明になる。

出力: sigma_ratio_result_v1/sigmaratio_N{n}_delta{delta}_seed{seed}.csv （列: tau,f,sigma_ratio）
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

RESULT_DIR = os.path.join(BASE_DIR, "sigma_ratio_result_v1")


def sigma_ratio_from_spectrum(sys_lr):
    """現在の生成子の σ_2/σ_1（第2σ/第1σ）。読み取り専用（Z・rng・wp を触らない）。
    平面が1枚以下なら 0.0。"""
    sig = sys_lr.sigma_spectrum()  # 降順の正σ
    if len(sig) >= 2 and sig[0] > 0.0:
        return float(sig[1] / sig[0])
    return 0.0


def run(n, delta, seed, cap, after, tol, sub):
    # --- 原本 run() を厳密複製（親構成・種・参照平面）---
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
    sr_hist = []   # (tau, sigma_ratio) をサブサンプルで記録
    max_ztz = 0.0
    crossed_at = None
    wp = rng.normal(size=sys_lr.m)
    t0 = time.time()
    for t in range(cap + 1):
        # --- 原本と同一の f 計算 ---
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
        # --- 追加: σ_ratio 記録（読み取り専用・軌道に影響しない）---
        if t % sub == 0:
            sr_hist.append((t, sigma_ratio_from_spectrum(sys_lr)))
        # --- 原本と同一の更新（wp・Z を進める）---
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    t_run = time.time() - t0

    # σ_ratio を tau→値 の辞書化（f 全ステップに対し、サンプル点のみ値・他は空）
    sr_map = {int(tau): sr for tau, sr in sr_hist}

    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = f"N{n:05d}_delta{delta:.0e}_seed{seed}"
    out_csv = os.path.join(RESULT_DIR, f"sigmaratio_{tag}.csv")
    with open(out_csv, "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["tau", "f", "sigma_ratio"])
        for t, f in enumerate(f_hist):
            wtr.writerow([t, f, sr_map.get(t, "")])
    summary = {
        "n": n, "m": sys_lr.m, "delta": delta, "seed": seed, "sub": sub,
        "parent_residual": residual,
        "f_initial": float(f_hist[0]), "f_final": float(f_hist[-1]),
        "crossing_tau": crossed_at, "steps_run": len(f_hist) - 1,
        "n_sigma_samples": len(sr_hist),
        "sigma_ratio_final": sr_hist[-1][1] if sr_hist else None,
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
    sub = int(flags.get("--sub", 10))
    run(n, delta, seed, cap, after, tol, sub)


if __name__ == "__main__":
    main()
