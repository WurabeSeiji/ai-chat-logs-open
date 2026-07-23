#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
大N自発的分裂走行 + σスペクトル（σ_1, σ_2, …）を「そのまま」記録（別プログラム・原本改変なし）

    python3 run_sigma_spectrum_largeN_v1.py N [delta] [--cap=N] [--seed=N] [--after=N] [--tol=X] [--sub=K] [--top=K]

比（σ_2/σ_1）に潰すと2つの値の個別変化が見えなくなる。ここでは σ を降順で並べた上位
top 個（既定 top=4）を「そのまま」各列に記録する（σ_1, σ_2, σ_3, σ_4 …）。何が正しい量かは
未知であり、それを見出すことが目的なので、加工せず生値を保存する。

run_spontaneous_splitting_largeN_v1.py の run() ループを厳密に複製し、休眠フラクション
f(τ)=1-h_plane1/h_total を毎ステップ記録すると同時に、σスペクトルを K ステップごとに記録する。
σスペクトルは読み取り専用 sigma_spectrum() で計算し Z・rng・wp を一切触らないため、軌道
（f も含む）は原本とビット単位で同一になる。

【内蔵テスト】出力の f 列が原本 fcurve と bit-exact 一致すれば、軌道が同一＝σ値が正しい
軌道上の値であることの証明になる。

出力: sigma_spectrum_result_v1/sigmaspec_N{n}_delta{delta}_seed{seed}.csv
      列: tau, f, sigma_1, sigma_2, …, sigma_{top}, n_active
      （n_active = その時刻の正σの本数＝活性平面数）
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

RESULT_DIR = os.path.join(BASE_DIR, "sigma_spectrum_result_v1")


def sigma_topk_from_spectrum(sys_lr, top):
    """現在の生成子の σ を降順に並べた上位 top 個と正σの本数。読み取り専用
    （Z・rng・wp を触らない）。本数が足りなければ 0.0 で右詰めパディング。"""
    sig = sys_lr.sigma_spectrum()  # 降順の正σ
    n_active = int(len(sig))
    vals = [float(sig[i]) if i < n_active else 0.0 for i in range(top)]
    return vals, n_active


def run(n, delta, seed, cap, after, tol, sub, top):
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
    spec_hist = []   # (tau, [σ_1..σ_top], n_active) をサブサンプルで記録
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
        # --- 追加: σスペクトル記録（読み取り専用・軌道に影響しない）---
        if t % sub == 0:
            vals, n_active = sigma_topk_from_spectrum(sys_lr, top)
            spec_hist.append((t, vals, n_active))
        # --- 原本と同一の更新（wp・Z を進める）---
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    t_run = time.time() - t0

    # tau→(σ列, n_active) の辞書化（f 全ステップに対し、サンプル点のみ値・他は空）
    spec_map = {int(tau): (vals, na) for tau, vals, na in spec_hist}

    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = f"N{n:05d}_delta{delta:.0e}_seed{seed}"
    out_csv = os.path.join(RESULT_DIR, f"sigmaspec_{tag}.csv")
    sigma_cols = [f"sigma_{i+1}" for i in range(top)]
    with open(out_csv, "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["tau", "f"] + sigma_cols + ["n_active"])
        for t, f in enumerate(f_hist):
            if t in spec_map:
                vals, na = spec_map[t]
                wtr.writerow([t, f] + vals + [na])
            else:
                wtr.writerow([t, f] + [""] * top + [""])
    last_vals, last_na = (spec_hist[-1][1], spec_hist[-1][2]) if spec_hist else (None, None)
    summary = {
        "n": n, "m": sys_lr.m, "delta": delta, "seed": seed, "sub": sub, "top": top,
        "parent_residual": residual,
        "f_initial": float(f_hist[0]), "f_final": float(f_hist[-1]),
        "crossing_tau": crossed_at, "steps_run": len(f_hist) - 1,
        "n_sigma_samples": len(spec_hist),
        "sigma_final": last_vals, "n_active_final": last_na,
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
    top = int(flags.get("--top", 4))
    run(n, delta, seed, cap, after, tol, sub, top)


if __name__ == "__main__":
    main()
