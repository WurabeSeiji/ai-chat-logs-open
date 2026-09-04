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


# ==== 相互作用スワップ: 以下4関数は run_and_plot_N3_N33_legacyparent_20260903.py から
# 逐語コピー（振幅込み Hermitian H・固定 Δτ=2π/den のユニタリ写像）。他は一切不変。 ====
def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A
def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
def one_step(z,A,den):
    H=H_of(z,A); w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)


def run(n, delta, seed, cap, after, tol=1e-8, den=None):
    if den is None:
        den = n
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + seed)
    t0 = time.time()
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=tol)
    t_parent = time.time() - t0
    progress(f"N={n} 親構成完了 残差={residual:.2e} 平面数={len(sig)} ({t_parent:.1f}s)")

    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    Z0_state = Z.copy()  # 状態保存版の追記: step0 状態
    ztz0 = abs(complex(Z @ Z))

    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    f_hist = []
    max_ztz = 0.0
    crossed_at = None
    wp = rng.normal(size=sys_lr.m)
    A_int = adjacency(n)  # 相互作用スワップ: 新力学の隣接行列（辺順序は build_edges と同一を検証済み）
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
        # 相互作用スワップ: 旧 set_theta/sigma_max_power/cayley_step の3行を one_step に置換
        Z = one_step(Z, A_int, den)
    t_run = time.time() - t0

    f_arr = np.array(f_hist)
    lo, hi = max(10.0 * f_arr[0], 1e-300), 1e-3
    mask = (f_arr > lo) & (f_arr < hi)
    idx = np.where(mask)[0]
    rate = float(np.polyfit(idx, np.log(f_arr[idx]), 1)[0]) if len(idx) >= 5 else None

    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = f"N{n:05d}_delta{delta:.0e}_seed{seed}_den{den}"
    with open(os.path.join(RESULT_DIR, f"fcurve_{tag}.csv"), "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["tau", "f"])
        for t, f in enumerate(f_hist):
            wtr.writerow([t, f])
    summary = {
        "n": n, "m": sys_lr.m, "delta": delta, "seed": seed,
        "den": den, "interaction": "one_step_exp_minus_i_2pi_over_den_H",
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
    # 状態保存版の追記: step0 と最終 step の全複素状態（力学には一切影響しない）
    np.savez_compressed(os.path.join(RESULT_DIR, f"states_{tag}.npz"),
                        Z0=Z0_state, Zfinal=Z, tau_final=np.int64(len(f_hist) - 1),
                        p=p, q=q)
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
    den = int(flags.get("--den", n))
    run(n, delta, seed, cap, after, tol=tol, den=den)


if __name__ == "__main__":
    main()
