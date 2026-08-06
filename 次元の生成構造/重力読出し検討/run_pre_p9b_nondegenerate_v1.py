#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P9b: 縮退破り版の万能空間射影——生きた力学でのゲージ歪み測定

P9の重大発見（事前記録）: a=b 初期条件は φ=2r·Im(b̄a)=0 で非弾性が永久に
発火せず、弾性回転も b∝a を実係数で保つため、全構成が実スカラー倍で静止
（τ_t≡0・δτ_x全T同一で実証）。P8系の全結果は凍結初期条件だった。
Z₂シート反転の正体＝累積実スカラーの符号。

P9bの初期条件: a = 海 + 塊、b = 海·e^{iπ/4}（相対位相で縮退を破る。
塊はaのみ）。対照 = 同構成で塊なし。測定はP9と同一（τ_x, τ_t, σ_x）。
判定: P9節の事前固定に同じ。追加判定: τ_t(対照)≠0（時計が進む＝力学が
生きている確認）を最初に確認する。
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p9b", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

K_SEA = 3

def gauges(a_prev, a_now, n, ne):
    P0 = a_prev.reshape(n, ne); P1 = a_now.reshape(n, ne)
    ip_x = np.einsum("xe,xe->x", np.conj(P1), np.roll(P1, -1, axis=0))
    nr = np.linalg.norm(P1, axis=1)
    tau_x = np.angle(ip_x)
    sig_x = np.arccos(np.clip(np.abs(ip_x) / (nr * np.roll(nr, -1) + 1e-300), 0, 1))
    tau_t = np.angle(np.einsum("xe,xe->x", np.conj(P0), P1))
    return tau_x, sig_x, tau_t

def run(sp, n, ne, x, dxv, mass, T):
    sea_a = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    sea_b = sea_a * np.exp(1j * np.pi / 4)
    lump = mass * np.exp(-0.5 * (dxv / 3.0) ** 2)[:, None] * np.ones((1, ne))
    a = (sea_a + lump).reshape(-1).astype(complex)
    b = sea_b.reshape(-1).astype(complex)
    for j in range(T):
        a_prev = a.copy()
        a, b, _ = ex.collision_step_exact(a, b, sp)
    return gauges(a_prev, a, n, ne)

def prof(v, dxv, dds=(0, 6, 12, 20, 40, 100)):
    return {dd: float(np.mean(v[(dxv >= max(dd - 2, 0)) & (dxv <= dd + 2)])) for dd in dds}

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    far = dxv > 60
    results = {}
    for T in (50, 200, 400):
        tx0, sx0, tt0 = run(sp, n, ne, x, dxv, 0.0, T)
        print(f"\n=== T={T} ===")
        print(f"[対照] τ_x: 平均={np.mean(tx0):+.6f} std={np.std(tx0):.2e} | "
              f"τ_t: 平均={np.mean(tt0):+.6f} std={np.std(tt0):.2e}  ←時計が進むか")
        for MASS in (0.05, 0.1, 0.2, 0.4):
            tx, sx, tt = run(sp, n, ne, x, dxv, MASS, T)
            dtx, dtt, dsx = tx - tx0, tt - tt0, sx - sx0
            key = f"T{T}_M{MASS}"
            results[key] = {"d_tau_x": prof(dtx, dxv), "d_tau_t": prof(dtt, dxv),
                            "d_sigma_x": prof(dsx, dxv),
                            "far_std_tau_t": float(np.std(dtt[far])),
                            "far_std_tau_x": float(np.std(dtx[far])),
                            "shapiro_sum": float(np.sum(dtx))}
            r = results[key]
            print(f"[M={MASS}] δτ_t: 塊上={r['d_tau_t'][0]:+.3e} dx6={r['d_tau_t'][6]:+.3e} "
                  f"dx12={r['d_tau_t'][12]:+.3e} dx20={r['d_tau_t'][20]:+.3e} "
                  f"dx40={r['d_tau_t'][40]:+.3e} 遠方std={r['far_std_tau_t']:.1e}")
            print(f"         δτ_x: 塊上={r['d_tau_x'][0]:+.3e} dx6={r['d_tau_x'][6]:+.3e} "
                  f"dx12={r['d_tau_x'][12]:+.3e} dx20={r['d_tau_x'][20]:+.3e} "
                  f"遠方std={r['far_std_tau_x']:.1e} Σ={r['shapiro_sum']:+.3e}")
            if abs(r['d_tau_t'][0]) > 1e-12:
                print(f"         γ_model(塊上)={r['d_tau_x'][0]/r['d_tau_t'][0]:+.4f}")
    out = HERE / "result_pre_p9b_nondegenerate_v1.json"
    out.write_text(json.dumps(results, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"\n保存: {out.name}  完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
