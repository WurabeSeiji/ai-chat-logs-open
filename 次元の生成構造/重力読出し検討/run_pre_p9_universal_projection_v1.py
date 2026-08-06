#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P9: 万能空間射影——関係量からのゲージ創生と局所歪みの読出し

設計（木原・事前記録）: 波の万能演算の空間版。空間を創生している射影関数に
万能性を持たせ、ゲージの局所歪みを関係量から作り込む。空間は自発的に
生じているため、無名性を破らず（力学への注入ゼロで）演算を持ち込める。

定義: 各セル x のηファイバー ψ_x ∈ C^16（aチャネル）について
  空間ゲージ目盛 τ_x(x) = arg⟨ψ_x, ψ_{x+1}⟩   （隣接ファイバー相対位相）
  時間ゲージ目盛 τ_t(x) = arg⟨ψ_x(T), ψ_x(T+1)⟩ （局所時計レート・純力学量）
  形状成分     σ_x(x) = arccos(|⟨ψ_x,ψ_{x+1}⟩|/‖ψ_x‖‖ψ_{x+1}‖)（FS距離）
全て純関係量・大域位相不変・種非依存。

判定（事前固定）: 分析ノート P9 節のとおり（対照一様性／δτ_t=時間ゲージ歪みの
初検出判定／δτ_x のT=0分離／γ_model＝空間/時間歪み比／累積Shapiro）。
使い方: python3 run_pre_p9_universal_projection_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p9", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

K_SEA = 3

def gauges(a_prev, a_now, n, ne):
    """状態(1ステップ前後)からゲージ目盛3種を読む。純関係量のみ。"""
    P0 = a_prev.reshape(n, ne)
    P1 = a_now.reshape(n, ne)
    ip_x = np.einsum("xe,xe->x", np.conj(P1), np.roll(P1, -1, axis=0))
    nr = np.linalg.norm(P1, axis=1)
    tau_x = np.angle(ip_x)
    sig_x = np.arccos(np.clip(np.abs(ip_x) / (nr * np.roll(nr, -1) + 1e-300), 0, 1))
    ip_t = np.einsum("xe,xe->x", np.conj(P0), P1)
    tau_t = np.angle(ip_t)
    return tau_x, sig_x, tau_t

def run(sp, n, ne, x, dxv, mass, T):
    sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    lump = mass * np.exp(-0.5 * (dxv / 3.0) ** 2)[:, None] * np.ones((1, ne))
    a = (sea + lump).reshape(-1).astype(complex)
    b = a.copy()
    a_prev = a.copy()
    for j in range(T):
        a_prev = a.copy()
        a, b, _ = ex.collision_step_exact(a, b, sp)
    return gauges(a_prev, a, n, ne)

def prof(v, dxv, dds=(0, 6, 12, 20, 40, 100)):
    out = {}
    for dd in dds:
        sel = (dxv >= max(dd - 2, 0)) & (dxv <= dd + 2)
        out[dd] = float(np.mean(v[sel]))
    return out

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n)
    x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    far = dxv > 60
    results = {}

    for T in (50, 200, 400):
        tx0, sx0, tt0 = run(sp, n, ne, x, dxv, 0.0, T)
        print(f"\n=== T={T} ===")
        print(f"[対照] τ_x: 平均={np.mean(tx0):+.6f} std={np.std(tx0):.2e} | "
              f"τ_t: 平均={np.mean(tt0):+.6f} std={np.std(tt0):.2e} | "
              f"σ_x: 平均={np.mean(sx0):.3e}")
        for MASS in (0.05, 0.1, 0.2, 0.4):
            tx, sx, tt = run(sp, n, ne, x, dxv, MASS, T)
            dtx, dtt, dsx = tx - tx0, tt - tt0, sx - sx0
            key = f"T{T}_M{MASS}"
            results[key] = {
                "d_tau_x": prof(dtx, dxv), "d_tau_t": prof(dtt, dxv),
                "d_sigma_x": prof(dsx, dxv),
                "far_std_tau_t": float(np.std(dtt[far])),
                "shapiro_sum": float(np.sum(dtx)),
            }
            r = results[key]
            print(f"[M={MASS}] δτ_t(時間ゲージ): 塊上={r['d_tau_t'][0]:+.3e} "
                  f"dx6={r['d_tau_t'][6]:+.3e} dx12={r['d_tau_t'][12]:+.3e} "
                  f"dx20={r['d_tau_t'][20]:+.3e} 遠方std={r['far_std_tau_t']:.1e}")
            print(f"         δτ_x(空間ゲージ): 塊上={r['d_tau_x'][0]:+.3e} "
                  f"dx6={r['d_tau_x'][6]:+.3e} dx12={r['d_tau_x'][12]:+.3e} "
                  f"環一周Σδτ_x={r['shapiro_sum']:+.3e}")
            print(f"         δσ_x(形状): 塊上={r['d_sigma_x'][0]:+.3e} dx12={r['d_sigma_x'][12]:+.3e}")
            if abs(r['d_tau_t'][0]) > 1e-12:
                print(f"         γ_model(塊上 空間/時間)={r['d_tau_x'][0]/r['d_tau_t'][0]:+.4f}")

    # T=0分離（空間目盛のみ・時間目盛はT=0で定義されない=純力学量）
    print("\n=== T=0（初期条件成分・空間目盛のみ） ===")
    sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    for MASS in (0.05, 0.4):
        lump = MASS * np.exp(-0.5 * (dxv / 3.0) ** 2)[:, None] * np.ones((1, ne))
        a_m = (sea + lump).reshape(-1).astype(complex)
        a_r = sea.reshape(-1).astype(complex)
        txm, _, _ = gauges(a_m, a_m, n, ne)
        txr, _, _ = gauges(a_r, a_r, n, ne)
        d0 = txm - txr
        print(f"[M={MASS}] δτ_x(T=0): 塊上={np.mean(d0[dxv <= 2]):+.3e} "
              f"dx12={np.mean(d0[(dxv >= 10) & (dxv <= 14)]):+.3e}")

    out = HERE / "result_pre_p9_universal_projection_v1.json"
    out.write_text(json.dumps(results, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"\n保存: {out.name}  完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
