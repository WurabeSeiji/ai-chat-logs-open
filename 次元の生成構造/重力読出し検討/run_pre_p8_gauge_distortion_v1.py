#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P8: ゲージ間隔の歪み測定——重力＝ゲージが等間隔でなくなること（木原の枠組み）

原理（事前記録）: 落ちるのではない。空間・時間のゲージ間隔が場所依存になる
＝それが加速度に見える。測るのは運動でなく、**海の位相刻み（1セルあたりの
位相の進み）の空間プロファイル** Δφ(x) = arg(ψ(x+1)ψ̄(x))。
質量塊のまわりでこれが場所依存になれば、ゲージ歪み＝重力場の読出し。
純関係量（隣接セル位相差）・時計不要・運動不要・注入なし（塊と海は初期条件）。

系: 二体(χ,η)格子。海=弱い低k平面波（a,b両チャネル・m=0的）。
塊=重い狭ガウス束（χ位置 x_L に局在）。対照=塊なし同一走行。
判定（記述的）: 差分プロファイル δ(x)=Δφ_塊(x)−Δφ_対照(x) が
(i) 塊の台の外でゼロ→接触のみ（長距離場なし） (ii) 台の外で距離依存の
裾を持つ→**長距離ゲージ歪み＝重力場様作用の初検出**。
使い方: python3 run_pre_p8_gauge_distortion_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p8", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    x = np.arange(n)

    # 海: 弱い低k平面波（両チャネル・η一様=m0的）
    K_SEA = 3
    sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    # 塊: 重い狭ガウス（χ位置 x_L・η一様）
    x_L = n // 2
    for MASS in (0.0, 4.0, 8.0):
        dx = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
        lump = MASS * np.exp(-0.5 * (dx / 3.0) ** 2)[:, None] * np.ones((1, ne))
        a = (sea + lump).reshape(-1).astype(complex)
        b = (sea + lump).reshape(-1).astype(complex)
        for j in range(400):
            a, b, _ = ex.collision_step_exact(a, b, sp)
        A2 = a.reshape(shape).mean(axis=1)   # η平均の場
        dphi = np.angle(A2[np.r_[1:n, 0]] * np.conj(A2))   # 隣接セル位相刻み
        if MASS == 0.0:
            dphi_ref = dphi.copy()
            print(f"[対照 MASS=0] 位相刻み: 平均={np.mean(dphi):.5f} 標準偏差={np.std(dphi):.2e}")
            continue
        delta = dphi - dphi_ref
        # 台の内外
        support = dx <= 9
        out_near = (dx > 9) & (dx <= 30)
        out_far = dx > 60
        print(f"[MASS={MASS}] δ(位相刻み): 台内 max|δ|={np.max(np.abs(delta[support])):.4f}  "
              f"近傍(9<dx≤30) max|δ|={np.max(np.abs(delta[out_near])):.2e}  "
              f"遠方(dx>60) max|δ|={np.max(np.abs(delta[out_far])):.2e}")
        # 距離プロファイル（近傍の減衰形）
        for dd in (12, 16, 24, 40, 80):
            sel = (dx >= dd - 2) & (dx <= dd + 2)
            print(f"    dx≈{dd}: ⟨|δ|⟩={np.mean(np.abs(delta[sel])):.3e}")
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
