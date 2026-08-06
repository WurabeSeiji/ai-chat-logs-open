#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P8c: シート補正付きδk抽出——質量によるゲージ間隔シフトの定量

P8bの発見（事前記録）:
(A) π汚染の正体＝**大域Z₂シート反転**。(MASS,T)依存で場全体が−1倍される
    （時計二重被覆のシート）。0.2856 rad量子化はシート混合の窓平均だった疑い。
(B) 反転なし走行では、キャリアoffが距離に比例成長・⟨off⟩≈0（符号対称）
    → off ≈ δk·(x−x_L) 型＝キャリア波数シフト＝**ゲージ間隔の一様変化**。
    MASS 0.05→0.4 で約8倍＝質量に線形。

P8cの測定:
(1) シート補正: 遠方の中央値位相で大域符号を判定し C_m を ±1 補正。
(2) off(x) を符号付き距離 s=x−x_L（周期折返し）に線形回帰 → δk。
(3) スケーリング判定（事前固定）: δk ∝ MASS^p の指数 p、δk の T 依存。
    ゲージ歪みが静的場なら T 非依存、蓄積量なら T 比例。
使い方: python3 run_pre_p8c_sheet_corrected_dk_v1.py
"""
from __future__ import annotations
import importlib.util, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p8c", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

K_SEA = 3

def run_field(sp, n, ne, x, mass, x_L, T):
    sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    dx = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    lump = mass * np.exp(-0.5 * (dx / 3.0) ** 2)[:, None] * np.ones((1, ne))
    a = (sea + lump).reshape(-1).astype(complex)
    b = (sea + lump).reshape(-1).astype(complex)
    for _ in range(T):
        a, b, _ = ex.collision_step_exact(a, b, sp)
    return a.reshape(n, ne).mean(axis=1)

def carrier_project(A, n):
    F = np.fft.fft(A)
    mask = np.zeros(n, dtype=bool)
    for k in (K_SEA - 1, K_SEA, K_SEA + 1):
        mask[k % n] = True
    return np.fft.ifft(np.where(mask, F, 0.0))

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n)
    x_L = n // 2
    # 符号付き距離（周期折返し）
    s = ((x - x_L + n // 2) % n) - n // 2
    dxv = np.abs(s)
    far = (dxv > 15) & (dxv < n // 2 - 5)

    for T in (100, 200, 400):
        A_ref = run_field(sp, n, ne, x, 0.0, x_L, T)
        C_ref = carrier_project(A_ref, n)
        print(f"\n=== T={T} ===")
        for MASS in (0.05, 0.1, 0.2, 0.4):
            C_m = carrier_project(run_field(sp, n, ne, x, MASS, x_L, T), n)
            off0 = np.angle(C_m * np.conj(C_ref))
            flipped = np.median(np.abs(off0[far])) > np.pi / 2
            if flipped:
                C_m = -C_m
            off = np.angle(C_m * np.conj(C_ref))
            # 線形回帰 off = δk·s + c （遠方のみ）
            dk, c = np.polyfit(s[far], off[far], 1)
            resid = np.std(off[far] - (dk * s[far] + c))
            print(f"[MASS={MASS}] シート={'−1(反転)' if flipped else '+1'}  "
                  f"δk={dk:+.3e} rad/cell  切片={c:+.2e}  残差std={resid:.2e}")
    print(f"\n完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
