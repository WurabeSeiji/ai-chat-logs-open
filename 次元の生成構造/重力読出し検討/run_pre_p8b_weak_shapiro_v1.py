#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P8b: 弱結合Shapiro測定——海キャリアの蓄積位相オフセット（時間ゲージ歪みの積分）

背景（事前記録）: P8で位相刻み（微分）の質量依存歪みを検出したが接触支配。
蓄積オフセット（積分＝Shapiro型）は強結合域（塊4-8 vs 海0.2）で場がドメイン
再構成し |off|≈π の汚染（平均が0.2856 radの整数倍に量子化）で読めなかった。

P8bの設計:
(1) 弱結合: 塊振幅 MASS ∈ {0.05, 0.1, 0.2, 0.4}（海0.2と同程度以下）。
(2) 汚染対策: 生のη平均場でなく**海キャリア帯域 k∈[K_SEA−1, K_SEA+1] に
    射影した場** A_c(x) の位相で読む。散乱生成物（他k）を帯域分離で除去。
(3) 診断: 生場のオフセット分布のπ近傍分率も併記（P8汚染の正体確認）。
判定（事前固定）:
  (i) 質量スケーリング: off が MASS（または MASS²）に比例 → 場的読出し。
  (ii) 時間蓄積: off が衝突数 T に対し増加 → Shapiro型（積分量）。
  (iii) どちらも欠ければ「蓄積時間ゲージ信号なし（この次数）」と記録。
使い方: python3 run_pre_p8b_weak_shapiro_v1.py
"""
from __future__ import annotations
import importlib.util, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p8b", UIM / "run_ignition_fate_exact_v3.py")
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
    return a.reshape(n, ne).mean(axis=1), dx

def carrier_project(A, n):
    """海キャリア帯域 k∈{K_SEA-1, K_SEA, K_SEA+1} のみ残す帯域射影。"""
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

    for T in (100, 400):
        A_ref, _ = run_field(sp, n, ne, x, 0.0, x_L, T)
        C_ref = carrier_project(A_ref, n)
        print(f"\n=== T={T} 衝突 ===")
        for MASS in (0.05, 0.1, 0.2, 0.4):
            A_m, dxv = run_field(sp, n, ne, x, MASS, x_L, T)
            C_m = carrier_project(A_m, n)
            # キャリア位相オフセット（帯域分離済み）
            off_c = np.angle(C_m * np.conj(C_ref))
            # 生場オフセット（P8汚染の診断用）
            off_raw = np.angle(A_m * np.conj(A_ref))
            far = dxv > 30
            frac_pi = np.mean(np.abs(np.abs(off_raw[far]) - np.pi) < 0.3)
            print(f"[MASS={MASS}] キャリアoff: "
                  f"台内max|off|={np.max(np.abs(off_c[dxv <= 9])):.3e}  "
                  f"遠方(dx>30)⟨off⟩={np.mean(off_c[far]):+.3e} "
                  f"std={np.std(off_c[far]):.2e}  |  "
                  f"生場π分率(遠方)={frac_pi:.2f}")
            for dd in (12, 20, 40, 80):
                sel = (dxv >= dd - 2) & (dxv <= dd + 2)
                print(f"    dx≈{dd}: キャリア⟨off⟩={np.mean(off_c[sel]):+.3e} "
                      f"⟨|off|⟩={np.mean(np.abs(off_c[sel])):.3e}")
    print(f"\n完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
