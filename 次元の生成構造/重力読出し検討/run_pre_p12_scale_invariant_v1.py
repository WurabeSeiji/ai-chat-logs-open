#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P12: スケール不変局所射影——滑らかな核で「裾は核か力学か」を判別する

設計（事前記録）: P11診断（急峻二帯マスクの固有スケールが普遍性を破る）を受け、
射影を次で置換する。
  パリティ分割: 偶k/奇k は半環並進 x→x+n/2 の厳密対称性＝スケールフリー。
  IRロールオフ: 基本波帯の除外を急峻エッジでなく L(k)=exp(−(|k|/3)⁴) で行う
    （モデル唯一の物理IRスケール。滑らか核→Dirichlet 1/x振動が消える）。
  W_f(k) = [k偶]·(1−L(k)),  W_b = 1−W_f
  f_loc(x) = Σ_η |ifft(Fa·W_f)|²+|ifft(Fb·W_f)|²（b_locも同様）
  θ(x) = atan2(√f_loc, √b_loc)。以下P10と同一ステップ。

事前予測と判別（最重要）:
  P10の1/x裾は急峻核の人工物の疑い。滑らか核では射影の直接影響は
  ガウス減衰＝短距離のはず。**それでも裾がTとともに外向きに育つなら、
  核でなく力学（伝播）が作る本物の遠距離場**。
判定（事前固定）:
 (i) 閉塞ドリフト ≤1e-12  (ii) 真空固定点（海k=3は奇→W_f=0→θ=0厳密）
 (iii) δτ_t プロファイルの T 依存: 前線が外に動くか（伝播速度）
 (iv) 裾が育つ場合のみ M×σ 普遍性回帰（包絡線・自己相似形）
使い方: python3 run_pre_p12_scale_invariant_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p12", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

K_SEA = 3

def make_weights(n):
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    L = np.exp(-((np.abs(k) / 3.0) ** 4))
    Wf = ((k % 2) == 0).astype(float) * (1.0 - L)
    return Wf, 1.0 - Wf

def step(a2, b2, Wf, Wb):
    Fa = np.fft.fft(a2, axis=0); Fb = np.fft.fft(b2, axis=0)
    f_loc = (np.sum(np.abs(np.fft.ifft(Fa * Wf[:, None], axis=0)) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * Wf[:, None], axis=0)) ** 2, axis=1))
    b_loc = (np.sum(np.abs(np.fft.ifft(Fa * Wb[:, None], axis=0)) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * Wb[:, None], axis=0)) ** 2, axis=1))
    th = np.arctan2(np.sqrt(f_loc), np.sqrt(b_loc + 1e-300))
    c, s_ = np.cos(th)[:, None], np.sin(th)[:, None]
    a2, b2 = c * a2 - s_ * b2, s_ * a2 + c * b2
    phi = 2.0 * (np.sin(th) ** 2)[:, None] * np.imag(np.conj(b2) * a2)
    cp, sp_ = np.cos(phi), np.sin(phi)
    return cp * a2 - sp_ * b2, sp_ * a2 + cp * b2

def run_tau(n, ne, x, dxv, mass, sig, T, Wf, Wb):
    sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    lump = mass * np.exp(-0.5 * (dxv / sig) ** 2)[:, None] * np.ones((1, ne))
    a2 = (sea + lump).astype(complex)
    b2 = (sea * np.exp(1j * np.pi / 4)).astype(complex)
    C0 = complex(np.sum(a2 * a2) + np.sum(b2 * b2))
    for j in range(T):
        ap = a2.copy()
        a2, b2 = step(a2, b2, Wf, Wb)
    dC = abs(complex(np.sum(a2 * a2) + np.sum(b2 * b2)) - C0)
    return np.angle(np.einsum("xe,xe->x", np.conj(ap), a2)), dC

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    Wf, Wb = make_weights(n)

    # (i)(ii) 対照と閉塞
    tt0_by_T = {}
    for T in (50, 100, 200, 400):
        tt0, dC0 = run_tau(n, ne, x, dxv, 0.0, 3.0, T, Wf, Wb)
        tt0_by_T[T] = tt0
        if T == 50:
            print(f"[対照] τ_t std={np.std(tt0):.1e} 閉塞ドリフト={dC0:.2e}")

    # (iii) 前線のT依存（M=0.1, σ=3）
    print("\n== 前線プロファイル δτ_t（M=0.1, σ=3） ==")
    for T in (50, 100, 200, 400):
        tt, dC = run_tau(n, ne, x, dxv, 0.1, 3.0, T, Wf, Wb)
        dtt = tt - tt0_by_T[T]
        prof = []
        for dd in (0, 8, 15, 25, 40, 60, 90, 130, 180, 240):
            sel = (dxv >= max(dd - 3, 0)) & (dxv <= dd + 3)
            prof.append(f"dx{dd}={np.mean(dtt[sel]):+.2e}")
        # 前線位置: |δτ_t| が閾値1e-8を超える最遠 dx
        above = np.abs(dtt) > 1e-8
        front = int(np.max(dxv[above])) if np.any(above) else 0
        print(f"T={T:>3}: 前線(|δτ|>1e-8)dx={front:>3} 閉塞={dC:.1e}\n      " + " ".join(prof))
    print(f"\n完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
