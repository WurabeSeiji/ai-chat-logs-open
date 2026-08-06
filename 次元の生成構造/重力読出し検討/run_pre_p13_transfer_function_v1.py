#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P13: 伝達関数検定——重力場は凍結源密度の普遍的線形汎関数か（ポアソン法則）

構造事実（事前記録）: 本力学の回転は各(x,η)で |a|²+|b|² を厳密保存
→ 源のエネルギー密度 Δρ(x)=Σ_η(|a|²+|b|²)−対照 は**永久凍結**（B段で検証）。
設計: 平衡平均場 v(x)=⟨δτ_t(x)⟩ について、線形応答なら v̂(k)=T(k)·ρ̂(k)。
**普遍性の判定＝伝達関数 T(k) が全(M,σ)走行で単一関数に潰れるか**（CV<10%）。
潰れれば結合定数と力の法則（ニュートン型なら T(k)∝−1/k²）を同時に読める。
潰れない場合の帰属: 非線形応答か実装不適（枠組みは未判定）。

段階: A) 定常性: 窓[800,1000] vs [1800,2000] の平均場ドリフト
      B) Δρ 凍結検証（T=2000後 vs T=0）
      C) 12走行 (M∈{0.05,0.08,0.12}×σ∈{1,3,5,8}) の T(k) 崩壊検定
         （帯域: |ρ̂(k)|>1e-2·max のkのみ・k≤16）
使い方: python3 run_pre_p13_transfer_function_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p13", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

K_SEA = 3

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    L = np.exp(-((np.abs(k) / 3.0) ** 4))
    Wf = ((k % 2) == 0).astype(float) * (1.0 - L); Wb = 1.0 - Wf

    def step(a2, b2):
        Fa = np.fft.fft(a2, axis=0); Fb = np.fft.fft(b2, axis=0)
        f = (np.sum(np.abs(np.fft.ifft(Fa * Wf[:, None], axis=0)) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * Wf[:, None], axis=0)) ** 2, axis=1))
        bo = (np.sum(np.abs(np.fft.ifft(Fa * Wb[:, None], axis=0)) ** 2, axis=1)
              + np.sum(np.abs(np.fft.ifft(Fb * Wb[:, None], axis=0)) ** 2, axis=1))
        th = np.arctan2(np.sqrt(f), np.sqrt(bo + 1e-300))
        c, s_ = np.cos(th)[:, None], np.sin(th)[:, None]
        a2, b2 = c * a2 - s_ * b2, s_ * a2 + c * b2
        phi = 2.0 * (np.sin(th) ** 2)[:, None] * np.imag(np.conj(b2) * a2)
        cp, sp_ = np.cos(phi), np.sin(phi)
        return cp * a2 - sp_ * b2, sp_ * a2 + cp * b2

    def density(a2, b2):
        return np.sum(np.abs(a2) ** 2 + np.abs(b2) ** 2, axis=1)

    def run_full(mass, sig, Tburn=1000, Tavg=1000):
        sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
        lump = mass * np.exp(-0.5 * (dxv / sig) ** 2)[:, None] * np.ones((1, ne))
        a2 = (sea + lump).astype(complex)
        b2 = (sea * np.exp(1j * np.pi / 4)).astype(complex)
        rho0 = density(a2, b2)
        acc1 = np.zeros(n); acc2 = np.zeros(n); h = Tavg // 2
        for j in range(Tburn + Tavg):
            ap = a2.copy()
            a2, b2 = step(a2, b2)
            if j >= Tburn:
                tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
                if j < Tburn + h: acc1 += tt
                else: acc2 += tt
        rho_end = density(a2, b2)
        return acc1 / h, acc2 / h, rho0, float(np.max(np.abs(rho_end - rho0)))

    # 対照
    v1c, v2c, rho_c, frz_c = run_full(0.0, 3.0)
    print(f"[対照] 平均場 std={np.std((v1c+v2c)/2):.1e}  ρ凍結逸脱={frz_c:.2e}")

    rows = []
    print(f"\n{'M':>5} {'σ':>2} | {'定常ドリフト':>12} {'ρ凍結':>9}")
    for M in (0.05, 0.08, 0.12):
        for S in (1.0, 3.0, 5.0, 8.0):
            v1, v2, rho0, frz = run_full(M, S)
            v = 0.5 * (v1 + v2) - 0.5 * (v1c + v2c)
            drift = float(np.max(np.abs(v1 - v2)) / (np.max(np.abs(v)) + 1e-300))
            drho = rho0 - rho_c
            rows.append({"M": M, "S": S, "v": v, "drho": drho})
            print(f"{M:>5} {S:>2.0f} | {drift:>12.3f} {frz:>9.2e}")

    # C) 伝達関数
    print("\n== 伝達関数 T(k)=v̂(k)/ρ̂(k)（実部・対称） ==")
    ks = np.arange(1, 17)
    Tmat = np.full((len(rows), len(ks)), np.nan)
    for i, r in enumerate(rows):
        vh = np.fft.fft(r["v"]); rh = np.fft.fft(r["drho"])
        thr = 1e-2 * np.max(np.abs(rh))
        for jj, kk in enumerate(ks):
            if abs(rh[kk]) > thr:
                Tmat[i, jj] = (vh[kk] / rh[kk]).real
    hdr = "  M   σ | " + " ".join(f"k={kk:<2}" for kk in ks[:8])
    print(hdr)
    for i, r in enumerate(rows):
        cells = " ".join(f"{Tmat[i,jj]:+.3f}" if np.isfinite(Tmat[i,jj]) else "  ---" for jj in range(8))
        print(f"{r['M']:>4} {r['S']:>3.0f} | {cells}")
    print("\nkごとの崩壊度 CV=std/|mean|（有効走行数）:")
    for jj, kk in enumerate(ks):
        col = Tmat[:, jj]; col = col[np.isfinite(col)]
        if len(col) >= 6:
            print(f"  k={kk:>2}: T平均={np.mean(col):+.4f} CV={np.std(col)/abs(np.mean(col)):.1%} (n={len(col)})")
    np.save(HERE / "result_pre_p13_Tmat_v1.npy", Tmat)
    print(f"\n完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
