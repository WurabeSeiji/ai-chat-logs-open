#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P14: 帯純度倍音パケット源——局在定常粒子によるリンギング根絶と構成則検定

設計（木原指摘の実装・事前記録）: 局在した定常粒子は倍音を増やせば簡単に作れる。
源をフェルミオン帯（偶|k|≥4）の倍音だけで組む:
  lump(x) = M·Σ_{k偶≥4} W(k)·cos(2πk(x−x_L)/n),  W=exp(−k²/2σ_k²)（ピーク規格化）
帯純度により、ガウス塊で起きた帯間スロッシング（奇k・低k漏れ⇄フェルミオン帯の
唸り）が構造的に消えるはず。偶kのみ→n/2周期＝対蹠に厳密な複製（既知の二源幾何）。

判定（事前固定）:
 (i) リンギング指標（窓[800,1000] vs [1800,2000] の平均場ドリフト比）が
     ガウス塊（~1.0）から大幅低下するか。
 (ii) 低下したら構成則 |v|=K·ρ_even^γ を M×σ_k 掃引で再検定（γ自由+固定）。
     K の CV<10% で普遍性成立。
使い方: python3 run_pre_p14_bandpure_source_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p14", UIM / "run_ignition_fate_exact_v3.py")
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
    carrier = np.exp(2j * np.pi * K_SEA * x / n)

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

    def bandpure_lump(mass, sig_k):
        ks = np.arange(4, n // 2, 2)
        Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
        prof = np.zeros(n)
        for kk, w in zip(ks, Wk):
            prof += w * np.cos(2 * np.pi * kk * (x - x_L) / n)
        prof /= np.max(np.abs(prof))
        return mass * prof

    def run_mean(lump_prof, Tburn=1000, Tavg=1000):
        sea = 0.2 * carrier[:, None] * np.ones((1, ne))
        a2 = (sea + lump_prof[:, None] * np.ones((1, ne))).astype(complex)
        b2 = (sea * np.exp(1j * np.pi / 4)).astype(complex)
        rho0 = np.sum(np.abs(a2) ** 2 + np.abs(b2) ** 2, axis=1)
        acc1 = np.zeros(n); acc2 = np.zeros(n); h = Tavg // 2
        for j in range(Tburn + Tavg):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
                if j < Tburn + h: acc1 += tt
                else: acc2 += tt
        return acc1 / h, acc2 / h, rho0

    # 対照
    z = np.zeros(n)
    v1c, v2c, rho_c = run_mean(z)
    vc = 0.5 * (v1c + v2c)
    print(f"[対照] std={np.std(vc):.1e}")

    # (i) リンギング比較: ガウス塊 vs 帯純度パケット（同ピークM=0.05・類似幅）
    g = 0.05 * np.exp(-0.5 * (dxv / 3.0) ** 2)
    for name, prof in (("ガウス塊(σ=3)", g), ("帯純度パケット(σ_k=32)", bandpure_lump(0.05, 32.0))):
        v1, v2, rho0 = run_mean(prof)
        v = 0.5 * (v1 + v2) - vc
        drift = float(np.max(np.abs(v1 - v2)) / (np.max(np.abs(v)) + 1e-300))
        print(f"[{name}] リンギング指標={drift:.3f}  max|v|={np.max(np.abs(v)):.3e}")

    # (ii) 掃引: M×σ_k
    print("\n== 構成則検定（帯純度源・M×σ_k） ==")
    X, Y, labels = [], [], []
    Klist = {}
    for M in (0.02, 0.05, 0.1):
        for SK in (8.0, 16.0, 32.0):
            v1, v2, rho0 = run_mean(bandpure_lump(M, SK))
            v = 0.5 * (v1 + v2) - vc
            drift = float(np.max(np.abs(v1 - v2)) / (np.max(np.abs(v)) + 1e-300))
            drho = rho0 - rho_c
            rho_even = 0.5 * (drho + drho[(x + n // 2) % n])
            sel = (rho_even > 1e-6) & (np.abs(v) > 1e-12)
            Kfix = float(np.exp(np.mean(np.log(np.abs(v[sel])) - 0.5 * np.log(rho_even[sel]))))
            Klist[(M, SK)] = Kfix
            sgn = float(np.mean(np.sign(v[sel])))
            print(f"[M={M} σ_k={SK:.0f}] リンギング={drift:.3f} K(γ=½)={Kfix:.4f} 符号={sgn:+.2f} 点数={int(np.sum(sel))}")
            for ri, vi in zip(rho_even[sel], np.abs(v[sel])):
                X.append([1.0, np.log(ri)]); Y.append(np.log(vi)); labels.append((M, SK))
    X = np.array(X); Y = np.array(Y)
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ coef
    R2 = 1 - np.sum((Y - pred) ** 2) / np.sum((Y - np.mean(Y)) ** 2)
    Ks = np.array(list(Klist.values()))
    print(f"\n自由フィット: γ={coef[1]:+.4f} K={np.exp(coef[0]):.4f} R²={R2:.4f}")
    print(f"γ=1/2固定: K平均={np.mean(Ks):.4f} CV={np.std(Ks)/np.mean(Ks):.1%}")
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
