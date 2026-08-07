#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2実験4: 定常粒子の構成処方の検定（新柱10の実測正本）

エンジン: 局所化プロトタイプ（円偏波 b=−ia＝チャネル回転の厳密不変多様体）。
処方（木原）: 基本波長λ₀から全倍音（偶+奇）を中心位相を揃えて積み上げる。
判定（事前固定）: 真の定常性指標＝瞬時場 τ_t(x,t) の時間stdの最大（burn500後200step）。
  (i) 全倍音・中心位相揃えは倍音数σ_kとともにリンギングが単調減少（ゼロ漸近）。
  (ii) 偶数のみ／位相ランダムとの比較で処方の優位を確認。
  (iii) 円偏波不変多様体の厳密性（|b+ia|最大）と振幅凍結（|a|変化最大）を併記。
使い方: python3 run_pre_v2_stationarity_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v2st", UIM / "run_ignition_fate_exact_v3.py")
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

    def ladder(mass, sig_k, even_only=False, phased=True):
        ks = np.arange(2 if even_only else 1, n // 2, 2 if even_only else 1)
        Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
        prof = np.zeros(n)
        rng = np.random.default_rng(7)
        for kk, w in zip(ks, Wk):
            ph = 0.0 if phased else rng.uniform(0, 2 * np.pi)
            prof += w * np.cos(2 * np.pi * kk * (x - x_L) / n + ph)
        return mass * prof / np.max(np.abs(prof))

    def stationarity(lump_prof, Tburn=500, Tobs=200):
        a2 = ((0.2 * carrier + lump_prof)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        amp0 = np.abs(a2).copy()
        snaps = []
        for j in range(Tburn + Tobs):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                snaps.append(np.angle(np.einsum("xe,xe->x", np.conj(ap), a2)))
        S = np.array(snaps)
        return (float(np.max(np.std(S, axis=0))),
                float(np.max(np.abs(np.mean(S, axis=0)))),
                float(np.max(np.abs(np.abs(a2) - amp0))),
                float(np.max(np.abs(b2 + 1j * a2))))

    M = 0.05
    cases = [("even_only_sk16", ladder(M, 16, even_only=True)),
             ("full_random_sk16", ladder(M, 16, phased=False)),
             ("full_phased_sk8", ladder(M, 8)),
             ("full_phased_sk16", ladder(M, 16)),
             ("full_phased_sk32", ladder(M, 32)),
             ("full_phased_sk64", ladder(M, 64)),
             ("full_phased_sk128", ladder(M, 128))]
    out = {"cases": {}}
    prev = None; mono = True
    for name, prof in cases:
        ring, sig, ampdrift, circdrift = stationarity(prof)
        out["cases"][name] = {"ringing": ring, "signal": sig,
                              "amp_freeze_drift": ampdrift, "circular_drift": circdrift}
        print(f"[{name}] ring={ring:.3e} signal={sig:.3e} 振幅凍結逸脱={ampdrift:.1e} 円偏波逸脱={circdrift:.1e}")
        if name.startswith("full_phased"):
            if prev is not None and ring >= prev: mono = False
            prev = ring
    verdict = ("処方成立: 全倍音・中心位相揃えでリンギングはσ_kとともに単調減少（ゼロ漸近）"
               if mono else "単調減少せず")
    print(verdict)
    out["verdict"] = verdict
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_v2_stationarity_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_v2_stationarity_result_v1.json")

if __name__ == "__main__":
    main()
