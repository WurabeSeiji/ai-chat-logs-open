#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G4正本: k̇=∇ω 運動量則と γ(t)=t/L のパラメータフリー導出

判定（事前固定）:
(i) γ(t)=⟨|δτ_x|⟩/⟨|δτ_t|⟩|_band が時間線形（R²>0.9）——γは定数でない
    （P21b解釈の訂正を正本化）。
(ii) 核幾何予言 γ(t)=t/L（L=⟨|∇ω|⟩/⟨|ω|⟩の逆数・解析ω場から・力学なし）が
    実測とCV<20%——空間ゲージの成長はδτ_x(t)=t·∇ω＝海の運動量が時計勾配
    レートで増える（k̇=∇ω）。
付録出力: ω(x)実測プロファイル vs 解析予言（柱G2の図用）。
使い方: python3 run_g4_momentum_law_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g4", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base
K_SEA = 3

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n)
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

    def ladder_real(fsrc, amp, center, sig_k=32.0):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    c = 100; f_src, amp = 0.6, 0.07
    a2 = ((0.2 * carrier + ladder_real(f_src, amp, c))[:, None] * np.ones((1, ne))).astype(complex)
    b2 = -1j * a2
    dxv = np.minimum(np.abs(x - c), n - np.abs(x - c))
    band = (dxv > 20) & (dxv < 60)
    tau_x0 = 2 * np.pi * K_SEA / n
    gam = {}; om_meas = None
    for j in range(1, 901):
        ap = a2.copy(); a2, b2 = step(a2, b2)
        if j in (200, 400, 600, 800):
            tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
            tx = np.angle(np.einsum("xe,xe->x", np.conj(a2), np.roll(a2, -1, axis=0)))
            gam[j] = float(np.mean(np.abs(tx[band] - tau_x0)) / np.mean(np.abs(tt[band])))
            if j == 200:
                om_meas = tt.copy()
    ts = np.array(sorted(gam)); gs = np.array([gam[t] for t in ts])
    lin = np.polyfit(ts, gs, 1)
    R2 = 1 - np.sum((gs - np.polyval(lin, ts)) ** 2) / np.sum((gs - gs.mean()) ** 2)
    # 解析ω場から L
    a1 = 0.2 * carrier + ladder_real(f_src, amp, c)
    F = np.fft.fft(a1)
    f_loc = 2 * np.abs(np.fft.ifft(F * Wf)) ** 2
    b_loc = 2 * np.abs(np.fft.ifft(F * Wb)) ** 2
    th = np.arctan2(np.sqrt(f_loc), np.sqrt(b_loc + 1e-300))
    om_pred = th + 2 * (np.sin(th) ** 2) * np.abs(a1) ** 2
    grad = np.abs(np.diff(np.r_[om_pred, om_pred[0]]))
    invL = float(np.mean(grad[band]) / np.mean(np.abs(om_pred[band])))
    ratio = gs / (ts * invL)
    ok = R2 > 0.9 and abs(float(np.mean(ratio)) - 1) < 0.2
    print(f"γ(t)実測: " + " ".join(f"t={t}:{g:.2f}" for t, g in gam.items()))
    print(f"線形性 R²={R2:.4f}  1/L={invL:.4e}（L={1/invL:.1f}セル）")
    print(f"予言γ=t/L との比: 平均={np.mean(ratio):.3f} CV={np.std(ratio)/np.mean(ratio):.1%}")
    verdict = "柱G4成立（k̇=∇ω・γ(t)=t/L）" if ok else "要精査"
    print(verdict)
    out = {"gamma_t": {int(t): float(g) for t, g in gam.items()},
           "linearity_R2": float(R2), "invL": invL, "L_cells": 1 / invL,
           "ratio_mean": float(np.mean(ratio)), "ratio_cv": float(np.std(ratio) / np.mean(ratio)),
           "omega_profile_meas": om_meas.tolist(), "omega_profile_pred": om_pred.tolist(),
           "center": c, "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g4_momentum_law_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
