#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G3正本: 質量生成の因子化と解析導出（ヒッグス節）

判定（事前固定）:
(i) 因子化: g(v)=m(f,v)/m(f,v0) の f間CV が全vで <5% なら m=y(f)·g(v) 成立。
(ii) 解析導出: 読出し式に初期構成を代入した予言 m_pred（力学なし・自由
    パラメータなし）が動力学実測15点と比CV<5%なら、質量法則は核の解析的定理。
使い方: python3 run_g3_higgs_factorization_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g3", UIM / "run_ignition_fate_exact_v3.py")
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

    def ladder_at(fsrc, center, sig_k=32.0, amp=0.05):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof * np.exp(2j * np.pi * K_SEA * center / n)

    def mass_meas(fsrc, v, c=100, T=150):
        a2 = ((v * carrier + ladder_at(fsrc, c))[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        dA = np.minimum(np.abs(x - c), n - np.abs(x - c)) <= 3
        Phi = 0.0
        for j in range(T):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
            Phi += float(np.mean(tt[dA]))
        return Phi / T

    def mass_pred(fsrc, v, c=100):
        a1 = v * carrier + ladder_at(fsrc, c)
        F = np.fft.fft(a1)
        f_loc = 2 * np.abs(np.fft.ifft(F * Wf)) ** 2
        b_loc = 2 * np.abs(np.fft.ifft(F * Wb)) ** 2
        th = np.arctan2(np.sqrt(f_loc), np.sqrt(b_loc + 1e-300))
        om = th + 2 * (np.sin(th) ** 2) * np.abs(a1) ** 2
        dA = np.minimum(np.abs(x - c), n - np.abs(x - c)) <= 3
        return float(np.mean(om[dA]))

    vs = [0.10, 0.14, 0.20, 0.28, 0.40]; fs = [0.3, 0.6, 0.9]
    M = np.array([[mass_meas(f, v) for v in vs] for f in fs])
    P = np.array([[mass_pred(f, v) for v in vs] for f in fs])
    G = M / M[:, 2:3]
    cv_fact = [float(np.std(G[:, j]) / np.mean(G[:, j])) for j in range(len(vs))]
    p_exp = float(np.polyfit(np.log(vs), np.log(G.mean(axis=0)), 1)[0])
    y = (M / M[1:2, :]).mean(axis=1)
    ratio = M / P
    cv_pred = float(ratio.std() / ratio.mean())
    ok = max(cv_fact) < 0.05 and cv_pred < 0.05
    print(f"因子化: f間CV最大={max(cv_fact):.3%}  g(v)∝v^{p_exp:.3f}  "
          f"y(f)/y(0.6)={y[0]:.4f},{y[1]:.4f},{y[2]:.4f}")
    print(f"解析導出: 実測/予言 平均={ratio.mean():.4f} CV={cv_pred:.3%}")
    verdict = "柱G3成立（因子化＋解析的定理）" if ok else "要精査"
    print(verdict)
    out = {"vs": vs, "fs": fs, "M_meas": M.tolist(), "M_pred": P.tolist(),
           "factorization_cv_max": max(cv_fact), "g_exponent": p_exp,
           "y_ratios": y.tolist(), "pred_ratio_mean": float(ratio.mean()),
           "pred_ratio_cv": cv_pred, "verdict": verdict,
           "runtime_sec": time.time() - t0}
    (HERE / "result_g3_higgs_factorization_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
