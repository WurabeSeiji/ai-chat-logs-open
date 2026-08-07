#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S1: 開度の閉形式・一体候補の反証実験（§18-2 の記録）

候補（事前記録）: 動的チャネル開度 = solo源のウォーク生成スペクトル比
√(P_2m/P_m)。成立すれば頂点代数の強度法則が一体量で閉じる。

判定（事前固定）:
 (S1a) 単元自己同型のスペクトル版: √(P_6/P_3)(m=3) = √(P_2/P_1)(m=1)（<1%）。
 (S1b) 一体候補: √(P_2m/P_m) が実測開度 (1,2)=0.69 と (2,4)=0.22 を
       ともに相対30%内で再現。

結果（実行済・本文§18-2）: S1a 通過（3クラス厳密同一 0.494）・S1b 反証——
一体量はクラス分裂 0.69/0.22 を持たない。開度は二体動力学（対共鳴）の性質。
使い方: python3 run_s1_opening_onebody_refutation_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_s1", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); eta = np.arange(ne)
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    L = np.exp(-((np.abs(k) / 3.0) ** 4))
    Wf = ((k % 2) == 0).astype(float) * (1.0 - L); Wb = 1.0 - Wf
    GOLD = 0.6180339887498949
    sea = np.zeros(n, complex)
    for kk in (1, 3, 5, 7, 9, 11):
        sea += (0.2 / np.sqrt(6)) * np.exp(2j * np.pi * kk * x / n
                                           + 2j * np.pi * ((kk * GOLD) % 1.0))

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

    def ladder_prof(center, fsrc=0.6, amp=0.05, sig_k=32.0):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    def charged_lump(m, center):
        return ladder_prof(center)[:, None] * np.exp(2j * np.pi * m * eta / ne)[None, :]

    def run_tau(lump2, Tburn=500, Tavg=200):
        a2 = (sea[:, None] * np.ones((1, ne)) + lump2).astype(complex)
        b2 = -1j * a2
        acc = np.zeros(n)
        for j in range(Tburn + Tavg):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                acc += np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
        return acc / Tavg


    cA = 100
    mask = (np.abs((x - cA + n // 2) % n - n // 2) <= 24)

    def eta_spectrum(m, T=700, snaps=(500, 600, 700)):
        a2 = (sea[:, None] * np.ones((1, ne)) + charged_lump(m, cA)).astype(complex)
        b2 = -1j * a2
        a0 = (sea[:, None] * np.ones((1, ne))).astype(complex)
        b0 = -1j * a0  # 海のみ対照
        P = np.zeros(ne); cnt = 0
        for j in range(T + 1):
            if j in snaps:
                da = a2 - a0; db = b2 - b0
                Fa = np.fft.fft(da[mask, :], axis=1); Fb = np.fft.fft(db[mask, :], axis=1)
                P += np.sum(np.abs(Fa) ** 2 + np.abs(Fb) ** 2, axis=0); cnt += 1
            if j < T:
                a2, b2 = step(a2, b2); a0, b0 = step(a0, b0)
        return P / cnt

    ratios = {}
    for m in (1, 2, 3):
        P = eta_spectrum(m)
        r1 = float(np.sqrt(P[(2 * m) % ne] / P[m % ne]))
        r2 = float(np.sqrt(P[(4 * m) % ne] / P[m % ne]))
        ratios[m] = {"sqrt_P2m_over_Pm": r1, "sqrt_P4m_over_Pm": r2}
        print(f"m={m}: √(P_2m/P_m)={r1:.3f}  √(P_4m/P_m)={r2:.3f}")
    MEAS = {"(1,2)": 0.69, "(2,4)": 0.22, "(1,4)": 0.03}  # G9c実測（/静的coh）
    s1a = abs(ratios[3]["sqrt_P2m_over_Pm"] - ratios[1]["sqrt_P2m_over_Pm"]) \
        / ratios[1]["sqrt_P2m_over_Pm"] < 0.01
    d12 = abs(ratios[1]["sqrt_P2m_over_Pm"] - MEAS["(1,2)"]) / MEAS["(1,2)"]
    d24 = abs(ratios[2]["sqrt_P2m_over_Pm"] - MEAS["(2,4)"]) / MEAS["(2,4)"]
    s1b = d12 < 0.30 and d24 < 0.30
    print(f"(S1a) 単元自己同型スペクトル版: {'通過' if s1a else '不成立'}")
    print(f"(S1b) 一体候補 相対差 (1,2):{d12:.2f} (2,4):{d24:.2f}（判定<0.30両方）: "
          f"{'通過' if s1b else '反証'}")
    verdict = ("一体候補成立" if s1b else
               "一体候補の反証: soloスペクトルはクラス厳密同一だが開度はクラス分裂"
               "——開度は二体動力学（対共鳴）の性質（§18-2）")
    print(verdict)
    out = {"ratios": ratios, "measured_openings": MEAS,
           "S1a": bool(s1a), "S1b": bool(s1b), "verdict": verdict,
           "runtime_sec": time.time() - t0}
    (HERE / "result_s1_opening_onebody_refutation_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
