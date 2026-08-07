#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AC-0: 場の読出し万能関数へのスピン読出し追加と検定（電子の反跳実験）

拡張（Gの規約——純関係量・パラメータなし・種別分岐なし——に従う）:
  二チャネル (a,b) の局所 Stokes 射影族をスピン読出しとして追加する。
    S0(x)=Σ_η(|a|²+|b|²),  S1(x)=Σ_η(|a|²−|b|²),
    S2(x)=2ReΣ_η(ā·b),     S3(x)=2ImΣ_η(ā·b)
    σ_i(x)=S_i/S0（i=1,2,3）・測定角 α のスピン読出し σ_α=σ1cosα+σ2sinα
  円偏波 b=−ia（確定電荷の厳密不変多様体）は σ1=σ2=0, σ3=−1。
  σ1σ2 面が未確定＝スピン灰色。EPRの測定角 α は σ1σ2 面の回転として実装できる。

判定（事前固定）:
 (A0a) 読出し無害性: スピン読出しは状態に非接触（読出し前後で状態ビット同一・
       閉塞ドリフトは素の力学と同一）。
 (A0b) 灰色誕生: 電子lump（b=−ia）の台上 max|σ1|,|σ2| < 1e-12・σ3=−1。
 (A0c) 灰色不変性: 素の力学 T=2000 の間、台上 max|σ1|,|σ2| < 1e-10 を維持
       ——スピンは相互作用なしには決まらない（不変多様体の実測）。
 (A0d) 光子C: 単色ボゾン1モードは時を刻まない（純ボゾン海の固定点・
       τ_t std < 1e-12）——観測装置C自身は時計を持たない。
使い方: python3 run_ac0_spin_readout_extension_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_ac0", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

M_ELECTRON = -3
FSRC = 0.7


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

    def prof(center, fsrc=FSRC):
        p = np.zeros(n); AMP = 0.05 * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / 32.0) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            p += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * p

    def lump(m, center):
        return prof(center)[:, None] * np.exp(2j * np.pi * m * eta / ne)[None, :]

    # --- G拡張: スピン読出し（Stokes射影族・純関係量・パラメータなし） ---
    def spin_readout(a2, b2):
        S0 = np.sum(np.abs(a2) ** 2 + np.abs(b2) ** 2, axis=1)
        S1 = np.sum(np.abs(a2) ** 2 - np.abs(b2) ** 2, axis=1)
        Z = np.sum(np.conj(a2) * b2, axis=1)
        S2 = 2.0 * np.real(Z); S3 = 2.0 * np.imag(Z)
        d = np.maximum(S0, 1e-300)
        return S1 / d, S2 / d, S3 / d

    cA = 100
    mask = (np.abs((x - cA + n // 2) % n - n // 2) <= 24)

    # (A0a) 読出し無害性: 読出しは状態に非接触（同一初期から読出し有無で状態一致）
    a1 = (sea[:, None] * np.ones((1, ne)) + lump(M_ELECTRON, cA)).astype(complex)
    b1 = -1j * a1
    a2, b2 = a1.copy(), b1.copy()
    for j in range(200):
        a1, b1 = step(a1, b1)                    # 読出しなし
        a2, b2 = step(a2, b2); spin_readout(a2, b2)  # 毎步読出し
    same = np.array_equal(a1, a2) and np.array_equal(b1, b2)
    C1 = float(np.abs(np.sum(a1 ** 2 + b1 ** 2)))
    print(f"(A0a) 読出し無害性: 状態ビット同一={same}・|Σz²|={C1:.3e} → "
          f"{'通過' if same else '不成立'}")

    # (A0b) 灰色誕生
    a0 = (sea[:, None] * np.ones((1, ne)) + lump(M_ELECTRON, cA)).astype(complex)
    b0 = -1j * a0
    s1, s2, s3 = spin_readout(a0, b0)
    m12_0 = float(max(np.max(np.abs(s1[mask])), np.max(np.abs(s2[mask]))))
    s3_0 = float(np.mean(s3[mask]))
    okb = (m12_0 < 1e-12) and abs(s3_0 + 1) < 1e-12
    print(f"(A0b) 灰色誕生: max|σ1,σ2|={m12_0:.2e}（<1e-12）・σ3={s3_0:+.6f} → "
          f"{'通過' if okb else '不成立'}")

    # (A0c) 灰色不変性（素の力学 T=2000）
    aa, bb = a0.copy(), b0.copy()
    worst = 0.0
    for j in range(2000):
        aa, bb = step(aa, bb)
        if j % 100 == 99:
            s1, s2, _ = spin_readout(aa, bb)
            worst = max(worst, float(max(np.max(np.abs(s1[mask])),
                                         np.max(np.abs(s2[mask])))))
    okc = worst < 1e-10
    print(f"(A0c) 灰色不変性(T=2000): max|σ1,σ2|={worst:.2e}（<1e-10） → "
          f"{'通過' if okc else '不成立'}")

    # (A0d) 光子C: 単色ボゾン1モード（奇k・海と同帯）は時を刻まない
    kC = 7
    C = 0.2 * np.exp(2j * np.pi * kC * x / n)
    aC = (C[:, None] * np.ones((1, ne))).astype(complex); bC = -1j * aC
    taus = []
    for j in range(300):
        ap = aC.copy(); aC, bC = step(aC, bC)
        taus.append(np.angle(np.einsum("xe,xe->x", np.conj(ap), aC)))
    tau_std = float(np.std(np.array(taus)))
    okd = tau_std < 1e-12
    print(f"(A0d) 単色光子C: τ_t std={tau_std:.2e}（<1e-12・時を刻まない） → "
          f"{'通過' if okd else '不成立'}")

    ok = same and okb and okc and okd
    verdict = ("スピン読出し拡張成立: Stokes射影族はGの規約を満たし、"
               "確定電荷=スピン灰色（σ1=σ2=0厳密）は素の力学の不変多様体、"
               "光子Cは時計を持たない" if ok else "要精査")
    print(verdict)
    out = {"A0a_state_identical": bool(same), "closure_abs": C1,
           "A0b_gray_birth": {"max_s12": m12_0, "s3": s3_0, "ok": bool(okb)},
           "A0c_gray_invariance": {"worst_s12": worst, "ok": bool(okc)},
           "A0d_photon_clockless": {"tau_std": tau_std, "ok": bool(okd)},
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_ac0_spin_readout_extension_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
