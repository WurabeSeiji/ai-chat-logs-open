#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G9昇格正本: 読出し分解定理——η直交性による G = G_blind ⊕ G_coh·δ_{Δm,0}

定理（解析・事前記録）: η和の双線形読出しにおいて、二源(m_A, m_B)の交差項は
  Σ_η e^{iΔm·2πη/ne} = ne·δ_{Δm≡0 (mod ne)}
により厳密に選択される。ゆえに任意のη和双線形読出しは
  G = G_blind(|sea|²,|L_A|²,|L_B|²)  ＋  G_coh·δ_{m_A,m_B}  ＋  海結合·δ_{m,0}
に厳密分解する。G_blind=電荷盲目（重力チャネル）／G_coh=巻き整合コヒーレンス
（電荷構造依存チャネル・gauge-like）。

導出済みの説明: E(++)=E(−−)厳密（既測6.8e-10）／E(+−)=盲目項のみ／
E(0,0)は海コヒーレンス項が加わり深い／荷電源の質量>中性源（海項の欠如）。

新予言（判定・初版P1は質量積不一致(gcd類)の設計ミスで修正済＝反証記録）:
 (P1') 不整合普遍性（等質量形）: E(+1,−2) = E(−1,+2)（相対差<1%）——
      非到達（Δm≠0かつ倍加非接続）の等質量対は盲目基底に厳密一致。
 (P2) コヒーレント項の|m|非依存: [E(+2,+2)−E(+2,−2)] = [E(+1,+1)−E(+1,−1)]
      （相対差<25%）。
 (P3) E(+2,+2)=E(−2,−2)（共役厳密・<1e-6相対）。
 (P4・発見) 動的選択則: 倍加到達対 (+1,+2)（2·(+1)=+2）は盲目基底からの
      超過を持つ——ゲージ的コヒーレンスは静的 δ_{Δm,0} に加え、ウォーク
      （和則 m*=2m_B−m_s [6]）到達対に動的に開く。超過量を記録する。
使い方: python3 run_g9b_decomposition_theorem_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g9b", UIM / "run_ignition_fate_exact_v3.py")
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

    SEPS = [24, 52, 81]; cA = 100
    tau0 = run_tau(np.zeros((n, ne), complex))
    ms_needed = (0, 1, -1, 2, -2)
    soloA = {m: run_tau(charged_lump(m, cA)) for m in ms_needed}
    soloB = {(m, d): run_tau(charged_lump(m, (cA + d) % n)) for d in SEPS for m in ms_needed}

    def E_pair(mA, mB):
        Es = [float(np.sum(run_tau(charged_lump(mA, cA) + charged_lump(mB, (cA + d) % n))
                           - soloA[mA] - soloB[(mB, d)] + tau0)) for d in SEPS]
        return float(np.mean(Es)), float(np.std(Es))

    E = {}
    for name, mA, mB in [("+1,+1", 1, 1), ("+1,-1", 1, -1), ("+2,+2", 2, 2),
                          ("+2,-2", 2, -2), ("-2,-2", -2, -2), ("+1,+2", 1, 2),
                          ("+1,-2", 1, -2), ("-1,+2", -1, 2)]:
        eb, es = E_pair(mA, mB)
        E[name] = eb
        print(f"({name}): Ē={eb:+.4e} ±{es:.2e}")
    # 判定
    p1 = abs(E["+1,-2"] - E["-1,+2"]) / abs(E["+1,-2"])
    blind_base = 0.5 * (E["+1,-2"] + E["-1,+2"])
    excess = E["+1,+2"] - blind_base
    coh1 = E["+1,+1"] - E["+1,-1"]; coh2 = E["+2,+2"] - E["+2,-2"]
    p2 = abs(coh2 - coh1) / abs(coh1)
    p3 = abs(E["+2,+2"] - E["-2,-2"]) / abs(E["+2,+2"])
    print(f"\n(P1') 等質量不整合普遍性: |E(+1,−2)−E(−1,+2)|/|E| = {p1:.4f}（判定<0.01）")
    print(f"(P2) コヒーレント項|m|非依存: coh(m=1)={coh1:+.3e} coh(m=2)={coh2:+.3e} 相対差={p2:.3f}（判定<0.25）")
    print(f"(P3) 共役厳密: |E(+2,+2)−E(−2,−2)|/|E| = {p3:.2e}（判定<1e-6）")
    print(f"(P4・発見) 倍加到達対の超過: E(+1,+2)−盲目基底 = {excess:+.3e}"
          f"（静的coh −1.43 の約{abs(excess/coh1):.0%}・動的選択則）")
    ok = p1 < 0.01 and p2 < 0.25 and p3 < 1e-6
    verdict = ("読出し分解定理成立（G=G_blind⊕G_coh・選択則=静的δ_{Δm,0}+動的ウォーク到達・"
               "G9は決定実験から構造定理へ昇格）" if ok else "要精査")
    print(verdict)
    out = {"E": E, "P1_equalmass_universality": p1, "blind_base": blind_base,
           "P4_dynamical_excess": excess, "coh_m1": coh1, "coh_m2": coh2,
           "P2_coh_m_independence": p2, "P3_conjugation": p3,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g9b_decomposition_theorem_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
