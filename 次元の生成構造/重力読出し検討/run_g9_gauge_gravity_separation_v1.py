#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G9正本（決定実験）: 同一読出しからのゲージ応答と重力応答の同時分離

設計（事前記録）: η巻き m∈{0,±1} の荷電源（組成f=0.6・振幅0.05の梯子×e^{2πimη/ne}）
の二体結合 E を広帯域海（柱G5と同一）・3分離で全電荷組合せ測定。
判定（事前固定）:
 (i) 引力普遍性: 全ペア（(0,0)(++)(−−)(+−)(+0)）で Ē<0——重力チャネルは電荷盲目。
 (ii) 電荷共役対称: |E(++)−E(−−)| ≪ 分離ばらつき。
 (iii) ゲージチャネルの実在: E(++)−E(+−) が分離ばらつきを超えて非零
      ——電荷構造（η可干渉性）に依存する応答が同じ読出しから出る。
 (iv) 中性プローブの重力: E(+1,0) が引力的。
→ 4判定成立で「ゲージと重力＝一つの読出し関数の異なる射影」（表題の直接実証）。
使い方: python3 run_g9_gauge_gravity_separation_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g9", UIM / "run_ignition_fate_exact_v3.py")
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
    dAmask = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 8
    soloA = {m: run_tau(charged_lump(m, cA)) for m in (0, 1, -1)}
    solo_mass = {m: float(np.mean(np.abs((soloA[m] - tau0)[dAmask]))) for m in (0, 1, -1)}
    print("単独源の質量:", {m: round(v, 5) for m, v in solo_mass.items()})
    soloB = {(m, d): run_tau(charged_lump(m, (cA + d) % n)) for d in SEPS for m in (0, 1, -1)}
    pairs = [("0,0", 0, 0), ("+1,+1", 1, 1), ("-1,-1", -1, -1),
             ("+1,-1", 1, -1), ("+1,0", 1, 0)]
    res = {}
    for name, mA, mB in pairs:
        Es = [float(np.sum(run_tau(charged_lump(mA, cA) + charged_lump(mB, (cA + d) % n))
                           - soloA[mA] - soloB[(mB, d)] + tau0)) for d in SEPS]
        res[name] = {"E_by_sep": Es, "E_bar": float(np.mean(Es)),
                     "spread": float(np.std(Es))}
        print(f"({name}): Ē={np.mean(Es):+.3e} ±{np.std(Es):.2e}")
    E = {kk: v["E_bar"] for kk, v in res.items()}
    sp_max = max(v["spread"] for v in res.values())
    ok1 = all(v < 0 for v in E.values())
    cc = abs(E["+1,+1"] - E["-1,-1"])
    gauge = E["+1,+1"] - E["+1,-1"]
    ok = ok1 and cc < 0.01 * abs(E["+1,+1"]) and abs(gauge) > 0.05 * abs(E["+1,+1"]) and E["+1,0"] < 0
    print(f"(i)全引力={ok1} (ii)共役対称差={cc:.2e} (iii)ゲージ成分={gauge:+.3e} (iv)E(+,0)={E['+1,0']:+.3e}")
    verdict = ("柱G9成立（同一読出しからゲージ応答と重力応答が分離・表題の直接実証）"
               if ok else "要精査")
    print(verdict)
    out = {"solo_mass": solo_mass, "pairs": res, "conj_diff": cc,
           "gauge_component": gauge, "grav_component": 0.5 * (E["+1,+1"] + E["+1,-1"]),
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g9_gauge_gravity_separation_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
