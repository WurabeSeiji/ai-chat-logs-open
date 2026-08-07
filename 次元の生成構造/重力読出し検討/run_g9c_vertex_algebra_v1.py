#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G9完成: 頂点代数——二層選択則の閉じた規則（検証プログラム2番の実施）

倍加軌道（mod ne=16・海 m_s=0 で和則 m*=2m_B−m_s=2m_B）:
  1→2→4→8→0 ／ 3→6→12→8→0 ／ 5→10→4→8→0
  注: m=8（半Nyquist）は全ての奇数巻きの軌道の吸引点、0（海）は終点。

頂点規則（事前記録・検定対象）:
  チャネル(m_A,m_B)が開く ⟺ m_B≡2^k·m_A または m_A≡2^k·m_B (mod ne)
  （k=0 が静的コヒーレンス、k≥1 が動的ウォーク到達）。開度は k とともに減衰。
  符号（共役）不整合な対は閉じたまま（G9bで実証済の静的直交性）。

判定（事前固定）:
 (V1) 1段到達の開放普遍性: (1,2)(2,4)(3,6) の超過（等質量盲目基底との差）が
      すべて負（チャネル開）かつ |超過| > 0.3·|静的coh|。
 (V2) 非到達対の盲目性: (1,3)(1,5)(1,6)(2,3)(2,6)(3,4) は巻き符号反転に
      不変（|E(a,b)−E(a,−b)|/|E| < 0.01）——盲目項は符号を見ない。
 (V3) 到達次数減衰: |超過(1,2)|（1段） > |超過(1,4)|（2段）。
使い方: python3 run_g9c_vertex_algebra_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g9c", UIM / "run_ignition_fate_exact_v3.py")
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
    ms_needed = (1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6)
    soloA = {m: run_tau(charged_lump(m, cA)) for m in ms_needed}
    soloB = {(m, d): run_tau(charged_lump(m, (cA + d) % n)) for d in SEPS for m in ms_needed}

    def E_pair(mA, mB):
        Es = [float(np.sum(run_tau(charged_lump(mA, cA) + charged_lump(mB, (cA + d) % n))
                           - soloA[mA] - soloB[(mB, d)] + tau0)) for d in SEPS]
        return float(np.mean(Es))

    PAIRS = [  # (mA, mB)
        (1, 1), (3, 3),                       # 静的 k=0（cohスケール）
        (1, -1), (3, -3),                     # 静的の盲目基底
        (1, 2), (2, 4), (3, 6),               # 1段到達
        (1, -2), (-1, 2), (2, -4), (-2, 4), (3, -6), (-3, 6),  # 等質量盲目基底
        (1, 4), (1, -4), (-1, 4),             # 2段到達＋基底
        (1, 3), (1, -3), (1, 5), (-1, 5), (1, 6), (1, -6),     # 非到達＋符号対
        (2, 3), (2, -3), (2, 6), (2, -6), (3, 4), (-3, 4),
    ]
    E = {}
    for mA, mB in PAIRS:
        key = f"{mA:+d},{mB:+d}"
        E[key] = E_pair(mA, mB)
        print(f"({key}): Ē={E[key]:+.4e}")

    coh_static = 0.5 * ((E["+1,+1"] - E["+1,-1"]) + (E["+3,+3"] - E["+3,-3"]))
    exc = {
        "1,2": E["+1,+2"] - 0.5 * (E["+1,-2"] + E["-1,+2"]),
        "2,4": E["+2,+4"] - 0.5 * (E["+2,-4"] + E["-2,+4"]),
        "3,6": E["+3,+6"] - 0.5 * (E["+3,-6"] + E["-3,+6"]),
        "1,4": E["+1,+4"] - 0.5 * (E["+1,-4"] + E["-1,+4"]),
    }
    print(f"\n静的coh（(1,1)/(3,3)平均）= {coh_static:+.3e}")
    for kk, v in exc.items():
        print(f"超過({kk}) = {v:+.3e}（/|coh| = {abs(v/coh_static):.2f}）")
    v1 = all(exc[p] < 0 and abs(exc[p]) > 0.3 * abs(coh_static)
             for p in ("1,2", "2,4", "3,6"))
    blind_checks = {p: abs(E[a] - E[b]) / abs(E[a]) for p, a, b in [
        ("1,3", "+1,+3", "+1,-3"), ("1,5", "+1,+5", "-1,+5"),
        ("1,6", "+1,+6", "+1,-6"), ("2,3", "+2,+3", "+2,-3"),
        ("2,6", "+2,+6", "+2,-6"), ("3,4", "+3,+4", "-3,+4")]}
    v2 = all(v < 0.01 for v in blind_checks.values())
    v3 = abs(exc["1,2"]) > abs(exc["1,4"])
    print(f"\n(V1) 1段到達の開放普遍性（3対とも開・>0.3coh）: {'通過' if v1 else '不成立'}")
    print("(V2) 非到達対の符号盲目性: " + "  ".join(
        f"({p})={v:.4f}" for p, v in blind_checks.items())
        + f" → {'通過（全て<0.01）' if v2 else '不成立'}")
    print(f"(V3) 到達次数減衰 |超過(1,2)|>|超過(1,4)|: "
          f"{abs(exc['1,2']):.3f} vs {abs(exc['1,4']):.3f} → {'通過' if v3 else '不成立'}")
    verdict = ("頂点規則成立: チャネル(m_A,m_B)開 ⟺ 倍加軌道到達 m_B≡2^k·m_A (mod ne)・"
               "開度はkで減衰・符号不整合は閉" if (v1 and v2 and v3) else "要精査")
    print(verdict)
    out = {"E": E, "coh_static": coh_static, "excess": exc,
           "blind_sign_checks": blind_checks, "V1": v1, "V2": v2, "V3": v3,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g9c_vertex_algebra_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
