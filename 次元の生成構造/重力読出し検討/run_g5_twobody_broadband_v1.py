#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G5正本: 広帯域海での二体時計結合——引力の普遍性と質量積対応

背景（事前記録）: 単色海では二体結合は定在波干渉が支配（P21c・変調>平均・
P15のE∝M₁M₂とd^{-1/2}は汚染込みと訂正）。実海に近い広帯域海で再測する。
海=奇k6モード{1,3,5,7,9,11}・黄金比位相（決定論）・全モードW_f=0（固定点維持）。
判定（事前固定）:
(i) 引力の普遍性: 全ペア×全分離で E<0。
(ii) |Ē| の回帰: 不変量積 √(M_t²+σ²) が M_t積と同等以上に潰す（R²比較）。
(iii) 直交対（S3/S4: 同M_t・異σ）: 共通の重い相手に対する E比>1 なら
     σも重力荷に寄与（等価原理=閉塞恒等式の方向）。
使い方: python3 run_g5_twobody_broadband_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g5", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n)
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

    def ladder_at(fsrc, amp, center, sig_k=32.0):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    def run_tau(lump, Tburn=500, Tavg=200):
        a2 = ((sea + lump)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        acc = np.zeros(n)
        for j in range(Tburn + Tavg):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            if j >= Tburn:
                acc += np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
        return acc / Tavg

    SRCS = {"S1": (0.3, 0.05), "S2": (0.9, 0.05), "S3": (0.3, 0.10),
            "S4": (0.9, 0.10), "S5": (0.6, 0.07)}
    SEPS = [24, 52, 81]; cA = 100
    tau0 = run_tau(np.zeros(n))
    soloA = {nm: run_tau(ladder_at(f, a, cA)) for nm, (f, a) in SRCS.items()}
    soloB = {(nm, d): run_tau(ladder_at(*SRCS[nm], (cA + d) % n))
             for d in SEPS for nm in SRCS}
    results = {}; all_neg = True
    for (nA, nB) in combinations(SRCS.keys(), 2):
        Es = []
        for d in SEPS:
            t12 = run_tau(ladder_at(*SRCS[nA], cA) + ladder_at(*SRCS[nB], (cA + d) % n))
            Es.append(float(np.sum(t12 - soloA[nA] - soloB[(nB, d)] + tau0)))
        if any(e > 0 for e in Es): all_neg = False
        results[f"{nA}+{nB}"] = {"E_by_sep": Es, "E_bar": float(np.mean(Es))}
        print(f"{nA}+{nB}: Ē={np.mean(Es):+.3e} 符号列=" +
              "".join("+" if e > 0 else "-" for e in Es))
    dAmask = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 8
    comp = {}
    for nm in SRCS:
        dtt = soloA[nm] - tau0
        comp[nm] = {"Mt": float(np.mean(np.abs(dtt[dAmask]))),
                    "sig": float(np.std(dtt[dAmask]))}
    Ebars = np.array([abs(v["E_bar"]) for v in results.values()])
    keys = list(results.keys())
    regs = {}
    for name, P in (("Mt_product", np.array([comp[a[:2]]["Mt"] * comp[a[3:]]["Mt"] for a in keys])),
                    ("invariant_product",
                     np.array([np.sqrt(comp[a[:2]]["Mt"]**2 + comp[a[:2]]["sig"]**2)
                               * np.sqrt(comp[a[3:]]["Mt"]**2 + comp[a[3:]]["sig"]**2) for a in keys]))):
        ok = (Ebars > 1e-15) & (P > 1e-15)
        pc = np.polyfit(np.log(P[ok]), np.log(Ebars[ok]), 1)
        pred = np.polyval(pc, np.log(P[ok]))
        R2 = 1 - np.sum((np.log(Ebars[ok]) - pred) ** 2) / np.sum((np.log(Ebars[ok]) - np.log(Ebars[ok]).mean()) ** 2)
        regs[name] = {"exponent": float(pc[0]), "R2": float(R2)}
        print(f"回帰 {name}: p={pc[0]:+.3f} R²={R2:.4f}")
    r_o1 = abs(results["S3+S5"]["E_bar"] / results["S4+S5"]["E_bar"])
    r_o2 = abs(results["S2+S3"]["E_bar"] / results["S2+S4"]["E_bar"])
    print(f"直交対比: E(S3+S5)/E(S4+S5)={r_o1:.3f}  E(S2+S3)/E(S2+S4)={r_o2:.3f}（>1ならσ寄与）")
    verdict = ("柱G5成立（普遍的引力・不変量方向）"
               if all_neg and r_o1 > 1 and r_o2 > 1 else "部分成立")
    print(verdict)
    out = {"pairs": results, "components": comp, "regressions": regs,
           "orthogonal_ratios": [r_o1, r_o2], "all_attractive": all_neg,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g5_twobody_broadband_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
