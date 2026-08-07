#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2検証プログラム第7項の実施: 分裂読出しの決定実験（広帯域海版）

背景（事前記録）: v1実験(#31)は単色キャリア海で同組成対照にも偽Δωが立った
（v2 §12.2 の既知設計限界・反証条件9）。その後の重力プログラム（P21c/P22）で
汚染源が単色キャリアとの定在波干渉そのものと同定され、n=512では奇数Kの
整数波長が不在のため「整列」は不可能、正しい回収は**広帯域海**と確定した。

設計: 海=奇k6モード{1,3,5,7,9,11}・黄金比位相（決定論）・全モード W_f=0
（真空固定点は厳密維持）。分離 sep∈{24,52,81}。
判定（事前固定）:
 (i) 対照（同組成 fA=fB=0.6）: Δω_floor を定量。t_split が観測窓 T=3000 内に
     ないこと、または組成差対より系統的に遅いこと。
 (ii) 組成差対（Δω > 3×floor のもの）: t_split·Δω/π ∈ [0.5, 2] かつ
     CV が v1 の 67% から改善すれば π基準は決定実験通過（新柱9確定）。
 (iii) 対照が組成差対と同等の Δω・t_split を示せば π基準は棄却（反証条件9）。
使い方: python3 run_pre_v2_splitting_readout_v2.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v2s2", UIM / "run_ignition_fate_exact_v3.py")
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

    def ladder_at(fsrc, center, sig_k=32.0, amp=0.05):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof

    def run_pair(fA, fB, sep, cA=100, T=3000):
        cB = (cA + sep) % n
        lump = ladder_at(fA, cA) + ladder_at(fB, cB)
        a2 = ((sea + lump)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        dA = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 3
        dB = np.minimum(np.abs(x - cB), n - np.abs(x - cB)) <= 3
        PhiA = PhiB = 0.0; t_split = None
        omA, omB = [], []
        for j in range(T):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
            wA, wB = float(np.mean(tt[dA])), float(np.mean(tt[dB]))
            PhiA += wA; PhiB += wB
            omA.append(wA); omB.append(wB)
            if t_split is None and abs(PhiA - PhiB) > np.pi:
                t_split = j + 1
        dOm = abs(float(np.mean(omA)) - float(np.mean(omB)))
        return t_split, dOm

    SEPS = [24, 52, 81]
    print("== 対照（同組成 fA=fB=0.6）——Δω床とt_splitの定量 ==")
    floor_dOms = []
    ctrl_rows = []
    for sep in SEPS:
        ts, dOm = run_pair(0.6, 0.6, sep)
        floor_dOms.append(dOm)
        ctrl_rows.append({"sep": sep, "t_split": ts, "delta_omega": dOm})
        print(f"  sep={sep}: Δω={dOm:.3e} t_split={ts if ts else '>3000'}")
    floor = float(np.mean(floor_dOms))
    print(f"  Δω床 = {floor:.3e}")

    print("\n== 組成差対 ==")
    pairs = [(0.4, 0.9), (0.4, 0.7), (0.4, 0.5), (0.5, 0.6), (0.3, 0.9)]
    rows = []
    for fA, fB in pairs:
        for sep in SEPS:
            ts, dOm = run_pair(fA, fB, sep)
            ratio = (ts * dOm / np.pi) if ts else None
            above = dOm > 3 * floor
            rows.append({"fA": fA, "fB": fB, "sep": sep, "t_split": ts,
                         "delta_omega": dOm, "ratio": ratio, "above_floor": above})
            print(f"  fA={fA} fB={fB} sep={sep}: Δω={dOm:.3e}"
                  f"（床の{dOm/floor:.1f}倍） t_split={ts if ts else '>3000'}"
                  f" 比={f'{ratio:.3f}' if ratio else '—'}")
    ratios = np.array([r["ratio"] for r in rows if r["ratio"] and r["above_floor"]])
    if len(ratios):
        cv = float(ratios.std() / ratios.mean())
        in_range = np.all((ratios > 0.5) & (ratios < 2.0))
        ctrl_ok = all((r["t_split"] is None) or
                      (r["delta_omega"] < floor * 1.5 and (r["t_split"] or 1e9) >
                       np.median([q["t_split"] or 1e9 for q in rows if q["above_floor"]]))
                      for r in ctrl_rows)
        print(f"\n判定: 床上ペア n={len(ratios)}  t·Δω/π 平均={ratios.mean():.3f} CV={cv:.1%}"
              f"（v1: 1.29/67%）")
        verdict = ("π基準は決定実験通過（新柱9確定）"
                   if in_range and cv < 0.67 and ctrl_ok else
                   "部分通過/要精査")
        print(f"対照の健全性={'OK' if ctrl_ok else 'NG'} → {verdict}")
    else:
        verdict = "床上ペアなし（判定不能）"
        print(verdict)
    out = {"floor": floor, "controls": ctrl_rows, "pairs": rows,
           "ratio_mean": float(ratios.mean()) if len(ratios) else None,
           "ratio_cv": cv if len(ratios) else None,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "pre_v2_splitting_readout_result_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s → pre_v2_splitting_readout_result_v2.json")

if __name__ == "__main__":
    main()
