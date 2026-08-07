#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2実験3: 分裂の読出し検定（新柱9の実測正本）

エンジン: 局所化プロトタイプ（v2実験2と同一・円偏波）。
定義: 二パケット（組成 fA≠fB・間隔24）の中心3セル平均時計位相の累積差 |ΦA−ΦB|。
判定基準（事前固定・πは自然基準で調整パラメータなし）:
  |ΦA−ΦB| < π のあいだ「1粒子」、超えたら「2粒子」。予言 t_split=π/Δω。
  検定: t_split·Δω/π が全ペアでオーダー1（0.1〜10）なら読出し成立。
既知の設計限界（正直な登録）: 海キャリア位相（λ=n/3≈170.7）がパケット位置で
異なるため、同組成対照でも偽Δωが生じる（キャリア位相整列配置は検証プログラム
(f)に登録済み・本実験は非整列のまま実測値を記録する）。
使い方: python3 run_pre_v2_splitting_readout_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v2s", UIM / "run_ignition_fate_exact_v3.py")
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
        return AMP * prof

    def run_pair(fA, fB, sep=24, T=3000):
        cA, cB = n // 2 - sep // 2, n // 2 + sep // 2
        lump = ladder_at(fA, cA) + ladder_at(fB, cB)
        a2 = ((0.2 * carrier + lump)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        dA = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 3
        dB = np.minimum(np.abs(x - cB), n - np.abs(x - cB)) <= 3
        PhiA = PhiB = 0.0; dphi_series = []; t_split = None
        omA, omB = [], []
        for j in range(T):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
            wA, wB = float(np.mean(tt[dA])), float(np.mean(tt[dB]))
            PhiA += wA; PhiB += wB
            omA.append(wA); omB.append(wB)
            if j % 5 == 0:
                dphi_series.append(abs(PhiA - PhiB))
            if t_split is None and abs(PhiA - PhiB) > np.pi:
                t_split = j + 1
        dOm = abs(float(np.mean(omA)) - float(np.mean(omB)))
        return t_split, dOm, dphi_series

    pairs = [(0.4, 0.9), (0.4, 0.7), (0.4, 0.5), (0.5, 0.6), (0.45, 0.5), (0.3, 0.9)]
    rows = []
    for fA, fB in pairs:
        ts, dOm, series = run_pair(fA, fB)
        ratio = (ts * dOm / np.pi) if ts else None
        rows.append({"fA": fA, "fB": fB, "t_split": ts, "delta_omega": dOm,
                     "ratio": ratio, "dphi_decimated": series[:200], "decimation": 5})
        print(f"fA={fA} fB={fB}: Δω={dOm:.3e} t_split={ts} 比={ratio if ratio else '—'}")
    ratios = np.array([r["ratio"] for r in rows if r["ratio"]])
    ok = np.all((ratios > 0.1) & (ratios < 10))
    verdict = (f"読出し成立（t·Δω/π 平均={ratios.mean():.2f} CV={ratios.std()/ratios.mean():.1%}・"
               f"全ペアがオーダー1）" if ok else "不成立")
    # 対照（同組成・キャリア位相非整列の限界を記録）
    ts0, dOm0, _ = run_pair(0.6, 0.6)
    print(f"対照(fA=fB): Δω={dOm0:.2e} t_split={ts0}（キャリア位相非整列の偽Δω・設計限界として記録）")
    print(verdict)
    out = {"rows": rows, "ratio_mean": float(ratios.mean()), "ratio_cv": float(ratios.std() / ratios.mean()),
           "control_same_comp": {"t_split": ts0, "delta_omega": dOm0,
                                  "note": "海キャリア位相の非対称による偽Δω（検証プログラム(f)=位相整列で除去予定）"},
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "pre_v2_splitting_readout_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_v2_splitting_readout_result_v1.json")

if __name__ == "__main__":
    main()
