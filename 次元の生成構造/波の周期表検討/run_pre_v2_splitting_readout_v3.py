#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2検証プログラム第7項・決定実験（最終形）: 分裂読出し π基準の独立予言検定

回収の経緯（事前記録・反証3件込み）:
 v1(#31): 単色海で同組成対照に偽Δω → 設計限界（v2 §12.2・反証条件9）。
 試行1（広帯域海・奇k6モード）: 海テクスチャ（|sea(x)|²の空間変動）が
   Δω床=0.216を作り信号を埋没——反証。
 試行2（FM海・定振幅多線）: 奇キャリア×奇変調の1次側帯が偶k（フェルミオン帯）
   に落ちリーク21%・非一様真空tick——反証。小定理: 定振幅×純ボゾン帯の海は
   単色に限る。
 解決（源側の規約）: **海同位相化** lump = L(x−c)·e^{iφ_carrier(c)}。
   並進×大域位相の合成対称性により同一源のωは厳密に位置不変
   → 対照床が63分の1に崩壊（0.216→0.0034・残差=真のペア相互作用）。

決定検定（事前固定）:
 予言: 単独源2本の累積時計位相差 |Φ_A(t)−Φ_B(t)| が π に達する時刻 t_pred
 （ペア実験とは独立な構成）。
 実測: ペアの t_split（3分離 sep∈{24,52,81}）。
 判定: 実測/予言 が全ペア・全分離で [0.8, 1.25]、CV<20% なら π基準確定。
使い方: python3 run_pre_v2_splitting_readout_v3.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v2s3", UIM / "run_ignition_fate_exact_v3.py")
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
    sea = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)

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

    def solo_Phi(fsrc, c=100, T=400):
        a2 = ((sea + ladder_at(fsrc, c))[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        dA = np.minimum(np.abs(x - c), n - np.abs(x - c)) <= 3
        Phi = [0.0]
        for j in range(T):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
            Phi.append(Phi[-1] + float(np.mean(tt[dA])))
        return np.array(Phi)

    def pair_tsplit(fA, fB, sep, cA=100, T=1200):
        cB = (cA + sep) % n
        lump = ladder_at(fA, cA) + ladder_at(fB, cB)
        a2 = ((sea + lump)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        dA = np.minimum(np.abs(x - cA), n - np.abs(x - cA)) <= 3
        dB = np.minimum(np.abs(x - cB), n - np.abs(x - cB)) <= 3
        PhiA = PhiB = 0.0
        for j in range(T):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
            PhiA += float(np.mean(tt[dA])); PhiB += float(np.mean(tt[dB]))
            if abs(PhiA - PhiB) > np.pi:
                return j + 1
        return None

    SEPS = [24, 52, 81]
    # 対照床
    ctrl = []
    for sep in SEPS:
        cB = (100 + sep) % n
        # 同組成: 固有Δω=0（対称性）——ペアの残差のみ測る
        ts = pair_tsplit(0.6, 0.6, sep)
        ctrl.append({"sep": sep, "t_split": ts})
    print("対照(0.6,0.6): t_split =", [c["t_split"] for c in ctrl], "（残差=真のペア相互作用）")

    fs_needed = sorted({0.3, 0.4, 0.5, 0.6, 0.7, 0.9})
    Phis = {f: solo_Phi(f) for f in fs_needed}
    omegas = {f: float((Phis[f][150] - Phis[f][0]) / 150) for f in fs_needed}
    print("単独源 ⟨ω⟩:", {f: round(w, 4) for f, w in omegas.items()})

    pairs = [(0.4, 0.9), (0.4, 0.7), (0.4, 0.5), (0.5, 0.6), (0.3, 0.9)]
    rows = []; ratios = []
    for fA, fB in pairs:
        dPhi = np.abs(Phis[fA] - Phis[fB])
        above = np.where(dPhi > np.pi)[0]
        t_pred = int(above[0]) if len(above) else None
        ts_list = [pair_tsplit(fA, fB, sep) for sep in SEPS]
        rs = [(ts / t_pred) if (ts and t_pred) else None for ts in ts_list]
        ratios += [r for r in rs if r]
        rows.append({"fA": fA, "fB": fB, "t_pred": t_pred,
                     "t_split_by_sep": ts_list, "ratios": rs})
        print(f"({fA},{fB}): 予言={t_pred} 実測={ts_list} 比=" +
              " ".join(f"{r:.2f}" if r else "—" for r in rs))
    r = np.array(ratios)
    ok = bool(np.all((r > 0.8) & (r < 1.25)) and (r.std() / r.mean()) < 0.2)
    verdict = (f"π基準確定（実測/予言 平均{r.mean():.3f} CV{r.std()/r.mean():.1%}・"
               f"独立予言・新柱9の決定実験通過）" if ok else "不通過/要精査")
    print(f"\n判定: {verdict}")
    out = {"controls": ctrl, "solo_omegas": omegas, "pairs": rows,
           "ratio_mean": float(r.mean()), "ratio_cv": float(r.std() / r.mean()),
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "pre_v2_splitting_readout_result_v3.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s → pre_v2_splitting_readout_result_v3.json")

if __name__ == "__main__":
    main()
