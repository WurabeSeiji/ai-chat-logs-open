#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v14: 被覆度×χパリティ×約数類のクロス表——モデル版スピン統計定理の検定

背景: v13bで帯電census(χ偶・m混在)=Z₂量子化、中性束(χ奇・m=0)=連続と判明したが
両軸が交絡。本実験は χパリティ{偶,奇} × 巻き数{0,1,2,4} の8セルを独立構成し、
較正済み判定器で被覆度を測って「二値性がどの軸に付くか」を決める。

判定器（v13b較正の教訓を反映・事前固定）:
  回帰ピーク列（C局所最大・実効閾）で位相 Φ_k=arg A(J_k) を収集。
  量子化スコア Qz2 = ピークのうち min(|Φ|,|Φ−π|,|Φ+π|)<0.1π の割合。
  Qz2>0.9 かつ 両符号出現 → 被覆度2（Z₂）／Qz2>0.9 かつ 全正 → 被覆度1／
  それ以外 → 連続（被覆なし）。

スピン統計仮説（事前記録）:
  H_par: 被覆度2 ⟺ χ偶（フェルミオン分類）——パリティ軸に付く
  H_wind: 被覆度2 ⟺ m≠0（帯電）——巻き数軸に付く
  クロス表がどちらか（または積・別構造）を決める。

使い方: python3 run_pre_spin_statistics_cross_v14.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v14", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_MAX = 6000

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    def species(chi_bins, m_t):
        a = single_winding(v1.make_bundle(sp, chi_bins, "A", scale=1.0)) * S
        b = single_winding(v1.make_bundle(sp, chi_bins, "B", scale=1.0)) * S
        a = project_eta(a, {m_t}); b = project_eta(b, {m_t})
        pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
        if pw < 1e-9:
            return None
        sc = np.sqrt(64.0 / pw)
        return a * sc, b * sc

    def measure(a, b):
        a0_, b0_ = a.copy(), b.copy()
        s0 = np.imag(np.conj(b0_) * a0_); s0c = s0 - s0.mean()
        nrm = float(np.sum(np.abs(a0_) ** 2) + np.sum(np.abs(b0_) ** 2))
        A, C = [], []
        for j in range(J_MAX):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            A.append(complex((np.vdot(a0_, a) + np.vdot(b0_, b)) / nrm))
            sj = np.imag(np.conj(b) * a); sjc = sj - sj.mean()
            C.append(float(np.dot(s0c, sjc) /
                            max(np.linalg.norm(s0c) * np.linalg.norm(sjc), 1e-300)))
        A = np.array(A); C = np.array(C)
        th = 0.85 * float(C[10:].max())
        peaks = []
        for j in range(11, J_MAX - 1):
            if C[j] > th and C[j] >= C[j - 1] and C[j] >= C[j + 1]:
                if not peaks or j - peaks[-1] > 5:
                    peaks.append(j)
            if len(peaks) >= 15:
                break
        phis = np.array([np.angle(A[j]) for j in peaks])
        res = np.array([abs(A[j]) for j in peaks])
        if len(phis) < 4:
            return {"n_peaks": len(phis), "verdict": "判定不能", "Qz2": None,
                     "C_max": float(C[10:].max())}
        d0 = np.abs(phis)
        dpi = np.minimum(np.abs(phis - np.pi), np.abs(phis + np.pi))
        near = np.minimum(d0, dpi) < 0.1 * np.pi
        Qz2 = float(np.mean(near))
        signs = np.sign(np.cos(phis[near])) if near.any() else np.array([])
        if Qz2 > 0.9 and (signs < 0).any() and (signs > 0).any():
            verdict = "被覆度2(Z₂)"
        elif Qz2 > 0.9 and (signs > 0).all():
            verdict = "被覆度1"
        else:
            verdict = "連続(被覆なし)"
        return {"n_peaks": len(phis), "Qz2": Qz2, "verdict": verdict,
                 "C_max": float(C[10:].max()),
                 "phis_over_pi": [float(x / np.pi) for x in phis],
                 "absA": [float(x) for x in res]}

    out = {"J_MAX": J_MAX, "cells": {}}
    print(f"{'χパリティ':>8} {'m':>3} {'ピーク数':>5} {'Qz2':>6} {'判定':>14}")
    for pname, bins in (("偶(F分類)", (30, 32, 34)), ("奇(B分類)", (29, 31, 33))):
        for m_t in (0, 1, 2, 4):
            spc = species(bins, m_t)
            if spc is None:
                print(f"{pname:>8} {m_t:>3}  構成不能"); continue
            r = measure(*spc)
            out["cells"][f"{pname}|m={m_t}"] = r
            q = f"{r['Qz2']:.2f}" if r["Qz2"] is not None else "  —"
            print(f"{pname:>8} {m_t:>3} {r['n_peaks']:>5} {q:>6} {r['verdict']:>14}")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_spin_statistics_cross_result_v14.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_spin_statistics_cross_result_v14.json")

if __name__ == "__main__":
    main()
