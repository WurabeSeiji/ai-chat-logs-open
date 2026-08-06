#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v17: ν行の白黒——m=0×フェルミオン帯（偶χ周波数）種の被覆度測定

背景: v14でSW済み束のm=0射影は空。非SW束(29,31,33=偶χ周波数=フェルミオン帯)の
m=0射影は存在する(v7海)が、v13bで連続位相だった。SWの有無が被覆度を決めるのか、
m=0が決めるのかを分離する。
ケース（事前固定）: A=非SW偶周波数束m=1／B=非SW偶周波数束m=0（ν候補）／
C=SW偶周波数束m=1（v14でZ₂確認済の対照）。
判定 H_nu: BがZ₂（Qz2>0.9・両符号）→ ν行成立（m=0フェルミオンあり）。
Bが連続かつAがZ₂ → m=0が被覆を壊す＝ν行不成立（中性フェルミオン禁止）。
AもBも連続 → SWが被覆の必要条件（構成の問題、ν判定は保留）。
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v17", UIM / "run_ignition_fate_exact_v3.py")
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

    eta = 2 * np.pi * np.arange(ne) / ne
    def shift_eta(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)

    def sp_case(sw, m_t):
        # 束は構成上 m=+1 の毛を持つ（実測）→ 目標巻き m_t へは射影でなくシフトで移す
        a = v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S
        b = v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S
        if sw:
            a = single_winding(a); b = single_winding(b)
        a = shift_eta(a, m_t - 1); b = shift_eta(b, m_t - 1)
        pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
        if pw < 1e-9: return None
        sc = np.sqrt(64.0 / pw)
        return a * sc, b * sc

    def covering(a, b):
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
                if not peaks or j - peaks[-1] > 5: peaks.append(j)
            if len(peaks) >= 12: break
        phis = np.array([np.angle(A[j]) for j in peaks])
        if len(phis) < 4: return "判定不能", None, len(phis)
        d0 = np.abs(phis); dpi = np.minimum(np.abs(phis - np.pi), np.abs(phis + np.pi))
        near = np.minimum(d0, dpi) < 0.1 * np.pi
        Qz2 = float(np.mean(near))
        signs = np.sign(np.cos(phis[near])) if near.any() else np.array([])
        if Qz2 > 0.9 and (signs < 0).any() and (signs > 0).any(): return "被覆度2(Z₂)", Qz2, len(phis)
        if Qz2 > 0.9 and (signs > 0).all(): return "被覆度1", Qz2, len(phis)
        return "連続", Qz2, len(phis)

    out = {"J_MAX": J_MAX, "cases": {}}
    for name, sw, m_t in (("A_非SW偶周波m=1", False, 1), ("B_非SW偶周波m=0(ν候補)", False, 0),
                            ("C_SW偶周波m=1(対照)", True, 1)):
        st = sp_case(sw, m_t)
        if st is None:
            print(f"{name}: 構成不能"); out["cases"][name] = {"verdict": "構成不能"}; continue
        verdict, qz2, npk = covering(*st)
        print(f"{name}: {verdict}（Qz2={qz2}, ピーク{npk}）")
        out["cases"][name] = {"verdict": verdict, "Qz2": qz2, "n_peaks": npk}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_neutrino_row_result_v17.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
