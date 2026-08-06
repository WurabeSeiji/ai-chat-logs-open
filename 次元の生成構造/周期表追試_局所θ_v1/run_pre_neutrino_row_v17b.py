#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v17b正式版: ν行の決定検定——純m=0・フェルミオン帯種の被覆度
構成: PJ{1}(SW(bundle(29,31,33)))を正規化後 e^{-iη} シフト（純m=0・偶χ周波数帯）。
（v17の教訓: 束は毛m=+1を厳密1単位保有、m=0は射影でなくシフトで作る。
本結果は当初ヒアドキュメント実行——正式スクリプト化・再現一致検証済み）
判定: Qz2>0.9かつ両符号→被覆度2＝ν行成立（中性フェルミオン許容）。"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v17b", UIM / "run_ignition_fate_exact_v3.py")
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
    eta = 2 * np.pi * np.arange(ne) / ne
    def SW(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho"); f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def PJ(v, mset):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(mset)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)
    def SH(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)
    a = PJ(SW(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0)) * S, {1})
    b = PJ(SW(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0)) * S, {1})
    pw = float(np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2))
    a *= (64.0 / pw) ** 0.5; b *= (64.0 / pw) ** 0.5
    a = SH(a, -1); b = SH(b, -1)
    f = np.fft.fft(a.reshape(shape), axis=1, norm="ortho")
    Pm = np.sum(np.abs(f) ** 2, axis=0)
    occ = [(int(mm[i]), float(Pm[i] / Pm.sum())) for i in range(ne) if Pm[i] / Pm.sum() > 1e-6]
    print("ν候補のη占有:", occ)
    a0_, b0_ = a.copy(), b.copy()
    s0 = np.imag(np.conj(b0_) * a0_); s0c = s0 - s0.mean()
    nrm = float(np.sum(np.abs(a0_) ** 2) + np.sum(np.abs(b0_) ** 2))
    A, C = [], []
    for j in range(J_MAX):
        a, b, _ = ex.collision_step_exact(a, b, sp)
        A.append(complex((np.vdot(a0_, a) + np.vdot(b0_, b)) / nrm))
        sj = np.imag(np.conj(b) * a); sjc = sj - sj.mean()
        C.append(float(np.dot(s0c, sjc) / max(np.linalg.norm(s0c) * np.linalg.norm(sjc), 1e-300)))
    A = np.array(A); C = np.array(C)
    th = 0.85 * float(C[10:].max())
    peaks = []
    for j in range(11, J_MAX - 1):
        if C[j] > th and C[j] >= C[j - 1] and C[j] >= C[j + 1]:
            if not peaks or j - peaks[-1] > 5: peaks.append(j)
        if len(peaks) >= 15: break
    phis = np.array([np.angle(A[j]) for j in peaks])
    d0 = np.abs(phis); dpi = np.minimum(np.abs(phis - np.pi), np.abs(phis + np.pi))
    near = np.minimum(d0, dpi) < 0.1 * np.pi
    Qz2 = float(np.mean(near))
    signs = np.sign(np.cos(phis[near])) if near.any() else np.array([])
    if Qz2 > 0.9 and (signs < 0).any() and (signs > 0).any(): v = "被覆度2(Z₂)"
    elif Qz2 > 0.9 and (signs > 0).all(): v = "被覆度1"
    else: v = "連続"
    print(f"ν候補（純m=0・フェルミオン帯）: ピーク{len(phis)} Qz2={Qz2:.2f} → {v}")
    h_nu = bool(v == "被覆度2(Z₂)")
    print(f"H_nu（ν行成立）= {h_nu}")
    out = {"J_MAX": J_MAX, "eta_occupancy": occ, "n_peaks": len(phis), "Qz2": Qz2,
           "phis_over_pi": [float(x / np.pi) for x in phis],
           "peaks_J": [int(j + 1) for j in peaks],
           "verdict": v, "H_nu": h_nu, "runtime_sec": time.time() - t0}
    (HERE / "pre_neutrino_row_result_v17b.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
