#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G6正本: 時計読出しの双線形署名（グラビトン節・double copy）

設計（事前記録）: リンギング源（σ_k=8梯子・外部駆動なし=無名性保持）で
場（線形量: 摂動形状への複素オーバーラップ）と時計場 τ_t の分光を比較。
判定（事前固定）:
(i) 場の主要2線 ν₁, ν₂ に対し、時計場の主線が |ν₁−ν₂|（差周波数）に一致
    （相対誤差<10%）。
(ii) 時計スペクトルにおいて場自身の線（ν₁, ν₂）が主線より抑制されている。
→ 双方成立で「時計読出し=場の双線形（重力振幅=場×場・double copy型）」の
  直接実証。素朴な単線2ν検定は多線源では不適（設計注記）。
使い方: python3 run_g6_bilinear_spectroscopy_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g6", UIM / "run_ignition_fate_exact_v3.py")
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

    def ladder_at(fsrc, center, sig_k, amp=0.05):
        prof = np.zeros(n); AMP = amp * np.sqrt(n) * 0.1
        for parity, wgt in (("even", np.sqrt(fsrc)), ("odd", np.sqrt(1 - fsrc))):
            ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
            Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
            sub = np.zeros(n)
            for kk, w in zip(ks, Wk):
                sub += w * np.cos(2 * np.pi * kk * (x - center) / n)
            prof += wgt * sub / np.sqrt(np.sum(sub ** 2))
        return AMP * prof * np.exp(2j * np.pi * K_SEA * center / n)

    c = 100
    lump = ladder_at(0.6, c, sig_k=8.0)
    a2 = ((sea + lump)[:, None] * np.ones((1, ne))).astype(complex)
    b2 = -1j * a2
    dxv = np.minimum(np.abs(x - c), n - np.abs(x - c))
    supp = dxv <= 12
    shape = lump.copy(); nrm = float(np.sum(np.abs(shape) ** 2)) * ne
    T = 1024
    Of = []; cn = []
    for j in range(T):
        ap = a2.copy(); a2, b2 = step(a2, b2)
        Of.append(complex(np.sum(np.conj(shape[:, None]) * (a2 - sea[:, None]))) / nrm)
        tt = np.angle(np.einsum("xe,xe->x", np.conj(ap), a2))
        cn.append(float(np.mean(tt[supp])))
    Of = np.array(Of); cn = np.array(cn)

    def spectrum(sig, nburn=200):
        s = sig[nburn:] - np.mean(sig[nburn:])
        S = np.abs(np.fft.rfft(s * np.hanning(len(s))))
        fr = np.arange(len(S)) * 2 * np.pi / len(s)
        return fr, S

    fr_f, S_f = spectrum(np.real(Of * np.exp(-1j * np.angle(Of[300]))))
    fr_c, S_c = spectrum(cn)
    # 場の速い2線（ν>0.05 帯で上位2）
    fast = fr_f > 0.05
    idx = np.argsort(S_f[fast])[::-1][:2]
    nus = sorted(fr_f[fast][idx])
    dnu = abs(nus[1] - nus[0])
    # 時計の主線
    ipk = int(np.argmax(S_c[1:])) + 1
    nu_clock = float(fr_c[ipk])
    err = abs(nu_clock - dnu) / dnu if dnu > 0 else 1.0
    # 場の線の時計側での抑制
    i1 = int(np.argmin(np.abs(fr_c - nus[0]))); i2 = int(np.argmin(np.abs(fr_c - nus[1])))
    supp_ratio = float(max(S_c[i1], S_c[i2]) / S_c[ipk])
    ok = err < 0.10 and supp_ratio < 0.5
    print(f"場の速い2線: ν₁={nus[0]:.4f} ν₂={nus[1]:.4f}  差={dnu:.4f}")
    print(f"時計主線: {nu_clock:.4f}（差との相対誤差 {err:.1%}）")
    print(f"場の線の時計側抑制: S_clock(ν_field)/S_clock(主線)={supp_ratio:.3f}")
    verdict = "柱G6成立（時計読出し=場の双線形・差周波数櫛）" if ok else "要精査"
    print(verdict)
    out = {"nu_field": nus, "delta_nu": dnu, "nu_clock_main": nu_clock,
           "rel_err": err, "field_line_suppression": supp_ratio,
           "spec_field": {"fr": fr_f[::2].tolist(), "S": S_f[::2].tolist()},
           "spec_clock": {"fr": fr_c[::2].tolist(), "S": S_c[::2].tolist()},
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g6_bilinear_spectroscopy_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
