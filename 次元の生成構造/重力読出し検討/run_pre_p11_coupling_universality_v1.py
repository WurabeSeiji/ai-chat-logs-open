#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P11: 結合普遍性検定——時間ゲージ歪みの裾係数は単一の質量測度に潰れるか

設計（事前記録・帰属規約）: P10の空間版万能演算で、塊の幅σと振幅Mを独立に
掃引し、δτ_t の裾係数 A を測る。候補の質量測度
  m1 = M²σ（パワー総和ノルム）  m2 = Mσ（振幅総和）  m3 = M（ピーク）
のどれかで A/m が (M,σ) 全格子で単一定数に潰れれば、それがこのモデルの
「質量」の定義であり、比 K=A/m が結合定数。
**帰属規約（木原訂正の反映）**: どの測度でも潰れない場合、反証されるのは
「θ局所化のこの実装」であって理論枠組みではない。別の局所化（マスク形・
射影核）を探す指針として散らばりのパターンを記録する。

測定: T=50（伝播変形前の最清プロファイル）。裾窓 dx∈[20,170]（対蹠跳ねと
台を除外）で |δτ_t|=A/x^p をフィット。参照点 dx=30,60 の生値も併記。
弱結合域に限定: M ∈ {0.05,0.08,0.12} × σ ∈ {2,3,5,8}。
判定: 各測度の K=A/m の変動係数 CV=std/mean を比較、CV<10%なら潰れたと判定。
使い方: python3 run_pre_p11_coupling_universality_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p11", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base

K_SEA = 3

def make_masks(n):
    k = np.rint(np.fft.fftfreq(n, d=1.0 / n)).astype(int)
    mf = (np.abs(k) >= 4) & (np.abs(k) % 2 == 0)
    return mf, ~mf

def local_theta(a2, b2, mf, mb):
    Fa = np.fft.fft(a2, axis=0, norm="ortho"); Fb = np.fft.fft(b2, axis=0, norm="ortho")
    f_loc = (np.sum(np.abs(np.fft.ifft(Fa * mf[:, None], axis=0, norm="ortho")) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * mf[:, None], axis=0, norm="ortho")) ** 2, axis=1))
    b_loc = (np.sum(np.abs(np.fft.ifft(Fa * mb[:, None], axis=0, norm="ortho")) ** 2, axis=1)
             + np.sum(np.abs(np.fft.ifft(Fb * mb[:, None], axis=0, norm="ortho")) ** 2, axis=1))
    return np.arctan2(np.sqrt(f_loc), np.sqrt(b_loc + 1e-300))

def step_local(a2, b2, mf, mb):
    th = local_theta(a2, b2, mf, mb)
    c, s_ = np.cos(th)[:, None], np.sin(th)[:, None]
    a2, b2 = c * a2 - s_ * b2, s_ * a2 + c * b2
    phi = 2.0 * (np.sin(th) ** 2)[:, None] * np.imag(np.conj(b2) * a2)
    cp, sp_ = np.cos(phi), np.sin(phi)
    return cp * a2 - sp_ * b2, sp_ * a2 + cp * b2

def run_tau_t(n, ne, x, dxv, mass, sigma, T, mf, mb):
    sea_a = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    lump = mass * np.exp(-0.5 * (dxv / sigma) ** 2)[:, None] * np.ones((1, ne))
    a2 = (sea_a + lump).astype(complex)
    b2 = (sea_a * np.exp(1j * np.pi / 4)).astype(complex)
    for j in range(T):
        a_prev = a2.copy()
        a2, b2 = step_local(a2, b2, mf, mb)
    return np.angle(np.einsum("xe,xe->x", np.conj(a_prev), a2))

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    mf, mb = make_masks(n)
    T = 50
    tt0 = run_tau_t(n, ne, x, dxv, 0.0, 3.0, T, mf, mb)

    fitw = (dxv >= 20) & (dxv <= 170)
    rows = []
    print(f"{'M':>5} {'σ':>3} | {'A(fit)':>10} {'p':>6} | {'|δτ_t|@30':>10} {'@60':>10}")
    for MASS in (0.05, 0.08, 0.12):
        for SIG in (2.0, 3.0, 5.0, 8.0):
            dtt = run_tau_t(n, ne, x, dxv, MASS, SIG, T, mf, mb) - tt0
            v = np.abs(dtt[fitw]); d = dxv[fitw]
            ok = v > 1e-14
            slope, intc = np.polyfit(np.log(d[ok]), np.log(v[ok]), 1)
            A = float(np.exp(intc)); p = float(-slope)
            v30 = float(np.mean(np.abs(dtt[(dxv >= 27) & (dxv <= 33)])))
            v60 = float(np.mean(np.abs(dtt[(dxv >= 55) & (dxv <= 65)])))
            rows.append({"M": MASS, "sig": SIG, "A": A, "p": p, "v30": v30, "v60": v60})
            print(f"{MASS:>5} {SIG:>3.0f} | {A:>10.3e} {p:>6.3f} | {v30:>10.3e} {v60:>10.3e}")

    print("\n== 規格化検定（K=測定量/質量測度, CV=std/mean） ==")
    for name, mfun in (("m1=M²σ", lambda r: r["M"]**2 * r["sig"]),
                       ("m2=Mσ", lambda r: r["M"] * r["sig"]),
                       ("m3=M", lambda r: r["M"]),
                       ("m4=M²σ²", lambda r: (r["M"] * r["sig"])**2 / r["sig"] * r["sig"])):
        for qname in ("v30", "v60", "A"):
            Ks = np.array([r[qname] / mfun(r) for r in rows])
            cv = np.std(Ks) / np.mean(Ks)
            flag = " ←潰れ" if cv < 0.10 else ""
            print(f"  {name:>9} / {qname:>3}: K平均={np.mean(Ks):.3e} CV={cv:.1%}{flag}")

    (HERE / "result_pre_p11_coupling_universality_v1.json").write_text(
        json.dumps(rows, indent=1), encoding="utf-8")
    print(f"\n完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
