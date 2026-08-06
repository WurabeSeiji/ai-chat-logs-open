#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P10: 空間版万能演算——θ読出しの局所射影化によるゲージ局所歪みの創生

設計（木原指示の実装・事前記録）: 波の万能演算の空間版。空間を創生する
射影関数（θ読出し）に万能性を持たせ、ゲージの局所歪みを関係量から作り込む。
空間は自発的に生じているため無名性を破らない。

実装: 大域θ（全格子スペクトル和の1スカラー）を局所場 θ(x) に昇格する。
  フェルミオン帯 = 偶数かつ|k|≥4（元のマスクと同一・追加パラメータなし）
  a_f = バンドパス(a), a_bos = 補帯域(a)（bも同様）
  f_loc(x) = Σ_η |a_f|²+|b_f|²,  b_loc(x) = Σ_η |a_bos|²+|b_bos|²
  θ(x) = atan2(√f_loc, √b_loc),  r(x) = sin²θ(x)
  ステップ: セルごとSO(2)回転（θ(x)）→ 点ごと非弾性 φ=2r(x)·Im(b̄a)。
閉塞 C=Σ(a²+b²) はセルごと実直交回転で厳密保存（構成的）。

判定（事前固定）:
(i) 閉塞ドリフト ≤ 1e-12（ゲージ側操作の資格）。
(ii) 対照（純海）: f_loc=0 → θ=0 → 固定点維持（真空は時を刻まない）。
(iii) 塊あり: δτ_t(x) が塊の台（ガウスσ=3, dx≳15で初期寄与消滅）の外で
     距離依存の裾を持つか。持てば log-log 勾配で減衰指数 p を測る。
     射影核（急峻マスク→Dirichlet型1/x裾）から p≈2 (θ²経由) が理論候補。
(iv) 質量スケーリング: 裾振幅 vs MASS。
使い方: python3 run_pre_p10_local_universal_op_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_p10", UIM / "run_ignition_fate_exact_v3.py")
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
    af = np.fft.ifft(Fa * mf[:, None], axis=0, norm="ortho")
    bf = np.fft.ifft(Fb * mf[:, None], axis=0, norm="ortho")
    ab_ = np.fft.ifft(Fa * mb[:, None], axis=0, norm="ortho")
    bb_ = np.fft.ifft(Fb * mb[:, None], axis=0, norm="ortho")
    f_loc = np.sum(np.abs(af) ** 2 + np.abs(bf) ** 2, axis=1)
    b_loc = np.sum(np.abs(ab_) ** 2 + np.abs(bb_) ** 2, axis=1)
    return np.arctan2(np.sqrt(f_loc), np.sqrt(b_loc + 1e-300))

def step_local(a2, b2, mf, mb):
    th = local_theta(a2, b2, mf, mb)
    c, s_ = np.cos(th)[:, None], np.sin(th)[:, None]
    a2, b2 = c * a2 - s_ * b2, s_ * a2 + c * b2
    r = (np.sin(th) ** 2)[:, None]
    phi = 2.0 * r * np.imag(np.conj(b2) * a2)
    cp, sp_ = np.cos(phi), np.sin(phi)
    a2, b2 = cp * a2 - sp_ * b2, sp_ * a2 + cp * b2
    return a2, b2

def gauges(P0, P1):
    ip_x = np.einsum("xe,xe->x", np.conj(P1), np.roll(P1, -1, axis=0))
    tau_x = np.angle(ip_x)
    tau_t = np.angle(np.einsum("xe,xe->x", np.conj(P0), P1))
    return tau_x, tau_t

def run(n, ne, x, dxv, mass, T, mf, mb):
    sea_a = 0.2 * np.exp(2j * np.pi * K_SEA * x / n)[:, None] * np.ones((1, ne))
    lump = mass * np.exp(-0.5 * (dxv / 3.0) ** 2)[:, None] * np.ones((1, ne))
    a2 = (sea_a + lump).astype(complex)
    b2 = (sea_a * np.exp(1j * np.pi / 4)).astype(complex)
    C0 = complex(np.sum(a2 * a2) + np.sum(b2 * b2))
    for j in range(T):
        a_prev = a2.copy()
        a2, b2 = step_local(a2, b2, mf, mb)
    dC = abs(complex(np.sum(a2 * a2) + np.sum(b2 * b2)) - C0)
    return (*gauges(a_prev, a2), dC)

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
    mf, mb = make_masks(n)
    results = {}

    for T in (50, 200):
        tx0, tt0, dC0 = run(n, ne, x, dxv, 0.0, T, mf, mb)
        print(f"\n=== T={T} ===")
        print(f"[対照] τ_t: 平均={np.mean(tt0):+.2e} std={np.std(tt0):.1e} 閉塞ドリフト={dC0:.2e}")
        for MASS in (0.05, 0.1, 0.2, 0.4):
            tx, tt, dC = run(n, ne, x, dxv, MASS, T, mf, mb)
            dtt = tt - tt0; dtx = tx - tx0
            key = f"T{T}_M{MASS}"
            pr = {}
            for dd in (0, 6, 12, 20, 30, 45, 70, 110, 170, 240):
                sel = (dxv >= max(dd - 3, 0)) & (dxv <= dd + 3)
                pr[dd] = float(np.mean(dtt[sel]))
            results[key] = {"d_tau_t": pr, "dC": dC,
                            "d_tau_x_lump": float(np.mean(dtx[dxv <= 3])),
                            "d_tau_x_far": float(np.mean(np.abs(dtx[dxv > 60])))}
            print(f"[M={MASS}] 閉塞ドリフト={dC:.2e}")
            print("   δτ_t: " + "  ".join(f"dx{dd}={pr[dd]:+.2e}" for dd in (0, 12, 20, 30, 45, 70, 110, 170, 240)))
            # 減衰指数（台の外 20..240 で |δτ_t| の log-log 勾配）
            dds = np.array([20, 30, 45, 70, 110, 170, 240])
            vals = np.array([abs(pr[d]) for d in dds])
            if np.all(vals > 1e-14):
                p_exp = np.polyfit(np.log(dds), np.log(vals), 1)[0]
                print(f"   減衰指数 p(20-240) = {p_exp:+.3f}   δτ_x: 塊上={results[key]['d_tau_x_lump']:+.2e} 遠方|平均|={results[key]['d_tau_x_far']:.2e}")
            else:
                print(f"   台外で機械零（裾なし）  δτ_x: 塊上={results[key]['d_tau_x_lump']:+.2e}")
    out = HERE / "result_pre_p10_local_universal_op_v1.json"
    out.write_text(json.dumps(results, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"\n保存: {out.name}  完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
