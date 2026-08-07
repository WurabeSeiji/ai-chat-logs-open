#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2実験2: 線幅・寿命関係の検定（新柱8の実測正本）

エンジン: 局所化プロトタイプ（滑らかIRロールオフ×パリティ分割・円偏波 b=−ia）。
整合性バッテリー（周期表追試_局所θ_v1）で整数・位相構造の頑健性を確認済みの系。
定義: σ_ω=台上の時間平均時計レートの空間std（線幅）／τ_coh=自己重なりC=1/2到達時間
（寿命）／質量=台上|⟨τ_t⟩|平均。
判定（事前固定）: τ_coh·σ_ω が構成横断でオーダー1定数（CV<50%）なら線幅・寿命
関係成立。σ_ωと質量の相関符号も記録。打ち切り（C>0.5のままTobs終了）は除外を明示。
使い方: python3 run_pre_v2_linewidth_lifetime_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v2l", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base
K_SEA = 3

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n); x_L = n // 2
    dxv = np.minimum(np.abs(x - x_L), n - np.abs(x - x_L))
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

    def parity_ladder(sig_k, parity):
        ks = np.arange(4 if parity == "even" else 5, n // 2, 2)
        Wk = np.exp(-0.5 * (ks / sig_k) ** 2)
        prof = np.zeros(n)
        for kk, w in zip(ks, Wk):
            prof += w * np.cos(2 * np.pi * kk * (x - x_L) / n)
        return prof / np.sqrt(np.sum(prof ** 2))

    def run_cfg(fsrc, amp_scale, sig_k, Tburn=200, Tobs=800):
        Lf = parity_ladder(sig_k, "even"); Lb = parity_ladder(sig_k, "odd")
        AMP = amp_scale * np.sqrt(n) * 0.1
        lump = AMP * (np.sqrt(fsrc) * Lf + np.sqrt(1 - fsrc) * Lb)
        supp = dxv <= max(6, int(3 * n / (2 * np.pi * sig_k)))
        a2 = ((0.2 * carrier + lump)[:, None] * np.ones((1, ne))).astype(complex)
        b2 = -1j * a2
        ctrl = (0.2 * carrier[:, None] * np.ones((1, ne))).astype(complex)
        cb = -1j * ctrl
        snaps = []; ref = None; C = []
        for j in range(Tburn + Tobs):
            ap = a2.copy(); a2, b2 = step(a2, b2)
            ctrl, cb = step(ctrl, cb)
            if j >= Tburn:
                snaps.append(np.angle(np.einsum("xe,xe->x", np.conj(ap), a2)))
                da = (a2 - ctrl)[supp, :]
                if ref is None:
                    ref = da.copy(); nr = float(np.sum(np.abs(ref) ** 2))
                C.append(abs(complex(np.sum(np.conj(ref) * da))) /
                         max(np.sqrt(nr * float(np.sum(np.abs(da) ** 2))), 1e-300))
        S = np.array(snaps); mean_tau = np.mean(S, axis=0)
        sigma_om = float(np.std(mean_tau[supp]))
        C = np.array(C); below = np.where(C < 0.5)[0]
        tau_coh = int(below[0]) if len(below) else None
        mass = float(np.mean(np.abs(mean_tau[supp])))
        return sigma_om, tau_coh, mass, float(C[-1])

    cfgs = [(0.3, 0.05, 32), (0.5, 0.05, 32), (0.7, 0.05, 32), (0.9, 0.05, 32),
            (1.0, 0.05, 32), (0.7, 0.10, 32), (0.7, 0.05, 64), (0.5, 0.10, 64)]
    rows = []
    for fsrc, amp, sk in cfgs:
        so, tc, m, cl = run_cfg(fsrc, amp, sk)
        rows.append({"f_src": fsrc, "amp": amp, "sig_k": sk, "sigma_omega": so,
                     "tau_coh": tc, "mass": m, "C_final": cl,
                     "product": (tc * so) if tc else None})
        print(f"f={fsrc} amp={amp} σk={sk:.0f}: σ_ω={so:.3e} τ={tc} 積={tc*so if tc else '—'} 質量={m:.3e}")
    prods = np.array([r["product"] for r in rows if r["product"] is not None])
    cv = float(prods.std() / prods.mean())
    verdict = f"線幅・寿命関係成立（積平均={prods.mean():.1f} CV={cv:.1%}<50%）" if cv < 0.5 else "不成立"
    sos = [r["sigma_omega"] for r in rows]; ms = [r["mass"] for r in rows]
    corr = float(np.corrcoef(sos, ms)[0, 1])
    print(f"積 τ·σ_ω: 平均={prods.mean():.2f} CV={cv:.1%} → {verdict}")
    print(f"質量⟨ω⟩とσ_ωの相関: {corr:+.3f}")
    out = {"rows": rows, "product_mean": float(prods.mean()), "product_cv": cv,
           "n_censored": sum(1 for r in rows if r["product"] is None),
           "mass_sigma_corr": corr, "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "pre_v2_linewidth_lifetime_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_v2_linewidth_lifetime_result_v1.json")

if __name__ == "__main__":
    main()
