#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2実験1: 偶奇非対称性の検定（新柱6・7の実測正本）

判定（事前固定）:
(A) 等振幅初期条件: フェルミオン分率 f0=0.500000 に二分法で厳密調整した束を
    公開系（大域θ）で3000衝突進化。f=0.5 が不動点なら偶奇結合は対称（新柱6棄却）。
    不動点でなければ非対称の実証。末500平均±stdを記録。
(B) 半並進射影: 全倍音・中心位相揃えの梯子（σ_k=32）に P±=(1±T_{n/2})/2 を適用。
    偶成分の奇kリーク／奇成分の偶kリークが機械精度なら「等振幅でも厳密分離可能」
    （区別子は振幅でなく変換応答）の実証。
使い方: python3 run_pre_v2_parity_asymmetry_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v2p", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, base = ex.v1, ex.base

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    out = {}

    # (A) 等振幅流出
    a_e = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=1.0)
    a_o = v1.make_bundle(sp, v1.ODD_KS, "A", scale=1.0)
    b_e = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=1.0)
    b_o = v1.make_bundle(sp, v1.ODD_KS, "B", scale=1.0)
    def f_of(w):
        a = a_e + w * a_o; b = b_e + w * b_o
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        return v1.fermionic_power_raw(a, b, sp) / tot
    lo, hi = 0.0, 5.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f_of(mid) > 0.5: hi = mid
        else: lo = mid
    w = 0.5 * (lo + hi)
    a = a_e + w * a_o; b = b_e + w * b_o
    f0 = f_of(w)
    fs = []
    for j in range(3000):
        a, b, _ = ex.collision_step_exact(a, b, sp)
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        fs.append(v1.fermionic_power_raw(a, b, sp) / tot)
    fs = np.array(fs)
    f_star = float(fs[-500:].mean()); f_std = float(fs[-500:].std())
    verdict_A = "非対称（f=0.5は不動点でない）" if abs(f_star - 0.5) > 3 * f_std / np.sqrt(50) else "判定不能"
    print(f"(A) f0={f0:.6f} (w={w:.4f}) → f*={f_star:.4f}±{f_std:.4f}  {verdict_A}")
    out["A"] = {"w": w, "f0": f0, "f_star": f_star, "f_std": f_std,
                "fs_decimated": fs[::5].tolist(), "decimation": 5, "verdict": verdict_A}

    # (B) 半並進射影の厳密分離
    x = np.arange(n); x_L = n // 2
    ks = np.arange(1, n // 2)
    Wk = np.exp(-0.5 * (ks / 32.0) ** 2)
    prof = np.zeros(n)
    for kk, wgt in zip(ks, Wk):
        prof += wgt * np.cos(2 * np.pi * kk * (x - x_L) / n)
    prof /= np.max(np.abs(prof))
    T = np.roll(prof, n // 2)
    even = 0.5 * (prof + T); odd = 0.5 * (prof - T)
    Fe = np.fft.fft(even); Fo = np.fft.fft(odd)
    leak_e = float(np.max(np.abs(Fe[1::2])) / np.max(np.abs(Fe)))
    leak_o = float(np.max(np.abs(Fo[0::2][1:])) / np.max(np.abs(Fo)))
    Pe = float(np.sum(even ** 2) / np.sum(prof ** 2))
    Po = float(np.sum(odd ** 2) / np.sum(prof ** 2))
    print(f"(B) パワー分割 偶={Pe:.4f} 奇={Po:.4f}  リーク e={leak_e:.1e} o={leak_o:.1e}")
    out["B"] = {"P_even": Pe, "P_odd": Po, "leak_even": leak_e, "leak_odd": leak_o,
                "verdict": "厳密分離（等振幅でも変換応答で区別可能）" if max(leak_e, leak_o) < 1e-12 else "リークあり"}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_v2_parity_asymmetry_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_v2_parity_asymmetry_result_v1.json")

if __name__ == "__main__":
    main()
