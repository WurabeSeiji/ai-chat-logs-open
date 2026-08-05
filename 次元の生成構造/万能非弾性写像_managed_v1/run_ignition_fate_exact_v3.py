#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""点火後の運命 厳密解版 v3——頂点の閉形式解の発見と検証

発見（2026-08-05、査読対応の検証中）: 頂点流
    da/dτ = i(|b|²a − b²ā) = −2s·b,  db/dτ = +2s·a,  s = Im(b̄a)
は s を厳密に保存する（ds/dτ = 2s·Im(|a|²) − 2s·Im(|b|²) = 0）。
ゆえに流れは各格子点で角速度 2s の厳密な回転であり、閉形式
    a' = cosφ·a − sinφ·b,  b' = sinφ·a + cosφ·b,  φ(x) = 2R·Im(b̄a)(x)
で厳密可解。閉塞 a²+b² とパワー |a|²+|b|² は点ごとに恒等保存（積分器不要）。

判定（実行前固定）:
    E1 厳密性: 3000衝突の零閉塞ドリフト・ノルムドリフトが ≤1e-12（機械精度）。
    E2 物理の同一性: f*（後半500平均）が RK4 値 0.4690 の ±0.005。
    E3 点火則の同一性: C が RK4 値 11.45 の ±2%。

使い方: python3 run_ignition_fate_exact_v3.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("v3_exact", HERE / "universal_inelastic_map_v3.py")
v3 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = v3
spec.loader.exec_module(v3)
v1, toy, base = v3.v1, v3.toy, v3.base

J = 3000
S = 8.0
SEED_AMP = 0.1


def collision_step_exact(a, b, sp):
    """弾性部（無変更）＋非弾性部の閉形式厳密解。"""
    ro = toy.theta_from_ab(a, b, sp)
    a, b = toy.rotate_ab(a, b, ro.theta)
    r = float(ro.reflection_rate)
    if r > 0.0:
        phi = 2.0 * r * np.imag(np.conj(b) * a)
        c, s_ = np.cos(phi), np.sin(phi)
        a, b = c * a - s_ * b, s_ * a + c * b
    return a, b, ro


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)

    # ---- E1+E2: 運命 ----
    a = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    a = a + v1.make_bundle(sp, v1.ODD_KS, "A", scale=SEED_AMP * S)
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    c0 = abs(complex(np.sum(a * a) + np.sum(b * b)))
    fs, norms = [], []
    for _ in range(J):
        a, b, _ = collision_step_exact(a, b, sp)
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        fs.append(v1.fermionic_power_raw(a, b, sp) / tot)
        norms.append(tot / tot0)
    fs = np.array(fs)
    c1 = abs(complex(np.sum(a * a) + np.sum(b * b)))
    closure_drift = abs(c1 - c0)
    norm_drift = float(max(abs(n - 1) for n in norms))
    f_eq = float(fs[-500:].mean())
    imax = int(fs.argmax())
    e1 = bool(closure_drift <= 1e-12 * max(tot0, 1.0) and norm_drift <= 1e-12)
    e2 = bool(abs(f_eq - 0.4690) <= 0.005)
    print(f"運命(厳密解): f最大={fs.max():.4f}(j={imax}) f*={f_eq:.4f}")
    print(f"  閉塞ドリフト={closure_drift:.2e} ノルムドリフト={norm_drift:.2e} → E1={e1}")
    print(f"  f*: 厳密解 {f_eq:.4f} vs RK4 0.4690 → E2={e2}")

    # ---- E3: C 掃引 ----
    rows = []
    for seed_amp in (1e-3, 1e-2, 1e-1):
        a2 = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
        a2 = a2 + v1.make_bundle(sp, v1.ODD_KS, "A", scale=seed_amp * S)
        b2 = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
        f_series = []
        for _ in range(200):
            a2, b2, _ = collision_step_exact(a2, b2, sp)
            tot = float(np.vdot(a2, a2).real + np.vdot(b2, b2).real)
            f_series.append(v1.fermionic_power_raw(a2, b2, sp) / tot)
        f_series = np.array(f_series)
        f0 = f_series[0]
        lnP = np.log(f_series)
        rate0 = float((lnP[20] - lnP[0]) / 20.0)
        C = rate0 / f0 ** 2
        rows.append({"seed_amp": seed_amp, "f0": float(f0), "C": float(C)})
        print(f"  掃引 seed_amp={seed_amp:.0e}: f0={f0:.3e} C={C:.3f}")
    Cs = [r_["C"] for r_ in rows]
    e3 = bool(all(abs(c_ / 11.45 - 1) <= 0.02 for c_ in Cs))
    print(f"  C = {[round(c_,3) for c_ in Cs]}（RK4: 11.45） → E3={e3}")

    out = {"criteria": {"E1": "drift <= machine (1e-12)", "E2": "f* = 0.4690±0.005",
                         "E3": "C = 11.45±2%"},
           "closure_drift": float(closure_drift), "norm_drift": norm_drift,
           "f_max": float(fs.max()), "f_equilibrium": f_eq, "C_sweep": rows,
           "E1": e1, "E2": e2, "E3": e3, "all_pass": bool(e1 and e2 and e3),
           "f_series_every10": [float(x) for x in fs[::10]],
           "runtime_sec": time.time() - t0}
    print(f"\n判定: {'ALL PASS——頂点は閉形式で厳密可解' if out['all_pass'] else '不成立あり'}")
    (HERE / "ignition_fate_exact_result_v3.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
