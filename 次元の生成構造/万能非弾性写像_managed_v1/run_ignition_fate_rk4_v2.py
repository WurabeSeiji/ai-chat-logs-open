#!/usr/bin/env python3
"""点火後の運命 RK4版 v2——二体論文修正1: 積分器ドリフトの解消

修正の経緯（2026-08-05）: v1（run_ignition_fate_v1.py）は中点法（RK2）の
部分刻みで、高R域の零閉塞ドリフト 2.9e-3/3000衝突 が残っていた。
多体エンジン（段階2）で RK4 部分刻みが点ごと閉塞 ~1e-11 を実証したため、
同じ RK4 を二体の運命実験に適用して再走行する。
初期条件・シード・パラメータ（S=8, seed_amp=0.1, J=3000, h_max=0.01）は
v1 と同一。変更は積分法のみ。

判定（実行前固定）:
    R1 ドリフト改善: 閉塞ドリフト |ΔC| が v1 の 2.9e-3 より2桁以上改善。
    R2 物理の保持: 統計的平衡 f*（後半500平均）が v1 値 0.469 の ±0.02、
       往復型の判定が保持される。
    R3 点火則の保持: seed割合掃引（1e-3/1e-2/1e-1）の C=rate/f² が
       v1系の 10.4 の ±10% で一定。

使い方: python3 run_ignition_fate_rk4_v2.py
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("v3_rk4", HERE / "universal_inelastic_map_v3.py")
v3 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = v3
spec.loader.exec_module(v3)
v1, toy, base = v3.v1, v3.toy, v3.base

H_MAX = 0.01
J = 3000
S = 8.0
SEED_AMP = 0.1


def collision_step_rk4(a, b, sp, h_max=H_MAX):
    """弾性部（無変更）＋非弾性流の RK4 部分刻み。"""
    ro = toy.theta_from_ab(a, b, sp)
    a, b = toy.rotate_ab(a, b, ro.theta)
    r = float(ro.reflection_rate)
    if r > 0.0:
        peak = float(max(np.max(np.abs(a)) ** 2, np.max(np.abs(b)) ** 2))
        n_sub = max(1, int(math.ceil(r * peak / h_max)))
        rs = r / n_sub
        for _ in range(n_sub):
            k1a, k1b = v3.vertex(a, b, 1.0)
            k2a, k2b = v3.vertex(a + 0.5 * rs * k1a, b + 0.5 * rs * k1b, 1.0)
            k3a, k3b = v3.vertex(a + 0.5 * rs * k2a, b + 0.5 * rs * k2b, 1.0)
            k4a, k4b = v3.vertex(a + rs * k3a, b + rs * k3b, 1.0)
            a = a + (rs / 6.0) * (k1a + 2 * k2a + 2 * k3a + k4a)
            b = b + (rs / 6.0) * (k1b + 2 * k2b + 2 * k3b + k4b)
    return a, b, ro


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)

    # ---- R1+R2: 運命再走行 ----
    a = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    a = a + v1.make_bundle(sp, v1.ODD_KS, "A", scale=SEED_AMP * S)
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    c0 = abs(complex(np.sum(a * a) + np.sum(b * b)))
    fs, norms = [], []
    for _ in range(J):
        a, b, _ = collision_step_rk4(a, b, sp)
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        fs.append(v1.fermionic_power_raw(a, b, sp) / tot)
        norms.append(tot / tot0)
    fs = np.array(fs)
    c1 = abs(complex(np.sum(a * a) + np.sum(b * b)))
    closure_drift = abs(c1 - c0)
    f_eq = float(fs[-500:].mean())
    imax = int(fs.argmax())
    fate = ("往復型（逆変換あり）→ 統計的平衡"
            if imax < J - 100 and fs[-1] < 0.9 * fs.max() else "その他")
    r1 = bool(closure_drift < 2.9e-3 / 100)
    r2 = bool(abs(f_eq - 0.469) <= 0.02 and fate.startswith("往復型"))
    print(f"運命RK4: f最大={fs.max():.4f}(j={imax}) 終値={fs[-1]:.4f} f*={f_eq:.4f}")
    print(f"  閉塞ドリフト={closure_drift:.2e}（v1: 2.9e-3） → R1={r1}")
    print(f"  ノルムドリフト最大={max(abs(n-1) for n in norms):.2e}")
    print(f"  判定={fate}（v1: 往復型 f*≈0.469） → R2={r2}")

    # ---- R3: C一定性の再確認（RK4） ----
    rows = []
    for seed_amp in (1e-3, 1e-2, 1e-1):
        a2 = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
        a2 = a2 + v1.make_bundle(sp, v1.ODD_KS, "A", scale=seed_amp * S)
        b2 = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
        f_series = []
        for _ in range(200):
            a2, b2, _ = collision_step_rk4(a2, b2, sp)
            tot = float(np.vdot(a2, a2).real + np.vdot(b2, b2).real)
            f_series.append(v1.fermionic_power_raw(a2, b2, sp) / tot)
        f_series = np.array(f_series)
        f0 = f_series[0]
        lnP = np.log(f_series)
        rate0 = float((lnP[20] - lnP[0]) / 20.0)
        C = rate0 / f0 ** 2
        rows.append({"seed_amp": seed_amp, "f0": float(f0), "rate0": rate0, "C": float(C)})
        print(f"  掃引 seed_amp={seed_amp:.0e}: f0={f0:.3e} C={C:.3f}")
    Cs = [r_["C"] for r_ in rows]
    spread = (max(Cs) - min(Cs)) / np.mean(Cs)
    r3 = bool(all(abs(c_ / 10.4 - 1) <= 0.10 for c_ in Cs))
    print(f"  C = {[round(c_,3) for c_ in Cs]}（v1系: 10.4） 広がり{spread*100:.1f}% → R3={r3}")

    out = {"criteria": {"R1": "closure drift < 2.9e-5", "R2": "f* within 0.469±0.02, oscillatory",
                         "R3": "C within 10.4±10%"},
           "closure_drift": float(closure_drift),
           "norm_drift_max": float(max(abs(n - 1) for n in norms)),
           "f_max": float(fs.max()), "f_final": float(fs[-1]), "f_equilibrium": f_eq,
           "fate": fate, "C_sweep": rows, "C_spread": float(spread),
           "R1": r1, "R2": r2, "R3": r3, "all_pass": bool(r1 and r2 and r3),
           "f_series_every10": [float(x) for x in fs[::10]],
           "runtime_sec": time.time() - t0}
    print(f"\n判定: {'ALL PASS——修正1完了' if out['all_pass'] else '不成立あり'}")
    (HERE / "ignition_fate_rk4_result_v2.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
