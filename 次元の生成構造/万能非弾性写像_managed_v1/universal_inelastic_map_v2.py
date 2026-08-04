#!/usr/bin/env python3
"""万能非弾性写像 v2：閉塞保存による頂点の一意化

導出（2026-08-04）:
    二波・点ごと・共通U(1)不変・自己散乱なしの三次頂点の一般形は
        da = i(g₁|b|²a + g₂b²ā),  db = i(g₁|a|²b + g₂a²b̄)
    零閉塞 C = Σ(a²+b²) の一次変化は
        dC = 2i(g₁+g₂)Σ(|b|²a² + |a|²b²)
    ⇒ 全状態で dC=0 ⟺ **g₂ = −g₁**（閉塞保存が相対係数を一意化）。
    v2 の頂点:
        da = i·g·(|b|²a − b²ā) = −2g·b·Im(b̄a)
        db = i·g·(|a|²b − a²b̄) = −2g·a·Im(āb)
    ——局所位相差 Im(b̄a) が無いところでは力が働かない（関係主義形）。
    総ノルムも一次で厳密保存（dN = 2Re(ig·実数) = 0）。
    積分器: 中点法（RK2）を採用し、離散化誤差を O(g³)/step に落とす。

テスト:
    T1 g=0 で原本と bitwise 一致（対応原理、v1 と同じ）
    T2 閉塞・ノルム保存の v1/v2 比較（同一 g・同一初期値で桁の改善を実測）
    T3 パリティ定理（純偶×純偶 → P_f 恒等0）
    T4 種付き利得の再測定（v2 頂点でも凝縮度依存利得が残るか——開いた問い）
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("uim_v1", HERE / "universal_inelastic_map_v1.py")
v1 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = v1
spec.loader.exec_module(v1)
toy = v1.toy
base = v1.base

J_MAX = 200
EVEN_KS = v1.EVEN_KS
ODD_KS = v1.ODD_KS


def vertex_v2(a: np.ndarray, b: np.ndarray, g: float):
    """閉塞保存頂点（g₂=−g₁）。(da, db) を返す。"""
    da = 1j * g * (np.abs(b) ** 2 * a - b * b * np.conj(a))
    db = 1j * g * (np.abs(a) ** 2 * b - a * a * np.conj(b))
    return da, db


def inelastic_pass_v2(a: np.ndarray, b: np.ndarray, g: float):
    """中点法一段（O(g³)/step）。g=0 で恒等。"""
    if g == 0.0:
        return a, b
    da1, db1 = vertex_v2(a, b, g)
    am, bm = a + 0.5 * da1, b + 0.5 * db1
    da2, db2 = vertex_v2(am, bm, g)
    return a + da2, b + db2


def collision_step_v2(a, b, sp, g: float):
    readout = toy.theta_from_ab(a, b, sp)
    a, b = toy.rotate_ab(a, b, readout.theta)
    a, b = inelastic_pass_v2(a, b, g)
    return a, b, readout


def closure_norm(a, b):
    return (abs(complex(np.sum(a * a) + np.sum(b * b))),
            float(np.vdot(a, a).real + np.vdot(b, b).real))


def main() -> None:
    t0 = time.time()
    print("万能非弾性写像 v2（閉塞保存・一意化頂点）単体テスト")
    params = base.Params(high_n=63, recursive_collision_count=J_MAX)
    sp = base.build_source_params(params)
    results = {"derivation": "閉塞保存 ⇔ g2 = -g1（一意化）。中点法積分"}

    # ---- T1 対応原理（g=0 bitwise）----
    a1 = v1.make_bundle(sp, EVEN_KS, "A"); b1 = v1.make_bundle(sp, ODD_KS, "B")
    a2, b2 = a1.copy(), b1.copy()
    for _ in range(50):
        ro = toy.theta_from_ab(a1, b1, sp)
        a1, b1 = toy.rotate_ab(a1, b1, ro.theta)
        a2, b2, _ = collision_step_v2(a2, b2, sp, 0.0)
    t1 = np.array_equal(a1, a2) and np.array_equal(b1, b2)
    print(f"\n[T1] g=0 bitwise 対応: {'PASS' if t1 else 'FAIL'}")
    results["T1"] = bool(t1)

    # ---- T2 保存の v1/v2 比較（同一条件）----
    print("\n[T2] 閉塞・ノルム保存: v1 vs v2（g=1e-3, 200衝突, 同一初期値）")
    comp = {}
    for tag, step in (("v1", lambda a, b: v1.collision_step(a, b, sp, 1e-3)),
                       ("v2", lambda a, b: collision_step_v2(a, b, sp, 1e-3))):
        a = v1.make_bundle(sp, EVEN_KS, "A"); b = v1.make_bundle(sp, ODD_KS, "B")
        c0, n0 = closure_norm(a, b)
        for _ in range(J_MAX):
            a, b, _ = step(a, b)
        c1, n1 = closure_norm(a, b)
        comp[tag] = {"closure_drift": c1 - c0, "norm_drift": abs(n1 - n0)}
        print(f"  {tag}: 閉塞ドリフト={c1 - c0:.3e}  ノルムドリフト={abs(n1 - n0):.3e}")
    improve = comp["v1"]["closure_drift"] / max(comp["v2"]["closure_drift"], 1e-300)
    print(f"  閉塞ドリフト改善率: v1/v2 = {improve:.1e}")
    results["T2"] = {**comp, "improvement_ratio": improve}

    # ---- T3 パリティ定理 ----
    a = v1.make_bundle(sp, EVEN_KS, "A"); b = v1.make_bundle(sp, EVEN_KS, "B")
    pf_max = 0.0
    for _ in range(J_MAX):
        a, b, ro = collision_step_v2(a, b, sp, 0.5)
        pf_max = max(pf_max, float(ro.fermionic_relation_power))
    t3 = pf_max < 1e-25
    print(f"\n[T3] パリティ定理（白猫×白猫, g=0.5）: P_f max={pf_max:.3e} → "
          f"{'PASS' if t3 else 'FAIL'}")
    results["T3"] = {"pass": bool(t3), "pf_max": pf_max}

    # ---- T4 種付き利得（v2 頂点）----
    print("\n[T4] 種付き利得（v2, 偶ポンプ＋奇シード1e-5, g=0.05, 占有率掃引）")
    gains = {}
    for s in (0.5, 1.0, 2.0, 4.0, 8.0):
        a = v1.make_bundle(sp, EVEN_KS, "A", scale=s)
        a = a + v1.make_bundle(sp, ODD_KS, "A", scale=1e-5 * s)
        b = v1.make_bundle(sp, EVEN_KS, "B", scale=s)
        pfs = [v1.fermionic_power_raw(a, b, sp)]
        cds = []
        c0, _ = closure_norm(a, b)
        for _ in range(J_MAX):
            a, b, _ = collision_step_v2(a, b, sp, 0.05)
            pfs.append(v1.fermionic_power_raw(a, b, sp))
        c1, _ = closure_norm(a, b)
        pfs = np.array(pfs)
        growth = float(pfs[-1] / max(pfs[0], 1e-300))
        with np.errstate(divide="ignore"):
            rate = float(np.polyfit(np.arange(len(pfs)),
                                     np.log(np.maximum(pfs, 1e-300)), 1)[0])
        gains[s] = {"pf_initial": float(pfs[0]), "pf_final": float(pfs[-1]),
                     "growth_ratio": growth, "log_rate_per_collision": rate,
                     "closure_drift": c1 - c0}
        print(f"  s={s:4.1f}（占有∝{s*s:5.2f}）: P_f {pfs[0]:.3e}→{pfs[-1]:.3e} "
              f"増幅率={growth:.3e} 利得/衝突={rate:+.4f} 閉塞ドリフト={c1-c0:.2e}")
    results["T4"] = gains

    results["runtime_sec"] = time.time() - t0
    (HERE / "universal_inelastic_map_test_result_v2.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
