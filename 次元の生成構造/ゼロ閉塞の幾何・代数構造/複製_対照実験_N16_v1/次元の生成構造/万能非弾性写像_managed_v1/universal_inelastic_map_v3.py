#!/usr/bin/env python3
"""万能非弾性写像 v3：結合定数の廃止——強さも万能読出し（反射率）から

v2 からの変更（2026-08-04 木原氏指示「g を持ち込む必要はない、反射率をそのまま使え」）:
    非弾性の強さ = readout.reflection_rate（R = sin²θ、theta_from_ab の万能読出し）。
    新設定数ゼロ。写像は完全に既存の万能関数＋閉塞保存の一意化頂点だけで決まる:

        衝突 = 弾性部（θ読出し→rotate_ab、無変更）
             ∘ 非弾性部: da = i·R·(|b|²a − b²ā),  db = i·R·(|a|²b − a²b̄)（中点法）

    帰結:
    - 白猫×白猫: R≡0 → 非弾性は厳密オフ → 原本と bitwise 一致（任意の占有率で）
      ——対応原理が近似でなく恒等になる
    - 自己触媒: R ∝ フェルミオン関係量割合 → フェルミオン内容が増えるほど結合が強まる
      （正帰還＝点火構造）。閾値は力学から出る
    - α 接続の内蔵: R と素電荷の橋（1−R = sin²θ = √(4πα)）は系列で確立済み

テスト:
    T1 恒等対応: 白猫×白猫（占有スケール s=8 でも）原本パイプラインと bitwise 一致・P_f≡0
    T2 保存: odd-odd（単位ノルム・R=61/64 で非弾性常時オン）200衝突の閉塞・ノルムドリフト
    T3 希薄還元: odd-odd の v3 と弾性のみの軌道差が s→0 で消えること（スケーリング実測）
    T4 種付き点火: 偶ポンプ＋奇シード・占有率掃引——自己触媒成長の形（点火の有無）を記録
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("uim_v1_for3", HERE / "universal_inelastic_map_v1.py")
v1 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = v1
spec.loader.exec_module(v1)
toy = v1.toy
base = v1.base

J_MAX = 200
EVEN_KS = v1.EVEN_KS
ODD_KS = v1.ODD_KS


def vertex(a, b, r):
    da = 1j * r * (np.abs(b) ** 2 * a - b * b * np.conj(a))
    db = 1j * r * (np.abs(a) ** 2 * b - a * a * np.conj(b))
    return da, db


def collision_step_v3(a, b, sp):
    """衝突一回。強さは readout.reflection_rate（新設定数なし）。"""
    readout = toy.theta_from_ab(a, b, sp)
    a, b = toy.rotate_ab(a, b, readout.theta)
    r = float(readout.reflection_rate)
    if r > 0.0:
        da1, db1 = vertex(a, b, r)
        am, bm = a + 0.5 * da1, b + 0.5 * db1
        da2, db2 = vertex(am, bm, r)
        a, b = a + da2, b + db2
    return a, b, readout


def closure_norm(a, b):
    return (abs(complex(np.sum(a * a) + np.sum(b * b))),
            float(np.vdot(a, a).real + np.vdot(b, b).real))


def seed_fraction_sweep() -> None:
    """T5 自己触媒則の検証（再現用・s=8固定・シード割合4桁掃引）。
    発見済みの法則: d(lnP_f)/d衝突 = C·f²（f=フェルミオン割合）、C(s=8)≈10.4。"""
    import json as _json
    params = base.Params(high_n=63, recursive_collision_count=J_MAX)
    sp = base.build_source_params(params)
    s = 8.0
    rows = []
    print("[T5] シード割合掃引（s=8, 200衝突）: 成長率/f₀² の一定性")
    for seed_amp in (1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
        a = v1.make_bundle(sp, EVEN_KS, "A", scale=s)
        a = a + v1.make_bundle(sp, ODD_KS, "A", scale=seed_amp * s)
        b = v1.make_bundle(sp, EVEN_KS, "B", scale=s)
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        pfs = [v1.fermionic_power_raw(a, b, sp)]
        for _ in range(J_MAX):
            a, b, _ = collision_step_v3(a, b, sp)
            pfs.append(v1.fermionic_power_raw(a, b, sp))
        pfs = np.array(pfs)
        f0 = float(pfs[0] / tot)
        lp = np.log(np.maximum(pfs[:21], 1e-300))
        rate0 = float(np.polyfit(np.arange(21), lp, 1)[0])
        C = rate0 / f0 ** 2
        rows.append({"seed_amp": seed_amp, "f0": f0, "rate0": rate0, "C": C,
                     "growth_200": float(pfs[-1] / pfs[0])})
        print(f"  seed={seed_amp:7.0e}: f0={f0:.3e} rate0={rate0:+.3e} C=rate0/f0²={C:.3f}")
    Cs = [r["C"] for r in rows]
    print(f"  C の変動: {min(Cs):.3f}〜{max(Cs):.3f}（比 {max(Cs)/min(Cs):.3f}——一定なら法則成立）")
    Path(__file__).parent.joinpath("seed_fraction_sweep_result_v3.json").write_text(
        _json.dumps({"s": s, "rows": rows, "law": "d(lnPf)/dj = C*f^2"},
                    ensure_ascii=False, indent=2), encoding="utf-8")
    print("saved: seed_fraction_sweep_result_v3.json")


def main() -> None:
    t0 = time.time()
    print("万能非弾性写像 v3（結合定数なし・強さ=反射率）単体テスト")
    params = base.Params(high_n=63, recursive_collision_count=J_MAX)
    sp = base.build_source_params(params)
    results = {"design": "coupling = reflection_rate (no new constants)"}

    # ---- T1 恒等対応（白猫×白猫, s=8）----
    s = 8.0
    a1 = v1.make_bundle(sp, EVEN_KS, "A", scale=s); b1 = v1.make_bundle(sp, EVEN_KS, "B", scale=s)
    a2, b2 = a1.copy(), b1.copy()
    pf_max = 0.0
    for _ in range(J_MAX):
        ro = toy.theta_from_ab(a1, b1, sp)
        a1, b1 = toy.rotate_ab(a1, b1, ro.theta)
        a2, b2, ro2 = collision_step_v3(a2, b2, sp)
        pf_max = max(pf_max, float(ro2.fermionic_relation_power))
    t1 = np.array_equal(a1, a2) and np.array_equal(b1, b2) and pf_max == 0.0
    print(f"\n[T1] 白猫×白猫 s=8: bitwise一致={np.array_equal(a1, a2) and np.array_equal(b1, b2)} "
          f"P_f max={pf_max:.1e} → {'PASS（R=0で非弾性は厳密オフ）' if t1 else 'FAIL'}")
    results["T1"] = bool(t1)

    # ---- T2 保存（odd-odd 単位ノルム: R=61/64 で常時オン）----
    a = v1.make_bundle(sp, ODD_KS, "A"); b = v1.make_bundle(sp, ODD_KS, "B")
    c0, n0 = closure_norm(a, b)
    r_first = None
    for _ in range(J_MAX):
        a, b, ro = collision_step_v3(a, b, sp)
        if r_first is None:
            r_first = float(ro.reflection_rate)
    c1, n1 = closure_norm(a, b)
    print(f"\n[T2] odd-odd（R₀={r_first:.6f}）200衝突: 閉塞 {c0:.2e}→{c1:.2e}"
          f"（ドリフト {c1-c0:.2e}） ノルムドリフト {abs(n1-n0):.2e}")
    results["T2"] = {"R_initial": r_first, "closure_drift": c1 - c0,
                      "norm_drift": abs(n1 - n0)}

    # ---- T3 希薄還元（v3 と弾性のみの軌道差のスケーリング）----
    print("\n[T3] 希薄還元スケーリング（odd-odd, 50衝突, 軌道差 vs s）")
    devs = {}
    for sc in (1.0, 0.5, 0.25):
        ae = v1.make_bundle(sp, ODD_KS, "A", scale=sc); be = v1.make_bundle(sp, ODD_KS, "B", scale=sc)
        av, bv = ae.copy(), be.copy()
        for _ in range(50):
            ro = toy.theta_from_ab(ae, be, sp)
            ae, be = toy.rotate_ab(ae, be, ro.theta)
            av, bv, _ = collision_step_v3(av, bv, sp)
        dev = float(np.linalg.norm(av - ae) + np.linalg.norm(bv - be)) / sc
        devs[sc] = dev
        print(f"  s={sc:4.2f}: 相対軌道差 = {dev:.3e}")
    ratio1 = devs[1.0] / max(devs[0.5], 1e-300)
    ratio2 = devs[0.5] / max(devs[0.25], 1e-300)
    t3 = devs[1.0] > devs[0.5] > devs[0.25]
    print(f"  縮小比: {ratio1:.2f}, {ratio2:.2f}（振幅²∝s²のスケーリングなら ≈4）→ "
          f"{'PASS（希薄極限で弾性へ還元）' if t3 else 'FAIL'}")
    results["T3"] = {"devs": devs, "ratios": [ratio1, ratio2], "pass": bool(t3)}

    # ---- T4 種付き点火（自己触媒成長）----
    print("\n[T4] 種付き点火（偶ポンプ＋奇シード1e-5, 占有率掃引・自己触媒の形を記録）")
    ign = {}
    for sc in (1.0, 2.0, 4.0, 8.0, 16.0):
        a = v1.make_bundle(sp, EVEN_KS, "A", scale=sc)
        a = a + v1.make_bundle(sp, ODD_KS, "A", scale=1e-5 * sc)
        b = v1.make_bundle(sp, EVEN_KS, "B", scale=sc)
        pfs = [v1.fermionic_power_raw(a, b, sp)]
        for _ in range(J_MAX):
            a, b, _ = collision_step_v3(a, b, sp)
            pfs.append(v1.fermionic_power_raw(a, b, sp))
        pfs = np.array(pfs)
        total = float(np.vdot(a, a).real + np.vdot(b, b).real)
        frac_final = float(pfs[-1] / total)
        growth = float(pfs[-1] / max(pfs[0], 1e-300))
        ign[sc] = {"pf_initial": float(pfs[0]), "pf_final": float(pfs[-1]),
                    "growth_ratio": growth, "fermionic_fraction_final": frac_final}
        print(f"  s={sc:4.1f}（占有∝{sc*sc:6.1f}）: P_f {pfs[0]:.3e}→{pfs[-1]:.3e} "
              f"増幅率={growth:.3e}  最終フェルミオン割合={frac_final:.3e}")
    results["T4"] = ign

    results["runtime_sec"] = time.time() - t0
    (HERE / "universal_inelastic_map_test_result_v3.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
