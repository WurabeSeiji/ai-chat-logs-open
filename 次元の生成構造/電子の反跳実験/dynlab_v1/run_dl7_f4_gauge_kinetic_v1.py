#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F4 走行 v1 — ゲージ運動項つき±則：速度レジスタの蓄積符号（実行前固定）

木原指摘: 運動はゲージ層（速度レジスタ v の carry）・位置は関係から毎瞬読む
（レジスタ不要）——DL4-0 の帳簿どおり。F3 系列は運動項なしの状態層を見ていた誤り
（CR0: 衝突は位置を巡らせない）。

ループ（CR1 正本 run() の忠実コピー・状態のみ差し替え可能に）:
  χ = 相対位置（circle_position・状態読み）
  r = theta_from_ab（状態読み）、κ = 1−r（透過率・素電荷側）
  acc = −4sin²(ω/2)·χ、v += acc、ω += κ·acc
  a ← K.k_translate_flat(a, −v)（万能運動関数＝ゲージ並進・片側適用が正本）
  (a,b) ← collision_step_exact

判定:
  F4-0 アンカー（前提）: CR1 初期状態で本ループが CR1 実測（周期 5.456・有界）を再現
       ——通らなければ電荷判定に進まない
  F4-1 電荷応答: フェルミオン型 φ∈{0,π} 対で v(τ)・χ(τ) の軌道が分岐するか
       （分岐量が丸め増幅を超えるか）
  F4-2 蓄積符号: v の長時間平均・χ の平均の φ 依存の符号（引力/斥力側の読み——
       評価出力として記録。5幾何の符号検定つき）

出力: result_dl7_f4_v1.json・dl7_f4_series_v1.npz
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
ODD = tuple(range(1, 18, 2))
T = 4000
OMEGA_CLOCK = np.pi / 72.0
GEOMS = [(-30.0, 20.0), (-30.0, 30.0), (-30.0, 45.0), (-20.0, 40.0), (-45.0, 30.0)]

def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    sys.modules[n] = m; s.loader.exec_module(m); return m

def main():
    t0 = time.time()
    uni = _load("uni_f4", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f4", EXP / "run_cr0_control_no_theta_v2.py")
    K = _load("kin_f4", UNI / "unified_kinetic_v1.py")
    toy = uni.two_body_v1.toy
    base = uni.two_body_base; step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    def loop(a, b, T_):
        """CR1 run() の力学部の忠実コピー（κ=transmission・v0=0）。"""
        omega, v = OMEGA_CLOCK, 0.0
        chi_s = np.empty(T_); v_s = np.empty(T_)
        for t in range(T_):
            pa, _ = cr0.circle_position(a, nc, ne)
            pb, _ = cr0.circle_position(b, nc, ne)
            chi = float(np.angle(np.exp(1j * (pa - pb))))
            r_now = float(toy.theta_from_ab(a, b, sp).reflection_rate)
            kappa = 1.0 - r_now
            acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
            v = v + acc
            omega = omega + kappa * acc
            a = K.k_translate_flat(a, -v, nc, ne)
            a, b, _ = step(a, b, sp)
            chi_s[t] = chi; v_s[t] = v
        return chi_s, v_s

    # ---- F4-0 アンカー: CR1 初期状態で周期・有界性を再現 ----
    aC, bC = cr0.make_pair(sp, cr0.shift_for_deg(-30.0, slope, icept),
                           cr0.shift_for_deg(+30.0, slope, icept))
    chiC, vC = loop(aC.copy(), bC.copy(), 3000)
    deg = np.degrees(chiC)
    sgn = np.sign(deg - deg.mean())
    zc = np.where(sgn[1:] * sgn[:-1] < 0)[0]
    period = 2.0 * float(np.mean(np.diff(zc))) if len(zc) > 3 else float("nan")
    canon = json.loads((EXP / "result_cr1_kinetic_feedback_v1.json").read_text())
    per_c = float(canon["metrics"]["period_chi"])
    F40 = {"period": period, "canonical": per_c, "ratio": period / per_c,
           "bounded": bool(np.max(np.abs(deg)) < 61.0),
           "pass": bool(abs(period / per_c - 1) < 0.05 and np.max(np.abs(deg)) < 61.0)}

    # ---- F4-1/F4-2: フェルミオン型 φ 対・5幾何 ----
    def fpair(dega, degb, phi0):
        c = base.explicit_packet_case(
            mode=f"f4_{dega}_{degb}", packet_a=ODD, packet_b=ODD,
            packet_a_shift=cr0.shift_for_deg(dega, slope, icept),
            packet_b_shift=cr0.shift_for_deg(degb, slope, icept))
        a = base.make_case_state(sp, c, "A", hair_enabled=True)
        b = base.make_case_state(sp, c, "B", hair_enabled=True)
        a /= np.sqrt(np.vdot(a, a).real); b /= np.sqrt(np.vdot(b, b).real)
        return a, np.exp(1j * phi0) * b

    rows = []
    store = {}
    for g in GEOMS:
        a0, b0 = fpair(g[0], g[1], 0.0)
        aP, bP = fpair(g[0], g[1], np.pi)
        chi0, v0 = loop(a0, b0, T)
        chiP, vP = loop(aP, bP, T)
        d_abschi = float(np.mean(np.abs(np.degrees(chi0)))
                         - np.mean(np.abs(np.degrees(chiP))))
        d_v = float(np.mean(np.abs(v0)) - np.mean(np.abs(vP)))
        dev200 = float(np.max(np.abs(np.degrees(chi0[:200] - chiP[:200]))))
        rows.append({"geom": g, "mean_abs_chi_phi0": float(np.mean(np.abs(np.degrees(chi0)))),
                     "mean_abs_chi_phiPi": float(np.mean(np.abs(np.degrees(chiP)))),
                     "diff_abs_chi": d_abschi, "sign": int(np.sign(d_abschi)),
                     "diff_mean_abs_v": d_v, "dev200_deg": dev200})
        store[f"{g}"] = (chi0, chiP, v0, vP)
        print(f"  geom={g}: |χ|平均 φ=0 {rows[-1]['mean_abs_chi_phi0']:7.2f}° "
              f"π {rows[-1]['mean_abs_chi_phiPi']:7.2f}°  差 {d_abschi:+8.3f}°  "
              f"dev200={dev200:.2e}°")
    signs = [r["sign"] for r in rows]
    unanimous = bool(abs(sum(signs)) == len(signs))
    res = {"config": {"T": T, "geoms": GEOMS, "packet_odd": list(ODD),
                      "loop": "CR1忠実コピー（κ=transmission・v0=0・片側translate）"},
           "F4_0_anchor": F40, "rows": rows,
           "F4_1_signal": bool(max(r["dev200_deg"] for r in rows) > 1e-6),
           "F4_2_sign_test": {"signs": signs, "unanimous": unanimous,
                              "verdict": ("同シート=|χ|大（斥力側）全会一致" if unanimous and signs[0] > 0
                                          else "同シート=|χ|小（引力側）全会一致" if unanimous
                                          else "配置で割れる——未決着")},
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f4_series_v1.npz",
                        chiC=chiC, vC=vC,
                        **{f"chi0_{i}": store[f"{g}"][0] for i, g in enumerate(GEOMS)},
                        **{f"chiP_{i}": store[f"{g}"][1] for i, g in enumerate(GEOMS)},
                        **{f"v0_{i}": store[f"{g}"][2] for i, g in enumerate(GEOMS)},
                        **{f"vP_{i}": store[f"{g}"][3] for i, g in enumerate(GEOMS)})
    (HERE / "result_dl7_f4_v1.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"F4-0 アンカー: period={period:.3f}（正本 {per_c:.3f}・比 {period/per_c:.3f}）"
          f" bounded={F40['bounded']} pass={F40['pass']}")
    print(f"F4-1 signal={res['F4_1_signal']}  F4-2: {res['F4_2_sign_test']['verdict']}")
    print(f"({res['elapsed_sec']:.0f}s)")

if __name__ == "__main__":
    main()
