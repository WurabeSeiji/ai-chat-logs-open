#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H4 裁定プローブ v1 — 巻き符号 (+,+) vs (+,−) の対で、誘導される相対運動の符号を測る

目的（DLdyn 導出ノート §6-5・DL6 判定 Q5 の理論側）:
  クーロン符号則の対決 H4 を、二体正本 collision_step_exact の上で直接裁定する。
  - W11 §4.3（実測正本）: 同符号＝斥力、逆符号＝引力
  - G9（実測）: E(++) が引力側
  - CR1: 符号非依存（χ 一変数縮約の人工物と診断済み）

方法:
  等振幅倍音 1..17 のパケットを A=−30°, B=+30° に置き（CR 系列と同一構成）、
  毛（η）方向の巻き e^{i m η} を刻印する。ケース: (mA,mB) = (0,0), (+3,+3), (+3,−3)。
  力学は二体正本 collision_step_exact（万能非弾性写像・厳密解）のみ。
  相対距離＝円周上の分離角 |Δθ|（重心第1モーメント、CR0 の circle_position）。
  判定: |Δθ|(τ) の初期勾配と時系列。(++) と (+−) で符号が分かれるか。

出所: 恒久ルールに従い F は unified_interaction_v1 経由の正本を read-only import。
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent          # dynlab_v1
EXP = HERE.parent                                # 電子の反跳実験
UNI = EXP.parent / "統一万能関数_v1"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


_uni = _load("uni_h4", UNI / "unified_interaction_v1.py")
_cr0 = _load("cr0_h4", EXP / "run_cr0_control_no_theta_v2.py")

base = _uni.two_body_base
collision_step_exact = _uni.collision_step_exact

PACKET = tuple(range(1, 18))
DEG_A, DEG_B = -30.0, +30.0
T_STEPS = 400
CASES = {"m0_0": (0, 0), "m+3_+3": (3, 3), "m+3_-3": (3, -3)}


def imprint_winding(psi, m, n_chi, n_eta):
    """毛（η）方向の巻き e^{i m η} を刻印する。m=0 は恒等。"""
    if m == 0:
        return psi
    A2 = psi.reshape(n_chi, n_eta)
    eta = 2.0 * np.pi * np.arange(n_eta) / n_eta
    return (A2 * np.exp(1j * m * eta)[None, :]).reshape(-1)


def sep_deg(a, b, n_chi, n_eta):
    """分離角 |Δθ|［度］（各チャネルの重心角の差、円周でラップ）。"""
    ta, _ = _cr0.circle_position(a, n_chi, n_eta)
    tb, _ = _cr0.circle_position(b, n_chi, n_eta)
    d = np.degrees(np.angle(np.exp(1j * (ta - tb))))
    return abs(float(d))


def run_case(sp, n_chi, n_eta, slope, icept, mA, mB):
    case = base.explicit_packet_case(
        mode=f"h4_{mA}_{mB}", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = base.make_case_state(sp, case, "A", hair_enabled=True)
    b = base.make_case_state(sp, case, "B", hair_enabled=True)
    a = imprint_winding(a, mA, n_chi, n_eta)
    b = imprint_winding(b, mB, n_chi, n_eta)

    # 巻きの検算（符号付き平均巻き数）
    _, _, wA = _cr0.winding_spectrum(a, n_chi, n_eta)
    _, _, wB = _cr0.winding_spectrum(b, n_chi, n_eta)

    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    seps, rs = [], []
    for _ in range(T_STEPS):
        a, b, ro = collision_step_exact(a, b, sp)
        seps.append(sep_deg(a, b, n_chi, n_eta))
        rs.append(float(ro.reflection_rate))
    tot1 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    seps = np.array(seps)

    # 初期勾配（最初の 20 步の最小二乗）と全区間平均勾配
    k = 20
    t = np.arange(k)
    slope0 = float(np.polyfit(t, seps[:k], 1)[0])
    slope_all = float(np.polyfit(np.arange(len(seps)), seps, 1)[0])
    return {
        "winding_A": wA, "winding_B": wB,
        "sep0_deg": float(seps[0]), "sep_end_deg": float(seps[-1]),
        "sep_min": float(seps.min()), "sep_max": float(seps.max()),
        "slope_first20_deg_per_step": slope0,
        "slope_all_deg_per_step": slope_all,
        "r_mean": float(np.mean(rs)),
        "norm_drift": abs(tot1 / tot0 - 1.0),
        "sep_series_head": [float(x) for x in seps[:10]],
    }


def main():
    t0 = time.time()
    sp = base.build_source_params(
        base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    print(f"格子 n_chi={n_chi} n_eta={n_eta}  T={T_STEPS}  packets=1..17 at ∓30°")

    out = {"T": T_STEPS, "packet": list(PACKET), "deg": [DEG_A, DEG_B], "cases": {}}
    for name, (mA, mB) in CASES.items():
        res = run_case(sp, n_chi, n_eta, slope, icept, mA, mB)
        out["cases"][name] = res
        print(f"[{name}]  巻き実測 A={res['winding_A']:+.3f} B={res['winding_B']:+.3f}  "
              f"r̄={res['r_mean']:.4f}")
        print(f"   |Δθ|: {res['sep0_deg']:.3f}° → {res['sep_end_deg']:.3f}°  "
              f"(min {res['sep_min']:.3f} / max {res['sep_max']:.3f})")
        print(f"   勾配: 初期20步 {res['slope_first20_deg_per_step']:+.5f}°/步  "
              f"全区間 {res['slope_all_deg_per_step']:+.5f}°/步  "
              f"ノルムドリフト {res['norm_drift']:.2e}")

    s_pp = out["cases"]["m+3_+3"]["slope_first20_deg_per_step"]
    s_pm = out["cases"]["m+3_-3"]["slope_first20_deg_per_step"]
    print("\n判定（初期勾配。+＝離れる＝斥力側 / −＝近づく＝引力側）:")
    print(f"  (+,+): {s_pp:+.5f}°/步   (+,−): {s_pm:+.5f}°/步")
    if s_pp * s_pm < 0:
        w11 = (s_pp > 0 and s_pm < 0)
        print("  → 符号が分かれた。" + ("W11 側（同符号=斥力・逆符号=引力）" if w11
                                       else "W11 と逆（同符号=引力・逆符号=斥力）"))
    else:
        print("  → 符号は分かれない（この構成では巻き符号に非依存）")

    out["elapsed_sec"] = time.time() - t0
    (HERE / "result_h4_sign_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"\n保存: result_h4_sign_v1.json  ({out['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
