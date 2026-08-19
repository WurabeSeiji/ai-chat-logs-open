#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H4 第三段 v1 — 毛ゲージ不変性の証明が予言する例外：巻き一致対では不変性が破れるか

証明の骨子（DLdyn ノートに書き下し）:
  純巻き状態 a=α(χ)e^{i w_A η}, b=β(χ)e^{i w_B η} に対し、二体力学のチャネル間
  交差項はすべて位相因子 e^{i(w_A−w_B)η} を持つ。w_A≠w_B (mod n_η) なら η 和
  Σ_η e^{i(w_A−w_B)η}=0 により交差項が消え、χ周辺化観測量は巻きに厳密不変（第一段の実測）。
  **例外**: w_A=w_B のとき交差項が生き残り、不変性の証明が成立しない。

予言（本プローブが検定）:
  P-1 不等巻き類（w_A≠w_B）: 全ケース互いに厳密一致（第一段の再確認を含む）
  P-2 等巻き類（w_A=w_B）: 共通巻きシフトの下で互いに厳密一致（共通巻き＝ゲージ）
  P-3 等巻き類は不等巻き類と**異なる力学**を示す（交差項の活性化）

ケース（imprint (mA,mB)。素の巻きは A:+1, B:+2）:
  不等巻き類: (0,0)→w=(1,2) / (+3,−3)→w=(4,−1)
  等巻き類:   (+1,0)→w=(2,2) / (0,−1)→w=(1,1) / (+2,+1)→w=(3,3)

力学は二体正本 collision_step_exact のみ。判定は分離角軌道の最大差。
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


_uni = _load("uni_h43", UNI / "unified_interaction_v1.py")
_cr0 = _load("cr0_h43", EXP / "run_cr0_control_no_theta_v2.py")
base = _uni.two_body_base
step = _uni.collision_step_exact

PACKET = tuple(range(1, 18))
T_STEPS = 200
CASES = {
    "uneq_0_0": (0, 0), "uneq_p3_m3": (3, -3),
    "eq_p1_0": (1, 0), "eq_0_m1": (0, -1), "eq_p2_p1": (2, 1),
}


def make_ab(sp, slope, icept, mA, mB, nc, ne):
    case = base.explicit_packet_case(
        mode=f"h4s3_{mA}_{mB}", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=_cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(+30.0, slope, icept))
    a = base.make_case_state(sp, case, "A", hair_enabled=True)
    b = base.make_case_state(sp, case, "B", hair_enabled=True)
    eta = 2.0 * np.pi * np.arange(ne) / ne
    a = (a.reshape(nc, ne) * np.exp(1j * mA * eta)[None, :]).reshape(-1)
    b = (b.reshape(nc, ne) * np.exp(1j * mB * eta)[None, :]).reshape(-1)
    return a, b


def sep_series(sp, slope, icept, mA, mB, nc, ne):
    a, b = make_ab(sp, slope, icept, mA, mB, nc, ne)
    _, _, wA = _cr0.winding_spectrum(a, nc, ne)
    _, _, wB = _cr0.winding_spectrum(b, nc, ne)
    seps = []
    for _ in range(T_STEPS):
        a, b, _ = step(a, b, sp)
        ta, _ = _cr0.circle_position(a, nc, ne)
        tb, _ = _cr0.circle_position(b, nc, ne)
        seps.append(abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb)))))))
    return np.array(seps), wA, wB


def main():
    t0 = time.time()
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, nc, ne)
    print(f"格子 n_chi={nc} n_eta={ne}  T={T_STEPS}")

    series, winds = {}, {}
    for name, (mA, mB) in CASES.items():
        s, wA, wB = sep_series(sp, slope, icept, mA, mB, nc, ne)
        series[name] = s
        winds[name] = (wA, wB)
        print(f"[{name:10s}] 巻き実測 w=({wA:+.2f},{wB:+.2f})  "
              f"sep: {s[0]:.3f}°→{s[-1]:.3f}°  min={s.min():.3f}")

    def d(x, y):
        return float(np.max(np.abs(series[x] - series[y])))

    print("\n--- 予言の検定（分離角軌道の最大差 [deg]）---")
    p1 = d("uneq_0_0", "uneq_p3_m3")
    print(f"P-1 不等巻き類の内部一致: |uneq(0,0) − uneq(+3,−3)| = {p1:.3e}")
    p2a = d("eq_p1_0", "eq_p2_p1")
    p2b = d("eq_p1_0", "eq_0_m1")
    print(f"P-2 等巻き類の内部一致:   |eq(2,2) − eq(3,3)| = {p2a:.3e}   "
          f"|eq(2,2) − eq(1,1)| = {p2b:.3e}")
    p3 = d("uneq_0_0", "eq_p1_0")
    print(f"P-3 類の間の差:           |uneq(0,0) − eq(2,2)| = {p3:.3e}")
    verdict = {
        "P1_unequal_class_internal": p1,
        "P2_equal_class_internal_22_33": p2a,
        "P2_equal_class_internal_22_11": p2b,
        "P3_between_classes": p3,
    }
    if p3 > 1e-6 and max(p1, p2a, p2b) < 1e-9:
        print("→ 予言どおり: 等巻き対だけが異なる力学（巻き一致選択則）")
    elif p3 < 1e-9:
        print("→ 予言は不成立: 等巻きでも不変（交差項が観測量に現れない）")
    else:
        print("→ 混在: 類の内部一致が崩れている。診断が必要")

    out = {"T": T_STEPS, "cases": {k: {"m": list(v), "w": [winds[k][0], winds[k][1]],
                                        "sep_head": [float(x) for x in series[k][:8]],
                                        "sep_end": float(series[k][-1])}
                                    for k, v in CASES.items()},
           "verdict": verdict, "elapsed_sec": time.time() - t0}
    (HERE / "result_h4_sign_v2_matching_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"\n保存: result_h4_sign_v2_matching_v1.json ({out['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
