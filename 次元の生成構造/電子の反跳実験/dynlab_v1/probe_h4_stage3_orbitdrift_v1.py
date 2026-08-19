#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H4 第三段b v1 — 束縛状態まわりの微小ドリフト摂動：向きの符号は永年効果を持つか

第二段の教訓: 大きなドリフト（0.5°/步）では接近/離反の自由経路の幾何が異なり、
力の符号を分離できない。本プローブは束縛状態のまわりの**微小**摂動で幾何を揃える。

構成:
  対称形（A=B=1..17）は χ 反転＋チャネル交換の対称性により ±δ が同値になりうるため、
  **非対称パケット（A=1..17, B=1..3。CR1 と同一）**で対称性を破る。
  δ = ±0.02°/步 の相対ドリフト（approach / recede）＋衝突を T=600。

判定:
  後期窓（後半300步）の平均分離角と振動振幅を ±δ で比較。
  差が出れば「相対ドリフトの向きが束縛状態の永年量に結合する」（ゲージ層符号則の初証拠）。
  出なければ「この構成では向き結合なし」。どちらも二値で記録する。
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


_uni = _load("uni_h44", UNI / "unified_interaction_v1.py")
K = _load("kin_h44", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_h44", EXP / "run_cr0_control_no_theta_v2.py")
base = _uni.two_body_base
step = _uni.collision_step_exact

PACKET_A = tuple(range(1, 18))
PACKET_B = (1, 2, 3)
V_DEG = 0.02
T_STEPS = 600
CASES = {"static": (0.0, 0.0), "approach": (+V_DEG, -V_DEG), "recede": (-V_DEG, +V_DEG)}


def main():
    t0 = time.time()
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, nc, ne)

    # 並進較正（1° 指定の実移動から符号込みで決める）
    case0 = base.explicit_packet_case(
        mode="h4s3b_cal", packet_a=PACKET_A, packet_b=PACKET_B,
        packet_a_shift=_cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(+30.0, slope, icept))
    a0 = base.make_case_state(sp, case0, "A", hair_enabled=True)
    probe = np.radians(1.0)
    p0, _ = _cr0.circle_position(a0, nc, ne)
    p1, _ = _cr0.circle_position(K.k_translate_flat(a0, probe, nc, ne), nc, ne)
    moved = float(np.degrees(np.angle(np.exp(1j * (p1 - p0)))))
    om_per_deg = probe * (1.0 if moved > 0 else -1.0) / abs(moved)
    print(f"格子 {nc}x{ne}  δ={V_DEG}°/步  T={T_STEPS}  A=1..17 B=1..3（非対称）")

    out = {"v_deg": V_DEG, "T": T_STEPS, "cases": {}}
    late = {}
    for name, (vA, vB) in CASES.items():
        case = base.explicit_packet_case(
            mode=f"h4s3b_{name}", packet_a=PACKET_A, packet_b=PACKET_B,
            packet_a_shift=_cr0.shift_for_deg(-30.0, slope, icept),
            packet_b_shift=_cr0.shift_for_deg(+30.0, slope, icept))
        a = base.make_case_state(sp, case, "A", hair_enabled=True)
        b = base.make_case_state(sp, case, "B", hair_enabled=True)
        seps = []
        for _ in range(T_STEPS):
            if vA != 0.0:
                a = K.k_translate_flat(a, om_per_deg * vA, nc, ne)
            if vB != 0.0:
                b = K.k_translate_flat(b, om_per_deg * vB, nc, ne)
            a, b, _ = step(a, b, sp)
            ta, _ = _cr0.circle_position(a, nc, ne)
            tb, _ = _cr0.circle_position(b, nc, ne)
            seps.append(abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb)))))))
        seps = np.array(seps)
        lw = seps[T_STEPS // 2:]
        late[name] = lw
        out["cases"][name] = {
            "late_mean": float(lw.mean()), "late_std": float(lw.std()),
            "late_min": float(lw.min()), "late_max": float(lw.max()),
        }
        print(f"[{name:8s}] 後期窓: 平均={lw.mean():8.4f}°  振幅std={lw.std():7.4f}  "
              f"[{lw.min():.3f}, {lw.max():.3f}]")

    dm = out["cases"]["approach"]["late_mean"] - out["cases"]["recede"]["late_mean"]
    ds = out["cases"]["approach"]["late_std"] - out["cases"]["recede"]["late_std"]
    base_m = out["cases"]["static"]["late_mean"]
    print(f"\n判定: 後期平均分離の差 approach−recede = {dm:+.5f}°"
          f"（static 基準 {base_m:.4f}°）  振幅差 = {ds:+.5f}")
    out["approach_minus_recede_late_mean"] = dm
    out["approach_minus_recede_late_std"] = ds
    out["elapsed_sec"] = time.time() - t0
    (HERE / "result_h4_stage3_orbitdrift_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"保存: result_h4_stage3_orbitdrift_v1.json ({out['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
