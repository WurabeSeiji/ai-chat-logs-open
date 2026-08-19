#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G v5（二体対読出し）の資格審査 — 三部作第一部のアンカー5件（実行前固定）

A1 ヌル定理: 搬送波あり（hair_enabled=True）単一倍音全対 7×7 で交差 ≤3e-17
A2 梯子: 搬送波なし（hair_enabled=False）で Δk=±2 の全対＋特別対 (1,1)（負周波数枝＝
   正本 c0 アンカー）が c=0.500000（厳密）、他は零（結合地図の実測で確定した規則）
A3 流れの恒等式: 混合状態×プローブ回転×位相8点で dN_B 実測と g_pair_flow 予言が
   機械精度一致、かつ重なり項が φ 半回転で厳密反転（反対称）
A4 θ 整合: 標準対で toy.theta_from_ab の θ を渡した g_pair_flow が実装衝突の
   dN_B と一致（θ 正本の受け渡し検査）
A5 ゲージ構造: 共通位相変換で全対読出し不変（≤1e-14）、片側 π シフトで
   全 φ_ch が π 反転（電荷共役＝状態準備であることの確認）

出力: qualification_pair_v5_result.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    G = _load("g5q", HERE / "unified_readout_v5.py")
    uni = _load("uni_g5q", HERE / "unified_interaction_v1.py")
    base = uni.two_body_base
    toy = uni.two_body_v1.toy
    sp = base.build_source_params(base.Params(high_n=63,
                                              recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)

    def single(k, which, hair):
        case = base.explicit_packet_case(mode=f"g5q_{which}_{k}_{hair}",
                                         packet_a=(k,), packet_b=(k,))
        s = base.make_case_state(sp, case, which, hair_enabled=hair)
        return s / np.sqrt(np.vdot(s, s).real)

    # ---- A1 ヌル定理（搬送波あり） ----
    max_on = 0.0
    for j in range(1, 8):
        for k in range(1, 8):
            aj = single(j, "A", True)
            bk = single(k, "B", True)
            max_on = max(max_on, abs(complex(np.vdot(aj, bk))))
    A1 = {"max_cross_carrier_on": max_on, "pass": bool(max_on < 1e-14)}

    # ---- A2 梯子（搬送波なし・Δk=±2・c=0.5） ----
    dev_ladder = 0.0
    max_off_ladder = 0.0
    for j in range(1, 8):
        for k in range(1, 8):
            aj = single(j, "A", False)
            bk = single(k, "B", False)
            c = abs(complex(np.vdot(aj, bk)))
            if k == j + 2 or k == j - 2 or (j == 1 and k == 1):
                dev_ladder = max(dev_ladder, abs(c - 0.5))
            else:
                max_off_ladder = max(max_off_ladder, c)
    A2 = {"ladder_dev_from_half": dev_ladder,
          "max_off_ladder": max_off_ladder,
          "pass": bool(dev_ladder < 1e-12 and max_off_ladder < 1e-12)}

    # ---- A3 流れの恒等式＋符号反転 ----
    a5 = single(5, "A", False)
    b7_0 = single(7, "B", False)
    b1 = single(1, "B", False)
    dev_flow = 0.0
    ovl = {}
    for th in (0.05, 0.1, 0.2, 0.3):
        for i, phi in enumerate(np.linspace(0, 2 * np.pi, 8, endpoint=False)):
            b = 0.6 * b1 + 0.8 * np.exp(1j * phi) * b7_0
            a = 0.5 * single(1, "A", False) + np.sqrt(1 - 0.25) * a5
            fl = G.g_pair_flow(a, b, th)
            a2 = a * np.cos(th) - b * np.sin(th)
            b2 = a * np.sin(th) + b * np.cos(th)
            dNb = float(np.vdot(b2, b2).real - np.vdot(b, b).real)
            dev_flow = max(dev_flow, abs(dNb - fl["dN_B_pred"]))
            ovl[(th, i)] = fl["flow_overlap_term"]
    # 反対称の検査: φ と φ+π（i と i+4）の重なり項の和は φ 非依存成分の2倍で一定
    anti = 0.0
    for th in (0.05, 0.1, 0.2, 0.3):
        sums = [ovl[(th, i)] + ovl[(th, (i + 4) % 8)] for i in range(4)]
        anti = max(anti, max(sums) - min(sums))
    A3 = {"max_flow_identity_dev": dev_flow, "sign_antisymmetry_dev": anti,
          "pass": bool(dev_flow < 1e-13 and anti < 1e-13)}

    # ---- A4 θ 正本の受け渡し ----
    PACK = tuple(range(1, 18))
    case = base.explicit_packet_case(mode="g5q_std", packet_a=PACK, packet_b=PACK)
    a_std = base.make_case_state(sp, case, "A", hair_enabled=True)
    b_std = base.make_case_state(sp, case, "B", hair_enabled=True)
    ro = toy.theta_from_ab(a_std, b_std, sp)
    th_std = float(ro.theta)
    fl = G.g_pair_flow(a_std, b_std, th_std)
    a2 = a_std * np.cos(th_std) - b_std * np.sin(th_std)
    b2 = a_std * np.sin(th_std) + b_std * np.cos(th_std)
    dNb = float(np.vdot(b2, b2).real - np.vdot(b_std, b_std).real)
    A4 = {"theta": th_std, "dev": abs(dNb - fl["dN_B_pred"]),
          "pass": bool(abs(dNb - fl["dN_B_pred"]) < 1e-13)}

    # ---- A5 ゲージ構造 ----
    a = 0.6 * single(1, "A", False) + 0.8 * a5
    b = 0.6 * b1 + 0.8 * np.exp(1j * 0.9) * b7_0
    p0 = G.g_pair_panel(a, b, 0.2, nc, ne)
    g = np.exp(1j * 1.3)
    p1 = G.g_pair_panel(g * a, g * b, 0.2, nc, ne)
    keys = ("N_A", "N_B", "overlap_re", "overlap_im", "dN_B_pred")
    gauge_dev = max(abs(p1[k] - p0[k]) for k in keys)
    p2 = G.g_pair_charge_phase(a, -b, nc, ne)
    ph0 = G.g_pair_charge_phase(a, b, nc, ne)
    flip_dev = max(abs(np.angle(np.exp(1j * (p2["phi_ch"][ch]
                                             - ph0["phi_ch"][ch] - np.pi))))
                   for ch in ph0["phi_ch"])
    A5 = {"common_gauge_dev": float(gauge_dev),
          "one_side_pi_flip_dev": float(flip_dev),
          "pass": bool(gauge_dev < 1e-14 and flip_dev < 1e-12)}

    res = {"A1_null_theorem": A1, "A2_ladder": A2, "A3_flow_identity": A3,
           "A4_theta_canonical": A4, "A5_gauge_structure": A5,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("A1_null_theorem", "A2_ladder", "A3_flow_identity",
                           "A4_theta_canonical", "A5_gauge_structure"))
    (HERE / "qualification_pair_v5_result.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("A1_null_theorem", "A2_ladder", "A3_flow_identity",
              "A4_theta_canonical", "A5_gauge_structure"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  {res[k]}")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
