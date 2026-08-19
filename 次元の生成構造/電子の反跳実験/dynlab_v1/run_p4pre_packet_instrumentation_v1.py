#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第4編 予備走行 v1 — パケット計装検算（第3編 修正3 で登録した L4/L5 変種・実行前固定）

目的: 読出し計装が既存系列と同一に働くことの検算。等振幅倍音 1..17 パケットに対する
厳密有理数を、二体正本計器と N 体エンジン計器の両方で検定する。

【実行前固定の予言（各計器の定義から導出）】
  PRE1 二体 r（toy.theta_from_ab・搬送波±1シフト計器）: r = 15/34
       （規則: 倍音3→1/2・奇数≥5→1・他0、w=8.5/17 → (8.5/17+8.5/17)/2＝probe_harmonic_
        composition_v1 で確立済みの規則の単一ケース再現）
  PRE2 二体 |z|（circle_position・τ=0）: |z| = 1 − 1/17 = 16/17（[F1] §4.4 の固定点読出し）
  PRE3 エンジン r（UnifiedEngine._readout・帯レジスタ）: r = 9/17
       （エンジン計器には搬送波シフトがない。等振幅帯 1..17 の奇数帯 9 本／17 本、
        孤立辺では Sagg=0 なので comb=Pk → sin²(atan2(√9,√8)) = 9/17。
        15/34 との差は計器の定義差であり、両者とも定義から導出される厳密有理数）
  PRE4 エンジン |z|（帯スペクトルの χ 円一次モーメント）: |z| = 16/17

構成: 二体は probe_harmonic_composition_v1 と同一の正本経路（import・独自再実装なし）。
エンジン側は build_standard_universe(16, 0, Nn=36, Neta=8) の器に、辺 e0 の帯 1..17 を
等振幅で置いた計装検査状態を与える（計装検査のための状態構成であり、物理走行ではない）。

出力: result_p4pre_packet_instrumentation_v1.json
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

DEG_A, DEG_B = -30.0, +30.0
PACK = tuple(range(1, 18))


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_p4pre", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_p4pre", EXP / "run_cr0_control_no_theta_v2.py")
    toy = _load("cr1_p4pre", EXP / "run_cr1_kinetic_feedback_v1.py").toy

    # ---- A. 二体正本計器 ----
    sp = uni.two_body_base.build_source_params(
        uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, n_chi, n_eta)
    case = uni.two_body_base.explicit_packet_case(
        mode="p4pre", packet_a=PACK, packet_b=PACK,
        packet_a_shift=cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=cr0.shift_for_deg(DEG_B, slope, icept))
    a = uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)
    r2 = float(toy.theta_from_ab(a, b, sp).reflection_rate)
    _, za = cr0.circle_position(a, n_chi, n_eta)
    _, zb = cr0.circle_position(b, n_chi, n_eta)
    PRE1 = {"measured": r2, "predicted": 15 / 34,
            "abs_err": abs(r2 - 15 / 34), "pass": bool(abs(r2 - 15 / 34) < 1e-10)}
    PRE2 = {"measured_A": float(abs(za)), "measured_B": float(abs(zb)),
            "predicted": 16 / 17,
            "abs_err": max(abs(abs(za) - 16 / 17), abs(abs(zb) - 16 / 17)),
            "pass": bool(max(abs(abs(za) - 16 / 17), abs(abs(zb) - 16 / 17)) < 1e-10)}

    # ---- B. N 体エンジン計器 ----
    N, Nn, Neta = 16, 36, 8
    eng, _, _ = uni.build_standard_universe(N, 0.0, Nn=Nn, Neta=Neta)
    C2 = np.zeros((eng.m, Nn, Neta), complex)
    amp = 1.0 / np.sqrt(len(PACK))
    for n in PACK:
        C2[0, n, 0] = amp                      # 辺 e0・帯 1..17・毛 0
    eng.C = C2.reshape(eng.m, -1)
    R = eng._readout()
    r_eng = float(R[0] / eng.scale)
    PRE3 = {"measured": r_eng, "predicted": 9 / 17,
            "abs_err": abs(r_eng - 9 / 17), "pass": bool(abs(r_eng - 9 / 17) < 1e-12)}
    # χ 円一次モーメント（帯スペクトル → χ 像 |ψ|² の一次モーメント）
    psi = np.fft.ifft(C2[0, :, 0]) * Nn
    P = np.abs(psi) ** 2
    z = np.sum(P * np.exp(1j * 2 * np.pi * np.arange(Nn) / Nn)) / P.sum()
    PRE4 = {"measured": float(abs(z)), "predicted": 16 / 17,
            "abs_err": abs(abs(z) - 16 / 17), "pass": bool(abs(abs(z) - 16 / 17) < 1e-12)}

    res = {"config": {"packet": list(PACK), "deg": [DEG_A, DEG_B],
                      "two_body": {"high_n": 63}, "engine": {"N": N, "Nn": Nn, "Neta": Neta},
                      "note": "15/34（二体・搬送波±1シフト計器）と 9/17（エンジン帯計器・"
                              "シフトなし）は計器定義の差。いずれも定義から導出した厳密有理数"},
           "PRE1_twobody_r": PRE1, "PRE2_twobody_z": PRE2,
           "PRE3_engine_r": PRE3, "PRE4_engine_z": PRE4,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("PRE1_twobody_r", "PRE2_twobody_z",
                           "PRE3_engine_r", "PRE4_engine_z"))
    (HERE / "result_p4pre_packet_instrumentation_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("PRE1_twobody_r", "PRE2_twobody_z", "PRE3_engine_r", "PRE4_engine_z"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  {res[k]}")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
