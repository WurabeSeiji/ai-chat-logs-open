#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H4 裁定 第二段 v1 — 相対ドリフト方向の符号は、衝突の距離レートに効くか

第一段の結果: 巻き（Q型電荷）符号は距離動力学に厳密に無効（毛ゲージ不変性）。
更新仮説: 符号依存の力はゲージ層（相対ドリフト方向）に住む（W11 §4.3 のゲージ語版）。

本プローブはその前提条件を二体正本で直接測る:
  ドリフト（並進）と衝突を交互に適用し、衝突が距離レートに与える寄与が
  相対ドリフトの「向き」に依存するかを調べる。

ケース（1步あたりの並進 [deg]、A at −30°, B at +30°）:
  static   : (0, 0)          基準
  approach : (+v, −v)        接近（A→B, B→A）
  recede   : (−v, +v)        離反
  comove   : (+v, +v)        並走（相対速度 0——静止と一致すれば並進共変性）

各ケースで:
  sep_col(τ)  = 並進＋衝突 の分離角
  sep_free(τ) = 並進のみ（衝突なし）の分離角（運動学基準）
  dev(τ) = sep_col − sep_free   （衝突が距離に与えた正味の寄与）

判定:
  J-a 並進共変性: dev[comove] と dev[static] が一致するか
  J-b 方向符号性: dev[approach] と dev[recede] が（時間対称性を超えて）
      符号の異なる系統的寄与を持つか

力学は二体正本（collision_step_exact）と正本並進（k_translate_flat）のみ。
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


_uni = _load("uni_h42", UNI / "unified_interaction_v1.py")
K = _load("kin_h42", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_h42", EXP / "run_cr0_control_no_theta_v2.py")

base = _uni.two_body_base
collision_step_exact = _uni.collision_step_exact

PACKET = tuple(range(1, 18))
DEG_A, DEG_B = -30.0, +30.0
V_DEG = 0.5
T_STEPS = 200
CASES = {"static": (0.0, 0.0), "approach": (+V_DEG, -V_DEG),
         "recede": (-V_DEG, +V_DEG), "comove": (+V_DEG, +V_DEG)}


def sep_deg(a, b, n_chi, n_eta):
    ta, _ = _cr0.circle_position(a, n_chi, n_eta)
    tb, _ = _cr0.circle_position(b, n_chi, n_eta)
    return abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))


def make_ab(sp, slope, icept):
    case = base.explicit_packet_case(
        mode="h4s2", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = base.make_case_state(sp, case, "A", hair_enabled=True)
    b = base.make_case_state(sp, case, "B", hair_enabled=True)
    return a, b


def run_case(sp, n_chi, n_eta, slope, icept, om1_per_deg, vA, vB, with_collision):
    a, b = make_ab(sp, slope, icept)
    omA = om1_per_deg * vA
    omB = om1_per_deg * vB
    seps = []
    for _ in range(T_STEPS):
        if vA != 0.0:
            a = K.k_translate_flat(a, omA, n_chi, n_eta)
        if vB != 0.0:
            b = K.k_translate_flat(b, omB, n_chi, n_eta)
        if with_collision:
            a, b, _ = collision_step_exact(a, b, sp)
        seps.append(sep_deg(a, b, n_chi, n_eta))
    return np.array(seps)


def main():
    t0 = time.time()
    sp = base.build_source_params(
        base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)

    # 並進の較正（CR6 と同じ手順: 1° あたりの omega1 を実測で決める）
    a0, b0 = make_ab(sp, slope, icept)
    probe = np.radians(1.0)
    p0, _ = _cr0.circle_position(a0, n_chi, n_eta)
    p1, _ = _cr0.circle_position(K.k_translate_flat(a0, probe, n_chi, n_eta),
                                 n_chi, n_eta)
    moved = float(np.degrees(np.angle(np.exp(1j * (p1 - p0)))))
    om1_per_deg = probe * (1.0 if moved > 0 else -1.0) / abs(moved) * 1.0
    print(f"格子 n_chi={n_chi} n_eta={n_eta}  v={V_DEG}°/步  T={T_STEPS}")
    print(f"並進較正: 1°指定で {moved:+.6f}° 移動 → om1_per_deg={om1_per_deg:+.8f}")

    out = {"v_deg": V_DEG, "T": T_STEPS, "cases": {}}
    devs = {}
    for name, (vA, vB) in CASES.items():
        s_col = run_case(sp, n_chi, n_eta, slope, icept, om1_per_deg, vA, vB, True)
        s_free = run_case(sp, n_chi, n_eta, slope, icept, om1_per_deg, vA, vB, False)
        dev = s_col - s_free
        devs[name] = dev
        out["cases"][name] = {
            "sep_col_end": float(s_col[-1]), "sep_free_end": float(s_free[-1]),
            "dev_mean": float(dev.mean()), "dev_end": float(dev[-1]),
            "dev_min": float(dev.min()), "dev_max": float(dev.max()),
        }
        print(f"[{name:8s}] sep_col端={s_col[-1]:8.3f}°  sep_free端={s_free[-1]:8.3f}°  "
              f"dev: 平均{dev.mean():+8.4f}° 端{dev[-1]:+8.4f}°")

    # J-a 並進共変性: comove と static の dev 一致
    ja = float(np.max(np.abs(devs["comove"] - devs["static"])))
    # J-b 方向符号性: approach と recede の dev の差（系統的か）
    d_ap, d_re = devs["approach"], devs["recede"]
    jb_diff_mean = float((d_ap - d_re).mean())
    jb_corr = float(np.corrcoef(d_ap, -d_re)[0, 1]) if d_ap.std() > 0 else float("nan")
    print(f"\nJ-a 並進共変性: max|dev(comove)−dev(static)| = {ja:.3e}")
    print(f"J-b 方向符号性: mean[dev(approach)−dev(recede)] = {jb_diff_mean:+.5f}°")
    print(f"    （参考）dev(approach) と −dev(recede) の相関 = {jb_corr:+.4f}")

    out["J_a_comove_vs_static"] = ja
    out["J_b_approach_minus_recede_mean"] = jb_diff_mean
    out["J_b_corr_ap_vs_negre"] = jb_corr
    out["elapsed_sec"] = time.time() - t0
    (HERE / "result_h4_stage2_drift_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"\n保存: result_h4_stage2_drift_v1.json ({out['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
