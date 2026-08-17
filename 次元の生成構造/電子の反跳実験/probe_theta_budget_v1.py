#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""θ の収支——1 步を演算ごとに切って、何が θ と A のパワーを動かすか測る

読んで分かること（先に確認済み・本スクリプトはこれを実測で確かめる）
--------------------------------------------------------------------
θ は AB 合成の χ パワー |A_k|²+|B_k|² のうち |k|≥4 かつ偶数の占有率だけで
決まる（theta_from_ab）。各演算のこの量への作用は:

  並進 k_translate_flat : A_k → A_k·e^{ikω₁}。|A_k| 不変 → θ を動かさないはず。
  弾性 rotate_ab(θ)     : 大域スカラー回転。ビンごとに合成パワー不変
                          → θ を動かさないはず。
  非弾性 φ=2r·Im(b̄a)    : 位置依存の回転なので k ビンを混ぜる
                          → これだけが θ を動かすはず。

また R′² = T²（cone_components）は **A 単体の総パワー**の二乗。並進は位相
だけなので不変、A↔B のパワー移動でのみ動くはず。

測ること
--------
  T1 1 步を4点で切り、各点で θ・P_A・P_B・Re⟨a,b⟩ を測る。
       (0) 步頭  (1) 並進後  (2) 弾性回転後  (3) 非弾性回転後＝步末
     各演算がもたらす Δθ・ΔP_A を分離する。
  T2 CR0（並進なし・対照）と CR4（並進あり）を同じ形式で比べる。
  T3 3ケース（case17_3 / case3_3 / case17_17）全部。
     case3_3 は r の変動が最大（3.65%）なのに R′ が凍る。その内訳を見る。

推測はしない。上の「はず」が外れたらそれを記録する。

使い方: python3 probe_theta_budget_v1.py
出力  : result_theta_budget_v1.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent / "統一万能関数_v1"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


_uni = _load("uni_tb", UNI / "unified_interaction_v1.py")
K = _load("kin_tb", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_tb", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_tb", HERE / "run_cr1_kinetic_feedback_v1.py")
toy = _cr1.toy

DEG_A, DEG_B = -30.0, +30.0
OMEGA0 = np.pi / 72.0
CASES = {
    "case17_3":  (tuple(range(1, 18)), tuple(range(1, 4))),
    "case3_3":   (tuple(range(1, 4)),  tuple(range(1, 4))),
    "case17_17": (tuple(range(1, 18)), tuple(range(1, 18))),
}


def probe(a, b, sp, n_chi, n_eta):
    """その瞬間の θ・A のパワー・B のパワー・重なり・R′。"""
    ro = toy.theta_from_ab(a, b, sp)
    Rp, _, _, _ = _cr1.cone_components(a, n_chi, n_eta)
    return {"theta": float(ro.theta),
            "PA": float(np.vdot(a, a).real),
            "PB": float(np.vdot(b, b).real),
            "ovR": float(np.vdot(a, b).real),
            "ovI": float(np.vdot(a, b).imag),
            "Rp": float(Rp)}


def make_ab(sp, slope, icept, tag, pa, pb):
    case = _uni.two_body_base.explicit_packet_case(
        mode=tag, packet_a=pa, packet_b=pb,
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    return (_uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True),
            _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True))


def run_budget(sp, n_chi, n_eta, slope, icept, tag, translate, T):
    """1 步を4点で切って収支を取る。translate=False で CR0 対照。"""
    pa_, pb_ = CASES[tag]
    a, b = make_ab(sp, slope, icept, "tb_" + tag, pa_, pb_)
    omega, v = OMEGA0, 0.0
    acc = {k: [] for k in ("d_tr", "d_el", "d_in",      # θ の増分
                           "p_tr", "p_el", "p_in",      # P_A の増分
                           "theta", "PA", "PB", "ovR", "Rp")}
    for _ in range(T):
        s0 = probe(a, b, sp, n_chi, n_eta)

        # --- 並進（CR4 のみ）---
        if translate:
            chi = float(np.angle(np.exp(1j * (
                _cr0.circle_position(a, n_chi, n_eta)[0]
                - _cr0.circle_position(b, n_chi, n_eta)[0]))))
            r_now = float(toy.theta_from_ab(a, b, sp).reflection_rate)
            acc_ = -4.0 * np.sin(omega / 2.0) ** 2 * chi
            v += acc_
            omega += (1.0 - r_now) * acc_
            a = K.k_translate_flat(a, -v, n_chi, n_eta)
        s1 = probe(a, b, sp, n_chi, n_eta)

        # --- 弾性回転（collision_step_exact の前半と同一）---
        ro = toy.theta_from_ab(a, b, sp)
        a, b = toy.rotate_ab(a, b, ro.theta)
        s2 = probe(a, b, sp, n_chi, n_eta)

        # --- 非弾性回転（後半と同一）---
        r = float(ro.reflection_rate)
        if r > 0.0:
            phi = 2.0 * r * np.imag(np.conj(b) * a)
            c, s_ = np.cos(phi), np.sin(phi)
            a, b = c * a - s_ * b, s_ * a + c * b
        s3 = probe(a, b, sp, n_chi, n_eta)

        acc["d_tr"].append(s1["theta"] - s0["theta"])
        acc["d_el"].append(s2["theta"] - s1["theta"])
        acc["d_in"].append(s3["theta"] - s2["theta"])
        acc["p_tr"].append(s1["PA"] - s0["PA"])
        acc["p_el"].append(s2["PA"] - s1["PA"])
        acc["p_in"].append(s3["PA"] - s2["PA"])
        for k in ("theta", "PA", "PB", "ovR", "Rp"):
            acc[k].append(s3[k])
    return {k: np.asarray(v_) for k, v_ in acc.items()}


def main() -> None:
    t0 = time.time()
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    T = 400
    out = {}

    for translate, label in ((False, "CR0 対照（並進なし）"),
                             (True, "CR4 本走行（並進あり）")):
        print(f"\n{'='*96}\n{label}  T={T}\n{'='*96}")
        print("【θ を動かすのはどの演算か】1 步あたりの |Δθ| 平均［rad］")
        print(f"  {'ケース':11} {'並進':>12} {'弾性回転':>12} {'非弾性回転':>12} "
              f"{'合計':>12} {'θ全変動%':>9}")
        for tag in CASES:
            H = run_budget(sp, n_chi, n_eta, slope, icept, tag, translate, T)
            out[f"{'cr4' if translate else 'cr0'}_{tag}"] = {
                k: [float(x) for x in v[:50]] for k, v in H.items()}
            tot = H["theta"]
            print(f"  {tag:11} {np.abs(H['d_tr']).mean():12.4e} "
                  f"{np.abs(H['d_el']).mean():12.4e} "
                  f"{np.abs(H['d_in']).mean():12.4e} "
                  f"{np.abs(H['d_tr']+H['d_el']+H['d_in']).mean():12.4e} "
                  f"{100*tot.std()/tot.mean():9.4f}")
            out[f"{'cr4' if translate else 'cr0'}_{tag}_summary"] = {
                "d_tr": float(np.abs(H["d_tr"]).mean()),
                "d_el": float(np.abs(H["d_el"]).mean()),
                "d_in": float(np.abs(H["d_in"]).mean()),
                "theta_cv": float(100 * tot.std() / tot.mean()),
                "PA_cv": float(100 * H["PA"].std() / H["PA"].mean()),
                "Rp_cv": float(100 * H["Rp"].std() / H["Rp"].mean()),
                "ovR_absmax": float(np.abs(H["ovR"]).max()),
                "PA_first": float(H["PA"][0]), "PA_last": float(H["PA"][-1]),
                "PB_first": float(H["PB"][0]), "PB_last": float(H["PB"][-1]),
                "p_tr": float(np.abs(H["p_tr"]).mean()),
                "p_el": float(np.abs(H["p_el"]).mean()),
                "p_in": float(np.abs(H["p_in"]).mean())}

        print("\n【A のパワーを動かすのはどの演算か】1 步あたりの |ΔP_A| 平均")
        print(f"  {'ケース':11} {'並進':>12} {'弾性回転':>12} {'非弾性回転':>12} "
              f"{'P_A 変動%':>10} {'R′ 変動%':>10} {'|Re<a,b>|最大':>13}")
        for tag in CASES:
            s = out[f"{'cr4' if translate else 'cr0'}_{tag}_summary"]
            print(f"  {tag:11} {s['p_tr']:12.4e} {s['p_el']:12.4e} "
                  f"{s['p_in']:12.4e} {s['PA_cv']:10.5f} {s['Rp_cv']:10.5f} "
                  f"{s['ovR_absmax']:13.4e}")

        print("\n【P_A と P_B の推移】（合計が保存し、片方に偏るかを見る）")
        print(f"  {'ケース':11} {'P_A 初→終':>26} {'P_B 初→終':>26} {'合計 初→終':>22}")
        for tag in CASES:
            s = out[f"{'cr4' if translate else 'cr0'}_{tag}_summary"]
            print(f"  {tag:11} {s['PA_first']:12.8f}→{s['PA_last']:12.8f} "
                  f"{s['PB_first']:12.8f}→{s['PB_last']:12.8f} "
                  f"{s['PA_first']+s['PB_first']:10.7f}→"
                  f"{s['PA_last']+s['PB_last']:10.7f}")

    p = HERE / "result_theta_budget_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n保存 {p.name}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
