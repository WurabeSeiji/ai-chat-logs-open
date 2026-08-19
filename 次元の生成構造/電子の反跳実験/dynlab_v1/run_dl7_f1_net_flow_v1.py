#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F1 走行 v1 — 正味定常流：電荷項の時間平均 DC(φ)（実行前固定・G v5 使用）

導出: DL7_導出ノート.md §10。三部作第一部 §6.4(ii) の登録課題
「符号は各衝突で厳密だが Rabi 型に振動する——正味の定常流は時間平均の議論を要する
（Bjerknes 同型）」の決着。

構成（テストデータ・±電荷の表現）:
  部分共有対（§6.4(iii)）: a=√(1−w)A_on+√w A_off、b=√(1−w)B_on+√w e^{iφ}B_off
  A_on/B_on = 標準パケット（1..17・搬送波あり・∓30°）——対角θ読出しを担う
  A_off = A5（搬送波なし）、B_off = B7（搬送波なし）——共有チャネル c6=0.5 を担う
  電荷ラベル = φ（片側の共有成分位相。φシフト＝電荷共役＝状態準備）

力学: collision_step_exact（状態駆動θ・無変更）。読出し: G v5 g_pair_flow
（θ は各衝突の ro.theta＝正本読み）。

判定（実行前固定）:
  F1-1 反対称: DC の φ 奇成分が有意（DC(φ)−DC(φ+π) の半分＝電荷項 DC）で、
       DC_odd(φ+π) = −DC_odd(φ)（構成上厳密——検算）。有意性は cos φ 振幅 vs 残差
  F1-2 形: DC(φ) = C0 + C1 cos φ + S1 sin φ で高調波残差が小さい（記録）
  F1-3 対照: 大きさの項の時間平均は φ 非依存（重力文法の符号盲目）

出力: result_dl7_f1_net_flow_v1.json・dl7_f1_series_v1.npz
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
PACK = tuple(range(1, 18))
T_STEPS = 4000
W_LIST = (0.2, 0.5)
NPHI = 8


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    G = _load("g5_f1", UNI / "unified_readout_v5.py")
    uni = _load("uni_f1", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f1", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base
    step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63,
                                              recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    case_std = base.explicit_packet_case(
        mode="f1_std", packet_a=PACK, packet_b=PACK,
        packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
    A_on = base.make_case_state(sp, case_std, "A", hair_enabled=True)
    B_on = base.make_case_state(sp, case_std, "B", hair_enabled=True)
    A_on = A_on / np.sqrt(np.vdot(A_on, A_on).real)
    B_on = B_on / np.sqrt(np.vdot(B_on, B_on).real)

    def single(k, which):
        c = base.explicit_packet_case(mode=f"f1_{which}_{k}",
                                      packet_a=(k,), packet_b=(k,))
        s = base.make_case_state(sp, c, which, hair_enabled=False)
        return s / np.sqrt(np.vdot(s, s).real)

    A_off = single(5, "A")
    B_off = single(7, "B")

    phis = np.linspace(0, 2 * np.pi, NPHI, endpoint=False)
    out = {}
    for w in W_LIST:
        dc_ovl, dc_mag = [], []
        for phi in phis:
            a = np.sqrt(1 - w) * A_on + np.sqrt(w) * A_off
            b = np.sqrt(1 - w) * B_on + np.sqrt(w) * np.exp(1j * phi) * B_off
            s_ovl = 0.0
            s_mag = 0.0
            for _ in range(T_STEPS):
                a2, b2, ro = step(a, b, sp)
                fl = G.g_pair_flow(a, b, float(ro.theta))
                s_ovl += fl["flow_overlap_term"]
                s_mag += fl["flow_mag_term"]
                a, b = a2, b2
            dc_ovl.append(s_ovl / T_STEPS)
            dc_mag.append(s_mag / T_STEPS)
        out[w] = {"phi": phis.tolist(), "dc_ovl": dc_ovl, "dc_mag": dc_mag}
        print(f"  w={w}: DC_ovl(φ) = " +
              " ".join(f"{v:+.3e}" for v in dc_ovl))

    res = {"config": {"T": T_STEPS, "w_list": list(W_LIST), "n_phi": NPHI,
                      "shared_channel": "A5-B7 (c6=0.5)", "grid": [nc, ne]}}
    for w in W_LIST:
        dc = np.array(out[w]["dc_ovl"])
        dm = np.array(out[w]["dc_mag"])
        # フーリエ分解（8点）
        C0 = float(dc.mean())
        C1 = float((dc * np.cos(phis)).mean() * 2)
        S1 = float((dc * np.sin(phis)).mean() * 2)
        resid = dc - (C0 + C1 * np.cos(phis) + S1 * np.sin(phis))
        dc_odd = 0.5 * (dc - np.roll(dc, NPHI // 2))
        amp1 = float(np.hypot(C1, S1))
        res[f"w={w}"] = {
            "DC_ovl": dc.tolist(), "DC_mag": dm.tolist(),
            "C0": C0, "C1_cos": C1, "S1_sin": S1,
            "harmonic_residual_max": float(np.max(np.abs(resid))),
            "F1_1_odd_amplitude": float(np.max(np.abs(dc_odd))),
            "F1_1_signal": bool(amp1 > 10 * np.max(np.abs(resid))),
            "F1_3_mag_phi_dependence": float(dm.max() - dm.min()),
            "F1_3_mag_mean": float(dm.mean()),
        }
    res["elapsed_sec"] = time.time() - t0
    np.savez_compressed(HERE / "dl7_f1_series_v1.npz",
                        **{f"dc_ovl_w{w}": np.array(out[w]["dc_ovl"]) for w in W_LIST},
                        **{f"dc_mag_w{w}": np.array(out[w]["dc_mag"]) for w in W_LIST},
                        phi=phis)
    (HERE / "result_dl7_f1_net_flow_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for w in W_LIST:
        r = res[f"w={w}"]
        print(f"  w={w}: C0={r['C0']:+.3e} C1cos={r['C1_cos']:+.3e} "
              f"S1sin={r['S1_sin']:+.3e} 残差max={r['harmonic_residual_max']:.2e} "
              f"奇振幅={r['F1_1_odd_amplitude']:.3e} signal={r['F1_1_signal']}")
        print(f"        大きさ項のφ依存={r['F1_3_mag_phi_dependence']:.2e} "
              f"(平均 {r['F1_3_mag_mean']:+.3e})")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
