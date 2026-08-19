#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F1b — 共有チャネル位相 φ_ch(τ) の追跡診断（実行前固定・G v5 使用）

F1 の DC(φ) が cos φ 構造を示さなかった原因の切り分け:
 仮説 H-b1: 各瞬間の流れは cos φ_ch(τ) に厳密比例（恒等式）だが、φ_ch 自身が
   歳差・拡散して初期ラベルの記憶が平均から消える（＝電荷の持続にはロックが要る。
   三部作第三部の Z₂ 量子化 {0,π} が持続条件）。
検査:
 (a) φ_mean(τ)（g_pair_charge_phase）の系列——歳差か・ロックか・拡散か
 (b) 瞬間相関: flow_overlap(τ) と sin2θ·Σc·uv·cos(φ_ch) の一致（恒等式の走行内検算）
 (c) 初期 φ0=0 と π の φ_mean(τ) の差の系列——ラベル記憶の寿命
出力: result_dl7_f1b_v1.json・dl7_f1b_series_v1.npz
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
PACK = tuple(range(1, 18))
T = 4000
W = 0.5

def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    sys.modules[n] = m; s.loader.exec_module(m); return m

def main():
    t0 = time.time()
    G = _load("g5_f1b", UNI / "unified_readout_v5.py")
    uni = _load("uni_f1b", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f1b", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base; step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)
    case_std = base.explicit_packet_case(
        mode="f1b_std", packet_a=PACK, packet_b=PACK,
        packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
    A_on = base.make_case_state(sp, case_std, "A", hair_enabled=True)
    B_on = base.make_case_state(sp, case_std, "B", hair_enabled=True)
    A_on /= np.sqrt(np.vdot(A_on, A_on).real); B_on /= np.sqrt(np.vdot(B_on, B_on).real)
    def single(k, which):
        c = base.explicit_packet_case(mode=f"f1b_{which}_{k}", packet_a=(k,), packet_b=(k,))
        s = base.make_case_state(sp, c, which, hair_enabled=False)
        return s / np.sqrt(np.vdot(s, s).real)
    A_off, B_off = single(5, "A"), single(7, "B")

    def run(phi0):
        a = np.sqrt(1 - W) * A_on + np.sqrt(W) * A_off
        b = np.sqrt(1 - W) * B_on + np.sqrt(W) * np.exp(1j * phi0) * B_off
        phim = np.empty(T); ovl = np.empty(T); pred = np.empty(T)
        for t in range(T):
            ph = G.g_pair_charge_phase(a, b, nc, ne)
            fl_ov = G.g_pair_overlap(a, b)
            a2, b2, ro = step(a, b, sp)
            th = float(ro.theta)
            fl = G.g_pair_flow(a, b, th)
            phim[t] = ph["phi_mean"]
            ovl[t] = fl["flow_overlap_term"]
            pred[t] = float(np.sin(2 * th) * fl_ov["overlap_abs"]
                            * np.cos(fl_ov["overlap_phase"]))
            a, b = a2, b2
        return phim, ovl, pred

    ph0, ov0, pr0 = run(0.0)
    phP, ovP, _ = run(np.pi)
    ident = float(np.max(np.abs(ov0 - pr0)))
    dphi = np.angle(np.exp(1j * (phP - ph0)))
    uw0 = np.unwrap(ph0)
    prec_rate = float(np.polyfit(np.arange(T), uw0, 1)[0])
    res = {"config": {"T": T, "w": W},
           "b_identity_max_dev": ident,
           "a_phase_precession_rate_rad_per_step": prec_rate,
           "a_phase_std_late": float(np.std(np.angle(np.exp(1j * ph0[T//2:])))),
           "c_label_memory": {
               "dphi_at_0": float(dphi[0]), "dphi_at_100": float(dphi[100]),
               "dphi_at_1000": float(dphi[1000]), "dphi_at_end": float(dphi[-1]),
               "dphi_late_mean_abs": float(np.mean(np.abs(dphi[T//2:])))},
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f1b_series_v1.npz",
                        phi0=ph0, phiP=phP, ovl0=ov0, pred0=pr0)
    (HERE / "result_dl7_f1b_v1.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"(b) 恒等式の走行内検算 max|ovl−pred| = {ident:.3e}")
    print(f"(a) φ_mean 歳差率 = {prec_rate:+.4e} rad/步  後半std = {res['a_phase_std_late']:.3f}")
    print(f"(c) ラベル差 Δφ: τ=0: {np.degrees(dphi[0]):+.1f}° → 100: "
          f"{np.degrees(dphi[100]):+.1f}° → 1000: {np.degrees(dphi[1000]):+.1f}° → "
          f"end: {np.degrees(dphi[-1]):+.1f}°  後半平均|Δφ|={np.degrees(res['c_label_memory']['dphi_late_mean_abs']):.1f}°")
    print(f"({res['elapsed_sec']:.0f}s)")

if __name__ == "__main__":
    main()
