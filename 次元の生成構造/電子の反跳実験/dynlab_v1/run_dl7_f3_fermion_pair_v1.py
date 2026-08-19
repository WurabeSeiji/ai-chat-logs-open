#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F3 走行 v1 — フェルミオン型（帯電種型）対での±則（実行前固定・G v5 使用）

木原指摘（2026-08-19）: F1/F2 の状態は奇偶混合パケットの一成分に位相を貼っただけで
フェルミオン型の波形設計になっていない。帯電種の実測条件（周期表 §5.2）は
χパリティ偶＝**奇数倍音のみ**の種であり、そのとき位相が {0,π} に Z₂ ロックする。

設計（実行前固定）:
  a, b = 純奇数倍音パケット（1,3,5,...,17・9本・搬送波なし＝共有チャネル活性）
  電荷ラベル = b 全体の Z₂ シート位相 φ∈{0, π}（種全体・一成分ではない）
  ※搬送波なし純種のため θ は状態駆動でなくプローブ回転 θ_p=0.2（計器設定・
    三部作 §6.3 の方式）を第一走行とし、状態駆動（collision_step_exact）を対照に併走。

判定:
  F3-0 ロック（前提）: 回帰列位相（s=Im(b̄a) 型・G v5 overlap_phase の回帰点列）が
       {0,π} 近傍に量子化されるか（周期表 §5.2 の再現）。立たなければ力判定に進まない
  F3-1 正味流: プローブ回転反復下の累積移乗 ΣdN_B が φ=0 と π で符号反転するか
  F3-2 分離分岐（状態駆動対照）: 周辺化分離の平均が φ=0/π で分岐するか

出力: result_dl7_f3_v1.json・dl7_f3_series_v1.npz
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
ODD = tuple(range(1, 18, 2))          # 1,3,...,17（純フェルミオン帯）
T_PROBE = 2000
T_DYN = 12000
TH_P = 0.2

def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    sys.modules[n] = m; s.loader.exec_module(m); return m

def main():
    t0 = time.time()
    G = _load("g5_f3", UNI / "unified_readout_v5.py")
    uni = _load("uni_f3", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f3", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base; step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    def fermion(which, hair, shift_deg):
        c = base.explicit_packet_case(
            mode=f"f3_{which}_{hair}", packet_a=ODD, packet_b=ODD,
            packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
            packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
        s = base.make_case_state(sp, c, which, hair_enabled=hair)
        return s / np.sqrt(np.vdot(s, s).real)

    # ---- F3-1 プローブ回転走行（搬送波なし・純種・φ = Z₂ シート） ----
    a0 = fermion("A", False, -30.0)
    b0 = fermion("B", False, +30.0)
    ov0 = G.g_pair_overlap(a0, b0)
    print(f"  純フェルミオン対の共有: |<a|b>|={ov0['overlap_abs']:.6f} "
          f"phase={np.degrees(ov0['overlap_phase']):.1f}°")

    def probe_run(phi0):
        a = a0.copy(); b = np.exp(1j * phi0) * b0
        cum = np.empty(T_PROBE); phi_ser = np.empty(T_PROBE)
        s_cum = 0.0
        for t in range(T_PROBE):
            fl = G.g_pair_flow(a, b, TH_P)
            s_cum += fl["flow_overlap_term"]
            cum[t] = s_cum
            phi_ser[t] = fl["overlap_phase"]
            a, b = (a * np.cos(TH_P) - b * np.sin(TH_P),
                    a * np.sin(TH_P) + b * np.cos(TH_P))
        return cum, phi_ser

    cum0, phis0 = probe_run(0.0)
    cumP, phisP = probe_run(np.pi)
    F31_anti = float(np.max(np.abs(cum0 + cumP)))
    F31 = {"cum_end_phi0": float(cum0[-1]), "cum_end_phiPi": float(cumP[-1]),
           "antisymmetry_max": F31_anti,
           "pass": bool(F31_anti < 1e-10 * max(1.0, abs(cum0[-1])) or
                        (np.sign(cum0[-1]) != np.sign(cumP[-1]) and
                         abs(cum0[-1]) > 1e-6))}

    # ---- F3-0/F3-2 状態駆動走行（Z₂ ロックの検査＋分離分岐） ----
    def dyn_run(phi0):
        a = fermion("A", True, -30.0)
        b = np.exp(1j * phi0) * fermion("B", True, +30.0)
        sep = np.empty(T_DYN); ph = np.empty(T_DYN)
        for t in range(T_DYN):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            sep[t] = abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))
            ph[t] = G.g_pair_overlap(a, b)["overlap_phase"]
        return sep, ph

    s0, ph0d = dyn_run(0.0)
    sP, phPd = dyn_run(np.pi)
    # F3-0: 位相の {0,π} 量子化度（後半窓・cos² 距離: 1=完全量子化, 0.5=一様）
    def z2_score(ph):
        late = ph[T_DYN // 2:]
        return float(np.mean(np.cos(late) ** 2))
    F30 = {"z2_score_phi0": z2_score(ph0d), "z2_score_phiPi": z2_score(phPd),
           "note": "cos²(位相) の後半平均。{0,π} 量子化なら→1、一様なら0.5",
           "locked": bool(min(z2_score(ph0d), z2_score(phPd)) > 0.9)}
    NBLK = 12
    def blk(x):
        b_ = x.reshape(NBLK, -1).mean(axis=1)
        return float(b_.mean()), float(b_.std(ddof=1) / np.sqrt(NBLK))
    m0, e0 = blk(s0); mP, eP = blk(sP)
    d = m0 - mP; se = float(np.hypot(e0, eP))
    F32 = {"sep_phi0": [m0, e0], "sep_phiPi": [mP, eP],
           "diff": d, "sigma": abs(d) / max(se, 1e-300),
           "split": bool(abs(d) > 3 * se)}

    res = {"config": {"packet_odd": list(ODD), "theta_probe": TH_P,
                      "T_probe": T_PROBE, "T_dyn": T_DYN,
                      "shared_abs": ov0["overlap_abs"]},
           "F3_0_lock": F30, "F3_1_net_flow_probe": F31, "F3_2_separation": F32,
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f3_series_v1.npz",
                        cum0=cum0, cumP=cumP, phis0=phis0, phisP=phisP,
                        sep0=s0, sepP=sP, ph0=ph0d, phP=phPd)
    (HERE / "result_dl7_f3_v1.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"F3-0 Z₂ロック: score(φ=0)={F30['z2_score_phi0']:.3f} "
          f"score(π)={F30['z2_score_phiPi']:.3f} locked={F30['locked']}")
    print(f"F3-1 累積移乗: φ=0: {cum0[-1]:+.4e}  φ=π: {cumP[-1]:+.4e}  "
          f"反対称max={F31_anti:.2e}  pass={F31['pass']}")
    print(f"F3-2 分離: φ=0: {m0:.3f}±{e0:.3f}°  π: {mP:.3f}±{eP:.3f}°  "
          f"差={d:+.4f}° ({F32['sigma']:.2f}σ) split={F32['split']}")
    print(f"({res['elapsed_sec']:.0f}s)")

if __name__ == "__main__":
    main()
