#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F2 走行 v1 — 距離辞書接続：φ={0,π} 対の径方向運動の分岐（実行前固定・G v5）

導出: DL7_導出ノート.md §10。判定 F2-1（分離の統計的分岐）・F2-2（向き＝
Bjerknes 型か電磁気型かの自動判定）。F1b の診断（ラベルは約100衝突で崩壊）を
受け、初期過渡（ラベル生存窓 τ<100）と長時間平均の両方を測る。

出力: result_dl7_f2_v1.json・dl7_f2_series_v1.npz
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
PACK = tuple(range(1, 18))
T = 12000
W = 0.5
NBLK = 12

def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    sys.modules[n] = m; s.loader.exec_module(m); return m

def main():
    t0 = time.time()
    uni = _load("uni_f2", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f2", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base; step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)
    case_std = base.explicit_packet_case(
        mode="f2_std", packet_a=PACK, packet_b=PACK,
        packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
    A_on = base.make_case_state(sp, case_std, "A", hair_enabled=True)
    B_on = base.make_case_state(sp, case_std, "B", hair_enabled=True)
    A_on /= np.sqrt(np.vdot(A_on, A_on).real); B_on /= np.sqrt(np.vdot(B_on, B_on).real)
    def single(k, which):
        c = base.explicit_packet_case(mode=f"f2_{which}_{k}", packet_a=(k,), packet_b=(k,))
        s = base.make_case_state(sp, c, which, hair_enabled=False)
        return s / np.sqrt(np.vdot(s, s).real)
    A_off, B_off = single(5, "A"), single(7, "B")

    def run(phi0):
        a = np.sqrt(1 - W) * A_on + np.sqrt(W) * A_off
        b = np.sqrt(1 - W) * B_on + np.sqrt(W) * np.exp(1j * phi0) * B_off
        sep = np.empty(T)
        for t in range(T):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            sep[t] = abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))
        return sep

    s0 = run(0.0)
    sP = run(np.pi)

    def blk(x):
        b = x.reshape(NBLK, -1).mean(axis=1)
        return float(b.mean()), float(b.std(ddof=1) / np.sqrt(NBLK))

    m0, e0 = blk(s0); mP, eP = blk(sP)
    d = m0 - mP; se = float(np.hypot(e0, eP))
    # 初期過渡（ラベル生存窓）
    d_early = float(s0[:100].mean() - sP[:100].mean())
    res = {"config": {"T": T, "w": W, "n_blocks": NBLK},
           "F2_1": {"sep_mean_phi0": [m0, e0], "sep_mean_phiPi": [mP, eP],
                    "diff": d, "combined_se": se,
                    "sigma": abs(d) / max(se, 1e-300),
                    "split": bool(abs(d) > 3 * se)},
           "F2_early_transient": {"diff_first100_deg": d_early,
                                  "note": "ラベル生存窓（F1b: 崩壊~100衝突）内の分離差"},
           "F2_2": ("同位相の分離が小=Bjerknes型（同位相=引力）" if (d > 0) else
                    "同位相の分離が大=電磁気型側") if abs(d) > 3 * se
           else "未分岐（長時間平均ではラベル崩壊により判定不能——F1b と整合）",
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f2_series_v1.npz", sep_phi0=s0, sep_phiPi=sP)
    (HERE / "result_dl7_f2_v1.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"F2-1: φ=0: {m0:.3f}±{e0:.3f}°  φ=π: {mP:.3f}±{eP:.3f}°  "
          f"差={d:+.4f}° ({res['F2_1']['sigma']:.2f}σ) split={res['F2_1']['split']}")
    print(f"初期過渡(τ<100)の分離差 = {d_early:+.3f}°")
    print(f"F2-2: {res['F2_2']}")
    print(f"({res['elapsed_sec']:.0f}s)")

if __name__ == "__main__":
    main()
