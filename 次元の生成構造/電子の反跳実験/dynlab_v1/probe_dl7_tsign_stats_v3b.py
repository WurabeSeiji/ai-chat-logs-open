#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7 v3b — t 符号則の統計決着：長時間走行とブロック平均（実行前固定）

v3 の結果: t符号応答は実在（dev(200)=0.43°）、相対反転対称は機械精度（6e-13）、
しかし平均分離差 0.01° は統計未決着。大域 C は非共変（力学が時間の向きを持つ）。

本走行: T=12000、ブロック平均（12ブロック×1000步）で各ケースの平均分離と
標準誤差を出し、力水準の判定を決着させる。

ケース: base=(a,b) 同t符号 / tflip=(a,C(b)) 逆t符号 / cc=(C(a),C(b)) 同t符号・大域反転
判定（実行前固定）:
  U1 力の分岐: |mean(base) − mean(tflip)| > 3×(SE合成) なら t符号の力が統計的に確定。
     向きも記録（W11 予言: base（同符号）の方が分離大＝斥力側）
  U2 相対符号性: |mean(cc) − mean(base)| vs |mean(cc) − mean(tflip)|——
     相対符号だけが効くなら cc は base 側に来る。大域向きが効くなら別──記録
  U3 振動周期の比較: 各ケースの支配周期（束縛の強さの独立指標）

出力: result_dl7_tsign_v3b.json・dl7_tsign_series_v3b.npz
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
PACKET = tuple(range(1, 18))
T_STEPS = 12000
NBLK = 12


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_dl7tb", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_dl7tb", EXP / "run_cr0_control_no_theta_v2.py")
    bm = uni.two_body_base
    step = uni.collision_step_exact
    sp = bm.build_source_params(bm.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)
    case = bm.explicit_packet_case(
        mode="dl7tb", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
    a0 = bm.make_case_state(sp, case, "A", hair_enabled=True)
    b0 = bm.make_case_state(sp, case, "B", hair_enabled=True)

    def run_pair(a, b):
        a = a.copy(); b = b.copy()
        seps = np.empty(T_STEPS)
        for t in range(T_STEPS):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            seps[t] = abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))
        return seps

    S = {"base": run_pair(a0, b0),
         "tflip": run_pair(a0, np.conj(b0)),
         "cc": run_pair(np.conj(a0), np.conj(b0))}

    def block_stats(x):
        blocks = x.reshape(NBLK, -1).mean(axis=1)
        return float(blocks.mean()), float(blocks.std(ddof=1) / np.sqrt(NBLK))

    def period_of(x):
        y = x - x.mean()
        Y = np.abs(np.fft.rfft(y * np.hanning(len(y))))
        f = np.fft.rfftfreq(len(y)); Y[0] = 0
        return float(1.0 / f[int(np.argmax(Y))])

    stats = {k: block_stats(v) for k, v in S.items()}
    periods = {k: period_of(v) for k, v in S.items()}
    diff = stats["base"][0] - stats["tflip"][0]
    se = float(np.hypot(stats["base"][1], stats["tflip"][1]))
    U1 = {"mean": {k: stats[k][0] for k in stats},
          "se": {k: stats[k][1] for k in stats},
          "base_minus_tflip": diff, "combined_se": se,
          "significance_sigma": abs(diff) / max(se, 1e-300),
          "split": bool(abs(diff) > 3 * se),
          "direction": ("同t符号の分離大（W11: 斥力側）" if diff > 0
                        else "逆t符号の分離大（W11 と逆）") if abs(diff) > 3 * se
          else "統計的に未分岐"}
    d_cb = abs(stats["cc"][0] - stats["base"][0])
    d_ct = abs(stats["cc"][0] - stats["tflip"][0])
    U2 = {"cc_minus_base": stats["cc"][0] - stats["base"][0],
          "cc_minus_tflip": stats["cc"][0] - stats["tflip"][0],
          "closer_to": "base（相対符号のみが効く）" if d_cb < d_ct
          else "tflip（大域向きも効く）"}
    U3 = {"period_steps": periods}

    res = {"config": {"T": T_STEPS, "n_blocks": NBLK, "grid": [nc, ne]},
           "U1_force_split": U1, "U2_relative_vs_global": U2, "U3_periods": U3,
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_tsign_series_v3b.npz", **S)
    (HERE / "result_dl7_tsign_v3b.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"U1: base={stats['base'][0]:.4f}±{stats['base'][1]:.4f}° "
          f"tflip={stats['tflip'][0]:.4f}±{stats['tflip'][1]:.4f}° "
          f"cc={stats['cc'][0]:.4f}±{stats['cc'][1]:.4f}°")
    print(f"    差(base−tflip)={diff:+.4f}°  有意度={U1['significance_sigma']:.2f}σ  "
          f"→ {U1['direction']}")
    print(f"U2: cc−base={U2['cc_minus_base']:+.4f}° cc−tflip={U2['cc_minus_tflip']:+.4f}° "
          f"→ {U2['closer_to']}")
    print(f"U3: 周期 base={periods['base']:.2f} tflip={periods['tflip']:.2f} "
          f"cc={periods['cc']:.2f} 步")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
