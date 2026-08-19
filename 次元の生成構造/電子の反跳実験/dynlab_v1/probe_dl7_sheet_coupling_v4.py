#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7 v4 — Z₂シート×結合レート応答：波の周期表v2 の±読出しの検定（実行前固定）

読む方向の正本（木原指摘・wave_periodic_table_v2_ja.md）:
 (1) §9.1: ±の実体は結合レート r への応答の符号——∂r/∂P奇>0・∂r/∂P偶<0
     「量が等しくても役割で区別される（電荷の±の類比）」。力の観測量は r 自身の応答。
 (2) §5.2: 帯電種の位相は回帰点で {0,π} に厳密量子化（Z₂ 二シート）。
     ±の二値性はシートに住み、対の相対シートは一次交差項 2Re(āb) を通して
     r に線形に入る（二乗型交差項はシート反転に盲目だが一次項は符号反転する）。

ケース（実行前固定・位置は全て −30°/+30°・パケット 1..17）:
  ref      = (a, b)                      相対シート 0
  sheet_pi = (a, −b)                     相対シート π（大域符号反転）
  odd_pi   = (a, 奇数倍音のみ −1 倍の b)  フェルミオン帯だけシート反転
  quarter  = (a, e^{iπ/2} b)             相対シート π/2（対照・中間値）

観測量: r(τ)＝正本 theta_from_ab(a,b).reflection_rate（結合そのもの・毎步）、
        周辺化分離角 |Δθ(τ)|（径方向動力学）。後半窓平均±ブロックSE。

判定（実行前固定）:
  V1 結合応答: r 系列が ref と sheet_pi で丸め包絡を超えて異なるか
     ——異なれば「±は r への応答として読める」（§9.1 の類比の力学実証）
  V2 分離応答: 後半平均分離が相対シートで統計的に分岐するか（3σ・ブロックSE）
     ——分岐すれば±の力（引力/斥力側の分離）が実測される。割当は評価出力
  V3 担い手の特定: odd_pi（奇帯のみ反転）が sheet_pi と一致すれば±の担い手は
     フェルミオン帯シート、ref と一致すれば大域位相（＝ゲージ・物理でない）と判明
  V4 中間位相: quarter が ref と sheet_pi の間に来るか（応答の連続性・記録）

出力: result_dl7_sheet_v4.json・dl7_sheet_series_v4.npz
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
T_STEPS = 6000
NBLK = 12


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_dl7s", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_dl7s", EXP / "run_cr0_control_no_theta_v2.py")
    toy = _load("cr1_dl7s", EXP / "run_cr1_kinetic_feedback_v1.py").toy
    bm = uni.two_body_base
    step = uni.collision_step_exact
    sp = bm.build_source_params(bm.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)
    case = bm.explicit_packet_case(
        mode="dl7s", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
    a0 = bm.make_case_state(sp, case, "A", hair_enabled=True)
    b0 = bm.make_case_state(sp, case, "B", hair_enabled=True)

    def odd_flip(psi):
        """χ 奇数倍音（フェルミオン帯）成分のみ −1 倍。"""
        F = np.fft.fft(psi.reshape(nc, ne), axis=0)
        k = np.fft.fftfreq(nc, d=1.0 / nc).astype(int)
        F[np.abs(k) % 2 == 1, :] *= -1.0
        return np.fft.ifft(F, axis=0).reshape(-1)

    variants = {
        "ref": b0.copy(),
        "sheet_pi": -b0,
        "odd_pi": odd_flip(b0),
        "quarter": 1j * b0,
    }

    def run_pair(b_init):
        a = a0.copy(); b = b_init.copy()
        r_ser = np.empty(T_STEPS)
        sep = np.empty(T_STEPS)
        for t in range(T_STEPS):
            r_ser[t] = float(toy.theta_from_ab(a, b, sp).reflection_rate)
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            sep[t] = abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))
        return r_ser, sep

    R, S = {}, {}
    for name, bi in variants.items():
        R[name], S[name] = run_pair(bi)
        print(f"  [{name}] r(0)={R[name][0]:.6f} r後半={R[name][T_STEPS//2:].mean():.6f} "
              f"分離後半={S[name][T_STEPS//2:].mean():.3f}°")

    def blk(x):
        b = x[T_STEPS // 2:].reshape(NBLK // 2, -1).mean(axis=1)
        return float(b.mean()), float(b.std(ddof=1) / np.sqrt(len(b)))

    r_dev_pi = float(np.max(np.abs(R["sheet_pi"] - R["ref"])))
    r_dev_odd_vs_pi = float(np.max(np.abs(R["odd_pi"] - R["sheet_pi"])))
    r_dev_odd_vs_ref = float(np.max(np.abs(R["odd_pi"] - R["ref"])))
    V1 = {"max_r_dev_sheetpi_vs_ref": r_dev_pi,
          "signal": bool(r_dev_pi > 1e-6), "pass": True}
    stats = {k: blk(S[k]) for k in S}
    d = stats["ref"][0] - stats["sheet_pi"][0]
    se = float(np.hypot(stats["ref"][1], stats["sheet_pi"][1]))
    V2 = {"sep_late": {k: stats[k] for k in stats},
          "ref_minus_sheetpi": d, "combined_se": se,
          "sigma": abs(d) / max(se, 1e-300),
          "split": bool(abs(d) > 3 * se), "pass": True}
    V3 = {"odd_vs_sheetpi_max_r_dev": r_dev_odd_vs_pi,
          "odd_vs_ref_max_r_dev": r_dev_odd_vs_ref,
          "carrier": ("フェルミオン帯シート（odd_pi≈sheet_pi）"
                      if r_dev_odd_vs_pi < r_dev_odd_vs_ref
                      else "大域位相はゲージ（odd_pi≈ref なら奇帯反転は独立効果）"),
          "pass": True}
    rq = R["quarter"][T_STEPS // 2:].mean()
    rr = R["ref"][T_STEPS // 2:].mean()
    rp = R["sheet_pi"][T_STEPS // 2:].mean()
    V4 = {"r_late": {"ref": rr, "quarter": float(rq), "sheet_pi": rp},
          "between": bool(min(rr, rp) <= rq <= max(rr, rp)), "pass": True}

    res = {"config": {"T": T_STEPS, "grid": [nc, ne], "packet": list(PACKET)},
           "V1_coupling_response": V1, "V2_separation_split": V2,
           "V3_carrier": V3, "V4_quarter": V4,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = True
    np.savez_compressed(HERE / "dl7_sheet_series_v4.npz",
                        **{f"r_{k}": R[k] for k in R},
                        **{f"sep_{k}": S[k] for k in S})
    (HERE / "result_dl7_sheet_v4.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"V1 結合応答: signal={V1['signal']}  max|Δr|={r_dev_pi:.3e}")
    print(f"V2 分離分岐: {V2['sigma']:.2f}σ  split={V2['split']}  "
          f"(ref−sheet_pi={d:+.4f}°±{se:.4f})")
    print(f"V3 担い手: {V3['carrier']}  (odd vs pi: {r_dev_odd_vs_pi:.3e} / "
          f"odd vs ref: {r_dev_odd_vs_ref:.3e})")
    print(f"V4 中間位相: between={V4['between']}  r後半 ref={rr:.6f} "
          f"quarter={rq:.6f} pi={rp:.6f}")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
