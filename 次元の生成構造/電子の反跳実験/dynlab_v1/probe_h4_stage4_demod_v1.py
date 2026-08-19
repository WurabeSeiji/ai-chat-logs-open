#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H4 第四段 v1 — 復調セクター検定：交差チャネルは相対巻き Δw の符号に依存するか

第一〜三段の総括と訂正（2026-08-19）:
  誤り1: ±3 の巻き刻印は 3≡0 (mod 3) で電荷中性——電荷対の比較になっていなかった。
  誤り2: 測った観測量（η 総和の χ 分離）は連続時計復調＝重力セクターの読出しであり、
         巻きに厳密盲目（二定理）。これは重力普遍性の導出であって、±則の不在証明ではない。
  訂正: セクター＝復調の違い（三種の観測時計）。交差チャネルは η モード Δw=w_A−w_B に
         厳密に生き残っており、±則は巡回復調成分に住むはず。

本プローブ:
  A の巻きを w_A=4 に固定し、B の巻きを w_B=2（Δw=+2）と w_B=6（Δw=−2）で比較。
  電荷読み（mod 3）: w=4≡+1, w=2≡−1, w=6≡0。さらに (4,1)（Δw=+3, 同電荷+1）も置く。
  観測量:
    (G) 重力側: η 総和の χ 分離（既知: Δw≠0 なら全ケース一致のはず——対照）
    (Q) 復調側: η モード分解の χ パワー分布とその重心・および
        チャネル間干渉のモード内容（a の η スペクトルの時間変化）
  判定:
    D-1 (G) は Δw=±2 で一致（重力側盲目の再確認・対照）
    D-2 (Q) は Δw=+2 と Δw=−2 で異なるか——異なれば「符号を運ぶ復調チャネルの存在」が確定

力学は二体正本 collision_step_exact のみ。全て保存・決定的。
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


_uni = _load("uni_h45", UNI / "unified_interaction_v1.py")
_cr0 = _load("cr0_h45", EXP / "run_cr0_control_no_theta_v2.py")
base = _uni.two_body_base
step = _uni.collision_step_exact

PACKET = tuple(range(1, 18))
T_STEPS = 200
# (mA, mB): 素の巻き A:+1, B:+2 → w_A=1+mA, w_B=2+mB
CASES = {
    "dw_p2": (3, 0),   # w=(4,2)  Δw=+2  電荷(+1,−1)
    "dw_m2": (3, 4),   # w=(4,6)  Δw=−2  電荷(+1, 0)
    "dw_p3": (3, -1),  # w=(4,1)  Δw=+3  電荷(+1,+1)（同電荷・巻き不等）
}


def make_ab(sp, slope, icept, mA, mB, nc, ne):
    case = base.explicit_packet_case(
        mode=f"h4s4_{mA}_{mB}", packet_a=PACKET, packet_b=PACKET,
        packet_a_shift=_cr0.shift_for_deg(-30.0, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(+30.0, slope, icept))
    a = base.make_case_state(sp, case, "A", hair_enabled=True)
    b = base.make_case_state(sp, case, "B", hair_enabled=True)
    eta = 2.0 * np.pi * np.arange(ne) / ne
    a = (a.reshape(nc, ne) * np.exp(1j * mA * eta)[None, :]).reshape(-1)
    b = (b.reshape(nc, ne) * np.exp(1j * mB * eta)[None, :]).reshape(-1)
    return a, b


def eta_mode_powers(psi, nc, ne):
    """η モード分解: P[m] = Σ_χ |FFT_η(ψ)(χ,m)|²（符号付きモード順に並べ替え）。"""
    f = np.fft.fft(psi.reshape(nc, ne), axis=1)
    P = np.sum(np.abs(f) ** 2, axis=0)
    m = ((np.arange(ne) + ne // 2) % ne) - ne // 2
    order = np.argsort(m)
    return m[order], P[order] / P.sum()


def run_case(sp, slope, icept, mA, mB, nc, ne):
    a, b = make_ab(sp, slope, icept, mA, mB, nc, ne)
    seps, aspec = [], []
    for _ in range(T_STEPS):
        a, b, _ = step(a, b, sp)
        ta, _ = _cr0.circle_position(a, nc, ne)
        tb, _ = _cr0.circle_position(b, nc, ne)
        seps.append(abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb)))))))
        _, Pa = eta_mode_powers(a, nc, ne)
        aspec.append(Pa)
    return np.array(seps), np.array(aspec)


def main():
    t0 = time.time()
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, nc, ne)
    print(f"格子 {nc}x{ne}  T={T_STEPS}")

    res = {}
    for name, (mA, mB) in CASES.items():
        seps, aspec = run_case(sp, slope, icept, mA, mB, nc, ne)
        res[name] = {"seps": seps, "aspec": aspec}
        wA, wB = 1 + mA, 2 + mB
        print(f"[{name}] w=({wA},{wB}) 電荷(mod3)=({wA % 3},{wB % 3})  "
              f"sep端={seps[-1]:.3f}°  aのηスペクトル端(上位3modes): "
              f"{np.argsort(res[name]['aspec'][-1])[::-1][:3]}")

    # D-1 重力側（η総和の分離角）: Δw=+2 vs −2 で一致するか
    d1 = float(np.max(np.abs(res["dw_p2"]["seps"] - res["dw_m2"]["seps"])))
    # D-2 復調側（a の η モードスペクトルの時系列）: Δw=+2 vs −2 で異なるか
    # モードの相対ラベルを揃える: dw_p2 の a(w=4) と dw_m2 の a(w=4) は同じ巻き。
    d2 = float(np.max(np.abs(res["dw_p2"]["aspec"] - res["dw_m2"]["aspec"])))
    # 参考: 同電荷(+1,+1)・巻き不等(Δw=+3) の重力側が他と一致するか
    d3 = float(np.max(np.abs(res["dw_p2"]["seps"] - res["dw_p3"]["seps"])))

    print(f"\nD-1 重力側分離角 |Δw=+2 − Δw=−2| = {d1:.3e}  （盲目なら機械精度）")
    print(f"D-2 復調側ηスペクトル |Δw=+2 − Δw=−2| = {d2:.3e}  （符号チャネルが在れば有限）")
    print(f"D-3 参考 |Δw=+2 − Δw=+3| (重力側) = {d3:.3e}")
    if d1 < 1e-9 and d2 > 1e-6:
        print("→ 判定: 重力側は盲目のまま、復調側は Δw の符号を区別する"
              "——符号を運ぶチャネルは復調セクターに実在する")
    elif d2 < 1e-9:
        print("→ 判定: 復調側も盲目（この観測量では符号チャネル検出されず）")
    else:
        print("→ 判定: 混在。診断が必要")

    out = {"T": T_STEPS,
           "D1_grav_sep_p2_vs_m2": d1,
           "D2_demod_spec_p2_vs_m2": d2,
           "D3_grav_sep_p2_vs_p3": d3,
           "cases": {k: {"sep_end": float(v["seps"][-1]),
                         "sep_head": [float(x) for x in v["seps"][:6]],
                         "aspec_end": [float(x) for x in v["aspec"][-1]]}
                     for k, v in res.items()},
           "elapsed_sec": time.time() - t0}
    (HERE / "result_h4_stage4_demod_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"保存: result_h4_stage4_demod_v1.json ({out['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
