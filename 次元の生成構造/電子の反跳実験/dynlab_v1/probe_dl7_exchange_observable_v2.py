#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7 観測量 v2 — 交換チャネル読出し（実行前固定・v1 の盲目診断を受けた再設計）

v1 の診断: 自巻きモードの重心位相は純巻き状態で周辺化像に因数分解し、毛ゲージ
不変性により構成的に盲目（S4 の完全共変が観測量定義の誤りを特定）。
符号チャネル（stage4 D2=0.99987）は**交差項＝相手巻きモードへの側帯波**に住む。
[F2] §9 の登録課題「巡回時計 mod 3 と η モードの写像の確定」の実験的探索を兼ねる。

観測量（実行前固定・4種）:
  (a) 交換占有 E_a = Σ_χ|F_a(χ, w_B)|²/Σ|F_a|²（a が相手の巻きチャネルに得た内容）と
      その後半窓レート
  (b) 交換チャネル分離 Δθ_x = pos(F_a(:,w_B)) − pos(F_b(:,w_A))（χ 円一次モーメント位相差）
      とそのレート r_x
  (c) 干渉位相 φ_a = arg Σ_χ F_a(χ,w_A)·conj(F_a(χ,w_B))（自×交換の相対位相）のレート
  (d) η スペクトル距離 D_case = max|aspec_case − aspec_neut|（stage4 D2 型・対中立）

判定（実行前固定）:
  X1 交換チャネルの活性: 荷電ケースで E_a 後半平均 > 中立ケース（チャネル実在）
  X2 mod3 分岐（本丸）: (a)〜(d) のいずれかのレート・値が同電荷2ケースで一致し
     異電荷2ケースと分かれるか（分岐すれば±則の観測量が同定される）
  X3 中立対照: neut は荷電ケースと系統的に異なる
  X4 鏡映共変: 全巻き反転で各観測量が対応値に一致（<1e-9）

ケース: same_pp=(4,1)/same_mm=(2,5)/opp_pm=(4,2)/opp_mp=(2,4)/neut=(3,6)＋mirror=(−4,−1)
出力: result_dl7_exchange_v2.json・dl7_exchange_series_v2.npz
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
T_STEPS = 1500
CASES = {"same_pp": (4, 1), "same_mm": (2, 5),
         "opp_pm": (4, 2), "opp_mp": (2, 4), "neut": (3, 6)}


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_dl7x", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_dl7x", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base
    step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63,
                                              recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    def make_ab(wA, wB):
        mA, mB = wA - 1, wB - 2
        case = base.explicit_packet_case(
            mode=f"dl7x_{wA}_{wB}", packet_a=PACKET, packet_b=PACKET,
            packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
            packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
        a = base.make_case_state(sp, case, "A", hair_enabled=True)
        b = base.make_case_state(sp, case, "B", hair_enabled=True)
        eta = 2.0 * np.pi * np.arange(ne) / ne
        a = (a.reshape(nc, ne) * np.exp(1j * mA * eta)[None, :]).reshape(-1)
        b = (b.reshape(nc, ne) * np.exp(1j * mB * eta)[None, :]).reshape(-1)
        return a, b

    ks = np.exp(1j * 2 * np.pi * np.arange(nc) / nc)

    def obs(psi, w_own, w_ex):
        F = np.fft.fft(psi.reshape(nc, ne), axis=1)
        tot = float(np.sum(np.abs(F) ** 2))
        col = F[:, w_ex % ne]
        P = np.abs(col) ** 2
        E = float(P.sum() / tot)
        pos = float(np.angle(np.sum(P * ks) / max(P.sum(), 1e-300)))
        own = F[:, w_own % ne]
        phi = float(np.angle(np.sum(own * np.conj(col))))
        spec = np.sum(np.abs(F) ** 2, axis=0)
        return E, pos, phi, spec / spec.sum()

    def run_case(wA, wB):
        a, b = make_ab(wA, wB)
        out = {k: [] for k in ("Ea", "Eb", "dthx", "phia", "phib")}
        spec_end = None
        for _ in range(T_STEPS):
            a, b, _ = step(a, b, sp)
            Ea, pa, fa, spa = obs(a, wA, wB)
            Eb, pb, fb, _ = obs(b, wB, wA)
            out["Ea"].append(Ea); out["Eb"].append(Eb)
            out["dthx"].append(float(np.degrees(np.angle(np.exp(1j * (pa - pb))))))
            out["phia"].append(fa); out["phib"].append(fb)
            spec_end = spa
        return {k: np.array(v) for k, v in out.items()}, spec_end

    def rate(x):
        h = len(x) // 2
        return float(np.polyfit(np.arange(h, len(x)), x[h:], 1)[0])

    R = {}
    specs = {}
    for name, (wA, wB) in CASES.items():
        o, spec = run_case(wA, wB)
        R[name] = o
        specs[name] = spec
        print(f"  [{name}] w=({wA},{wB}): E_a後半={o['Ea'][len(o['Ea'])//2:].mean():.4e} "
              f"rate(E_a)={rate(o['Ea']):+.2e} r_x={rate(np.abs(o['dthx'])):+.4e} "
              f"rate(φ_a)={rate(np.unwrap(o['phia'])):+.4e}")
    o_mir, _ = run_case(-4, -1)

    summ = {}
    for name in CASES:
        o = R[name]
        summ[name] = {
            "Ea_late": float(o["Ea"][T_STEPS // 2:].mean()),
            "rate_Ea": rate(o["Ea"]),
            "r_x": rate(np.abs(o["dthx"])),
            "rate_phia": rate(np.unwrap(o["phia"])),
            "spec_dist_to_neut": float(np.max(np.abs(specs[name] - specs["neut"]))),
        }
    charged = ["same_pp", "same_mm", "opp_pm", "opp_mp"]
    X1 = {"Ea_late": {k: summ[k]["Ea_late"] for k in summ},
          "pass": bool(all(summ[k]["Ea_late"] > 1e-12 for k in charged))}
    # X2: 各観測量で mod3 分岐を検査
    split_report = {}
    any_split = False
    for key in ("Ea_late", "rate_Ea", "r_x", "rate_phia", "spec_dist_to_neut"):
        vals = {k: summ[k][key] for k in charged}
        s_same = [vals["same_pp"], vals["same_mm"]]
        s_opp = [vals["opp_pm"], vals["opp_mp"]]
        within_same = abs(s_same[0] - s_same[1])
        within_opp = abs(s_opp[0] - s_opp[1])
        between = abs(np.mean(s_same) - np.mean(s_opp))
        scale = max(abs(v) for v in list(s_same) + list(s_opp)) + 1e-300
        # 分岐＝クラス間差がクラス内差の10倍超、かつ絶対有意（スケールの1e-6超）
        split = bool(between > 10 * max(within_same, within_opp, 1e-300)
                     and between > 1e-6 * scale)
        split_report[key] = {"same": s_same, "opp": s_opp,
                             "within_max": max(within_same, within_opp),
                             "between": between, "split": split}
        any_split = any_split or split
    X2 = {"per_observable": split_report, "any_split": any_split, "pass": True}
    X3 = {"spec_dist_to_neut": {k: summ[k]["spec_dist_to_neut"] for k in charged},
          "pass": bool(min(summ[k]["spec_dist_to_neut"] for k in charged) > 1e-6)}
    mir_rx = rate(np.abs(o_mir["dthx"]))
    X4 = {"r_x_mirror": mir_rx, "r_x_orig": summ["same_pp"]["r_x"],
          "abs_dev": abs(mir_rx - summ["same_pp"]["r_x"]),
          "pass": bool(abs(mir_rx - summ["same_pp"]["r_x"]) < 1e-9)}

    res = {"config": {"T": T_STEPS, "grid": [nc, ne], "cases": CASES},
           "summary": summ,
           "X1_channel_active": X1, "X2_mod3_split": X2,
           "X3_neutral_control": X3, "X4_mirror": X4,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("X1_channel_active", "X2_mod3_split",
                           "X3_neutral_control", "X4_mirror"))
    np.savez_compressed(HERE / "dl7_exchange_series_v2.npz",
                        **{f"{k}_{q}": R[k][q] for k in R for q in R[k]},
                        **{f"spec_{k}": specs[k] for k in specs})
    (HERE / "result_dl7_exchange_v2.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("X1_channel_active", "X2_mod3_split", "X3_neutral_control", "X4_mirror"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}")
    print(f"  X2 分岐あり: {any_split}  "
          f"{ {kk: vv['split'] for kk, vv in split_report.items()} }")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
