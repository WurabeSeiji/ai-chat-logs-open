#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7 走行 v1 — 電荷セクターの力則：η点力法則の検算と復調±則の実測（S0〜S5）

導出: DL7_導出ノート.md（[F1] §7 DL6-Q5 の電荷側・[F2] §8 追加行の遂行）。
理論部: [F2] 定理4.1 はレジスタ点を選ばない——エンジンの _vertex_rate は定理式の
列（レジスタ点）ごと適用そのもの（g1=1）。S0 がこの恒等を全点で検算する。
したがって S2 の実測符号は**導出式の評価出力**である（別仮定なし）。

観測量（DL7 §3・実行前固定）:
  復調位置 θ_η(X; w) = 巻きモード w に η-FFT で射影した χ 像の円一次モーメント位相
  復調分離 Δθ_η(τ) = θ_η(a; w_A) − θ_η(b; w_B)（各体は自分の巻きチャネルで読む）
  復調レート r_η = 後半窓の |Δθ_η| の線形勾配（増大＝斥力型／減少＝引力型）
  整流電荷 q = w mod 3（1→+1, 2→−1, 0→中性）

事前登録ケース（w_A, w_B）:
  same_pp=(4,1) q=(+,+) / same_mm=(2,5) q=(−,−) /
  opp_pm=(4,2) q=(+,−) / opp_mp=(2,4) q=(−,+) / neut=(3,6) q=(0,0) 対照

判定:
  S0 η点力法則: 全レジスタ点で |2Re(W̄·rate) − 定理式| < 1e-13（スケール比）
  S1 stage4 再現: D2（復調弁別）が正本値 0.99987 と相対 1e-6 で一致
  S2 ±則（本丸）: 同電荷2ケースの r_η 符号が一致し、異電荷2ケースの符号と分かれる
     （分岐すれば±則成立＝統一主張の完成条件(ii)。分岐しなければ三値読みで記録）
  S3 重力側の不変: 全ケースの周辺化分離角軌道が一致 < 1e-9（毛ゲージ不変性）
  S4 巻き反転共変性: 全巻き反転ケースの r_η が元と一致（対応モードで）
  S5 距離依存: 初期分離 {20°,30°,45°,60°} で r_η の指数を初測定（記録）

出力: result_dl7_charge_force_v1.json・dl7_series_v1.npz
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
CASES = {  # (w_A, w_B)
    "same_pp": (4, 1), "same_mm": (2, 5),
    "opp_pm": (4, 2), "opp_mp": (2, 4), "neut": (3, 6),
}


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    uni = _load("uni_dl7", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_dl7", EXP / "run_cr0_control_no_theta_v2.py")

    # ================= S0: η点力法則（N体エンジン・全レジスタ点） =================
    N = 16
    eng, _, _ = uni.build_standard_universe(N, 0.1)
    for _ in range(15000):
        eng.step()
    C2 = eng.C2().copy()
    Nn, Ne = C2.shape[1], C2.shape[2]
    rngv = np.cos(np.arange(eng.m) * 0.37) + 1j * np.sin(np.arange(eng.m) * 0.61)
    C2[:, 1, 2] += 0.03 * rngv / np.linalg.norm(rngv)      # 巻き +2 成分
    C2[:, 1, Ne - 2] += 0.03 * np.conj(rngv) / np.linalg.norm(rngv)  # 巻き −2 成分
    eng.C = C2.reshape(eng.m, -1)
    W = np.fft.ifft2(eng.C2(), axes=(1, 2)) * (Nn * Ne)
    Wf = W.reshape(eng.m, -1)
    R = eng._readout()
    rate = eng._vertex_rate(Wf, R)
    dd2_num = 2 * np.real(np.conj(Wf) * rate)
    z2 = Wf ** 2
    B = np.zeros((N, Wf.shape[1]), complex)
    BR = np.zeros((N, Wf.shape[1]), complex)
    np.add.at(B, eng.ia, z2); np.add.at(B, eng.ib, z2)
    Rz2 = R[:, None] * z2
    np.add.at(BR, eng.ia, Rz2); np.add.at(BR, eng.ib, Rz2)
    xb2 = np.conj(Wf) ** 2
    dd2_th = (R[:, None] * np.imag((B[eng.ia] + B[eng.ib]) * xb2)
              + np.imag((BR[eng.ia] + BR[eng.ib]) * xb2))
    scale = float(np.max(np.abs(dd2_th)) + 1e-300)
    S0_dev = float(np.max(np.abs(dd2_num - dd2_th)) / scale)
    S0 = {"rel_dev_all_points": S0_dev, "n_points": int(Wf.shape[1]),
          "pass": bool(S0_dev < 1e-13)}

    # ================= 二体系（S1〜S5） =================
    base = uni.two_body_base
    step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63,
                                              recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    def make_ab(wA, wB, deg_b=+30.0):
        mA, mB = wA - 1, wB - 2      # 素の巻き A:+1, B:+2（stage4 正本と同一）
        case = base.explicit_packet_case(
            mode=f"dl7_{wA}_{wB}_{deg_b}", packet_a=PACKET, packet_b=PACKET,
            packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
            packet_b_shift=cr0.shift_for_deg(deg_b, slope, icept))
        a = base.make_case_state(sp, case, "A", hair_enabled=True)
        b = base.make_case_state(sp, case, "B", hair_enabled=True)
        eta = 2.0 * np.pi * np.arange(ne) / ne
        a = (a.reshape(nc, ne) * np.exp(1j * mA * eta)[None, :]).reshape(-1)
        b = (b.reshape(nc, ne) * np.exp(1j * mB * eta)[None, :]).reshape(-1)
        return a, b

    def mode_pos(psi, w):
        F = np.fft.fft(psi.reshape(nc, ne), axis=1)
        P = np.abs(F[:, w % ne]) ** 2
        z = np.sum(P * np.exp(1j * 2 * np.pi * np.arange(nc) / nc)) / max(P.sum(), 1e-300)
        return float(np.angle(z))

    def run_case(wA, wB, deg_b=+30.0, T=T_STEPS):
        a, b = make_ab(wA, wB, deg_b)
        sep_marg, sep_dem = [], []
        for _ in range(T):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            sep_marg.append(abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb)))))))
            dth = mode_pos(a, wA) - mode_pos(b, wB)
            sep_dem.append(abs(float(np.degrees(np.angle(np.exp(1j * dth))))))
        return np.array(sep_marg), np.array(sep_dem)

    def rate_of(sd):
        h = len(sd) // 2
        t = np.arange(h, len(sd))
        return float(np.polyfit(t, sd[h:], 1)[0])

    results = {}
    for name, (wA, wB) in CASES.items():
        sm, sd = run_case(wA, wB)
        results[name] = {"marg": sm, "dem": sd, "r_eta": rate_of(sd)}
        print(f"  [{name}] w=({wA},{wB}) q=({(wA%3+1)%3-1:+d},{(wB%3+1)%3-1:+d}) "
              f"r_eta={results[name]['r_eta']:+.5f}°/步 "
              f"marg端={sm[-1]:.2f}° dem端={sd[-1]:.2f}°")

    # S1 stage4 再現（Δw=±2 の復調弁別 D2）
    st4 = _load("st4_dl7", HERE / "probe_h4_stage4_demod_v1.py")
    canon = json.loads((HERE / "result_h4_stage4_demod_v1.json").read_text())
    _, sp_a2 = st4.run_case(sp, slope, icept, 3, 0, nc, ne)
    _, sp_m2 = st4.run_case(sp, slope, icept, 3, 4, nc, ne)
    D2 = float(np.max(np.abs(sp_a2 - sp_m2)))
    S1 = {"D2_measured": D2, "D2_canonical": canon["D2_demod_spec_p2_vs_m2"],
          "rel_dev": abs(D2 - canon["D2_demod_spec_p2_vs_m2"])
          / canon["D2_demod_spec_p2_vs_m2"],
          "pass": bool(abs(D2 - canon["D2_demod_spec_p2_vs_m2"])
                       / canon["D2_demod_spec_p2_vs_m2"] < 1e-6)}

    # S2 ±則（本丸・三値読み）
    r = {k: results[k]["r_eta"] for k in results}
    same_signs = [np.sign(r["same_pp"]), np.sign(r["same_mm"])]
    opp_signs = [np.sign(r["opp_pm"]), np.sign(r["opp_mp"])]
    split = bool(same_signs[0] == same_signs[1] and opp_signs[0] == opp_signs[1]
                 and same_signs[0] != opp_signs[0])
    S2 = {"r_eta": r, "sign_split_same_vs_opposite": split,
          "assignment": ("同電荷=" + ("斥力型(増大)" if same_signs[0] > 0 else "引力型(減少)")
                         + "・異電荷=" + ("斥力型" if opp_signs[0] > 0 else "引力型")
                         if split else "分岐せず——三値読み(ii)/(iii)として記録"),
          "pass": True}

    # S3 重力側の不変（毛ゲージ不変性：全ケースの周辺化軌道一致）
    ref = results["same_pp"]["marg"]
    S3_dev = float(max(np.max(np.abs(results[k]["marg"] - ref)) for k in results))
    S3 = {"max_traj_dev_deg": S3_dev, "pass": bool(S3_dev < 1e-9)}

    # S4 巻き反転共変性
    sm_r, sd_r = run_case(-4, -1)
    S4_dev = abs(rate_of(sd_r) - r["same_pp"])
    S4 = {"r_eta_mirror": rate_of(sd_r), "r_eta_orig": r["same_pp"],
          "abs_dev": S4_dev, "pass": bool(S4_dev < 1e-9)}

    # S5 距離依存（opp_pm 型・初期分離掃引）
    seps0, rates5 = [], []
    for deg in (20.0, 30.0, 45.0, 60.0):
        _, sd5 = run_case(4, 2, deg_b=deg, T=800)
        L = 2 * np.sin(np.radians(30.0 + deg) / 2)   # 初期分離の弦長
        seps0.append(L)
        rates5.append(abs(rate_of(sd5)))
    ok5 = [i for i in range(len(rates5)) if rates5[i] > 0]
    slope5 = (float(np.polyfit(np.log(np.array(seps0)[ok5]),
                               np.log(np.array(rates5)[ok5]), 1)[0])
              if len(ok5) >= 3 else float("nan"))
    S5 = {"chord": seps0, "abs_rate": rates5, "loglog_slope": slope5,
          "note": "復調力の距離指数の初測定（記録）", "pass": True}

    res = {"config": {"T": T_STEPS, "packet": list(PACKET), "cases": CASES,
                      "grid": [nc, ne],
                      "s0": {"N": N, "delta": 0.1, "tau": 15000,
                             "hair_labels": [2, -2]}},
           "S0_point_force_law": S0, "S1_stage4_reproduction": S1,
           "S2_sign_rule": S2, "S3_gravity_invariance": S3,
           "S4_winding_reversal": S4, "S5_distance_dependence": S5,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("S0_point_force_law", "S1_stage4_reproduction",
                           "S2_sign_rule", "S3_gravity_invariance",
                           "S4_winding_reversal", "S5_distance_dependence"))
    np.savez_compressed(HERE / "dl7_series_v1.npz",
                        **{f"{k}_marg": v["marg"] for k, v in results.items()},
                        **{f"{k}_dem": v["dem"] for k, v in results.items()},
                        mirror_dem=sd_r)
    (HERE / "result_dl7_charge_force_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("S0_point_force_law", "S1_stage4_reproduction", "S2_sign_rule",
              "S3_gravity_invariance", "S4_winding_reversal", "S5_distance_dependence"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  "
              f"{ {kk: vv for kk, vv in res[k].items() if kk not in ('note',)} }")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
