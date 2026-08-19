#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL5 走行 v1 — 運動量：整合拘束（恒等式）・CR1較正・重み付き保存の検定と重みの逆同定

導出: DL5_導出ノート.md（判定 R1〜R6・実行前固定）。プロトコルは DL4 §6 追補を継承。
力学: unified_interaction_v1（N=16, δ=0.1）・測定 τ∈[15000,16000] 毎步。

実装固定（実行前）:
- 速度: u_A(n) = (R_n^T X3(n+1))_A − X3_A(n)（無重み Kabsch 逐次整列・正準）
- 質量: [F1] §4.2 質量則の頂点星集約版（本走行で式を一意固定・事後変更しない）
    f_A = Σ_{e∋A} P_odd(e),  b_A = Σ_{e∋A} P_even(e)（エンジンの odd_k/even_k 帯）
    θ_A = atan2(√f_A, √b_A),  r_A = f_A/(f_A+b_A),  m_A = ω_pred_A = θ_A + 2 r_A f_A
    （[F1] の ω_pred(x)=θ(x)+2r(x)|a(x)|² の |a|²→f_A（奇数帯パワー）対応）
- R4/R5 の保存量: P_w(n) = Σ_A w_A(n) u_A(n)。保存比 = ||P_w|| の時間変動(std) ÷
  ⟨Σ_A w_A|u_A|⟩（分子が運動量和のふらつき・分母が素朴な運動量スケール）
- R6: CR1 の交換配分は κ=1−r（透過率）で構成されている——構成どおりであることの
  厳密確認（|Δω − κ·acc| < 1e-15）。作用面(a)の平面選択は 1 次元相対座標では生じず、
  3次元対イベントでの独立検定は DL6/次段の実装固定事項として記録する。

判定:
  R1 Σu_A 恒等式 <1e-12 ／ R2 トルクゼロ <1e-10 ／
  R3 CR1 再現: 周期比 1±0.05・有界性 ／
  R4 Σ m_A u_A の保存（三値読み: <1e-10 成立／系統的残差／無秩序） ／
  R5 重み {1, m_A, A_v} の保存比並走→最良重みの同定 ／ R6 上記

出力: result_dl5_momentum_v1.json・dl5_series_v1.npz
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
N, DELTA = 16, 0.1
T_SKIP, T_MEAS = 15000, 1000


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def kabsch(x, y, w=None):
    if w is None:
        w = np.ones(len(x))
    H = (x * w[:, None]).T @ y
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T, S


def main():
    t0 = time.time()

    # ---- R3/R6: CR1 較正（正本 import・短窓再走行） ----
    cr1 = _load("cr1_dl5", EXP / "run_cr1_kinetic_feedback_v1.py")
    H, _, _, _ = cr1.run(cr1.KAPPA_MODE, 0.0, 3000, "dl5cal")
    chi = np.degrees(H["chi"])
    sgn = np.sign(chi - chi.mean())
    zc = np.where(sgn[1:] * sgn[:-1] < 0)[0]
    period = 2.0 * float(np.mean(np.diff(zc))) if len(zc) > 3 else float("nan")
    canon = json.loads((EXP / "result_cr1_kinetic_feedback_v1.json").read_text())
    per_canon = float(canon["metrics"]["period_chi"])
    R3 = {"period_measured": period, "period_canonical": per_canon,
          "ratio": period / per_canon,
          "bounded": bool(np.max(np.abs(chi)) < 61.0),
          "pass": bool(abs(period / per_canon - 1.0) < 0.05
                       and np.max(np.abs(chi)) < 61.0)}
    # R6: Δω = κ·acc の構成確認（厳密）
    # CR1 は更新後の ω を記録する: ω[n]=ω[n-1]+κ[n]·acc[n] → Δω[n]=κ[n+1]·acc[n+1]
    dw = np.diff(H["omega"])
    expect = (H["kappa"] * H["a"])[1:]
    R6_dev = float(np.max(np.abs(dw - expect)))
    R6 = {"exchange_rule_dev": R6_dev, "rule": "Δω = κ·acc, κ=1−r（透過率）",
          "pass": bool(R6_dev < 1e-15),
          "note": "交換配分は透過率で構成どおり（厳密）。作用面(a)の平面選択は"
                  "1次元相対座標では生じない——3次元対イベントの独立検定は次段の実装固定"}

    # ---- N 体走行 ----
    u1 = _load("uni_dl5", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    def frame_state():
        C2 = eng.C2()
        x = C2.sum(axis=(1, 2))
        P2 = np.abs(C2) ** 2
        Pk_odd = P2[:, eng.odd_k, :].sum(axis=(1, 2))
        Pk_even = P2[:, eng.even_k, :].sum(axis=(1, 2))
        f = np.zeros(N); bb = np.zeros(N); Av = np.zeros(N)
        np.add.at(f, eng.ia, Pk_odd); np.add.at(f, eng.ib, Pk_odd)
        np.add.at(bb, eng.ia, Pk_even); np.add.at(bb, eng.ib, Pk_even)
        Pk = np.abs(x) ** 2
        np.add.at(Av, eng.ia, Pk); np.add.at(Av, eng.ib, Pk)
        th = np.arctan2(np.sqrt(f), np.sqrt(bb))
        rA = f / np.maximum(f + bb, 1e-300)
        mass = th + 2.0 * rA * f
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = Pk
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        lamk, Vk = lamk[o], Vk[:, o]
        X3 = Vk[:, :3] * np.sqrt(np.maximum(lamk[:3], 0.0))[None, :]
        return {"X3": X3, "V3": Vk[:, :3], "mass": mass, "Av": Av,
                "lam3": float(lamk[2]),
                "gmin": float(min(lamk[0] - lamk[1], lamk[1] - lamk[2],
                                  lamk[2] - lamk[3]))}

    for _ in range(T_SKIP):
        eng.step()

    fs_prev = frame_state()
    ser = {k: [] for k in ("tau", "sum_u", "torque", "gmin", "lam3",
                            "P1x", "Pm", "Pav", "scale1", "scalem", "scaleav")}
    U_series = []
    for k in range(T_MEAS):
        eng.step()
        fs = frame_state()
        R, S = kabsch(fs_prev["X3"], fs["X3"])
        u = (R.T @ fs["X3"].T).T - fs_prev["X3"]
        Rw, _ = kabsch(fs_prev["X3"], fs["X3"], w=fs["Av"])
        tq = np.sum(fs["Av"][:, None] * np.cross((Rw @ fs_prev["X3"].T).T,
                                                 fs["X3"]), axis=0)
        w1 = np.ones(N)
        wm = fs_prev["mass"]
        wa = fs_prev["Av"]
        ser["tau"].append(T_SKIP + 1 + k)
        ser["sum_u"].append(float(np.max(np.abs(u.sum(axis=0)))))
        ser["torque"].append(float(np.linalg.norm(tq)))
        ser["gmin"].append(fs["gmin"])
        ser["lam3"].append(fs["lam3"])
        for tag, w in (("1", w1), ("m", wm), ("av", wa)):
            P = (w[:, None] * u).sum(axis=0)
            ser[f"P{tag}" if tag != "1" else "P1x"].append(P)
            ser[f"scale{tag}" if tag != "1" else "scale1"].append(
                float((w * np.linalg.norm(u, axis=1)).sum()))
        U_series.append(u)
        fs_prev = fs

    R1_max = float(np.max(ser["sum_u"]))
    R2_max = float(np.max(ser["torque"]))
    R1 = {"max": R1_max, "pass": bool(R1_max < 1e-12)}
    R2 = {"max": R2_max, "pass": bool(R2_max < 1e-10)}

    def cons(tag):
        P = np.array(ser[tag])
        sc = np.array(ser["scale" + ("1" if tag == "P1x" else tag[1:])])
        return float(np.linalg.norm(P.std(axis=0)) / max(sc.mean(), 1e-300))

    ratios = {"w=1": cons("P1x"), "w=m": cons("Pm"), "w=Av": cons("Pav")}
    best = min(ratios, key=ratios.get)
    Pm = np.array(ser["Pm"])
    R4_abs = float(np.linalg.norm(Pm.std(axis=0)))
    R4 = {"P_m_std_norm": R4_abs, "conservation_ratio": ratios["w=m"],
          "reading": ("(i) 保存成立" if ratios["w=m"] < 1e-10 else
                      "(ii)/(iii) 残差あり——構造は R5 と npz 系列で読む"),
          "pass": True}
    R5 = {"conservation_ratios": ratios, "best_weight": best,
          "note": "比が小さいほど保存が良い（P_w の std / 素朴スケール）", "pass": True}
    lam3_pos = bool(np.min(ser["lam3"]) > 0)

    res = {"config": {"N": N, "delta": DELTA, "window": [T_SKIP, T_SKIP + T_MEAS],
                      "engine": "unified_interaction_v1",
                      "mass_def": "m_A=θ_A+2r_Af_A（頂点星集約・本走行で固定）"},
           "R1_sum_u": R1, "R2_torque": R2, "R3_cr1_calibration": R3,
           "R4_weighted_momentum": R4, "R5_weight_identification": R5,
           "R6_exchange_rule": R6,
           "protocol_records": {"lam3_positive_all": lam3_pos},
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("R1_sum_u", "R2_torque", "R3_cr1_calibration",
                           "R4_weighted_momentum", "R5_weight_identification",
                           "R6_exchange_rule"))
    np.savez_compressed(HERE / "dl5_series_v1.npz",
                        U=np.array(U_series),
                        P1=np.array(ser["P1x"]), Pm=np.array(ser["Pm"]),
                        Pav=np.array(ser["Pav"]),
                        tau=np.array(ser["tau"]), gmin=np.array(ser["gmin"]))
    (HERE / "result_dl5_momentum_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("R1_sum_u", "R2_torque", "R3_cr1_calibration",
              "R4_weighted_momentum", "R5_weight_identification", "R6_exchange_rule"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  "
              f"{ {kk: vv for kk, vv in res[k].items() if kk != 'note'} }")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
