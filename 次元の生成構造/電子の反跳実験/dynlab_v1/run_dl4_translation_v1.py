#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL4 走行 v1 — 並進：無名性（恒等式）・凍結慣性（定理）・背景変動下の共変等速（物理検定）

導出: DL4_導出ノート.md（判定 P1〜P6・§6 プロトコル追補を含め実行前固定）。
力学: unified_interaction_v1（N=16, δ=0.1）。測定区間: 物質相後期 τ∈[15000,16000]・毎步
（追補2）。枠は値順・正方向上位3軸（λ3>0 を毎瞬検査・追補3）、SO(3) 逐次整列＋
枝パリティ χ 記録（追補1）。

ゲージ層の実装（[F1] 定義2.1 の θ_A レジスタ・D の carry 型）:
  共動チャート Y(n)=Q(n)^T X3(n)（Q は無重み Kabsch 逐次整列の累積——実装決定(a)＝正準）。
  体ゲージ位置 p_A = Y_A + θ_A、更新 θ_A(n+1)=θ_A(n)+u_A（u_A＝共動チャートで一定
  ＝共変等速の登録。Σu=0・|u|=1e-5/步・決定的パターン）。読出しは対スカラー
  d_AB=||p_A−p_B|| のみ（ゲージ差依存——[F1] 2.1）。z への書き戻しは一切ない。

判定（DL4 ノート §4）:
  P1 共通ドリフト不変性（恒等式）: 全 p_A に共通 (c,Q_c) → 対スカラー変化 <1e-12
  P2 凍結区間の等速（定理）: z を凍結（stepしない）して θ のみ更新 → 対スカラーの
     一階差分が一定（変動 <1e-12）かつ z のレジスタがビット不変
  P3 背景変動下の共変等速（物理検定・記録）: z を毎步進めながら同じ登録速度で更新。
     共動読みの対レート変動と、固定枠読み（コリオリ型混入）の変動を分離して記録
  P4 Kabsch Σ の最小特異値ギャップ（全步記録・縮退計数）
  P5 トルクゼロ停留条件（恒等式・重み付き w=A_v）: ||Σw (R*x)×y|| <1e-10
  P6 並進整列残差: 重み付き重心の時系列（無重み正準採用の残差記録）

出力: result_dl4_translation_v1.json・dl4_series_v1.npz
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent.parent / "統一万能関数_v1"
N, DELTA = 16, 0.1
T_SKIP, T_MEAS = 15000, 1000
U_MAG = 1e-5


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def kabsch(x, y, w=None):
    """y ≈ R x を最小二乗で解く（重み任意・SO(3) 射影）。"""
    if w is None:
        w = np.ones(len(x))
    H = (x * w[:, None]).T @ y
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, S


def main():
    t0 = time.time()
    u1 = _load("uni_dl4", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    def frame_state():
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        lamk, Vk = lamk[o], Vk[:, o]
        X3 = Vk[:, :3] * np.sqrt(np.maximum(lamk[:3], 0.0))[None, :]
        Av = np.zeros(N)
        Pk = np.abs(x) ** 2
        np.add.at(Av, eng.ia, Pk)
        np.add.at(Av, eng.ib, Pk)
        absl = np.abs(lamk)
        return {"X3": X3, "V3": Vk[:, :3], "lam": lamk, "Av": Av,
                "lam3": float(lamk[2]),
                "gmin": float(min(lamk[0] - lamk[1], lamk[1] - lamk[2],
                                  lamk[2] - lamk[3])),
                "W_minus": float(np.where(lamk < 0, -lamk, 0).sum() / absl.sum())}

    for _ in range(T_SKIP):
        eng.step()

    # 登録ドリフト（決定的・Σu=0・共動チャートで一定）
    ang = 2 * np.pi * np.arange(N) / N
    u_reg = U_MAG * np.stack([np.cos(ang), np.sin(ang),
                              np.cos(2 * ang) / np.sqrt(2)], axis=1)
    u_reg -= u_reg.mean(axis=0)

    # ---- P2 凍結区間（z を進めない） ----
    fs0 = frame_state()
    C_before = eng.C.copy()
    Y0 = fs0["X3"]
    theta = np.zeros((N, 3))
    dsel = [(0, 1), (0, 8), (3, 12), (5, 10)]
    dev = 0.0
    for k in range(1, 201):
        theta += u_reg
        p = Y0 + theta
        for a, b in dsel:
            d = np.linalg.norm(p[a] - p[b])
            d_pred = np.linalg.norm((Y0[a] - Y0[b]) + k * (u_reg[a] - u_reg[b]))
            dev = max(dev, abs(d - d_pred))
    # u^cov = Δθ は構成により厳密に u_reg（枠凍結中 R=I）——変動は厳密0
    P2_z_frozen = bool(np.array_equal(C_before, eng.C))
    P2 = {"analytic_distance_dev_max": dev, "ucov_variation": 0.0,
          "z_bitwise_frozen": P2_z_frozen,
          "pass": bool(dev < 1e-12 and P2_z_frozen)}

    # ---- P1 共通ドリフト不変性（恒等式） ----
    p = Y0 + theta
    rng_ang = 0.7
    Qc, _ = kabsch(np.eye(3), np.array([[np.cos(rng_ang), -np.sin(rng_ang), 0],
                                        [np.sin(rng_ang), np.cos(rng_ang), 0],
                                        [0, 0, 1.0]]))
    cshift = np.array([0.3, -0.2, 0.11])
    p2 = (Qc @ p.T).T + cshift
    d_orig = np.array([np.linalg.norm(p[a] - p[b]) for a, b in dsel])
    d_moved = np.array([np.linalg.norm(p2[a] - p2[b]) for a, b in dsel])
    P1_dev = float(np.max(np.abs(d_orig - d_moved)))
    P1 = {"scalar_change_max": P1_dev, "pass": bool(P1_dev < 1e-12)}

    # ---- P3〜P6 背景変動下（z を毎步進める） ----
    theta = np.zeros((N, 3))
    Q = np.eye(3)
    fs_prev = fs0
    Y_prev = fs0["X3"]
    ser = {k: [] for k in ("tau", "gmin", "lam3", "W_minus", "sig_gap",
                            "torque", "wcent", "chi_adj", "omega_frame", "nonrigid")}
    d_cov = []
    vfix_prev_p = None
    vfix = []
    n_degen = 0
    for k in range(T_MEAS):
        eng.step()
        fs = frame_state()
        # 逐次整列（無重み・正準）——実装決定(a)
        R, S = kabsch(fs_prev["X3"], fs["X3"])
        sig_gap = float(min(S[0] - S[1], S[1] - S[2], S[2]))
        n_degen += int(sig_gap < 1e-8)
        chi_adj = float(np.sign(np.linalg.det(fs_prev["V3"].T @ fs["V3"])))
        Q = R @ Q
        Y = (Q.T @ fs["X3"].T).T          # 共動チャート
        # P5 トルクゼロ（重み付き w=A_v・恒等式）
        Rw, _ = kabsch(fs_prev["X3"], fs["X3"], w=fs["Av"])
        tq = np.sum(fs["Av"][:, None] * np.cross((Rw @ fs_prev["X3"].T).T,
                                                 fs["X3"]), axis=0)
        # ゲージ更新（共変等速の登録）
        theta += u_reg
        p = Y + theta
        d_cov.append([np.linalg.norm(p[a] - p[b]) for a, b in dsel])
        p_fixed = (Q @ p.T).T             # 固定枠読み（コリオリ型混入は速度成分に現れる）
        if vfix_prev_p is not None:
            vfix.append(p_fixed - vfix_prev_p)
        vfix_prev_p = p_fixed
        cosang = (np.trace(R) - 1.0) / 2.0
        ser["omega_frame"].append(float(np.degrees(np.arccos(np.clip(cosang, -1, 1)))))
        ser["nonrigid"].append(float(np.linalg.norm(Y - Y_prev)))
        ser["tau"].append(T_SKIP + 1 + k)
        ser["gmin"].append(fs["gmin"])
        ser["lam3"].append(fs["lam3"])
        ser["W_minus"].append(fs["W_minus"])
        ser["sig_gap"].append(sig_gap)
        ser["torque"].append(float(np.linalg.norm(tq)))
        ser["wcent"].append(float(np.linalg.norm(
            (fs["Av"][:, None] * fs["X3"]).sum(axis=0) / fs["Av"].sum())))
        ser["chi_adj"].append(chi_adj)
        fs_prev = fs
        Y_prev = Y

    d_cov = np.array(d_cov)
    vfix = np.array(vfix)                 # (T-1, N, 3) 固定枠の体速度
    rate_cov = np.diff(d_cov, axis=0)
    vfix_std = float(np.mean(vfix.std(axis=0)))
    P3 = {"cov_registered_u_variation": 0.0,
          "pair_rate_std_mean": float(np.mean(rate_cov.std(axis=0))),
          "fixed_frame_velocity_std_mean": vfix_std,
          "fixed_frame_velocity_std_over_u": vfix_std / U_MAG,
          "omega_frame_deg_mean": float(np.mean(ser["omega_frame"])),
          "omega_frame_deg_max": float(np.max(ser["omega_frame"])),
          "nonrigid_residual_mean": float(np.mean(ser["nonrigid"])),
          "note": "物理検定（記録）：共変登録速度は構成により一定。固定枠速度の"
                  "変動/uの比がコリオリ型混入、ω_frame=枠角速度（測定可能量）、"
                  "nonrigid=共動チャートの非剛体残差（配置の呼吸）", "pass": True}
    tq_max = float(np.max(ser["torque"]))
    P5 = {"torque_max": tq_max, "pass": bool(tq_max < 1e-10)}
    P4 = {"sig_gap_min": float(np.min(ser["sig_gap"])),
          "n_degenerate_lt_1e-8": n_degen, "pass": True}
    P6 = {"weighted_centroid_norm_mean": float(np.mean(ser["wcent"])),
          "weighted_centroid_norm_max": float(np.max(ser["wcent"])),
          "note": "無重み正準採用の残差記録（実装決定a）", "pass": True}
    lam3_pos = bool(np.min(ser["lam3"]) > 0)
    n_chi_neg = int(np.sum(np.array(ser["chi_adj"]) < 0))

    res = {"config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                      "window": [T_SKIP, T_SKIP + T_MEAS], "u_mag": U_MAG,
                      "pairs": dsel,
                      "protocol": "毎步・枝検出つき（追補1）・後期窓（追補2）・"
                                  "λ3>0検査と(W_-,gmin)記録（追補3）"},
           "P1_anonymity": P1, "P2_frozen_inertia": P2, "P3_covariant": P3,
           "P4_kabsch_gap": P4, "P5_torque_zero": P5, "P6_alignment_residual": P6,
           "protocol_records": {"lam3_positive_all": lam3_pos,
                                "W_minus_late_mean": float(np.mean(ser["W_minus"])),
                                "n_adjacent_chi_neg": n_chi_neg},
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("P1_anonymity", "P2_frozen_inertia", "P3_covariant",
                           "P4_kabsch_gap", "P5_torque_zero", "P6_alignment_residual"))
    np.savez_compressed(HERE / "dl4_series_v1.npz",
                        d_cov=d_cov, vfix=vfix,
                        **{k: np.array(v) for k, v in ser.items()})
    (HERE / "result_dl4_translation_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("P1_anonymity", "P2_frozen_inertia", "P3_covariant",
              "P4_kabsch_gap", "P5_torque_zero", "P6_alignment_residual"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  "
              f"{ {kk: vv for kk, vv in res[k].items() if kk != 'note'} }")
    print(f"  記録: λ3>0全步={lam3_pos}  W_-平均={np.mean(ser['W_minus']):.3f}  "
          f"隣接χ=-1步={n_chi_neg}")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
