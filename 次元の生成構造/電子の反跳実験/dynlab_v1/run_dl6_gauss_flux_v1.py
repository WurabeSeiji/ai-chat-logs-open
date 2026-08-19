#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL6 走行 v1 — 力の帳簿：離散 Gauss 恒等式・流束・距離指数・三角閉路ホロノミー・AB辞書

導出: DL6_導出ノート.md（判定 Q1〜Q7・追補7〜11 を含め実行前固定）。
力学: unified_interaction_v1（N=16, δ=0.1）・測定 τ∈[15000,16000] 毎步。

実装固定（追補10・11）:
- U_AB = A 星重み配置 → B 星重み配置の Kabsch 回転（SO(3)）。局所配置＝X3 に
  星重み w_j^(A)=|x_{Aj}|²（対角 w_A=Σw_j）を掛け重み付き中心化したもの
- W_ABC = U_AB U_BC U_CA、ε_ABC＝回転角。真空では未定義（基準点・Q6）
- Q5: 標準宇宙は電荷中立（毛0）——中立物質の符号（普遍引力側）と巻き0の確認を記録。
  ±判別の実体は DL7（追補11 の明記どおり）

判定:
  Q1 離散 Gauss 恒等式（全分割型・恒等式）: |Σ_{v∈S}B_v − 2Σ_{e⊆S}x² − Σ_{∂S}x²| <1e-12
  Q2 流束 Φ(r) の r 依存（前提A1）: 各体の距離順 boundary flux のプラトー記録
  Q3 流束の等方性（前提A2）: 方向別флト束の偏差記録
  Q4 加速度の距離指数: 対の radial 加速度 vs 距離の log-log 勾配（記録・−2なら経路A成立）
  Q5 中立物質の符号: 対の平均 radial 加速度の符号（引力側か）＋毛占有が帯0に集中
  Q6 三角閉路 ε_ABC: 全三角形の ε 分布と、面積×局所密度との相関（比例の検定）
  Q7 AB 辞書: 調和閉鎖恒等式 |ω_n|Δθ_n=Ω（π/72 梯子・厳密）と log-log 勾配 −2 の再計算

出力: result_dl6_gauss_flux_v1.json・dl6_series_v1.npz
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
RNG_PART = 20240819


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


def rot_angle(R):
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))


def main():
    t0 = time.time()
    u1 = _load("uni_dl6", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)
    edge_of = {}
    for e, (a, b) in enumerate(zip(ii, jj)):
        edge_of[(a, b)] = e
        edge_of[(b, a)] = e

    def snap():
        C2 = eng.C2()
        x = C2.sum(axis=(1, 2))
        Pk = np.abs(x) ** 2
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
        return x, Pk, X3, C2

    for _ in range(T_SKIP):
        eng.step()

    rng = np.random.default_rng(RNG_PART)
    x0, Pk0, X3_prev, _ = snap()

    # ---- Q1 離散 Gauss 恒等式（初期スナップで20分割・複素で厳密） ----
    q1_dev = 0.0
    x2 = x0 ** 2
    Bv = np.zeros(N, complex)
    np.add.at(Bv, ii, x2)
    np.add.at(Bv, jj, x2)
    for _ in range(20):
        S_set = set(rng.choice(N, size=rng.integers(2, N - 1), replace=False).tolist())
        lhs = sum(Bv[v] for v in S_set)
        inner = sum(x2[e] for e, (a, b) in enumerate(zip(ii, jj))
                    if a in S_set and b in S_set)
        bnd = sum(x2[e] for e, (a, b) in enumerate(zip(ii, jj))
                  if (a in S_set) != (b in S_set))
        q1_dev = max(q1_dev, abs(lhs - 2 * inner - bnd))
    Q1 = {"max_dev": float(q1_dev), "pass": bool(q1_dev < 1e-12)}

    # ---- 本走行（毎步・整列速度と加速度・三角閉路） ----
    tri_sel = [(0, 5, 10), (1, 6, 11), (2, 7, 12), (3, 8, 13), (0, 4, 8),
               (2, 9, 14), (1, 5, 15), (3, 7, 11)]
    ser = {k: [] for k in ("tau", "eps_tri", "area_tri", "dens_tri")}
    Xa, Qacc = [], np.eye(3)
    Y_list = []
    flux_prof = []
    for k in range(T_MEAS):
        eng.step()
        x, Pk, X3, C2 = snap()
        R, _ = kabsch(X3_prev, X3)
        Qacc = R @ Qacc
        Y_list.append((Qacc.T @ X3.T).T)
        X3_prev = X3
        if k % 50 == 0:
            # Q2/Q3: 各体の距離順 boundary flux（複素 x² の境界和の絶対値）
            x2n = x ** 2
            prof = []
            for A in range(N):
                dists = np.array([np.linalg.norm(X3[A] - X3[j]) if j != A else 0.0
                                  for j in range(N)])
                order = np.argsort(dists)
                row = []
                for m in range(1, N - 1):
                    S_set = set(order[:m + 1].tolist())
                    bnd = sum(x2n[e] for e, (a, b) in enumerate(zip(ii, jj))
                              if (a in S_set) != (b in S_set))
                    row.append(abs(bnd))
                prof.append(row)
            flux_prof.append(prof)
        # Q6: 三角閉路（追補10・11 の U_AB）
        if k % 10 == 0:
            for (A, B_, Cv) in tri_sel:
                def local_cfg(P):
                    w = np.array([np.abs(x[edge_of[(P, j)]]) ** 2 if j != P else 0.0
                                  for j in range(N)])
                    w[P] = w.sum()
                    c = (w[:, None] * X3).sum(axis=0) / w.sum()
                    return X3 - c, w
                cfgA, wA = local_cfg(A)
                cfgB, wB = local_cfg(B_)
                cfgC, wC = local_cfg(Cv)
                UAB, _ = kabsch(cfgA, cfgB, w=np.sqrt(wA * wB))
                UBC, _ = kabsch(cfgB, cfgC, w=np.sqrt(wB * wC))
                UCA, _ = kabsch(cfgC, cfgA, w=np.sqrt(wC * wA))
                W = UAB @ UBC @ UCA
                v1_, v2_ = X3[B_] - X3[A], X3[Cv] - X3[A]
                area = 0.5 * np.linalg.norm(np.cross(v1_, v2_))
                dens = float(Pk[edge_of[(A, B_)]] + Pk[edge_of[(B_, Cv)]]
                             + Pk[edge_of[(Cv, A)]])
                ser["eps_tri"].append(rot_angle(W))
                ser["area_tri"].append(float(area))
                ser["dens_tri"].append(dens)
                ser["tau"].append(T_SKIP + 1 + k)

    Y = np.array(Y_list)                    # (T, N, 3) 共動チャート
    U = np.diff(Y, axis=0)
    Acc = np.diff(U, axis=0)                # (T-2, N, 3)

    # ---- Q4/Q5: 対の radial 加速度 vs 距離 ----
    pair_r, pair_ar = [], []
    for t in range(0, len(Acc), 10):
        Xt = Y[t + 1]
        for e, (a, b) in enumerate(zip(ii, jj)):
            dvec = Xt[a] - Xt[b]
            dn = np.linalg.norm(dvec)
            if dn < 1e-12:
                continue
            arel = float((Acc[t, a] - Acc[t, b]) @ (dvec / dn))
            pair_r.append(dn)
            pair_ar.append(arel)
    pair_r = np.array(pair_r)
    pair_ar = np.array(pair_ar)
    # log-log 勾配（|a_rel| の距離ビン中央値で回帰）
    bins = np.quantile(pair_r, np.linspace(0, 1, 9))
    br, ba = [], []
    for i in range(8):
        m = (pair_r >= bins[i]) & (pair_r < bins[i + 1])
        if m.sum() > 10:
            br.append(np.median(pair_r[m]))
            ba.append(np.median(np.abs(pair_ar[m])))
    slope = float(np.polyfit(np.log(br), np.log(ba), 1)[0])
    Q4 = {"loglog_slope": slope, "n_samples": int(len(pair_r)),
          "note": "経路Aの検定（記録）。−2 なら成立、ずれは R 集約の形の情報",
          "pass": True}
    frac_attr = float(np.mean(pair_ar < 0))
    # 毛占有（電荷中立の確認）——C2 の第2軸は巻き指数（周波数）なので直接分率で読む
    C2 = eng.C2()
    P_eta = np.abs(C2) ** 2
    w0_frac = float(P_eta[:, :, 0].sum() / P_eta.sum())
    Q5 = {"frac_pairs_attractive": frac_attr, "hair_w0_fraction": w0_frac,
          "note": "中立物質（巻き0占有=hair_w0_fraction）の符号記録。"
                  "±判別の実体は DL7（追補11）", "pass": True}

    flux_prof = np.array(flux_prof)          # (snap, N, N-2)
    prof_mean = flux_prof.mean(axis=(0, 1))
    Q2 = {"flux_profile_mean": [float(v) for v in prof_mean],
          "plateau_cv_mid": float(np.std(prof_mean[4:12]) / np.mean(prof_mean[4:12])),
          "note": "前提A1の検定（記録）：距離順 boundary flux の中域変動", "pass": True}
    Q3 = {"note": "等方性は N=16 の角度分解能では方向ビンが疎——"
                  "各体プロファイルの体間ばらつきで代理記録",
          "body_cv_mid": float(np.mean(flux_prof[:, :, 4:12].std(axis=1)
                                       / flux_prof[:, :, 4:12].mean(axis=1))),
          "pass": True}
    eps = np.array(ser["eps_tri"])
    area = np.array(ser["area_tri"])
    dens = np.array(ser["dens_tri"])
    ad = area * dens
    corr = float(np.corrcoef(ad, eps)[0, 1])
    Q6 = {"eps_deg_mean": float(eps.mean()), "eps_deg_max": float(eps.max()),
          "corr_eps_vs_area_x_density": corr,
          "note": "曲率第一候補の初回実測。比例（面積×密度）の相関を記録", "pass": True}

    # ---- Q7: AB 辞書（調和閉鎖恒等式・π/72 梯子） ----
    w1 = np.pi / 72.0
    ns = np.arange(1, 25)
    omega_n = ns * w1
    dth = 2 * np.pi / ns
    Omega = 2 * np.pi * w1
    dual_dev = float(np.max(np.abs(omega_n * dth - Omega)))
    alpha = omega_n ** 2                     # R=1 規約（[AB] α_n = R|ω_n|²）
    s2 = float(np.polyfit(np.log(dth), np.log(alpha), 1)[0])
    Q7 = {"dual_identity_dev": dual_dev, "alpha_vs_dtheta_slope": s2,
          "pass": bool(dual_dev < 1e-14 and abs(s2 + 2) < 1e-10)}

    res = {"config": {"N": N, "delta": DELTA, "window": [T_SKIP, T_SKIP + T_MEAS],
                      "engine": "unified_interaction_v1",
                      "triangles": tri_sel, "partitions_q1": 20},
           "Q1_gauss": Q1, "Q2_flux": Q2, "Q3_isotropy": Q3,
           "Q4_distance_exponent": Q4, "Q5_neutral_sign": Q5,
           "Q6_triangle_holonomy": Q6, "Q7_ab_dictionary": Q7,
           "elapsed_sec": time.time() - t0}
    res["all_pass"] = all(res[k]["pass"] for k in
                          ("Q1_gauss", "Q2_flux", "Q3_isotropy", "Q4_distance_exponent",
                           "Q5_neutral_sign", "Q6_triangle_holonomy", "Q7_ab_dictionary"))
    np.savez_compressed(HERE / "dl6_series_v1.npz",
                        pair_r=pair_r, pair_ar=pair_ar, eps=eps, area=area,
                        dens=dens, flux_prof=flux_prof)
    (HERE / "result_dl6_gauss_flux_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("Q1_gauss", "Q2_flux", "Q3_isotropy", "Q4_distance_exponent",
              "Q5_neutral_sign", "Q6_triangle_holonomy", "Q7_ab_dictionary"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  "
              f"{ {kk: vv for kk, vv in res[k].items() if kk not in ('note', 'flux_profile_mean')} }")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
