#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 解析 v1 — τ=15581 イベントの切り分け（第5回査読対応・保存済み系列のみ使用）

probe_dl3_event18_wide_v1 の基底不変量 Δφ で唯一 |Δφ|≈0 となった τ=15581 を、
保存済み系列（dl3_branch_lineage_v1.npz: 毎步 λ・ex34・χ・マージン）で切り分ける。

検査（実行前固定）:
 (a) その近傍の χ=-1 窓の有無と個数（このイベントは χ イベントか、系譜のみの検出か）
 (b) 近傍の隣接ギャップ最小値（λ1−λ2, λ2−λ3, λ3−λ4, λ4−λ5）——2×2 縮約が壊れる
     多枝イベント（第2枝または第5枝の同時接近）かどうか
 (c) 隣接步の割当マージン mu_adj の近傍最小——系譜割当の信頼度

出力: result_dl3_event15581_v1.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TAU_EV = 15581
HALF = 60


def main():
    D = np.load(HERE / "dl3_branch_lineage_v1.npz")
    tau, lam, ex34, chi = D["tau"], D["lam"], D["ex34"], D["chi"]
    mu_adj = D["mu_adj"]
    W = 10
    i_ev = int(np.where(tau == TAU_EV)[0][0])
    sl = slice(max(0, i_ev - HALF), i_ev + HALF)

    # (a) χ=-1 窓（窓 k は [tau[k], tau[k+10]] を覆う）
    k_lo = max(0, i_ev - HALF - W)
    k_hi = min(len(chi), i_ev + HALF)
    chi_local = chi[k_lo:k_hi]
    n_chi_windows = int(np.sum(chi_local < 0))
    covering = chi[max(0, i_ev - W):min(len(chi), i_ev + 1)]
    n_chi_covering = int(np.sum(covering < 0))

    # (b) 近傍の隣接ギャップ
    gaps = {}
    names = ["l1-l2", "l2-l3", "l3-l4", "l4-l5"]
    for j, nm in enumerate(names):
        g = lam[sl, j] - lam[sl, j + 1]
        gaps[nm] = {"min": float(g.min()),
                    "argmin_tau": int(tau[sl][int(np.argmin(g))])}

    # (c) 割当マージン
    mu_local = mu_adj[max(0, i_ev - HALF):i_ev + HALF]
    ex_local = np.where(ex34[max(0, i_ev - HALF):i_ev + HALF])[0]

    res = {
        "config": {"tau_event": TAU_EV, "half_window": HALF,
                   "source": "dl3_branch_lineage_v1.npz（決定的系列の事後解析）"},
        "a_chi_windows_nearby": n_chi_windows,
        "a_chi_windows_covering_event": n_chi_covering,
        "b_adjacent_gap_minima": gaps,
        "c_mu_adj_min_nearby": float(mu_local.min()),
        "c_n_ex34_nearby": int(len(ex_local)),
    }
    (HERE / "result_dl3_event15581_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"τ={TAU_EV}±{HALF}:")
    print(f" (a) χ=-1 窓: 近傍{n_chi_windows}本 / イベントを覆う窓{n_chi_covering}本")
    for nm in names:
        print(f" (b) {nm}: min={gaps[nm]['min']:.3e} @τ={gaps[nm]['argmin_tau']}")
    print(f" (c) 割当マージン近傍min={res['c_mu_adj_min_nearby']:.3e}  "
          f"近傍の3↔4交換数={res['c_n_ex34_nearby']}")


if __name__ == "__main__":
    main()
