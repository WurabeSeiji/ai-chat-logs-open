#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 解析 v1 — 局在の選択性・W_- 先行・η_{3+}（第3回査読対応・保存済み系列のみ使用）

データ源: dl3_full_scan_v1.npz（probe_dl3_full_scan_v1.py の毎步 λ 系列・χ・窓内最小ギャップ）。
再走行はしない（決定的系列の事後解析。判定基準は本スクリプト作成時に固定）。

(1) 局在の選択性:
    P(Δ34=Δmin | χ=-1) vs P(Δ34=Δmin | χ=+1) —— 113/113 が強い選択か基底率かを決着。
    P(χ=-1 | Δ34窓内最小 ∈ 五分位) —— 枝イベント率がギャップ閉鎖へ向けて増えるか。
(2) η_{3+} = (λ1+λ2+λ3)/Σ_{λ>0}λ: 正方向内部の3軸支配（W_- と分離した二層量）。
(3) W_- の先行性: τ_W（W_- が真空基線の10倍を50步連続で超える初出）と
    τ_-, τ_g の順序。dW_-/dn の立ち上がりと g_min 成長の時間差。

出力: result_dl3_branch_stats_v1.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W = 10


def first_sustained(taus, mask, run=50):
    cnt = 0
    for t, m in zip(taus, mask):
        cnt = cnt + 1 if m else 0
        if cnt >= run:
            return int(t - run + 1)
    return None


def main():
    D = np.load(HERE / "dl3_full_scan_v1.npz")
    tau, lam, chi, W_minus = D["tau"], D["lam"], D["chi"], D["W_minus"]
    gmin = D["gmin"]
    scan0 = int(D["scan_tau0"][0])
    i0 = scan0 - int(tau[0])
    gaps = lam[:, :3] - lam[:, 1:4]
    n_win = len(chi)

    # ---- (1) 局在の選択性 ----
    which_min = np.zeros(n_win, int)
    d34_min = np.zeros(n_win)
    for k in range(n_win):
        g_in = gaps[i0 + k:i0 + k + W + 1]
        which_min[k] = int(np.unravel_index(np.argmin(g_in), g_in.shape)[1])
        d34_min[k] = float(g_in[:, 2].min())
    ev = chi < 0
    p_d34min_ev = float(np.mean(which_min[ev] == 2)) if ev.any() else None
    p_d34min_nev = float(np.mean(which_min[~ev] == 2))
    qs = np.quantile(d34_min, [0.2, 0.4, 0.6, 0.8])
    bins = np.digitize(d34_min, qs)
    p_ev_by_quintile = [float(np.mean(ev[bins == b])) for b in range(5)]
    # イベント統計（連続イベント窓を併合し、イベントごとの最小 Δ34 の五分位）
    ev_int = ev.astype(int)
    starts = np.where(np.diff(np.concatenate([[0], ev_int])) == 1)[0]
    ends = np.where(np.diff(np.concatenate([ev_int, [0]])) == -1)[0]
    ev_quintiles = []
    for s, e in zip(starts, ends):
        ev_quintiles.append(int(np.digitize(d34_min[s:e + 1].min(), qs)))
    ev_quintile_counts = [int(sum(q == b for q in ev_quintiles)) for b in range(5)]

    # ---- (2) η_{3+} ----
    pos = np.where(lam > 0, lam, 0.0)
    eta3p = pos[:, :3].sum(axis=1) / pos.sum(axis=1)
    late = tau >= scan0

    # ---- (3) W_- の先行性 ----
    base = tau <= 8500
    Wb = float(np.median(W_minus[base]))
    gb = float(np.median(gmin[base]))
    tau_W = first_sustained(tau, W_minus > 10 * max(Wb, 1e-15))
    nu = (lam < -1e-12 * lam.sum(axis=1, keepdims=True)).sum(axis=1)
    tau_m = first_sustained(tau, nu > 0)
    tau_g = first_sustained(tau, gmin > 10 * gb)

    res = {
        "config": {"source": "dl3_full_scan_v1.npz（決定的系列の事後解析）",
                   "W": W, "quintile_edges": [float(q) for q in qs]},
        "selectivity": {
            "P_d34_is_min_given_event": p_d34min_ev,
            "P_d34_is_min_given_no_event": p_d34min_nev,
            "P_event_by_d34_quintile_low_to_high": p_ev_by_quintile,
            "n_events_merged": len(ev_quintiles),
            "event_counts_by_quintile_low_to_high": ev_quintile_counts},
        "eta3plus": {"late_mean": float(eta3p[late].mean()),
                     "late_min": float(eta3p[late].min()),
                     "late_max": float(eta3p[late].max())},
        "W_minus_lead": {"W_baseline_median": Wb,
                         "tau_W": tau_W, "tau_minus": tau_m, "tau_g": tau_g,
                         "order_W_before_g": (tau_W is not None and tau_g is not None
                                              and tau_W < tau_g)},
    }
    (HERE / "result_dl3_branch_stats_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"(1) P(Δ34=min | χ=-1)={p_d34min_ev}  P(Δ34=min | χ=+1)={p_d34min_nev:.4f}")
    print(f"    P(χ=-1 | Δ34五分位 低→高)={['%.4f' % p for p in p_ev_by_quintile]}")
    print(f"    イベント統計: {len(ev_quintiles)}件の五分位内訳(低→高)={ev_quintile_counts}")
    print(f"(2) η_3+ 後期: 平均={eta3p[late].mean():.4f} "
          f"[{eta3p[late].min():.4f}, {eta3p[late].max():.4f}]")
    print(f"(3) τ_W={tau_W}  τ_-={tau_m}  τ_g={tau_g}  "
          f"W_-がg_minに先行: {res['W_minus_lead']['order_W_before_g']}")


if __name__ == "__main__":
    main()
