#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 全スペクトル毎步走査（再査読対応・実行前固定）

同一の決定的物質走行（N=16, δ=0.1, F=v1）を τ=8000〜20000 の毎步で走査し、
再査読の4検定を一括で決着させる：

(A) 点火順序のロバスト性: ν_-(ε) の初出 τ_-(ε)（ε∈{1e-12,1e-10,1e-8,1e-6}×tr(B)）と
    ギャップ開通 τ_g（しきい {5,10,20}×g_b × 持続 {20,50,100}步）の全組合せで
    τ_- < τ_g が保存されるか（順序不変性）。
(B) 枝イベントの全窓走査: 後期 15000〜19990 の全スライド窓 [n,n+10] で
    枝パリティ χ(n) = sgn det(Q_n^T Q_{n+10}) · Π_k sgn det(Q_k^T Q_{k+1})
    （固有ベクトル列符号の反転は両因子で相殺——符号規約不変）。
    χ=-1 の時刻集合と、その窓内の隣接ギャップ最小値（どのギャップが閉じたか）を照合。
(C) M5b（Eckart–Young 最適性）: min_{i≤3}|λ_i| ≥ max_{j≥4}|λ_j| の成立割合。
(D) λ3>0（正方向3本の存在）の成立割合と W_-（負固有値の重み比）系列。

出力: result_dl3_full_scan_v1.json・dl3_full_scan_v1.npz（λ全系列・W_-系列）
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
T0, T1 = 8000, 20000
SCAN0 = 15000        # 枝イベント全窓走査の開始
W = 10
EPSILONS = [1e-12, 1e-10, 1e-8, 1e-6]
G_FACTORS = [5, 10, 20]
RUNS = [20, 50, 100]


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def first_sustained(taus, mask, run):
    cnt = 0
    for t, m in zip(taus, mask):
        cnt = cnt + 1 if m else 0
        if cnt >= run:
            return int(t - run + 1)
    return None


def main():
    t0c = time.time()
    u1 = _load("uni_dl3fs", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    for _ in range(T0):
        eng.step()

    n_steps = T1 - T0
    lams = np.zeros((n_steps, N - 1))
    V3s = []
    for k in range(n_steps):
        eng.step()
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        lams[k] = lamk[o]
        V3s.append(Vk[:, o][:, :3].copy())
    taus = np.arange(T0 + 1, T1 + 1)

    trB = lams.sum(axis=1)
    gaps = lams[:, :3] - lams[:, 1:4]          # (λ1-λ2, λ2-λ3, λ3-λ4)
    gmin = gaps.min(axis=1)
    absl = np.abs(lams)
    W_minus = np.where(lams < 0, -lams, 0.0).sum(axis=1) / absl.sum(axis=1)

    # ---- (A) 点火順序のロバスト性 ----
    base = taus <= 8500
    g_b = float(np.median(gmin[base]))
    grid = []
    all_ordered = True
    for eps in EPSILONS:
        nu = (lams < -eps * trB[:, None]).sum(axis=1)
        for run in RUNS:
            tm = first_sustained(taus, nu > 0, run)
            for gf in G_FACTORS:
                tg = first_sustained(taus, gmin > gf * g_b, run)
                ordered = (tm is not None and tg is not None and tm < tg)
                all_ordered &= bool(ordered)
                grid.append({"eps": eps, "run": run, "g_factor": gf,
                             "tau_minus": tm, "tau_g": tg, "ordered": ordered})

    # ---- (B) 枝イベントの全窓走査 ----
    i0 = SCAN0 - (T0 + 1)
    sgn_adj = np.array([np.sign(np.linalg.det(V3s[k].T @ V3s[k + 1]))
                        for k in range(i0, n_steps - 1)])
    events = []
    n_win = n_steps - i0 - W
    chi = np.zeros(n_win)
    for k in range(n_win):
        s_dir = np.sign(np.linalg.det(V3s[i0 + k].T @ V3s[i0 + k + W]))
        chi[k] = s_dir * np.prod(sgn_adj[k:k + W])
        if chi[k] < 0:
            g_in = gaps[i0 + k:i0 + k + W + 1]
            j = int(np.unravel_index(np.argmin(g_in), g_in.shape)[1])
            events.append({"tau_start": int(taus[i0 + k]),
                           "min_gap_in_window": float(g_in.min()),
                           "which_gap": ["l1-l2", "l2-l3", "l3-l4"][j]})
    # イベント窓 vs 全窓の窓内最小ギャップ
    win_min_gap = np.array([gaps[i0 + k:i0 + k + W + 1].min()
                            for k in range(n_win)])
    ev_mask = chi < 0
    # 連続窓は同一イベントに併合して計数
    ev_starts = int(np.sum(ev_mask[1:] & ~ev_mask[:-1]) + (1 if ev_mask[0] else 0))

    # ---- (C)(D) M5b・λ3>0・W_- ----
    m5b_ok = absl[:, :3].min(axis=1) >= absl[:, 3:].max(axis=1)
    lam3_pos = lams[:, 2] > 0
    late = taus >= SCAN0

    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "window": [T0, T1], "scan0": SCAN0, "W": W,
                   "g_baseline_median": g_b,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行"},
        "A_order_robustness": {
            "all_ordered": bool(all_ordered),
            "n_combos": len(grid),
            "n_ordered": int(sum(g["ordered"] for g in grid)),
            "grid": grid},
        "B_branch_scan": {
            "n_windows": int(n_win),
            "n_event_windows": int(ev_mask.sum()),
            "n_events_merged": ev_starts,
            "frac_event_windows": float(ev_mask.mean()),
            "min_gap_event_median": (float(np.median(win_min_gap[ev_mask]))
                                     if ev_mask.any() else None),
            "min_gap_all_median": float(np.median(win_min_gap)),
            "events": events[:50]},
        "C_m5b_optimality": {
            "frac_all": float(m5b_ok.mean()),
            "frac_late": float(m5b_ok[late].mean())},
        "D_lam3_W": {
            "lam3_pos_frac_all": float(lam3_pos.mean()),
            "lam3_pos_frac_late": float(lam3_pos[late].mean()),
            "W_minus_late_mean": float(W_minus[late].mean()),
            "W_minus_max": float(W_minus.max())},
        "elapsed_sec": time.time() - t0c,
    }
    np.savez_compressed(HERE / "dl3_full_scan_v1.npz",
                        tau=taus, lam=lams, W_minus=W_minus, gmin=gmin,
                        chi=chi, win_min_gap=win_min_gap,
                        scan_tau0=np.array([SCAN0]))
    (HERE / "result_dl3_full_scan_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"(A) 順序不変性: {res['A_order_robustness']['n_ordered']}/"
          f"{res['A_order_robustness']['n_combos']} 組合せで τ_-<τ_g  "
          f"all={all_ordered}")
    print(f"(B) 全窓走査: 窓{n_win}本中 イベント窓{int(ev_mask.sum())}"
          f"（併合{ev_starts}イベント・率{ev_mask.mean():.4f}）")
    print(f"    イベント窓の窓内最小ギャップ中央値="
          f"{res['B_branch_scan']['min_gap_event_median']}  "
          f"全窓中央値={res['B_branch_scan']['min_gap_all_median']:.2e}")
    if events:
        from collections import Counter
        print(f"    閉じたギャップの内訳: "
              f"{Counter(e['which_gap'] for e in events)}")
    print(f"(C) M5b 最適性: 全区間 {m5b_ok.mean():.4f} / 後期 {m5b_ok[late].mean():.4f}")
    print(f"(D) λ3>0: 全区間 {lam3_pos.mean():.4f} / 後期 {lam3_pos[late].mean():.4f}  "
          f"W_- 後期平均={W_minus[late].mean():.4f} 最大={W_minus.max():.4f}")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
