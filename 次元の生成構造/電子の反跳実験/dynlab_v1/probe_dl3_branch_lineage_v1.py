#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 固有枝の系譜追跡：交差か回避交差か（第3回査読対応・実行前固定）

「枝イベントは λ3−λ4 の最小ギャップに局在する」（新事実2）に残る一点を決着させる：
値順表示では常に λ3≥λ4 なので、「最小」と「交差」は同じでない。
本プローブは同一の決定的物質走行（N=16, δ=0.1, F=v1）の τ=15000〜20000 を毎步で走査し、

(1) 固有枝の系譜: 隣接步の全基底重なり |V(n)^T V(n+1)| の貪欲最大割当で
    15本の枝にラベルを付けて追跡。順位3の枝と順位4の枝がラベル交換していれば
    「真の交差」、接近して離れるだけなら「回避交差」。
(2) χ=−1 窓（枝パリティ・再計算で自己完結）と 3↔4 交換の対応:
    窓内の 3↔4 交換回数の偶奇と χ の一致率（Z₂ ホロノミーの機構同定）。
(3) 行列式マージン: 隣接步 μ=|det(V3(n)^T V3(n+1))| と直接 |det(V3(n)^T V3(n+10))| の
    イベント窓内最小値。機械零から十分離れていれば χ=−1 は数値符号ノイズでない。

出力: result_dl3_branch_lineage_v1.json
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
T_SKIP, T_END, W = 15000, 20000, 10


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def greedy_assign(O):
    """重なり行列 O (k×k) の貪欲最大割当。perm[i]=行iに対応する列。"""
    k = O.shape[0]
    O = O.copy()
    perm = np.full(k, -1)
    for _ in range(k):
        i, j = np.unravel_index(np.argmax(O), O.shape)
        perm[i] = j
        O[i, :] = -1.0
        O[:, j] = -1.0
    return perm


def main():
    t0 = time.time()
    u1 = _load("uni_dl3bl", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    for _ in range(T_SKIP):
        eng.step()

    n_steps = T_END - T_SKIP
    Vs, lams = [], np.zeros((n_steps, N - 1))
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
        Vs.append(Vk[:, o].copy())
    taus = np.arange(T_SKIP + 1, T_END + 1)

    # ---- (1) 系譜追跡 ----
    K = N - 1
    labels = np.arange(K)          # 順位 r にいる枝ラベル（初期＝順位）
    label_at_rank = np.zeros((n_steps, K), int)
    label_at_rank[0] = labels
    ex34 = np.zeros(n_steps - 1, bool)   # 步 n→n+1 で順位3・4の枝が交換
    for k in range(n_steps - 1):
        O = np.abs(Vs[k].T @ Vs[k + 1])
        perm = greedy_assign(O)          # 旧順位 i → 新順位（列）ではなく列対応:
        # perm[i]=j: 旧順位 i の枝が新步の順位 j に居る
        new_labels = np.zeros(K, int)
        new_labels[perm] = labels
        ex34[k] = (new_labels[2] == labels[3]) and (new_labels[3] == labels[2])
        labels = new_labels
        label_at_rank[k + 1] = labels

    # 全交換の統計（順位3↔4 に限らない隣接交換）
    n_ex34 = int(ex34.sum())

    # ---- (2) χ 再計算と 3↔4 交換の対応 ----
    sgn_adj = np.array([np.sign(np.linalg.det(Vs[k][:, :3].T @ Vs[k + 1][:, :3]))
                        for k in range(n_steps - 1)])
    mu_adj = np.array([abs(np.linalg.det(Vs[k][:, :3].T @ Vs[k + 1][:, :3]))
                       for k in range(n_steps - 1)])
    n_win = n_steps - W
    chi = np.zeros(n_win)
    mu_dir = np.zeros(n_win)
    ex_count_in_win = np.zeros(n_win, int)
    for k in range(n_win):
        d = np.linalg.det(Vs[k][:, :3].T @ Vs[k + W][:, :3])
        mu_dir[k] = abs(d)
        chi[k] = np.sign(d) * np.prod(sgn_adj[k:k + W])
        ex_count_in_win[k] = int(ex34[k:k + W].sum())
    ev = chi < 0
    odd_ex = (ex_count_in_win % 2) == 1
    agree = float(np.mean(odd_ex == ev))
    n_ev = int(ev.sum())
    # イベント窓での交換内訳
    ev_with_ex = int(np.sum(ev & (ex_count_in_win > 0)))
    ev_odd = int(np.sum(ev & odd_ex))

    # ---- (3) マージン ----
    mu_adj_min_event = (float(min(mu_adj[k:k + W].min() for k in np.where(ev)[0]))
                        if n_ev else None)
    mu_dir_min_event = float(mu_dir[ev].min()) if n_ev else None
    mu_adj_min_all = float(mu_adj.min())

    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "window": [T_SKIP, T_END], "W": W,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行。"
                           "系譜割当は隣接步全基底重なりの貪欲最大割当"},
        "lineage": {
            "n_rank34_exchanges": n_ex34,
            "n_event_windows": n_ev,
            "n_event_windows_with_exchange": ev_with_ex,
            "n_event_windows_with_odd_exchange": ev_odd,
            "chi_vs_odd_exchange_agreement": agree,
        },
        "margins": {
            "mu_adjacent_min_all": mu_adj_min_all,
            "mu_adjacent_min_in_events": mu_adj_min_event,
            "mu_direct_min_in_events": mu_dir_min_event,
        },
        "elapsed_sec": time.time() - t0,
    }
    np.savez_compressed(HERE / "dl3_branch_lineage_v1.npz",
                        tau=taus, lam=lams, ex34=ex34, chi=chi,
                        mu_adj=mu_adj, mu_dir=mu_dir,
                        ex_count_in_win=ex_count_in_win)
    (HERE / "result_dl3_branch_lineage_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"3↔4 交換: {n_ex34}回（毎步系譜）")
    print(f"χ=-1 窓: {n_ev}  うち交換あり {ev_with_ex}・奇数回 {ev_odd}")
    print(f"χ ⟺ 交換回数の偶奇 一致率: {agree:.4f}")
    print(f"マージン: 隣接min(全体)={mu_adj_min_all:.3e}  "
          f"隣接min(イベント)={mu_adj_min_event}  直接min(イベント)={mu_dir_min_event}")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
