#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 第18イベント（τ=16985）の決着：広窓＋基底不変量（第5回査読対応）

probe_dl3_avoided_crossing_v1 の ±25步窓で唯一 a−b の符号反転が検出されなかった
イベント τ=16985 を、±150步の広窓で再解析する。

基底依存性の注記（実行前固定）: 平行輸送基底 P の初期回転により (a−b, 2c) は
2次元ベクトルとして回転するため、「a−b の符号反転」自体は基底依存である。
基底に依らない判別量は
  (i) 軌跡ノルム √((a−b)²+(2c)²) = 有効固有値差（最小値＝最接近ギャップ）
  (ii) 軌跡の回転角 Δφ = unwrap(atan2(2c, a−b)) の総変化
      回避交差（枝交換を伴う）⇔ |Δφ| ≈ π（軌跡が原点の片側を回り込む）
      単なる接近（枝交換なし）⇔ |Δφ| ≈ 0（近づいて戻る）
判定: |Δφ| > π/2 なら回避交差と分類（18/18 が回避交差で閉じる）。
     |Δφ| ≤ π/2 なら「有限ギャップの深い接近イベント」として別分類を記録。
全18イベントについて同じ広窓 Δφ を併記する（±25窓の符号反転判定の基底依存性を除去）。

出力: result_dl3_event18_wide_v1.json
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
T_SKIP, T_END = 15000, 20000
HALF = 150
EVENT_TAUS = [None]  # 実行時に系譜追跡で再同定


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def greedy_assign(O):
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
    u1 = _load("uni_dl3e18", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    for _ in range(T_SKIP):
        eng.step()
    n_steps = T_END - T_SKIP
    Vs, Bs, lams = [], [], np.zeros((n_steps, N - 1))
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
        Bs.append(B)

    # 系譜追跡で 3↔4 交換步を再同定
    K = N - 1
    labels = np.arange(K)
    ex_steps = []
    for k in range(n_steps - 1):
        O = np.abs(Vs[k].T @ Vs[k + 1])
        perm = greedy_assign(O)
        nl = np.zeros(K, int)
        nl[perm] = labels
        if (nl[2] == labels[3]) and (nl[3] == labels[2]):
            ex_steps.append(k)
        labels = nl

    results = []
    for ke in ex_steps:
        k0, k1 = max(0, ke - HALF), min(n_steps - 1, ke + HALF)
        P = Vs[k0][:, 2:4].copy()
        ab, cc = [], []
        for k in range(k0, k1 + 1):
            Q = Vs[k][:, 2:4]
            Pp = Q @ (Q.T @ P)
            Uq, _, Vtq = np.linalg.svd(Pp, full_matrices=False)
            P = Uq @ Vtq
            M = P.T @ Bs[k] @ P
            ab.append(M[0, 0] - M[1, 1])
            cc.append(2 * M[0, 1])
        ab, cc = np.array(ab), np.array(cc)
        norm = np.sqrt(ab ** 2 + cc ** 2)
        phi = np.unwrap(np.arctan2(cc, ab))
        dphi = float(phi[-1] - phi[0])
        sgn = np.sign(ab)
        flips = int(np.sum(sgn[:-1] * sgn[1:] < 0))
        d34 = lams[k0:k1 + 1, 2] - lams[k0:k1 + 1, 3]
        results.append({
            "tau": int(T_SKIP + 1 + ke),
            "n_ab_sign_flips_wide": flips,
            "delta_phi_rad": dphi,
            "delta_phi_over_pi": dphi / np.pi,
            "min_norm": float(norm.min()),
            "min_d34": float(d34.min()),
            "avoided_crossing": bool(abs(dphi) > np.pi / 2),
        })

    n_avoided = sum(r["avoided_crossing"] for r in results)
    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "window": [T_SKIP, T_END], "half_window": HALF,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行。"
                           "判別量は基底不変（軌跡ノルム・回転角）"},
        "n_events": len(results),
        "n_avoided_by_dphi": n_avoided,
        "min_gap_overall": float(min(r["min_norm"] for r in results)),
        "min_d34_overall": float(min(r["min_d34"] for r in results)),
        "events": results,
        "elapsed_sec": time.time() - t0,
    }
    (HERE / "result_dl3_event18_wide_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"イベント {len(results)} 件、|Δφ|>π/2（回避交差）: {n_avoided} 件")
    for r in results:
        print(f"  τ={r['tau']}: Δφ/π={r['delta_phi_over_pi']:+.3f} "
              f"flips={r['n_ab_sign_flips_wide']} minGap={r['min_norm']:.3e} "
              f"→ {'回避交差' if r['avoided_crossing'] else '接近のみ'}")
    print(f"最小有効ギャップ（全18）={res['min_gap_overall']:.3e}  "
          f"最小Δ34={res['min_d34_overall']:.3e} ({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
