#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 厳密交差か回避交差か：有効2×2行列の off-diagonal 結合（第4回査読対応）

系譜ラベル交換（18回）は「枝交換」の証明であり「厳密縮退 λ3=λ4」の証明ではない
（回避交差でも固有ベクトルの性格は入れ替わる）。本プローブは各イベント近傍で
輸送された非固有2次元基底 P（初期化＝イベント前の [v3,v4]、毎步で現行 span{v3,v4} へ
射影して正規直交化＝断熱平面内の平行輸送）を構成し、有効2×2行列

    M(τ) = P^T B(τ) P = [[a, c], [c, b]]

の軌跡 (a−b, 2c) を測る。実対称2×2の固有値差は √((a−b)²+(2c)²) なので、
  厳密交差 ⇔ 軌跡が原点を通る（a−b と c が同時に零）
  回避交差 ⇔ a−b は符号反転するが |2c| が有限に残る（最小ギャップ ≈ 2|c*|）

判定（実行前固定）: 各イベントで a−b の符号反転点（線形補間）における 2|c*| を測り、
  2|c*| が整数步での最小 Δ34 と同程度（比 > 0.5）→ 回避交差（有限ギャップが解像済み）
  2|c*| ≪ 最小 Δ34（比 < 0.1）→ 厳密交差候補（さらなる細分化が必要）
  中間 → 個別に記録

出力: result_dl3_avoided_crossing_v1.json
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
HALF = 25   # イベント前後の解析窓


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
    u1 = _load("uni_dl3ac", UNI / "unified_interaction_v1.py")
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

    # 系譜追跡で 3↔4 交換イベント步を再同定（probe_dl3_branch_lineage_v1 と同一手続き）
    K = N - 1
    labels = np.arange(K)
    ex_steps = []
    for k in range(n_steps - 1):
        O = np.abs(Vs[k].T @ Vs[k + 1])
        perm = greedy_assign(O)
        new_labels = np.zeros(K, int)
        new_labels[perm] = labels
        if (new_labels[2] == labels[3]) and (new_labels[3] == labels[2]):
            ex_steps.append(k)
        labels = new_labels

    # 各イベントの有効2×2解析
    events = []
    for ke in ex_steps:
        k0 = max(0, ke - HALF)
        k1 = min(n_steps - 1, ke + HALF)
        # 断熱2次元平面内の平行輸送基底
        P = Vs[k0][:, 2:4].copy()
        ab, cc, d34 = [], [], []
        for k in range(k0, k1 + 1):
            Q = Vs[k][:, 2:4]
            Pp = Q @ (Q.T @ P)                     # span{v3,v4} へ射影
            Uq, _, Vtq = np.linalg.svd(Pp, full_matrices=False)
            P = Uq @ Vtq                            # 正規直交化（平面内平行輸送）
            M = P.T @ Bs[k] @ P
            ab.append(M[0, 0] - M[1, 1])
            cc.append(2 * M[0, 1])
            d34.append(lams[k, 2] - lams[k, 3])
        ab, cc, d34 = map(np.array, (ab, cc, d34))
        min_d34 = float(d34.min())
        # a−b の符号反転点での |2c|（線形補間）
        sgn = np.sign(ab)
        cross_idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
        c_at_cross = []
        for i in cross_idx:
            w = abs(ab[i]) / (abs(ab[i]) + abs(ab[i + 1]))
            c_at_cross.append(abs((1 - w) * cc[i] + w * cc[i + 1]))
        gap_eff = float(min(c_at_cross)) if c_at_cross else None
        # 軌跡の原点最短距離（離散点の (a−b, 2c) ノルム最小）
        min_norm = float(np.sqrt(ab ** 2 + cc ** 2).min())
        ratio = (gap_eff / min_d34) if (gap_eff is not None and min_d34 > 0) else None
        events.append({"tau": int(T_SKIP + 1 + ke),
                       "n_sign_flips_ab": int(len(cross_idx)),
                       "min_d34_integer_steps": min_d34,
                       "gap_at_ab_zero_2c": gap_eff,
                       "min_traj_norm": min_norm,
                       "ratio_gap_to_min_d34": ratio})

    ratios = [e["ratio_gap_to_min_d34"] for e in events
              if e["ratio_gap_to_min_d34"] is not None]
    n_avoided = sum(r > 0.5 for r in ratios)
    n_exact_cand = sum(r < 0.1 for r in ratios)
    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "window": [T_SKIP, T_END], "half_window": HALF,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行。"
                           "基底＝断熱2平面内の平行輸送（非固有基底）"},
        "n_events": len(events),
        "n_with_ab_sign_flip": sum(e["n_sign_flips_ab"] > 0 for e in events),
        "n_avoided_ratio_gt_0.5": n_avoided,
        "n_exact_candidate_ratio_lt_0.1": n_exact_cand,
        "ratio_median": float(np.median(ratios)) if ratios else None,
        "ratio_min": float(np.min(ratios)) if ratios else None,
        "gap_eff_median": float(np.median([e["gap_at_ab_zero_2c"] for e in events
                                           if e["gap_at_ab_zero_2c"] is not None])),
        "events": events,
        "elapsed_sec": time.time() - t0,
    }
    (HERE / "result_dl3_avoided_crossing_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"イベント {len(events)} 件（a−b 符号反転あり "
          f"{res['n_with_ab_sign_flip']} 件）")
    print(f"比 2|c*|/minΔ34: 中央値={res['ratio_median']}  最小={res['ratio_min']}")
    print(f"回避交差(比>0.5): {n_avoided}  厳密交差候補(比<0.1): {n_exact_cand}")
    print(f"有効ギャップ 2|c*| 中央値={res['gap_eff_median']:.3e}")
    print(f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
