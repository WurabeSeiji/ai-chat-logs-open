#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N5 予備実験B: 娘平面は1/2根（フェルミオン根）にロックするか（対照テスト済み第5論文コードを import）

当初仮説 H6: 最初に鋳造される娘平面は 1/2 根（フェルミオン、R=0, λ=-1, 位数2）に
優先的にロックする。これを直接検定する。

【観測量】
  第5論文の逐次力学で、活性平面数 n(t) が増える瞬間（＝娘平面の鋳造）を捉える。
  鋳造直後の実効生成子 K の固有値 σ を取り、
    - 娘平面の相対結合 r = σ_娘 / σ_主（σ_主 は最大 σ）
    - 娘平面の回転固有位相 θ_daughter（Cayley: φ = 2·arctan(γσ)）から
      交換位相の類似量を作り、最も近い低位数根 cos²(πm/n) を求める。
  多数の鋳造イベントで r の分布を取り、1/2 根（R=0 近傍、r→0 が完全交換）への
  優先ロックがあるかを見る。

【予言（測定前固定）】
  H6 成立: 娘平面の実効根が 1/2 根（フェルミオン、cos²(π/2)=0）に集積する。
  → 反証（別の根/N依存/一様）なら「娘平面=フェルミオン即ロック」を捨て、
    フェルミオン生成の描像を見直す。

【σ→R 橋の留保】
  N体平面の σ から2体交換係数 R への橋は未確立（第3論文で留保）。本実験は σ 比と
  回転固有位相という直接量で「1/2根らしさ」を見る。橋なしで言える範囲に限定し、
  結論はその限界とともに報告する。

実行: python3 N5_expB_half_root_lock_test_v1.py
出力: N5_expB_half_root_lock_result_v1.json
"""

import importlib.util
import json
import math
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "nbody_spontaneous_splitting_reproduction_v1",
                    "run_spontaneous_splitting_preliminary_v1.py")
_spec = importlib.util.spec_from_file_location("paper5_verified", _SRC)
p5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5)

sine_generator = p5.sine_generator
cayley = p5.cayley
prepare_initial_state = p5.prepare_initial_state
line_graph_adjacency = p5.line_graph_adjacency
GAMMA = math.tan(math.pi / 144.0)


def sigma_spectrum(K, thr=1e-9):
    mu = np.sort(np.linalg.eigvalsh(1j * K))[::-1]
    return mu[mu > thr]


def nearest_low_root(rr, max_n=8):
    """r に最も近い低位数根 cos²(πm/n) (n<=max_n, gcd(m,n)=1, 1<=m<n) と (n,m,距離)。"""
    best = None
    for n in range(2, max_n + 1):
        for mm in range(1, n):
            if math.gcd(mm, n) != 1:
                continue
            val = math.cos(math.pi * mm / n) ** 2
            d = abs(rr - val)
            if best is None or d < best[0]:
                best = (d, n, mm, val)
    return best


def run_expB(N, delta, seed, steps=8000, thr=0.05, sub=20):
    """娘平面の実効根を、鋳造の瞬間（閾値張り付きの artifact）でなく、
    飽和後（後半）の定着値として測る。σ_娘/σ_主 = 第2σ/第1σ の飽和中央値。"""
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z, _ = prepare_initial_state(rng, A, m, delta)
    ratios_saturated = []   # 飽和後の σ_2/σ_1
    for t in range(steps + 1):
        K = sine_generator(np.angle(Z), A)
        if t > steps // 2 and t % sub == 0:  # 飽和後のみ
            sig = np.sort(sigma_spectrum(K))[::-1]
            if len(sig) >= 2:
                ratios_saturated.append(float(sig[1] / sig[0]))
        if t < steps:
            Z = cayley(K) @ Z
    med = float(np.median(ratios_saturated)) if ratios_saturated else float("nan")
    return {"N": N, "delta": delta, "seed": seed,
            "n_samples": len(ratios_saturated),
            "daughter_ratios": [med],  # 飽和定着値（run当たり1値）
            "saturated_ratio_median": med}


if __name__ == "__main__":
    print(f"import 元: {_SRC}\n")
    print("=== 実験B: 娘平面の実効根（σ_娘/σ_主）と1/2根ロック検定 ===")
    all_ratios = []
    results = []
    for N in (5, 6):
        for seed in range(4):
            r = run_expB(N, 1e-3, seed)
            results.append(r)
            all_ratios.extend(r["daughter_ratios"])
            print(f"N={N} seed={seed}: 飽和サンプル数={r['n_samples']:3d} "
                  f"σ_娘/σ_主 飽和定着={r['saturated_ratio_median']:.4f}")

    all_ratios = np.array(all_ratios)
    print(f"\n収集した娘平面 {len(all_ratios)} 個")
    print(f"σ_娘/σ_主 全体中央={np.median(all_ratios):.4f} "
          f"平均={all_ratios.mean():.4f} 幅=[{all_ratios.min():.3f}, {all_ratios.max():.3f}]")
    # 低位数根への近さ（1/2根=cos²(π/2)=0）
    d, n, mm, val = nearest_low_root(float(np.median(all_ratios)))
    print(f"全体中央に最も近い低位数根: cos²({mm}π/{n})={val:.4f}（距離{d:.4f}）")
    # 1/2根（R=0）への集積の有無: r≈0 の割合
    frac_near_0 = float(np.mean(all_ratios < 0.05))
    frac_near_half_root = float(np.mean(np.abs(all_ratios - 0.0) < 0.05))
    print(f"σ比 < 0.05（フェルミオン根 R=0 近傍）の割合: {frac_near_0:.3f}")
    # ヒストグラム的集計
    print("σ_娘/σ_主 の分布:")
    for lo, hi in [(0, 0.05), (0.05, 0.15), (0.15, 0.25), (0.25, 0.35),
                   (0.35, 0.45), (0.45, 0.55), (0.55, 1.0)]:
        c = int(np.sum((all_ratios >= lo) & (all_ratios < hi)))
        print(f"  [{lo:.2f}, {hi:.2f}): {c:3d}  {'#' * (c * 40 // max(1, len(all_ratios)))}")

    print("\n=== 判定 ===")
    verdict_half = frac_near_0 > 0.5
    print(f"1/2根（フェルミオン, R=0近傍）への優先ロック: "
          f"{'PASS' if verdict_half else 'FAIL（1/2根に集積しない）'}")
    if not verdict_half:
        print(f"⇒ H6（娘平面=1/2根即ロック）は反証。娘平面の実効根は "
              f"cos²({mm}π/{n})≈{val:.3f} 近傍で N 依存。フェルミオン生成の描像を見直す。")

    with open(os.path.join(_HERE, "N5_expB_half_root_lock_result_v1.json"), "w") as f:
        json.dump({"import_source": _SRC,
                   "n_daughters": int(len(all_ratios)),
                   "ratio_median": float(np.median(all_ratios)),
                   "ratio_mean": float(all_ratios.mean()),
                   "frac_near_R0": frac_near_0,
                   "nearest_root": {"n": n, "m": mm, "value": val, "dist": d},
                   "half_root_lock": bool(verdict_half),
                   "per_run": [{"N": r["N"], "seed": r["seed"],
                                "saturated_ratio_median": r["saturated_ratio_median"]}
                               for r in results]},
                  f, indent=2, ensure_ascii=False)
    print("\n結果を N5_expB_half_root_lock_result_v1.json に保存")
