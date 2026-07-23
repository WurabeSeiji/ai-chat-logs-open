#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N2 rank 測定拡張（検証済み第5論文コードの import による）

【再現性の規約】
  本スクリプトは独自再実装を一切行わない。対照テスト（control_test_v1.py）で
  公開ベースラインの完全一致（426フィールド, 最大差 0.0）を確認した
  run_spontaneous_splitting_preliminary_v1.py の力学関数
  （sine_generator / cayley / prepare_initial_state / line_graph_adjacency）を
  そのまま import して用いる。したがって本測定の力学は、公開済み第5論文の力学と
  ビット単位で同一である。

【測定内容（N2 の主張 = 平面鋳造・偶数階段）】
  第5論文の自己参照(再構成)正弦生成子力学を N=5,6 で走らせ、実効生成子 K(τ) の
  階数を時系列で測る:
    (1) 固有値 μ(iK) の ±対性残差 max_i|μ_i+μ_{-i}|（実反対称ゆえ階数は構造的偶数）。
    (2) 活性平面数 n(τ)=#{σ>閾}（各平面=σ対=rank 2）と rank(τ)=2n(τ)。
    (3) 増分イベントを検出し、全増分が +2（+1皆無）かを判定。
    (4) 閾値非依存性（{0.02,0.05,0.10}）。

実行: python3 rank_measurement_extension_v1.py
出力: rank_measurement_result_v1.json
"""

import importlib.util
import json
import os
import numpy as np

# --- 検証済み第5論文モジュールを import（別実装しない） ---
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "nbody_spontaneous_splitting_reproduction_v1",
                    "run_spontaneous_splitting_preliminary_v1.py")
_spec = importlib.util.spec_from_file_location("paper5_verified", _SRC)
p5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5)

# 検証済み関数のみを使う
sine_generator = p5.sine_generator
cayley = p5.cayley
prepare_initial_state = p5.prepare_initial_state
line_graph_adjacency = p5.line_graph_adjacency


def pairing_residual(K):
    """固有値 μ(iK) の ±対性残差。実反対称なら非零 rank は偶数（定理T2）。"""
    mu = np.sort(np.linalg.eigvalsh(1j * K))
    return float(np.max(np.abs(mu + mu[::-1])))


def active_planes(K, thr):
    return int(np.sum(np.linalg.eigvalsh(1j * K) > thr))


def run_rank(N, delta, seed, steps=4320, sub=8, thresholds=(0.02, 0.05, 0.10)):
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)  # 第5論文と同じ種規約
    Z, _residual = prepare_initial_state(rng, A, m, delta)  # 検証済み関数は (Z, residual) を返す
    series = {thr: [] for thr in thresholds}
    max_pair = 0.0
    for t in range(steps + 1):
        K = sine_generator(np.angle(Z), A)  # 検証済み関数
        if t % sub == 0:
            max_pair = max(max_pair, pairing_residual(K))
            for thr in thresholds:
                series[thr].append(active_planes(K, thr))
        if t < steps:
            Z = cayley(K) @ Z  # 検証済み関数
    base = np.array(series[0.05])
    rank = 2 * base
    inc = np.diff(rank)
    inc = inc[inc != 0]
    return {
        "N": N, "delta": delta, "seed": seed, "m": m,
        "max_rank": int(rank.max()),
        "max_pairing_residual": max_pair,
        "n_increment_events": int(inc.size),
        "odd_rank_increments": int(np.sum(np.abs(inc) % 2 != 0)),
        "increments_of_2": int(np.sum(np.abs(inc) == 2)),
        "even_not_2": int(np.sum((np.abs(inc) % 2 == 0) & (np.abs(inc) != 2))),
        "threshold_max_rank": {str(thr): int(2 * np.array(series[thr]).max())
                                for thr in thresholds},
    }


if __name__ == "__main__":
    print(f"import 元: {_SRC}")
    print(f"{'N':>2} {'delta':>7} {'seed':>4} {'maxRank':>7} {'pairRes':>10} "
          f"{'#evt':>4} {'odd':>3} {'+-2':>4} {'evOth':>5}")
    results = []
    for N in (5, 6):
        for d in (1e-3, 1e-4):
            for s in range(3):
                r = run_rank(N, d, s)
                results.append(r)
                print(f"{r['N']:>2} {r['delta']:>7.0e} {r['seed']:>4} "
                      f"{r['max_rank']:>7} {r['max_pairing_residual']:>10.2e} "
                      f"{r['n_increment_events']:>4} {r['odd_rank_increments']:>3} "
                      f"{r['increments_of_2']:>4} {r['even_not_2']:>5}")
    tot_odd = sum(r["odd_rank_increments"] for r in results)
    tot_oth = sum(r["even_not_2"] for r in results)
    max_pair = max(r["max_pairing_residual"] for r in results)
    tot_evt = sum(r["n_increment_events"] for r in results)
    print("-" * 62)
    print(f"全増分イベント数: {tot_evt}")
    print(f"奇数ランク増分: {tot_odd}（0 が予言）")
    print(f"±2以外の偶数増分: {tot_oth}")
    print(f"±対性残差の最大: {max_pair:.2e}（機械精度で0＝rank偶数）")
    with open(os.path.join(_HERE, "rank_measurement_result_v1.json"), "w") as f:
        json.dump({"import_source": _SRC, "runs": results,
                   "totals": {"events": tot_evt, "odd": tot_odd,
                              "even_not_2": tot_oth, "max_pairing_residual": max_pair}},
                  f, indent=2, ensure_ascii=False)
    print("結果を rank_measurement_result_v1.json に保存")
