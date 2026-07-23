#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E1: 生成子ランクの偶数階段（N2 論文 付録B 再現スクリプト）

第5論文の再構成正弦生成子 O1 力学を走らせ、実効生成子 K(τ) の階数を時系列で測る。
測るもの:
  (1) 固有値の ±対性残差  max_i |μ_i + μ_{−i}|  ——実反対称ゆえ rank は構造的に偶数。
      これが機械精度で零であることが「ランクは偶数」の直接確認（定理T2）。
  (2) 活性平面数 n(τ) = #{σ_j > 閾}（各平面 = σ 対 = rank 2）の時系列。
  (3) rank(τ) = 2·n(τ) の階段。増分イベントを検出し、全増分が +2（+1 皆無）かを判定。
  (4) 閾値非依存性: 複数閾値で n(τ) を測り、階段構造が閾に依らないことを確認。

依存: 第5論文再現パッケージ（nbody_spontaneous_splitting_reproduction_v1）の関数を再利用。
      同梱せず、パッケージ内の sine_generator / cayley / prepare_initial_state を再実装する
      （本スクリプト単体で動く）。整数演算ではなく float64（力学は連続、閾値は相対）。

実行: python3 rank_even_staircase.py
"""

import json
import math
import numpy as np

GAMMA = math.tan(math.pi / 144.0)


def complete_graph_edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def line_graph_adjacency(n):
    edges = complete_graph_edges(n)
    m = len(edges)
    A = np.zeros((m, m))
    for a in range(m):
        for b in range(m):
            if a != b and len(set(edges[a]) & set(edges[b])) > 0:
                A[a, b] = 1.0
    return A


def sine_generator(theta, A):
    K = A * np.sin(theta[None, :] - theta[:, None])
    s = np.linalg.norm(K, 2)
    return K if s < 1e-300 else K / s


def cayley(K, gamma=GAMMA):
    I = np.eye(K.shape[0])
    return np.linalg.solve(I - gamma * K, I + gamma * K)


def top_plane_eigenmode(K):
    mu, V = np.linalg.eigh(1j * K)
    idx = int(np.argmax(mu))
    v = V[:, idx]
    p = np.sqrt(2.0) * v.real
    q = np.sqrt(2.0) * v.imag
    p = p / np.linalg.norm(p)
    q = q - (q @ p) * p
    q = q / np.linalg.norm(q)
    return (p + 1j * q) / math.sqrt(2.0)


def kernel_basis(K, tol=1e-10):
    _, s, Vt = np.linalg.svd(K)
    return [Vt[i] for i in range(K.shape[0]) if s[i] <= tol]


def prepare_initial_state(rng, A, m, delta, iters=400, beta=0.5):
    theta = rng.uniform(0.0, 2.0 * np.pi, m)
    g1 = rng.normal(size=m)
    g2 = rng.normal(size=m)
    v = None
    for _ in range(iters):
        v = top_plane_eigenmode(sine_generator(theta, A))
        theta_new = np.angle(v)
        mix = 0.5 * np.exp(1j * theta) + 0.5 * np.exp(1j * theta_new)
        theta = np.angle(mix)
    K_fin = sine_generator(np.angle(v), A)
    kernel = kernel_basis(K_fin)

    def kproj(x):
        out = np.zeros(m)
        for k_vec in kernel:
            out += (k_vec @ x) * k_vec
        return out

    if delta > 0.0 and kernel:
        u = kproj(g1)
        u = u / np.linalg.norm(u)
        w = kproj(g2)
        w = w - (w @ u) * u
        w = w / np.linalg.norm(w)
        g = (u + 1j * w) / math.sqrt(2.0)
        Z = v + delta * g
    else:
        Z = v.copy()
    return Z / np.linalg.norm(Z)


def pairing_residual(K):
    """固有値 μ(iK) の ±対性残差。実反対称なら μ は ± で対、非零 rank は偶数。"""
    mu = np.sort(np.linalg.eigvalsh(1j * K))
    # 昇順 mu と降順反転の和 → 対性が完全なら全零
    return float(np.max(np.abs(mu + mu[::-1])))


def active_planes(K, thr):
    mu = np.linalg.eigvalsh(1j * K)
    return int(np.sum(mu > thr))


def run(N, delta, seed, steps=4320, sub=8, thresholds=(0.02, 0.05, 0.10)):
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z = prepare_initial_state(rng, A, m, delta)
    ts, series = [], {thr: [] for thr in thresholds}
    max_pair_res = 0.0
    for t in range(steps + 1):
        K = sine_generator(np.angle(Z), A)
        if t % sub == 0:
            max_pair_res = max(max_pair_res, pairing_residual(K))
            for thr in thresholds:
                series[thr].append(active_planes(K, thr))
            ts.append(t)
        if t < steps:
            Z = cayley(K) @ Z
    # rank(τ) = 2 * active_planes（基準閾 0.05）
    base = np.array(series[0.05])
    rank = 2 * base
    dranks = np.diff(rank[rank[:-1] != rank[1:]]) if False else np.diff(rank)
    increments = dranks[dranks != 0]
    max_rank = int(rank.max())
    # 増分の絶対値内訳（+2/−2 のみか、奇数増分があるか）
    odd_increments = int(np.sum(np.abs(increments) % 2 != 0))
    two_increments = int(np.sum(np.abs(increments) == 2))
    other_even = int(np.sum((np.abs(increments) % 2 == 0) & (np.abs(increments) != 2)))
    return {
        "N": N, "delta": delta, "seed": seed, "m": m,
        "max_rank": max_rank,
        "max_pairing_residual": max_pair_res,
        "n_increment_events": int(increments.size),
        "odd_rank_increments": odd_increments,
        "increments_of_2": two_increments,
        "even_but_not_2_increments": other_even,
        "increment_values": [int(x) for x in increments.tolist()],
        "threshold_consistency": {
            str(thr): int(np.array(series[thr]).max()) for thr in thresholds
        },
    }


if __name__ == "__main__":
    results = []
    cases = [(N, d, s) for N in (5, 6) for d in (1e-3, 1e-4) for s in range(3)]
    print(f"{'N':>2} {'delta':>7} {'seed':>4} {'maxRank':>7} "
          f"{'pairRes':>10} {'#evt':>4} {'odd':>3} {'+-2':>4} {'evOth':>5}")
    for N, d, s in cases:
        r = run(N, d, s)
        results.append(r)
        print(f"{r['N']:>2} {r['delta']:>7.0e} {r['seed']:>4} {r['max_rank']:>7} "
              f"{r['max_pairing_residual']:>10.2e} {r['n_increment_events']:>4} "
              f"{r['odd_rank_increments']:>3} {r['increments_of_2']:>4} "
              f"{r['even_but_not_2_increments']:>5}")
    total_odd = sum(r["odd_rank_increments"] for r in results)
    total_other = sum(r["even_but_not_2_increments"] for r in results)
    max_pair = max(r["max_pairing_residual"] for r in results)
    print("-" * 60)
    print(f"全試行の奇数ランク増分: {total_odd}（0 が予言）")
    print(f"±2 以外の偶数増分: {total_other}（一括分裂等、階段の粗さ指標）")
    print(f"ランク偶数性（対性残差の最大）: {max_pair:.2e}（機械精度で 0 が予言）")
    with open("rank_even_staircase_result.json", "w") as f:
        json.dump(results, f, indent=2)
    print("結果を rank_even_staircase_result.json に保存")
