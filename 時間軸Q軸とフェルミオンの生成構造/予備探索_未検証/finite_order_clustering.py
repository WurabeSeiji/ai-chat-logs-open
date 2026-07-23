#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E4 / N3: モノドロミーの well-definedness と有限位数根への集積の測定

二段構成:
  第一段（前提作業・モノドロミー well-definedness）:
    自己参照(再構成)力学の飽和後の準安定状態で、各参照平面の回転位相速度を
    複数のスライド窓で測り、窓をずらしても値が一定に収束するかを調べる。
    収束すれば「平面ごとの実効回転は時間窓に依らず定義できる」——E4以降の
    測定量が well-defined であることの確認。
  第二段（H1/H3・根への集積）:
    飽和後の平面間位相速度比 ω_j/ω_1 を全試行で集め、低位数有理数
    （m/n, n≤6）への集積を測る。一様分布なら仮説棄却、低位数へのピークなら
    有限位数根構造の証拠。

測定量の定義:
  各参照平面 j に状態 Z を射影して複素振幅 a_j = p_j·Z, b_j = q_j·Z を得る。
  双線形レジスタ c_j = a_j² + b_j² の偏角 arg(c_j) は平面内回転の 2 倍角で進む。
  位相速度 ω_j = d arg(c_j)/dt を飽和窓で最小二乗推定する。

依存: 第5論文再現パッケージの力学を単体再実装（前スクリプトと共通）。float64。

実行: python3 finite_order_clustering.py
"""

import json
import math
from fractions import Fraction
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


def plane_decomposition(K, tol=1e-10):
    mu, V = np.linalg.eigh(1j * K)
    order = np.argsort(-mu)
    planes = []
    for idx in order:
        m_val = float(mu[idx])
        if m_val <= tol:
            continue
        v = V[:, idx]
        p = np.sqrt(2.0) * v.real
        q = np.sqrt(2.0) * v.imag
        p = p / np.linalg.norm(p)
        q = q - (q @ p) * p
        q = q / np.linalg.norm(q)
        planes.append({"sigma": m_val, "p": p, "q": q})
    return planes


def kernel_basis(K, tol=1e-10):
    _, s, Vt = np.linalg.svd(K)
    return [Vt[i] for i in range(K.shape[0]) if s[i] <= tol]


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


def prepare_initial_state(rng, A, m, delta, iters=400):
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
        u = kproj(g1); u = u / np.linalg.norm(u)
        w = kproj(g2); w = w - (w @ u) * u; w = w / np.linalg.norm(w)
        g = (u + 1j * w) / math.sqrt(2.0)
        Z = v + delta * g
    else:
        Z = v.copy()
    return Z / np.linalg.norm(Z)


def run_and_track(N, delta, seed, steps=12000, sub=2):
    """力学を走らせ、参照平面ごとの複素双線形レジスタ c_j(τ) の時系列を返す。
    参照平面は初期の自己無撞着生成子の平面分解＋核から張った固定枠。"""
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z = prepare_initial_state(rng, A, m, delta)
    # 参照枠: 初期生成子の全固有平面（活性+休眠）を固定基底として使う
    K0 = sine_generator(np.angle(Z), A)
    mu, V = np.linalg.eigh(1j * K0)
    order = np.argsort(-mu)
    frames = []
    used = np.zeros(m, dtype=bool)
    for idx in order:
        if mu[idx] <= 1e-9:
            continue
        v = V[:, idx]
        p = np.sqrt(2.0) * v.real; p /= np.linalg.norm(p)
        q = np.sqrt(2.0) * v.imag; q = q - (q @ p) * p; q /= np.linalg.norm(q)
        frames.append((p, q))
    # 休眠平面も加える（核から直交対を作る）
    ker = kernel_basis(K0)
    for a in range(0, len(ker) - 1, 2):
        p = ker[a] / np.linalg.norm(ker[a])
        q = ker[a + 1] - (ker[a + 1] @ p) * p
        nq = np.linalg.norm(q)
        if nq > 1e-9:
            frames.append((p, q / nq))
    nfr = len(frames)
    c_hist = np.empty((steps // sub + 1, nfr), dtype=complex)
    h_hist = np.empty((steps // sub + 1, nfr))
    ti = 0
    for t in range(steps + 1):
        if t % sub == 0:
            for j, (p, q) in enumerate(frames):
                a = p @ Z; b = q @ Z
                c_hist[ti, j] = a * a + b * b
                h_hist[ti, j] = abs(a) ** 2 + abs(b) ** 2
            ti += 1
        if t < steps:
            Z = cayley(sine_generator(np.angle(Z), A)) @ Z
    return c_hist[:ti], h_hist[:ti], sub


def phase_velocity(c_series, dt):
    """arg(c) の連続化位相速度を最小二乗で推定。"""
    ph = np.unwrap(np.angle(c_series))
    tt = np.arange(len(ph)) * dt
    A = np.vstack([tt, np.ones_like(tt)]).T
    slope, _ = np.linalg.lstsq(A, ph, rcond=None)[0]
    return slope


def nearest_low_order(x, max_n=6):
    """x に最も近い低位数有理数 m/n (n<=max_n) と距離。"""
    best = None
    for n in range(1, max_n + 1):
        for mm in range(0, n * 3 + 1):
            val = mm / n
            d = abs(x - val)
            if best is None or d < best[0]:
                best = (d, Fraction(mm, n))
    return best


if __name__ == "__main__":
    # ---- 第一段: モノドロミー well-definedness ----
    print("=== 第一段: モノドロミー well-definedness（主平面回転率の窓非依存性）===")
    print(f"{'N':>2} {'delta':>7} {'seed':>4} {'w1':>10} {'w2':>10} {'w3':>10} {'rel_spread':>10}")
    mono_spreads = []
    track_cache = {}
    for N in (5, 6):
        for delta in (1e-3, 1e-4):
            for seed in range(3):
                c_hist, h_hist, sub = run_and_track(N, delta, seed)
                track_cache[(N, delta, seed)] = (c_hist, h_hist, sub)
                nt = c_hist.shape[0]
                sat0 = nt // 2  # 飽和後（後半）
                # 主平面 = 飽和窓で平均レジスタ h が最大の枠
                hbar = h_hist[sat0:].mean(axis=0)
                jmax = int(np.argmax(hbar))
                seg = c_hist[sat0:, jmax]
                L = len(seg)
                # 3つのスライド窓で位相速度を測る
                ws = []
                for k in range(3):
                    a = k * L // 4
                    b = a + L // 2
                    ws.append(phase_velocity(seg[a:b], sub))
                ws = np.array(ws)
                spread = float(np.std(ws) / (abs(np.mean(ws)) + 1e-30))
                mono_spreads.append(spread)
                print(f"{N:>2} {delta:>7.0e} {seed:>4} {ws[0]:>10.5f} {ws[1]:>10.5f} "
                      f"{ws[2]:>10.5f} {spread:>10.2e}")
    print(f"主平面回転率の窓間相対ばらつき 中央値: {np.median(mono_spreads):.2e} "
          f"最大: {np.max(mono_spreads):.2e}")

    # ---- 第二段: 平面間位相速度比の集積 ----
    print("\n=== 第二段: 飽和後の平面間位相速度比 ω_j/ω_1 と低位数根への近さ ===")
    ratios = []
    for (N, delta, seed), (c_hist, h_hist, sub) in track_cache.items():
        nt = c_hist.shape[0]
        sat0 = nt // 2
        hbar = h_hist[sat0:].mean(axis=0)
        # 活性枠（相対レジスタが閾以上）だけ
        active = np.where(hbar / hbar.max() > 0.01)[0]
        if len(active) < 2:
            continue
        vel = {j: phase_velocity(c_hist[sat0:, j], sub) for j in active}
        jref = max(active, key=lambda j: hbar[j])
        w1 = vel[jref]
        if abs(w1) < 1e-9:
            continue
        for j in active:
            if j == jref:
                continue
            r = vel[j] / w1
            ratios.append(r)
    ratios = np.array(ratios)
    print(f"収集した平面間比の数: {len(ratios)}")
    if len(ratios):
        dists = np.array([float(nearest_low_order(abs(r))[0]) for r in ratios])
        print(f"|比| の低位数根(n<=6)への距離 中央値: {np.median(dists):.3f} "
              f"平均: {dists.mean():.3f}")
        # 一様分布(0..3)なら期待距離 ~ 平均最近傍間隔/4
        print("比の値（絶対値、最初の20個）:",
              np.round(np.sort(np.abs(ratios))[:20], 3).tolist())
        # ヒストグラム的集計
        for lo, hi, label in [(0, 0.05, "≈0(n=2根)"), (0.20, 0.30, "≈1/4"),
                               (0.45, 0.55, "≈1/2"), (0.95, 1.05, "≈1"),
                               (0.70, 0.80, "≈3/4")]:
            cnt = int(np.sum((np.abs(ratios) >= lo) & (np.abs(ratios) < hi)))
            print(f"  |比| ∈ [{lo},{hi}) {label}: {cnt}")
    out = {
        "monodromy_spread_median": float(np.median(mono_spreads)),
        "monodromy_spread_max": float(np.max(mono_spreads)),
        "n_ratios": int(len(ratios)),
        "ratio_lowroot_dist_median": float(np.median(dists)) if len(ratios) else None,
        "ratios_abs_sorted": np.round(np.sort(np.abs(ratios)), 4).tolist() if len(ratios) else [],
    }
    with open("finite_order_clustering_result.json", "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("\n結果を finite_order_clustering_result.json に保存")
