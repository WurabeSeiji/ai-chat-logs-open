#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
検証[1]: de Broglie 論文の除外率 1.5〜4.0% の理論式化（論文 展望(b) への回答）

論文の判定手続き（debroglie_align_lambda.py と §2.5 より）:
  A(λ) = |Σ_{n odd ≤ N_max} e^{i2πnr/λ}|² / K²,  K = (N_max+1)/2（成分数）
  u ≡ 2πr/λ とおくと A(u) = sin²(Ku)/(K² sin²u)（ディリクレ核²、周期 π）
  幾何窓 λ₀′ ± λ₀′²/(4r) は u 座標でちょうど半セル幅 ±π/2
  窓内の A≥0.98 連結帯の数がちょうど1でなければ NG（除外）

理論:
  共鳴帯の半幅 ε*: A(ε*) = A_th → 小角展開 A ≈ 1 − (K²−1)ε²/3
      ε* = √(3(1−A_th)/(K²−1))
  NG は λ₀′（→u₀）がセル中間点の幅 w = 2ε* の帯に落ちたとき（隣接2共鳴が両方窓に入る）
  波長±1%の乱択はセル間隔（相対10⁻⁹）より遥かに広い → u₀ はセル内一様
      除外率 = w/π = (2/π)·√(3(1−A_th)/(K²−1))

検証:
  (i) 帯数え上げの直接シミュレーションと閉形式の一致
  (ii) 論文の観測値 1.5〜4.0%・n_OK=192〜197/200 との比較（N_max=9→K=5, N_max=17→K=9）
"""
import math

PI = math.pi
A_TH = 0.98

def A_of_u(u, K):
    s = math.sin(u)
    if abs(s) < 1e-14:
        return 1.0
    return (math.sin(K * u) / (K * s)) ** 2

def eps_star(K, a_th=A_TH):
    """共鳴帯の半幅（厳密には二分法、比較用に小角式も）"""
    lo, hi = 0.0, PI / (2 * K)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if A_of_u(mid, K) > a_th:
            lo = mid
        else:
            hi = mid
    exact = lo
    small_angle = math.sqrt(3 * (1 - a_th) / (K * K - 1))
    return exact, small_angle

def rate_formula(K, a_th=A_TH):
    exact, sa = eps_star(K, a_th)
    return 2 * exact / PI, 2 * sa / PI

def rate_simulation(K, a_th=A_TH, n_u0=4000, n_scan=6000):
    """帯数え上げの直接シミュレーション: u₀ をセル内一様に置き、窓 [u₀-π/2, u₀+π/2]
    内の A≥a_th 連結帯を数える。帯数≠1 → NG。"""
    ng = 0
    for i in range(n_u0):
        u0 = (i + 0.5) / n_u0 * PI  # セル (0, π) 内一様
        lo, hi = u0 - PI / 2, u0 + PI / 2
        prev = False
        bands = 0
        for j in range(n_scan + 1):
            u = lo + (hi - lo) * j / n_scan
            cur = A_of_u(u, K) >= a_th
            if cur and not prev:
                bands += 1
            prev = cur
        if bands != 1:
            ng += 1
    return ng / n_u0

def main():
    print("[1] 除外率の理論式化と検証")
    print(f"    閉形式: 除外率 = (2/π)·√(3(1−A_th)/(K²−1)),  A_th = {A_TH}")
    print()
    print(f"    {'N_max':>6} {'K':>3} {'閉形式(厳密)':>11} {'閉形式(小角)':>11} {'シミュレーション':>13} {'200試行の期待NG数':>16}")
    for N_max in [5, 9, 13, 17, 33]:
        K = (N_max + 1) // 2
        r_exact, r_sa = rate_formula(K)
        r_sim = rate_simulation(K)
        print(f"    {N_max:>6} {K:>3} {r_exact:>10.4%} {r_sa:>10.4%} {r_sim:>12.4%} {200*r_exact:>13.1f}")
    print()
    print("    論文の観測値との比較:")
    for N_max in [9, 17]:
        K = (N_max + 1) // 2
        r, _ = rate_formula(K)
        import statistics
        sd = math.sqrt(200 * r * (1 - r)) / 200
        print(f"    N_max={N_max} (K={K}): 予言 {r:.2%} ± {sd:.2%}(二項揺らぎ, n=200)"
              f" → n_OK 期待 {200*(1-r):.1f}, 1σ範囲 [{200*(1-r)-200*sd:.0f}, {200*(1-r)+200*sd:.0f}]")
    print(f"    論文の観測: 除外率 1.5〜4.0%, n_OK = 192〜197/200")

if __name__ == "__main__":
    main()
