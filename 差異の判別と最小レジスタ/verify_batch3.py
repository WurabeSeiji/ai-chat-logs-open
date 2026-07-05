#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
検証バッチ3（2026-07-05）
[3] Born cos² への収束レート：有限N核読出しの可視度定理 V(m,N) = 1 − m/N
[S] Sargent 則 Q⁵ の格子セル計数：3体終状態の離散計数が連続体則を再現するか＋小Q偏差

[3] の理論:
  |W_N(θ)|² = Σ_{j=-(N-1)}^{N-1} (N-|j|) e^{i2jθ}  （三角窓＝Fejér型係数）
  よって正規化核 F_N = |W_N|²/(πN) との畳み込みは、周波数 2m の縞の振幅を
  ちょうど (1 - m/N) 倍する（m ≥ N で零＝分解能カットオフ）。
  検証: (i) フーリエ係数の直接計算で三角窓を確認
        (ii) cos² 縞との数値畳み込みで可視度 = 1 - m/N を確認

[S] の設定:
  β崩壊型の3体終状態（核は運動量吸収のみ）。電子・ニュートリノの運動量を
  整数格子（間隔1）に置き、超相対論近似 E=|p| でエネルギー殻 E_e+E_ν=Q の
  セル数を数える。連続体則: dN/dQ ∝ Q⁵（Sargent）。
  検証: log-log 傾き ≈ 5（大Q）、小Qでの離散偏差の定量化。
"""
import cmath
import math

PI = math.pi

# ---------- [3] 可視度定理 ----------

def W(theta, N):
    return sum(cmath.exp(1j * (2 * k - 1) * theta) for k in range(1, N + 1))

def fourier_coeff_of_kernel(N, j, M=4096):
    """|W_N(θ)|² の e^{i2jθ} 係数（数値、π周期上）"""
    s = 0.0 + 0.0j
    for i in range(M):
        th = PI * i / M
        s += abs(W(th, N)) ** 2 * cmath.exp(-2j * j * th)
    return (s / M).real

def visibility_by_convolution(N, m, M=4096):
    """縞 1+cos(2mθ) を核 F_N と数値畳み込みし、出力可視度を測る"""
    F = [abs(W(PI * i / M, N)) ** 2 for i in range(M)]
    total = sum(F)
    out_max, out_min = -1e18, 1e18
    for i in range(0, M, 16):  # 出力を粗くサンプル
        acc = 0.0
        for j in range(0, M, 4):  # 積分も間引き（十分滑らか）
            th = PI * (i - j) / M
            acc += F[j] * (1.0 + math.cos(2 * m * th))
        val = acc
        out_max = max(out_max, val)
        out_min = min(out_min, val)
    return (out_max - out_min) / (out_max + out_min)

def task3():
    print("[3] Born cos² への収束：可視度定理 V(m,N) = 1 − m/N")
    N = 9
    print(f"    (i) |W_{N}|² のフーリエ係数（理論値 N−|j|）:")
    for j in [0, 1, 3, 6, 8, 9]:
        c = fourier_coeff_of_kernel(N, j)
        print(f"        j={j}: 数値 {c:9.4f}  理論 {max(N - abs(j), 0):3d}")
    print(f"    (ii) 畳み込みで測った可視度 vs 1 − m/N:")
    for (N_, m) in [(9, 1), (9, 3), (9, 6), (9, 8), (5, 2), (17, 4)]:
        v = visibility_by_convolution(N_, m)
        print(f"        N={N_:2d}, m={m}: 測定 {v:.4f}  理論 {1 - m / N_:.4f}")

# ---------- [S] Sargent 則の格子計数 ----------

def shell_counts(R):
    """3次元整数格子点を |p| のビン（幅1、中心 E=1,2,...）に計数"""
    bins = [0] * (R + 2)
    for x in range(-R, R + 1):
        for y in range(-R, R + 1):
            zmax = int(math.isqrt(max(R * R - x * x - y * y, 0)))
            for z in range(-zmax, zmax + 1):
                E = math.sqrt(x * x + y * y + z * z)
                b = int(E + 0.5)
                if 1 <= b <= R:
                    bins[b] += 1
    return bins

def task_sargent(R=40):
    print()
    print("[S] Sargent 則 Q⁵ の格子セル計数（3体終状態、E=|p| 超相対論）")
    n = shell_counts(R)
    # 崩壊率(Q) = Σ_E n_e(E)·n_ν(Q−E)
    rate = {}
    for Q in range(2, R + 1):
        rate[Q] = sum(n[E] * n[Q - E] for E in range(1, Q))
    # log-log 傾き（大Q側 Q=15..R）
    xs = [math.log(Q) for Q in range(15, R + 1)]
    ys = [math.log(rate[Q]) for Q in range(15, R + 1)]
    nfit = len(xs)
    xm = sum(xs) / nfit
    ym = sum(ys) / nfit
    slope = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sum((x - xm) ** 2 for x in xs)
    print(f"    大Q側 (Q=15..{R}) の log-log 傾き = {slope:.3f}（Sargent 則の理論値 5）")
    # 小Q偏差: 連続体則 rate_cont = C·Q⁵ を大Q側で規格化し、小Qの相対偏差を見る
    C = rate[R] / R ** 5
    print(f"    小Q領域の連続体則 C·Q⁵ からの相対偏差:")
    print(f"    {'Q':>4} {'格子計数':>12} {'C·Q⁵':>12} {'偏差%':>8}")
    for Q in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, R]:
        cont = C * Q ** 5
        dev = 100 * (rate[Q] / cont - 1)
        print(f"    {Q:>4} {rate[Q]:>12d} {cont:>12.0f} {dev:>+8.1f}")

if __name__ == "__main__":
    task3()
    task_sargent()
