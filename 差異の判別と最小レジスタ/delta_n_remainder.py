#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
規格化・次元の剰余 δ(N) の計算（v2、2026-07-05 木原氏の訂正反映）

訂正：局在ピークを作るのは等振幅の奇数倍音（再生核型）。
振幅 1/n 減衰（矩形波係数）はエッジ（階段）であり局在しない。
→ 物理モデル（粒子=局在波）は等振幅・有限N。1/n 減衰は対照モデルに降格。

主モデル（等振幅局在波・de Broglie 論文型）:
  W_N(θ) = Σ_{k=1}^{N} e^{i(2k-1)θ} = e^{iNθ} sin(Nθ)/sin(θ)
  |W_N|² = sin²(Nθ)/sin²θ、ピーク N²、第一零点 θ₀ = π/(N... 正確には主ローブ半幅 π/(2N)×2）
  ここでは第一零点 θ₀: sin(Nθ)=0 の最小正解 θ₀ = π/N
  Parseval: (1/2π)∫|W_N|²dθ = N

  剰余の住所（等振幅では規格化は厳密なので、剰余は別の場所に住む）:
  (a) 漏れ率 L(N): 主ローブ外に漏れる強度の割合
  (b) モード数の識別可能性: |⟨Ŵ_N|Ŵ_{N+1}⟩|² = N/(N+1)
      → 隣接 N の区別可能性 = 1/(N+1) ＝ 論文の「鋭さ 1/(N+1)」と同一
  (c) 小 N での不確定性積の欠損: σ_n·θ₀ vs 極限 π/√3

対照モデル（矩形波・振幅 1/n）: 局在せず、d_eff → 1.5（非整数）
  「振幅減衰は局在と整数次元の両方を壊す」ことの対照例として保持
"""
import math

PI = math.pi

# ---------- 主モデル：等振幅局在波 ----------

def dirichlet_sq(theta, N):
    s = math.sin(theta)
    if abs(s) < 1e-12:
        return float(N * N)
    return (math.sin(N * theta) / s) ** 2

def leakage(N, samples=200000):
    """主ローブ [-π/N, π/N] の外に漏れる強度割合（θ∈[-π/2, π/2] を1周期分として数値積分）
    注: |W_N|² は π 周期（奇数倍音のみ）なので基本領域を [-π/2, π/2] に取る"""
    a = PI / N          # 第一零点
    half = PI / 2
    def integrate(lo, hi):
        n = samples
        h = (hi - lo) / n
        s = 0.5 * (dirichlet_sq(lo + 1e-12, N) + dirichlet_sq(hi - 1e-12, N))
        for i in range(1, n):
            s += dirichlet_sq(lo + i * h, N)
        return s * h
    total = integrate(-half, half)
    main = integrate(-a, a) if a < half else total
    return 1.0 - main / total

def model_main(N):
    sigma_n = math.sqrt((N * N - 1) / 3.0) if N > 1 else 0.0
    theta0 = PI / N
    product = sigma_n * theta0
    fidelity = N / (N + 1)          # |⟨Ŵ_N|Ŵ_{N+1}⟩|²
    disting = 1.0 / (N + 1)         # 1 - fidelity ＝ 鋭さ
    L = leakage(N) if N >= 2 else 0.0
    return theta0, L, fidelity, disting, product

# ---------- 対照モデル：矩形波（1/n 減衰、局在しない）----------

TOTAL_A = PI * PI / 8.0

def model_contrast(N):
    S2 = sum(1.0 / (2 * k - 1) ** 2 for k in range(1, N + 1))
    S4 = sum(1.0 / (2 * k - 1) ** 4 for k in range(1, N + 1))
    delta = 1.0 - S2 / TOTAL_A
    d_eff = S2 * S2 / S4
    return delta, d_eff

def main():
    print("== 主モデル：等振幅局在波（局在する。規格化・次元は厳密、剰余は別の場所に住む）==")
    print(f"{'N':>4} {'θ₀=π/N':>9} {'漏れ率L(N)':>10} {'忠実度N/(N+1)':>13} {'識別1/(N+1)':>12} {'σ_n·θ₀':>9}")
    for N in [1, 2, 3, 4, 5, 7, 10, 20, 50, 100]:
        theta0, L, fid, dis, prod = model_main(N)
        print(f"{N:>4} {theta0:>9.4f} {L:>10.4f} {fid:>13.4f} {dis:>12.4f} {prod:>9.4f}")
    print(f"不確定性積の極限 π/√3 = {PI/math.sqrt(3):.6f}")
    print()
    print("== 対照モデル：矩形波（1/n 減衰。局在しない・エッジ波形）==")
    print("振幅減衰は局在と整数次元の両方を壊す、ことの対照例")
    print(f"{'N':>4} {'δ(N)':>10} {'d_eff':>8}")
    for N in [1, 3, 10, 100, 1000]:
        delta, d_eff = model_contrast(N)
        print(f"{N:>4} {delta:>10.6f} {d_eff:>8.4f}")
    print(f"d_eff 極限 = 96/64 = 1.5（非整数）、剰余係数 2/π² = {2/PI**2:.6f}")

if __name__ == "__main__":
    main()
