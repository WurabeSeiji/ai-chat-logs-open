#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
検証[4]: 中心投影の対蹠同一視 Z₂ と奇数倍音の反周期性

幾何的事実（論文1・radial projection ノートで確認済み）:
  チャート点 x ∈ Π_R ↔ 原点を通る直線 ↔ 球面上の対蹠対 {Φ(x), −Φ(x)}
  接点を通る大円上では、チャート座標 x = R·tanθ、対蹠 θ→θ+π は同一チャート点

検証する定理:
  (a) 反周期性 f(θ+π) = −f(θ) ⟺ フーリエ台が奇数倍音のみ（数値FFT検証）
  (b) チャート関数（π周期＝偶数セクター）と対蹠ビットを運ぶ波（反周期＝奇数セクター）の分解
  (c) 偶数核と奇数核はどちらも θ=0, π に双子ピークを持つが、
      偶数核は同符号（対蹠対に「同」）、奇数核は逆符号（対蹠対に「異」＝1ビット）
"""
import cmath
import math
import random

PI = math.pi
random.seed(20260705)

def fft_coeffs(samples):
    """素朴DFT（依存なし）: c_m, m = 0..M-1"""
    M = len(samples)
    return [sum(samples[j] * cmath.exp(-2j * PI * m * j / M) for j in range(M)) / M
            for m in range(M)]

def task_a():
    """反周期関数のフーリエ台が奇数のみか"""
    M = 256
    # ランダム関数を反周期射影: g(θ) = (f(θ) - f(θ+π))/2
    f = [random.uniform(-1, 1) + 1j * random.uniform(-1, 1) for _ in range(M)]
    g = [(f[j] - f[(j + M // 2) % M]) / 2 for j in range(M)]
    c = fft_coeffs(g)
    even_power = sum(abs(c[m]) ** 2 for m in range(0, M, 2))
    odd_power = sum(abs(c[m]) ** 2 for m in range(1, M, 2))
    print("[a] 反周期射影 g(θ)=(f(θ)−f(θ+π))/2 のスペクトル:")
    print(f"    偶数調和のパワー = {even_power:.3e}（ゼロであるべき）")
    print(f"    奇数調和のパワー = {odd_power:.3e}")
    # 逆向き: 奇数倍音だけの波は反周期か
    K = 9
    def W(theta):
        return sum(cmath.exp(1j * (2 * k - 1) * theta) for k in range(1, K + 1))
    max_dev = max(abs(W(th + PI) + W(th))
                  for th in [0.1 + 0.13 * i for i in range(40)])
    print(f"    奇数倍音波の反周期性 max|W(θ+π)+W(θ)| = {max_dev:.3e}（ゼロであるべき）")

def task_c():
    """偶数核・奇数核の双子ピークの符号"""
    K = 9
    def W_odd(theta):
        return sum(cmath.exp(1j * (2 * k - 1) * theta) for k in range(1, K + 1))
    def W_even(theta):
        return sum(cmath.exp(1j * (2 * k) * theta) for k in range(0, K))  # 0,2,4,...
    e = 1e-9
    print()
    print("[c] 双子ピーク（θ=0 と θ=π）の振幅:")
    for name, Wf in (("奇数核", W_odd), ("偶数核", W_even)):
        a0 = Wf(e)
        api = Wf(PI + e)
        print(f"    {name}: W(0)={a0.real:+.3f}{a0.imag:+.3f}i, "
              f"W(π)={api.real:+.3f}{api.imag:+.3f}i, "
              f"強度比 |W(π)|²/|W(0)|² = {abs(api)**2/abs(a0)**2:.6f}")
    print("    → 強度はどちらも双子（対蹠で同値）。振幅は奇数核のみ反転（対蹠に「異」を書く）")

def task_b():
    """チャート同一視: 大円上 θ と θ+π が同一チャート点か（x = R tanθ）"""
    print()
    print("[b] チャート座標の対蹠同一視: x(θ) = R·tanθ")
    R = 1.0
    devs = []
    for i in range(1, 20):
        th = -PI / 2 + PI * i / 20
        if abs(math.cos(th)) < 1e-9:
            continue
        devs.append(abs(R * math.tan(th) - R * math.tan(th + PI)))
    print(f"    max|x(θ) − x(θ+π)| = {max(devs):.3e}（ゼロ＝同一チャート点）")

if __name__ == "__main__":
    task_a()
    task_b()
    task_c()
