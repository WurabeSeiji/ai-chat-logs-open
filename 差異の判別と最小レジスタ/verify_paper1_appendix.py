#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
論文「共役複素対の二つの積と射影の階梯」付録A・Bの数値検証（同梱プログラム）

付録A（可視度定理）:
  (A-i)  |W_K(θ)|² のフーリエ係数が三角窓 K−|j| に一致すること
  (A-ii) 正規化核との数値畳み込みが縞 1+cos(2mθ) の可視度を 1−m/K にすること
付録B（隣接モード数の識別可能性）:
  (B)    |⟨Ŵ_K, Ŵ_{K+1}⟩|² = K/(K+1) の数値確認
補助（命題5.3 誘導測度）:
  (C)    θ 一様 → u=cosθ の経験分布が逆正弦則に従うこと（モーメント比較）
"""
import cmath
import math
import random

PI = math.pi

def W(theta, K):
    return sum(cmath.exp(1j * (2 * k - 1) * theta) for k in range(1, K + 1))

def check_A_coeffs(K=9, M=4096):
    print(f"(A-i) |W_{K}|² のフーリエ係数 vs 三角窓 K−|j|:")
    ok = True
    for jj in [0, 1, 3, 6, K - 1, K]:
        s = sum(abs(W(PI * i / M, K)) ** 2 * cmath.exp(-1j * 2 * jj * (PI * i / M))
                for i in range(M)) / M
        theory = max(K - abs(jj), 0)
        ok &= abs(s.real - theory) < 1e-6
        print(f"    j={jj}: 数値 {s.real:9.4f}  理論 {theory:3d}")
    print(f"    → {'一致' if ok else '不一致'}")

def check_A_visibility(cases=((9, 1), (9, 3), (9, 6), (5, 2), (17, 4)), M=2048):
    print("(A-ii) 畳み込み後の可視度 vs 1 − m/K:")
    for K, m in cases:
        F = [abs(W(PI * i / M, K)) ** 2 for i in range(M)]
        out_max, out_min = -1e18, 1e18
        for i in range(0, M, 16):
            acc = sum(F[j] * (1.0 + math.cos(2 * m * PI * (i - j) / M))
                      for j in range(0, M, 4))
            out_max, out_min = max(out_max, acc), min(out_min, acc)
        v = (out_max - out_min) / (out_max + out_min)
        print(f"    K={K:2d}, m={m}: 測定 {v:.4f}  理論 {1 - m / K:.4f}")

def check_B(Ks=(1, 2, 3, 5, 9, 17), M=4096):
    print("(B) |⟨Ŵ_K, Ŵ_{K+1}⟩|² vs K/(K+1):")
    for K in Ks:
        ip = sum(W(2 * PI * i / M, K).conjugate() * W(2 * PI * i / M, K + 1)
                 for i in range(M)) / M
        fid = abs(ip) ** 2 / (K * (K + 1))
        print(f"    K={K:2d}: 数値 {fid:.6f}  理論 {K / (K + 1):.6f}")

def check_C(n=200000):
    print("(C) 誘導測度（逆正弦則）: u=cosθ, θ一様。モーメント比較:")
    random.seed(1)
    m2 = m4 = 0.0
    for _ in range(n):
        u = math.cos(random.uniform(0, 2 * PI))
        m2 += u * u
        m4 += u ** 4
    print(f"    E[u²] 数値 {m2/n:.4f}  理論 1/2 = 0.5000")
    print(f"    E[u⁴] 数値 {m4/n:.4f}  理論 3/8 = 0.3750")

if __name__ == "__main__":
    check_A_coeffs()
    check_A_visibility()
    check_B()
    check_C()
