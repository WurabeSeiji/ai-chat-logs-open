#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""同じ弦の上の二つの複素パルスの、親のノルム R(Δ)・重なりの角 Θ・零閉塞を計算する確認用（2026-10-29）。

各パルス i：c^{(i)}_m = (2/π) sin(mφ_i/2)/m（m = 1..M_i、奇・偶で選ぶ）、直流 c₀ = φ_i/2π、中心 θ_i
和の係数   c_m = c^{(1)}_m e^{−imθ₁} + c^{(2)}_m e^{−imθ₂}
親のノルム R² = Σ|c_m|² = R₁² + R₂² + 2 Σ c^{(1)}_m c^{(2)}_m cos(mΔ)、Δ = θ₁ − θ₂
重なりの角 cos Θ = Σ c^{(1)}_m c^{(2)}_m cos(mΔ) / (R₁ R₂)
視野角     tan α = R₂ / R₁
零閉塞     Σ c_m²（基準 ψ = 0 から見た本当の二乗和）
"""
import math, cmath

def coeffs(phi_deg, M, odd=True, even=True):
    phi = math.radians(phi_deg)
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    return {m: (2 / math.pi) * math.sin(m * phi / 2) / m for m in ms}, phi / (2 * math.pi)

def pair(c1, c2, th1_deg, th2_deg):
    t1, t2 = math.radians(th1_deg), math.radians(th2_deg)
    ms = sorted(set(c1) | set(c2))
    c = {m: c1.get(m, 0) * cmath.exp(-1j * m * t1) + c2.get(m, 0) * cmath.exp(-1j * m * t2) for m in ms}
    R1 = math.sqrt(sum(v * v for v in c1.values())); R2 = math.sqrt(sum(v * v for v in c2.values()))
    R2tot = sum(abs(v) ** 2 for v in c.values()); R = math.sqrt(R2tot)
    cross = sum(c1.get(m, 0) * c2.get(m, 0) * math.cos(m * (t1 - t2)) for m in ms)
    cosT = cross / (R1 * R2); T = math.degrees(math.acos(max(-1, min(1, cosT))))
    S2 = sum(v * v for v in c.values())
    alpha = math.degrees(math.atan2(R2, R1))
    return dict(R1=R1, R2=R2, R=R, cross=cross, cosT=cosT, T=T, S2=S2, alpha=alpha, ratio=R / (R1 + R2))

def table(title, c1, c2, cases):
    print(title)
    print("  theta1, theta2, Delta | R1, R2, R, R/(R1+R2) | cross/R^2, cosT, Theta[deg] | S2 (ref psi=0) | alpha[deg]")
    for th1, th2 in cases:
        d = pair(c1, c2, th1, th2)
        print(f"  {th1:6.1f}, {th2:6.1f}, {th1-th2:6.1f} | {d['R1']:.4f}, {d['R2']:.4f}, {d['R']:.4f}, {d['ratio']:.4f} | "
              f"{d['cross']/d['R']**2:+.4f}, {d['cosT']:+.4f}, {d['T']:7.3f} | {d['S2'].real:+.4f}{d['S2'].imag:+.4f}i | {d['alpha']:.3f}")

if __name__ == "__main__":
    c1, _ = coeffs(10, 10); c2, _ = coeffs(10, 10)
    deltas = [0, 5, 10, 18, 36, 60, 90, 180]
    table("同じパルス二つ（φ=10°, M=10）、基準を真ん中（θ1 = 180+Δ/2, θ2 = 180−Δ/2）", c1, c2, [(180 + d / 2, 180 - d / 2) for d in deltas])
    table("同じパルス二つ、基準をパルス 1（θ1 = 180, θ2 = 180−Δ）", c1, c2, [(180, 180 - d) for d in deltas])
    c3, _ = coeffs(10, 100)
    table("違うパルス（M1=10, M2=100, φ=10°）、基準をパルス 1", c1, c3, [(180, 180 - d) for d in deltas])
    c4, _ = coeffs(10, 10, True, False); c5, _ = coeffs(10, 10, False, True)
    table("奇のみと偶のみ（φ=10°, M=10）、基準をパルス 1", c4, c5, [(180, 180 - d) for d in deltas])
