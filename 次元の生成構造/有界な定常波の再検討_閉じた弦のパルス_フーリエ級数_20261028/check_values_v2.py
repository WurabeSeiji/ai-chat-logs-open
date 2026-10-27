#!/usr/bin/env python3
"""index_v2.html と同じ式を独立に計算して max z / min z / z(θ) を出す（確認用）。

パルス  z_M(ψ) = φ/360 + (2/π) Σ_{m=1}^{M} [ sin(mπφ/360) / m ] · cos( mπ(ψ−θ)/180 )
矩形波  s_M(ψ) = (4/π) Σ_{k=0}^{M−1} [ (−1)^k / (2k+1) ] · cos( (2k+1) N π(ψ−θ)/180 ),  N = round(180/φ)

ψ の刻みは index_v2.html と同じ（最高波数の 1 周期あたり 8 点以上、最低 1441 点）。
"""
import math

def pulse(psi, theta, phi, M):
    s = phi / 360
    for m in range(1, M + 1):
        s += (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m * math.cos(m * math.pi * (psi - theta) / 180)
    return s

def square(psi, theta, N, M):
    x = N * math.pi * (psi - theta) / 180
    return (4 / math.pi) * sum((-1) ** k / (2 * k + 1) * math.cos((2 * k + 1) * x) for k in range(M))

if __name__ == "__main__":
    theta, phi = 180.0, 10.0
    N = max(1, round(180 / phi))
    print(f"theta = {theta}, phi = {phi}, N = {N}")
    print("mode, M, kmax, n_points, max z, min z, z(theta), z(theta+phi)")
    for M in (1, 10, 100):
        for mode in ("pulse", "square"):
            kmax = (2 * M - 1) * N if mode == "square" else M
            n = min(40001, max(1441, 8 * kmax + 1))
            f = (lambda p: pulse(p, theta, phi, M)) if mode == "pulse" else (lambda p: square(p, theta, N, M))
            vals = [f(i * 360 / (n - 1)) for i in range(n)]
            print(f"{mode}, {M}, {kmax}, {n}, {max(vals):.4f}, {min(vals):.4f}, {f(theta):.4f}, {f(theta + phi):.4f}")
