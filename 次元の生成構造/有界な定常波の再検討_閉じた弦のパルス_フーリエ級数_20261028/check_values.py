#!/usr/bin/env python3
"""index.html と同じ式を独立に計算して max z / min z / z(θ) を出す（確認用）。

    z_M(ψ) = φ/360 + (2/π) Σ_{m=1}^{M} [ sin(mπφ/360) / m ] · cos( mπ(ψ−θ)/180 )

ψ は 0.25° 刻み 1441 点（index.html と同じ）。
"""
import math

def z(psi, theta, phi, M):
    s = phi / 360
    for m in range(1, M + 1):
        s += (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m * math.cos(m * math.pi * (psi - theta) / 180)
    return s

if __name__ == "__main__":
    theta, phi = 180.0, 10.0
    print(f"theta = {theta}, phi = {phi}")
    print("M, max z, psi at max, min z, z(theta), z(0)")
    for M in (1, 10, 50, 100):
        vals = [z(i * 0.25, theta, phi, M) for i in range(1441)]
        imax = vals.index(max(vals))
        print(f"{M}, {max(vals):.4f}, {imax * 0.25}, {min(vals):.4f}, {z(theta, theta, phi, M):.4f}, {z(0, theta, phi, M):.4f}")
