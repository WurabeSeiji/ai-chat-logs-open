#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
検証バッチ1（2026-07-05）
[6] 結晶学的制限の全数検証：整数2×2行列の有限位数は {1,2,3,4,6} のみか
[8] 小N不確定性積欠損の頑健性：幅の測度（第一零点/FWHM/RMS）依存性
"""
import math

PI = math.pi

# ---------- [6] 結晶学的制限 ----------

def mat_mul(A, B):
    return (
        (A[0][0]*B[0][0] + A[0][1]*B[1][0], A[0][0]*B[0][1] + A[0][1]*B[1][1]),
        (A[1][0]*B[0][0] + A[1][1]*B[1][0], A[1][0]*B[0][1] + A[1][1]*B[1][1]),
    )

I2 = ((1, 0), (0, 1))

def order_of(M, max_order=100):
    P = M
    for k in range(1, max_order + 1):
        if P == I2:
            return k
        P = mat_mul(P, M)
    return None  # 有限位数でない（max_order まで）

def task6(entry_range=3):
    orders_found = {}
    examples = {}
    r = range(-entry_range, entry_range + 1)
    count_checked = 0
    for a in r:
        for b in r:
            for c in r:
                for d in r:
                    det = a * d - b * c
                    if det not in (1, -1):
                        continue
                    count_checked += 1
                    M = ((a, b), (c, d))
                    k = order_of(M)
                    if k:
                        orders_found[k] = orders_found.get(k, 0) + 1
                        if k not in examples:
                            examples[k] = M
    print(f"[6] 結晶学的制限の全数検証（成分∈[-{entry_range},{entry_range}], det=±1, {count_checked}個検査）")
    print(f"    発見された有限位数の集合: {sorted(orders_found.keys())}")
    for k in sorted(orders_found):
        print(f"    位数{k}: {orders_found[k]}個  例 {examples[k]}")
    # 名指しの検算
    D90 = ((0, 1), (-1, 0))      # (a,b)->(b,-a)
    E1 = ((0, -1), (1, -1))      # (a,b)->(-b, a-b)
    E2 = ((1, -1), (1, 0))       # (a,b)->(a-b, a)
    print(f"    検算: 90°段 (a,b)->(b,-a) の位数 = {order_of(D90)} （Z₄）")
    print(f"    検算: (a,b)->(-b,a-b) の位数 = {order_of(E1)} （＝120°、Z₃）")
    print(f"    検算: (a,b)->(a-b,a) の位数 = {order_of(E2)} （＝60°、Z₆）")

# ---------- [8] 幅測度の頑健性 ----------

def kernel(theta, N):
    s = math.sin(theta)
    if abs(s) < 1e-14:
        return float(N * N)
    return (math.sin(N * theta) / s) ** 2

def fwhm(N):
    """半値全幅：K(θ)=N²/2 の解（主ローブ内、二分法）×2"""
    lo, hi = 0.0, PI / N
    target = N * N / 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if kernel(mid, N) > target:
            lo = mid
        else:
            hi = mid
    return 2.0 * lo

def rms_width(N, samples=400000):
    """RMS幅：基本領域 [-π/2, π/2] 上の √(∫θ²K/∫K)"""
    half = PI / 2
    h = 2 * half / samples
    s0 = s2 = 0.0
    for i in range(samples + 1):
        th = -half + i * h
        w = 1.0 if 0 < i < samples else 0.5
        K = kernel(th, N)
        s0 += w * K
        s2 += w * th * th * K
    return math.sqrt(s2 / s0)

def task8():
    print()
    print("[8] 小N不確定性積の頑健性（幅測度3種）")
    limit_zero = PI / math.sqrt(3)
    limit_fwhm = 2 * 1.391557 / math.sqrt(3)  # sinc²の半値半幅1.391557
    print(f"    極限値: 第一零点測度 π/√3={limit_zero:.4f}, FWHM測度≈{limit_fwhm:.4f}, RMS測度=収束せず（下表参照）")
    print(f"    {'N':>3} {'σ_n·θ₀':>9} {'欠損%':>7} {'σ_n·FWHM':>10} {'欠損%':>7} {'σ_n·σθ':>9}")
    ref_fwhm = None
    for N in [2, 3, 4, 5, 7, 10, 20, 50, 100]:
        sigma_n = math.sqrt((N * N - 1) / 3.0)
        p0 = sigma_n * (PI / N)
        p1 = sigma_n * fwhm(N)
        p2 = sigma_n * rms_width(N, samples=200000)
        d0 = 100 * (p0 / limit_zero - 1)
        d1 = 100 * (p1 / limit_fwhm - 1)
        print(f"    {N:>3} {p0:>9.4f} {d0:>+7.1f} {p1:>10.4f} {d1:>+7.1f} {p2:>9.4f}")

if __name__ == "__main__":
    task6()
    task8()
