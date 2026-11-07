#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""右パネルのクリップ範囲の再計算（2026-11-07）。index_v7c3.html と同じ規則：パルス i の |Z_i|² を同じ格子（N = n − 1 = 1440 点）で取り、
中心 θ_i の格子点から両側へ、|Z_i|² ≥ 最大の 5% である限り広げる。返り値は θ_i からの相対角（度）。内部位相は |Z_i|² に効かない。"""
import math
def coeffs(phi, M=10): return [(m, (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m) for m in range(1, M + 1)], phi / 360
def clip(theta, phi, n=1441, thr=0.05):
    N = n - 1; cs, c0 = coeffs(phi); P = []
    for k in range(N):
        x = 2 * math.pi * k / N - math.radians(theta)
        re = c0 + sum(c * math.cos(m * x) for m, c in cs); im = sum(c * math.sin(m * x) for m, c in cs); P.append(re * re + im * im)
    at = lambda k: k % N; k0 = at(round(theta / 360 * N)); pk = max(P)
    a = b = k0
    while k0 - a < N / 2 and P[at(a - 1)] >= thr * pk: a -= 1
    while b - k0 < N / 2 and P[at(b + 1)] >= thr * pk: b += 1
    w = lambda d: ((d + 180) % 360 + 360) % 360 - 180
    return w(at(a) * 360 / N - theta), w(at(b) * 360 / N - theta)
if __name__ == "__main__":
    for phi in (10, 20, 40): print(f"φ = {phi}: 山の範囲 {clip(180, phi)[0]:.1f}°〜+{clip(180, phi)[1]:.1f}°")
