#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パワー座標の再計算（2026-11-07）。index_v8c3.html と同じ規則：各パルス P_i(ψ) = √|Z_i(ψ)|²（内部位相は効かない）、S(ψ) = Σ_i P_i(ψ)。max_ψ S を返す。同じ格子 n = 1441"""
import math
def coeffs(phi, M): return [(m, (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m) for m in range(1, M + 1)], phi / 360
def smax(pulses, n=1441):
    best = 0.0
    for k in range(n):
        psi = 2 * math.pi * k / (n - 1); S = 0.0
        for th, ph, M in pulses:
            cs, c0 = coeffs(ph, M); x = psi - math.radians(th)
            re = c0 + sum(c * math.cos(m * x) for m, c in cs); im = sum(c * math.sin(m * x) for m, c in cs); S += math.hypot(re, im)
        best = max(best, S)
    return best
CASES = {
  "default": [(180, 20, 10), (180, 10, 10), (180, 10, 10)],
  "tri": [(180, 10, 10), (150, 10, 10), (100, 10, 10)],
  "mirror": [(180, 10, 10), (200, 10, 10), (160, 10, 10)],
  "std": [(180, 20, 10), (180, 10, 10), (90, 10, 10)],
  "newdefault": [(97, 1, 100), (194, 1, 100), (58, 1, 100)],
}
if __name__ == "__main__":
    for name, ps in CASES.items(): print(f"{name:12s} max Σ P_i = {smax(ps):.4f}")
