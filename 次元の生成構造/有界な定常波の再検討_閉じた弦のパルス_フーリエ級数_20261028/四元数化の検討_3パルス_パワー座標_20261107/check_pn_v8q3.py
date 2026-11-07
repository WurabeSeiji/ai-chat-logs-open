#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パワー座標の再計算（2026-11-07）。index_v8q3.html と同じ規則：各パルス P_i(ψ) = √|Q_i(ψ)|² = |c₀ + Σ_m c_m e^{im(ψ−θ_i)}|（内部位相は効かない）、
S(ψ) = Σ_i P_i(ψ) n_i。max_ψ |S| を返す。同じ格子 n = 1441（両端を含む）。"""
import math
def coeffs(phi, M): return [(m, (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m) for m in range(1, M + 1)], phi / 360
def unit(v):
    L = math.sqrt(sum(x * x for x in v)); return [0, 0, 1] if L < 1e-12 else [x / L for x in v]
def smax(pulses, n=1441):
    """pulses: [(theta, phi, M, n_vec)]"""
    best = 0.0
    for k in range(n):
        psi = 2 * math.pi * k / (n - 1); S = [0.0, 0.0, 0.0]
        for th, ph, M, nv in pulses:
            cs, c0 = coeffs(ph, M); x = psi - math.radians(th)
            re = c0 + sum(c * math.cos(m * x) for m, c in cs); im = sum(c * math.sin(m * x) for m, c in cs); P = math.hypot(re, im); u = unit(nv)
            for a in range(3): S[a] += P * u[a]
        best = max(best, math.sqrt(sum(x * x for x in S)))
    return best
K, I, J = [0, 0, 1], [1, 0, 0], [0, 1, 0]
CASES = {
  "default": [(180, 20, 10, K), (180, 10, 10, K), (180, 10, 10, K)],
  "tri": [(180, 10, 10, K), (150, 10, 10, K), (100, 10, 10, K)],
  "tri_ijk": [(180, 10, 10, I), (150, 10, 10, J), (100, 10, 10, K)],
  "mirror": [(180, 10, 10, K), (200, 10, 10, K), (160, 10, 10, K)],
  "mirror_k_i_mi": [(180, 10, 10, K), (200, 10, 10, I), (160, 10, 10, [-1, 0, 0])],
  "std": [(180, 20, 10, K), (180, 10, 10, K), (90, 10, 10, K)],
  "newdefault": [(97, 1, 100, K), (194, 1, 100, K), (58, 1, 100, K)],
  "newdefault_ijk": [(97, 1, 100, I), (194, 1, 100, J), (58, 1, 100, K)],
}
if __name__ == "__main__":
    for name, ps in CASES.items(): print(f"{name:16s} max |Σ P_i n_i| = {smax(ps):.4f}")
