#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""右パネル「線分」の山の再計算（2026-11-07）。index_v6c3.html と同じ規則を Python で独立に実装し、DOM の値と突き合わせる。
  規則：和のパワー P(ψ) = |Z(ψ)|² を同じ格子（n = max(1441, 8 kmax + 1)、両端を除く N = n − 1 点）で取り、極大（P ≥ 最大の 1 割）を候補に、
  高い順に、既に選んだ山の半値全幅（P > 頂点/2 の連続区間）の中にある極大は捨てる。線分 = 原点→(P cos ψ*, P sin ψ*)、
  円の直径 = 半値全幅の両端（極座標の曲線上の 2 点）の距離。
  単独パルス（M = 10）の Gibbs の副山は主山の 4.3% 以下なので閾値 1 割で落ちる。"""
import math

def coeffs(phi, M=10):
    return [(m, (2 / math.pi) * math.sin(m * math.pi * phi / 360) / m) for m in range(1, M + 1)], phi / 360

def lobes(pulses, n=1441, thr=0.1):
    """pulses: [(theta, phi, phase0)]。返り値 [(psi*, P, width_deg, diameter)]"""
    N = n - 1; P = [0.0] * N
    for k in range(N):
        psi = 2 * math.pi * k / N; re = im = 0.0
        for th, ph, p0 in pulses:
            cs, c0 = coeffs(ph); f = math.radians(p0); x = psi - math.radians(th)
            a = c0 + sum(c * math.cos(m * x) for m, c in cs); b = sum(c * math.sin(m * x) for m, c in cs)
            re += a * math.cos(f) - b * math.sin(f); im += a * math.sin(f) + b * math.cos(f)
        P[k] = re * re + im * im
    top = max(P); at = lambda i: i % N
    cand = []
    for i in range(N):
        if P[i] >= P[at(i - 1)] and P[i] > P[at(i + 1)] and P[i] >= thr * top:
            a = b = i
            while P[at(a - 1)] > P[i] / 2 and i - a < N: a -= 1
            while P[at(b + 1)] > P[i] / 2 and b - i < N: b += 1
            cand.append((i, a, b))
    cand.sort(key=lambda c: -P[c[0]])
    keep = []
    for c in cand:
        if not any(((c[0] - k[1]) % N) <= k[2] - k[1] for k in keep): keep.append(c)
    out = []
    for i, a, b in keep:
        Pk = P[i]; ps = i * 360 / N; ia, ib = at(a), at(b)
        pa = (P[ia] * math.cos(2 * math.pi * ia / N), P[ia] * math.sin(2 * math.pi * ia / N)); pb = (P[ib] * math.cos(2 * math.pi * ib / N), P[ib] * math.sin(2 * math.pi * ib / N))
        out.append((ps, Pk, (b - a) * 360 / N, math.hypot(pa[0] - pb[0], pa[1] - pb[1])))
    return out

CASES = {
  "default": [(180, 20, 0), (180, 10, 0), (180, 10, 0)],
  "tri":     [(180, 10, 0), (150, 10, 0), (100, 10, 0)],
  "mirror":  [(180, 10, 0), (200, 10, 45), (160, 10, -45)],
  "std":     [(180, 20, 0), (180, 10, 0), (90, 10, 0)],
}
if __name__ == "__main__":
    for name, ps in CASES.items():
        print(f"== {name} ==")
        for l in lobes(ps): print(f"  ψ* {l[0]:.1f}°, P {l[1]:.4f}, 幅 {l[2]:.1f}°, 直径 {l[3]:.4f}")
