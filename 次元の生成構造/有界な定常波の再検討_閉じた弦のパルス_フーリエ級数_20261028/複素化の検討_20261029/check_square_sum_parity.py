#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""単位振幅の複素倍音 c_m = e^{-imθ} について、ノルム Σ|c_m|² と本当の二乗和 Σ c_m² を
奇数倍音だけ・偶数倍音だけ・全倍音で、中心位相 θ を変えて計算する（2026-10-29 の検討）。

閉じた形：
  全  m=1..M      Σ e^{-2imθ} = e^{-i(M+1)θ} sin(Mθ)/sin θ
  奇  m=1,3,..,2K-1 Σ = e^{-2iKθ}     sin(2Kθ)/sin(2θ)
  偶  m=2,4,..,2K   Σ = e^{-2i(K+1)θ} sin(2Kθ)/sin(2θ)
"""
import cmath, math

def sums(ms, th):
    c = [cmath.exp(-1j * m * th) for m in ms]
    return sum(abs(x) ** 2 for x in c), sum(x * x for x in c)

def closed(kind, K, th):
    s2 = math.sin(2 * th)
    if kind == 'all':
        M = 2 * K
        return K * 2 if abs(math.sin(th)) < 1e-12 else cmath.exp(-1j * (M + 1) * th) * math.sin(M * th) / math.sin(th)
    if abs(s2) < 1e-12:
        return None
    if kind == 'odd':
        return cmath.exp(-2j * K * th) * math.sin(2 * K * th) / s2
    return cmath.exp(-2j * (K + 1) * th) * math.sin(2 * K * th) / s2

if __name__ == '__main__':
    K = 4
    sets = {'odd': [2 * k - 1 for k in range(1, K + 1)], 'even': [2 * k for k in range(1, K + 1)], 'all': list(range(1, 2 * K + 1))}
    print(f"K = {K}（奇 {sets['odd']}、偶 {sets['even']}、全 1..{2*K}）")
    print("kind, theta[deg], norm, Re sum c^2, Im sum c^2, |sum|, closed form check")
    for kind, ms in sets.items():
        for deg in (0, 30, 45, 60, 90, 120, 135, 150, 180):
            th = math.radians(deg)
            n, s = sums(ms, th)
            cf = closed(kind, K, th)
            chk = '' if cf is None else f"{abs(cf - s):.1e}"
            print(f"{kind}, {deg}, {n:.0f}, {s.real:+.4f}, {s.imag:+.4f}, {abs(s):.4f}, {chk}")
    print()
    print("M 奇数（1..7、奇 4 本・偶 3 本）の θ=90°：", sums(list(range(1, 8)), math.pi / 2)[1])
