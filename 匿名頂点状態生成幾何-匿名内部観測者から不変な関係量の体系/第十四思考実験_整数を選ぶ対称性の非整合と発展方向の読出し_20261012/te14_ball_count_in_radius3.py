# -*- coding: utf-8 -*-
"""
第十四思考実験 追記 第21.2節：半径3の四次元球に収まる四次元球の数（2026-10-14）
使い方: python te14_ball_count_in_radius3.py [小球の半径 r]（既定 0.5）
  r=0.5（直径1）、r=(2/pi^2)^(1/4)=0.6709383（一辺1の立方体と同じ体積）など。
出力: Z4 格子（原点に格子点／原点に穴）と D4 格子（最短距離1）について
  (完全に収まる個数, 各小球の内側体積割合で重み付けした個数)
小球の内側体積割合は、軸に沿った切り口が同心の三次元球になることを使って数値積分する。
"""
import itertools, math, sys
from functools import lru_cache
R = 3.0; r = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5

@lru_cache(None)
def frac(d2):
    d = math.sqrt(d2)
    if d + r <= R: return 1.0
    if d - r >= R: return 0.0
    n = 4000; tot = 0
    for k in range(n):
        x = -r + (k+0.5)*2*r/n
        s2 = r*r - x*x; S2 = R*R - (d+x)**2
        m = min(s2, S2)
        if m > 0: tot += (4*math.pi/3)*m**1.5
    tot *= 2*r/n
    return tot/(math.pi**2/2*r**4)

def count(points):
    full = sum(1 for p in points if math.sqrt(p) + r <= R + 1e-12)
    w = sum(frac(round(p, 10)) for p in points)
    return full, w

Z = [sum(c*c for c in v) for v in itertools.product(range(-4, 5), repeat=4)]
print('Z4 原点に格子点: 完全に収まる個数, 重み付き', count(Z))
Zh = [sum((c+.5)**2 for c in v) for v in itertools.product(range(-4, 4), repeat=4)]
print('Z4 原点に穴:', count(Zh))
D = [sum(c*c for c in v)/2 for v in itertools.product(range(-6, 7), repeat=4)
     if sum(v) % 2 == 0 and sum(c*c for c in v)/2 <= 12.5]
print('D4（最短距離1）原点に格子点:', count(D))
print('137 になる半径の範囲:', 3-math.sqrt(6), '< r <=', 3-math.sqrt(5))
