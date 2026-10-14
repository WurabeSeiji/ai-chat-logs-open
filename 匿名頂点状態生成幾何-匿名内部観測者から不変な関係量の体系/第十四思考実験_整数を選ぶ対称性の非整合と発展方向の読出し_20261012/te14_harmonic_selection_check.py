# -*- coding: utf-8 -*-
"""
第十四思考実験 倍音の選択則の確認（2026-10-14）
23/124 にロックした周期124の軌道（四次元フルヴィッツ格子、Z4 軸、外部時計なし）について、
(1) i を掛けた状態が軌道上の何ステップ先にあるか（四分の一周期の対称性）
(2) z1 成分の周波数分解：エネルギーが m ≡ 23 (mod 4)、すなわち (23+4k)/124 の系列だけに入るか
(3) 側帯波 k ごとの平均パワー比
を調べる。
※ 同じフォルダに E2 の te14_E2_hurwitz_124_23.py と E5 の te14_E5_gauge_backreaction.py
  （E3 の te14_E3_ground_and_levels.py も）を置いて実行する。初期点：半径64で60点（乱数種 7）。
"""
import math, cmath, collections
from te14_E2_hurwitz_124_23 import qmul, round_h
from te14_E5_gauge_backreaction import rand_points
PHI = 2*math.pi*23/124; u = (math.cos(PHI), math.sin(PHI), 0, 0); P = 124

def cycle(p):
    seen = {}; path = []
    while p not in seen:
        seen[p] = len(path); path.append(p); p = round_h(qmul(u, tuple(v/2 for v in p)))[0]
    return path[seen[p]:]

def main():
    pts = rand_points(60, 64, 7)
    shifts = collections.Counter(); acc = collections.defaultdict(float); frac = []; n = 0
    for p0 in pts:
        c = cycle(p0)
        if len(c) != P: continue
        n += 1
        iq = tuple(int(v) for v in qmul((0, 1, 0, 0), c[0]))          # i*q（2倍座標のまま整数）
        shifts[c.index(iq) if iq in c else None] += 1
        z = [complex(x[0], x[1])/2 for x in c]
        spec = {m: abs(sum(z[k]*cmath.exp(-2j*math.pi*m*k/P) for k in range(P)))**2 for m in range(P)}
        tot = sum(spec.values())
        frac.append(sum(spec[m] for m in range(P) if m % 4 == 23 % 4)/tot)
        for kk in range(-15, 16): acc[kk] += spec[(23+4*kk) % P]/tot
    print('周期124の軌道数', n)
    print('i*q が現れる位置（ステップ）', dict(shifts))
    print('(23+4k)/124 系列に入るパワーの割合 最小', min(frac))
    for kk in range(-8, 9): print('k', kk, 'm', 23+4*kk, 'power', f"{acc[kk]/n:.3e}")

if __name__ == '__main__':
    main()
