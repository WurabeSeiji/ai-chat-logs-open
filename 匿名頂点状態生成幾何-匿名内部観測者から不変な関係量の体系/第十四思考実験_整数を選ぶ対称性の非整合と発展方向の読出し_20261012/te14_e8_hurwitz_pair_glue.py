# -*- coding: utf-8 -*-
"""
第十四思考実験 追記（続）：E8 を二つのフルヴィッツ四元数の接着として確認する（2026-10-15）
L = {(x, y) : x, y はフルヴィッツ整数, x - y ∈ (1+i)H}
(1) L のノルム2のベクトル（根）の数が 240、ノルム1が 0 であること（E8 格子）
(2) 相手 y = u x（u は24個のフルヴィッツ単位）が接着条件を満たすのはどの u か
座標は2倍した整数で扱う（半整数を厳密に扱うため）。
"""
import itertools

def hur(maxc):
    return [v for v in itertools.product(range(-maxc, maxc+1), repeat=4)
            if all(x % 2 == 0 for x in v) or all(x % 2 for x in v)]

def mul(p, q):   # 2倍座標どうしの積 → 2倍座標
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    r = (a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
         a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2)
    return tuple(x//2 for x in r)

def ish(v): return all(x % 2 == 0 for x in v) or all(x % 2 for x in v)
INV = (1, -1, 0, 0)            # (1-i)/2 の2倍座標
def in_ideal(v): return ish(mul(INV, v))   # v ∈ (1+i)H  ⇔  (1+i)^{-1} v ∈ H

H = hur(4)
small = [v for v in H if sum(x*x for x in v) <= 8]
n2 = sum(1 for x in small for y in small
         if sum(a*a for a in x) + sum(a*a for a in y) == 8 and in_ideal(tuple(a-b for a, b in zip(x, y))))
n1 = sum(1 for x in small for y in small
         if sum(a*a for a in x) + sum(a*a for a in y) == 4 and in_ideal(tuple(a-b for a, b in zip(x, y))))
print('ノルム2のベクトル数（E8 なら240）:', n2, ' ノルム1:', n1)

units = [v for v in H if sum(x*x for x in v) == 4]
for q in [(2, 4, -2, 6), (1, 3, 5, -1), (2, 0, 0, 0)]:
    ok = [u for u in units if in_ideal(tuple(a-b for a, b in zip(q, mul(u, q))))]
    print('x =', q, ' 接着条件を満たす相手 u·x の u（2倍座標）:', sorted(ok))
print('→ u は ±1, ±i, ±j, ±k の8個（Z4・180°の反対点を含む）。三相の16個（(±1±i±j±k)/2）は満たさない。')
