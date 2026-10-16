# -*- coding: utf-8 -*-
"""E4 追加（E4(A) の結果を見た後に設定）：分母が4の倍数の根と、そうでない根の自己ロック率（Z4 軸、半径64、各30点）"""
import math, random, json
from fractions import Fraction
from te14_E4_selection_pullin_recovery import orbit, rand_point, rho_of, AX
random.seed(20261015)
out = {}
for t in ['21/116', '19/104', '25/136', '27/148', '23/120', '21/110', '25/134', '19/102']:
    fr = Fraction(t); hit = n = 0
    for _ in range(30):
        r = orbit(rand_point(64), AX, 2*math.pi*float(fr))
        if r is None: continue
        rh = rho_of(r)
        if rh is None: continue
        n += 1; hit += (rh == fr)
    out[t] = (hit, n, fr.denominator % 4 == 0); print(t, out[t], flush=True)
json.dump(out, open('te14_E4_mod4_summary.json', 'w'), indent=1)
