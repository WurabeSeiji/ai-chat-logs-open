# -*- coding: utf-8 -*-
"""E5 対照：歪みなし（s=0）を同じ初期点で計算し、歪みの有無で軌道が変わったかを比べる"""
import json
from te14_E5_gauge_backreaction import run, rand_points, xread
pts = rand_points(24, 64, 20261014)
a0 = 1/xread(23/124)
res = {s: [run(p, a0, s) for p in pts] for s in (0, +1, -1)}
diff = {s: sum(1 for r0, r in zip(res[0], res[s]) if (r0 and r) and (r0['rho'] != r['rho'] or r0['P'] != r['P'])) for s in (+1, -1)}
same_cycle = {s: sum(1 for r0, r in zip(res[0], res[s]) if (r0 and r) and r0['P'] == r['P'] and r0['rho'] == r['rho']) for s in (+1, -1)}
locked0 = sum(1 for r in res[0] if r and str(r['rho']) == '23/124')
out = dict(locked_s0=locked0, changed_vs_s0=diff, unchanged_vs_s0=same_cycle,
           rhos_s0=sorted({str(r['rho']) for r in res[0] if r}))
print(out); json.dump(out, open('te14_E5_control_summary.json', 'w'), indent=1)
