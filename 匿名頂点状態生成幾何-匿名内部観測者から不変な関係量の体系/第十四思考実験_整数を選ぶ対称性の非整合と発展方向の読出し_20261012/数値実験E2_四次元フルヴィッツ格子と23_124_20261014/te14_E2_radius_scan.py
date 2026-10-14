# -*- coding: utf-8 -*-
"""E2 補助：半径依存（初期半径を変えて同じ更新則を適用）。乱数種固定。観測値は入力しない。"""
import math, random, csv, json
from te14_E2_hurwitz_124_23 import run, axis, alpha_read_inv, PHI_TEST
random.seed(20261014)
RADII = [3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
NPER = 12
rows = []
for axn in ('Z4', 'Z6'):
    ax = axis(axn)
    for R in RADII:
        for _ in range(NPER):
            v = [random.gauss(0, 1) for _ in range(4)]; s = math.sqrt(sum(x*x for x in v))
            half = random.random() < 0.5
            q = [R*x/s for x in v]
            p0 = tuple(2*math.floor(x)+1 for x in q) if half else tuple(2*round(x) for x in q)
            r = run(p0, ax, PHI_TEST)
            if r is None: rows.append([axn, R, *p0, 'NA']); continue
            rho = r['ang']/(2*math.pi*r['P'])
            rows.append([axn, R, *p0, r['P'], f"{rho:.15f}", f"{rho-23/124:+.6e}",
                         f"{alpha_read_inv(rho):.9f}", r['n_changes'], r['nmin'], r['nmax']])
        sub = [x for x in rows if x[0] == axn and x[1] == R and x[6] != 'NA']
        print(axn, R, 'P:', sorted({x[6] for x in sub})[:8], 'rho-23/124:', sorted({x[8] for x in sub})[:4], flush=True)
with open('te14_E2_radius_scan.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['axis','R','a2','b2','c2','d2','P','rho','rho_minus_23_124','alpha_read_inv','norm_changes','N_min','N_max'])
    w.writerows(rows)
