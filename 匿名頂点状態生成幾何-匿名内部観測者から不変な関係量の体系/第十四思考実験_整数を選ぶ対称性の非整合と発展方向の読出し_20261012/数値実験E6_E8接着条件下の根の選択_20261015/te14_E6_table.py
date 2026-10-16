# -*- coding: utf-8 -*-
"""E6 集計：free（接着なし）と glued（接着あり）の x のロック率を対象ごとに並べ、E4 part A の値も併記する。
出力：te14_E6_table.md、te14_E6_fig.png/.svg。走行中でも、ある行だけで表を作る。"""
import csv, math, os
from fractions import Fraction
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(mode):
    p = f'te14_E6_rows_{mode}.csv'
    if not os.path.exists(p): return {}
    return {(r['target'], int(r['R'])): (int(r['lock_x']), int(r['lock_y']), int(r['n'])) for r in csv.DictReader(open(p))}

def load_live(mode):
    p = f'run_E6_{mode}_stdout.log'
    out = {}
    if not os.path.exists(p): return out
    for line in open(p):
        parts = line.split()
        if len(parts) >= 4 and parts[0] == mode:
            t, R = parts[1], int(parts[2]); tup = line[line.index('('):line.index(')')+1]
            x, y, n = [int(v) for v in tup.strip('()').split(',')]
            out[(t, R)] = (x, y, n)
    return out

free = load('free') or load_live('free'); glued = load('glued') or load_live('glued')
e4 = {(r['key'], int(r['param'])): (int(r['v1']), int(r['v2'])) for r in csv.DictReader(open('control/te14_E4_rows.csv')) if r['part'] == 'A'}
groups = [('mod 4 の類', ['17/92', '19/104', '21/116', '23/124', '25/136', '27/148', '23/120']),
          ('窓の隣', ['23/122', '23/123', '23/125', '24/127', '22/119']),
          ('E4 の対照', ['23/126', '24/131', '18/97'])]
N_of = lambda t: 4*math.pi/math.sin(math.pi*float(Fraction(t)))**4
lines = ['| 群 | 根 | 分母 mod 4 | N(R) | E4（x 単独） R32/64/96 | free（対・接着なし） | glued（接着あり） |', '|---|---|---|---|---|---|---|']
def fmt(d, t):
    cells = []
    for R in (32, 64, 96):
        v = d.get((t, R))
        cells.append('—' if v is None else f"{v[0]}/{v[2]}")
    return ' '.join(cells)
for g, ts in groups:
    for t in ts:
        e4c = ' '.join('—' if e4.get((t, R)) is None else f"{e4[(t,R)][0]}/{e4[(t,R)][1]}" for R in (32, 64, 96))
        lines.append(f"| {g} | {t} | {Fraction(t).denominator % 4} | {N_of(t):.1f} | {e4c} | {fmt(free, t)} | {fmt(glued, t)} |")
open('te14_E6_table.md', 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
# figure: lock rate at R=64 and R=96, free vs glued
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
for ax, R in zip(axes, (64, 96)):
    labels = [t for _, ts in groups for t in ts]
    fr = [free.get((t, R), (0, 0, 0)) for t in labels]; gl = [glued.get((t, R), (0, 0, 0)) for t in labels]
    xs = range(len(labels)); w = 0.38
    ax.bar([x-w/2 for x in xs], [v[0]/v[2] if v[2] else 0 for v in fr], w, label='free (pair, no gluing)', color='tab:gray')
    ax.bar([x+w/2 for x in xs], [v[0]/v[2] if v[2] else 0 for v in gl], w, label='glued (E8)', color='tab:red')
    ax.set_xticks(list(xs)); ax.set_xticklabels(labels, rotation=60, fontsize=8); ax.set_ylim(0, 1.05)
    ax.set_ylabel('lock rate of x (rho_top = target)'); ax.set_title(f'R = {R}'); ax.axvline(6.5, color='0.7', lw=1); ax.axvline(11.5, color='0.7', lw=1)
    ax.legend(fontsize=8)
fig.suptitle('E6: self-lock rate on the Hurwitz lattice, with and without the E8 gluing condition')
fig.savefig('te14_E6_fig.png', dpi=150); fig.savefig('te14_E6_fig.svg'); print('figure written')
