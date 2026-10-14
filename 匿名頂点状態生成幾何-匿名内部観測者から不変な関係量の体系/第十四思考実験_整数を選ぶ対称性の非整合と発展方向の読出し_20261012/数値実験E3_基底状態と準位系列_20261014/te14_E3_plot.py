# -*- coding: utf-8 -*-
"""第十四思考実験 E3 図化：te14_E3_orbits.csv から te14_E3_fig.png / .svg"""
import csv
from fractions import Fraction
from collections import Counter
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows = [r for r in csv.DictReader(open('te14_E3_orbits.csv')) if r['rho_top'] != 'W_inconsistent']
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for axn, col, off in (('Z4', 'C0', -0.00004), ('Z6', 'C3', 0.00004)):
    c = Counter(r['rho_top'] for r in rows if r['axis'] == axn)
    xs = [float(Fraction(k)) for k in c]; ys = list(c.values())
    ax[0].bar([x+off for x in xs], ys, width=0.00008, color=col, label=f'axis {axn}')
for k in ('3/16', '5/27', '13/70', '18/97', '23/124', '29/156'):
    v = float(Fraction(k)); ax[0].axvline(v, color='gray', lw=.5, ls=':'); ax[0].text(v, ax[0].get_ylim()[1]*0.95, k, fontsize=7, rotation=90, va='top')
ax[0].set_xlim(0.1845, 0.1880); ax[0].set_xlabel('topological rotation number rho_top = W/P')
ax[0].set_ylabel('number of orbits'); ax[0].set_title('(a) occupied levels'); ax[0].legend(fontsize=8)
radii = sorted({int(r['R']) for r in rows})
for axn, col in (('Z4', 'C0'), ('Z6', 'C3')):
    fr = []; fn = []
    for R in radii:
        s = [r for r in rows if r['axis'] == axn and int(r['R']) == R]
        fr.append(sum(1 for r in s if r['rho_top'] == '23/124')/len(s))
        fn.append(sum(1 for r in s if r['rho_top'] == '23/124' and r['norm_changes'] == '0')/len(s))
    ax[1].plot(radii, fr, 'o-', color=col, label=f'{axn}: rho_top = 23/124')
    ax[1].plot(radii, fn, 'x--', color=col, label=f'{axn}: 23/124 and norm constant')
ax[1].set_xscale('log'); ax[1].set_xlabel('initial radius'); ax[1].set_ylabel('fraction')
ax[1].set_title('(b) 23/124 lock vs complete lock'); ax[1].legend(fontsize=7)
fig.suptitle('TE14-E3: levels on the Hurwitz lattice with R_{124,23} rotation (topological rotation number)')
fig.tight_layout()
for ext in ('png', 'svg'): fig.savefig(f'te14_E3_fig.{ext}', dpi=150)
