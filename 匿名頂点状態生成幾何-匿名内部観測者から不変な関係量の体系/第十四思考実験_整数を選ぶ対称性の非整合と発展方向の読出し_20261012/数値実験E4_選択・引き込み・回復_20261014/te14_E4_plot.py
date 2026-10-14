# -*- coding: utf-8 -*-
"""第十四思考実験 E4 図化：te14_E4_summary.json, te14_E4_mod4_summary.json から te14_E4_fig.png / .svg"""
import json, math
from fractions import Fraction
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
S = json.load(open('te14_E4_summary.json')); M = json.load(open('te14_E4_mod4_summary.json'))
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
# (a) self-lock rate at R=64
lab = []; val = []; col = []
for k, (h, n) in S['A'].items():
    t, R = k.split('@R')
    if R == '64': lab.append(t); val.append(h/n); col.append('C0' if Fraction(t).denominator % 4 == 0 else 'C7')
for t, (h, n, m4) in M.items():
    lab.append(t); val.append(h/n); col.append('C0' if m4 else 'C7')
order = sorted(range(len(lab)), key=lambda i: Fraction(lab[i]))
ax[0].bar(range(len(lab)), [val[i] for i in order], color=[col[i] for i in order])
ax[0].set_xticks(range(len(lab))); ax[0].set_xticklabels([lab[i] for i in order], rotation=60, fontsize=8)
ax[0].set_ylabel('self-lock fraction (R=64)'); ax[0].set_title('(a) blue: denominator divisible by 4')
# (b) pull-in
for base, c in (('23/124', 'C0'), ('23/122', 'C3')):
    pts = []
    for k, (h, n) in S['B'].items():
        if k.startswith(base): pts.append((float(k[len(base):]), h/n))
    pts.sort(); ax[1].plot([p[0] for p in pts], [p[1] for p in pts], 'o-', color=c, label=base)
ax[1].set_xscale('symlog', linthresh=1e-6); ax[1].set_xlabel('eps (turns/step)')
ax[1].set_ylabel('fraction locked to base'); ax[1].set_title('(b) pull-in window'); ax[1].legend()
# (c) recovery
C = S['C']; X = 4*math.pi/math.sin(math.pi*23/124)**4
ax[2].scatter([c['transient'] for c in C], [c['tr_alpha_inv'] for c in C],
              c=['C0' if c['relock'] == '23/124' else 'C3' for c in C])
ax[2].axhline(X, color='k', lw=.8, label='x_23/124'); ax[2].axhline(137.035999177, color='gray', ls='--', lw=.8, label='CODATA (reference)')
ax[2].set_xscale('log'); ax[2].set_xlabel('recovery length (updates)'); ax[2].set_ylabel('alpha_read^-1 during transient')
ax[2].set_title('(c) recovery after 1-unit kick'); ax[2].legend(fontsize=8)
fig.suptitle('TE14-E4: selection, pull-in and recovery on the 4D Hurwitz lattice (Z4 axis)')
fig.tight_layout()
for ext in ('png', 'svg'): fig.savefig(f'te14_E4_fig.{ext}', dpi=150)
