# -*- coding: utf-8 -*-
"""
第十四思考実験 E5 図化（2026-10-15 追加）
入力 : te14_E5_rows.csv（te14_E5_gauge_backreaction.py の出力）
       te14_E5_control_summary.json（te14_E5_control.py の出力）
出力 : te14_E5_fig.png, te14_E5_fig.svg
(a) 23/124 にロックした軌道の数：歪みなし（s=0）、伸張（s=+1）、収縮（s=-1）
(b) 各軌道の位相幾何的回転数 rho と 23/124 の差（反復0回目）
(c) 歪みなしと比べて周期軌道が変わった軌道の数
実行順：python te14_E5_gauge_backreaction.py → python te14_E5_control.py → python te14_E5_plot.py
"""
import csv, json
from fractions import Fraction
import matplotlib; matplotlib.use('Agg'); matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt

rows = [r for r in csv.DictReader(open('te14_E5_rows.csv')) if r['iter'] == '0']
ctl = json.load(open('te14_E5_control_summary.json'))
T = Fraction(23, 124)
lock = {'0': ctl['locked_s0']}
for s in ('1', '-1'):
    lock[s] = sum(1 for r in rows if r['s'] == s and r['rho'] == '23/124')
n = sum(1 for r in rows if r['s'] == '1')

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
labels = ['s=0\n(no distortion)', 's=+1\n(stretch)', 's=-1\n(shrink)']
ax[0].bar(labels, [lock['0'], lock['1'], lock['-1']], color=['C7', 'C0', 'C3'])
for i, v in enumerate([lock['0'], lock['1'], lock['-1']]):
    ax[0].text(i, v + 0.3, f'{v}/{n}', ha='center')
ax[0].set_ylabel('orbits locked to 23/124'); ax[0].set_title('(a) locking with self-distortion alpha/r')

for s, c, dx in (('1', 'C0', -0.08), ('-1', 'C3', 0.08)):
    d = [float(Fraction(r['rho']) - T) for r in rows if r['s'] == s]
    ax[1].scatter([dx + (0 if s == '1' else 1)] * len(d), d, color=c, alpha=.6)
ax[1].axhline(0, color='k', lw=.8)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(['s=+1', 's=-1'])
ax[1].set_ylabel('rho_top - 23/124 (turns/step)'); ax[1].set_title('(b) rotation number of each orbit (iter 0)')

ch = [ctl['changed_vs_s0']['1'], ctl['changed_vs_s0']['-1']]
un = [ctl['unchanged_vs_s0']['1'], ctl['unchanged_vs_s0']['-1']]
ax[2].bar(['s=+1', 's=-1'], un, color='C7', label='same cycle as s=0')
ax[2].bar(['s=+1', 's=-1'], ch, bottom=un, color='C2', label='different cycle')
ax[2].set_ylabel('orbits'); ax[2].set_title('(c) effect of distortion vs s=0'); ax[2].legend(fontsize=8)

fig.suptitle('TE14-E5: gauge self-distortion x <- x(1 + s*alpha/r) on the 4D Hurwitz lattice (Z4 axis, R=64, 24 orbits)')
fig.tight_layout()
for ext in ('png', 'svg'): fig.savefig(f'te14_E5_fig.{ext}', dpi=150)
