#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位相正規化版走行のインフレーション図（読み出しのみ・単独大判）。
入力: results_staticparent_stage2init/timeseries_64bit_with124_N3_N40.csv の N=40 行。
様式は 7月 make_largeN_figure_v1.py と同じ semilog（縦軸 H⊥/H、横軸 step）。
6分母の曲線を1パネルに重ね描きする。"""
import csv
import os

import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, 'results_staticparent_stage2init', 'timeseries_64bit_with124_N3_N40.csv')
OUT = os.path.join(BASE, 'results_staticparent_stage2init',
                   'fig_inflation_N40_staticparent_stage2init.png')

curves = {}
with open(SRC) as fh:
    for row in csv.DictReader(fh):
        if row['N'] != '40':
            continue
        den = int(row['denominator'])
        curves.setdefault(den, ([], []))
        curves[den][0].append(float(row['step']))
        curves[den][1].append(float(row['Hperp_frac']))

fig, ax = plt.subplots(figsize=(8.2, 5.2))
for den in (38, 39, 40, 41, 42, 124):
    t, f = curves[den]
    label = f'2pi/{den}' if den == 124 else f'den={den} (N{den-40:+d})' if den != 40 else 'den=40 (N)'
    ax.semilogy(t, f, lw=1.4, label=label)
ax.set_xlabel('step')
ax.set_ylabel('Hperp/H (complement projection)')
ax.set_title('N=40 static July parent, stage2 at init only (z0 equimodularized), loop generator i*Im(H(z)), fixed dtau=2pi/den')
ax.grid(alpha=0.3, which='both')
ax.legend(loc='lower right', fontsize=8)
ax.set_ylim(1e-34, 2.0)
fig.tight_layout()
fig.savefig(OUT, dpi=160)
plt.close(fig)
print('wrote', OUT)
