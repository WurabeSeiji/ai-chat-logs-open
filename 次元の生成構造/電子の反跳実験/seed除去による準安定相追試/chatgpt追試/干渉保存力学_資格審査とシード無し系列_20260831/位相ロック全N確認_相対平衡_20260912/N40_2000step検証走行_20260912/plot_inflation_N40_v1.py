#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40・den=40・2000歩のインフレーション図（単図・読み出しのみ）。
N=64版 plot_inflation_N64_v1.py の忠実コピー（N/入出力名/タイトルのみ変更）。
様式は正本俯瞰図の各パネルと同一（semilogy、H⊥/H、ylim 1e-34..3）。
入力: results/timeseries_64bit_with124_N3_N40.csv（正本が書く固定名）
出力: fig_Hperp_N40_den40_2000.png"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
xs, ys = [], []
with open(os.path.join(BASE, 'results', 'timeseries_64bit_with124_N3_N40.csv')) as f:
    for r in csv.DictReader(f):
        if int(r['N']) == 40 and r['series'] == 'N':
            xs.append(int(r['step'])); ys.append(float(r['Hperp_frac']))
assert len(xs) == 2001, len(xs)
fig, ax = plt.subplots(figsize=(9, 6))
ax.semilogy(xs, ys, linewidth=1.0, label='N (2pi/40)')
ax.set_xlim(0, 2000); ax.set_ylim(1e-34, 3)
ax.set_xlabel('step'); ax.set_ylabel('Hperp/H')
ax.grid(alpha=.25); ax.legend(loc='lower right')
ax.set_title('N=40 (M=780) make_parent parent: Hperp/H, den=N, 2000 steps')
fig.tight_layout()
out = os.path.join(BASE, 'fig_Hperp_N40_den40_2000.png')
fig.savefig(out, dpi=180)
print('saved:', out)
