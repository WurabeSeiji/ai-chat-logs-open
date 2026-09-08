#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R反射系（固定R124,23）N=3..16 の標準図化（読み出しのみ）。
1) インフレーション図: Test（本走行）と Control（既存理論床走行 ../対称親v2_500step走行_20260906/results）を重ね描き
2) step0 / step500 複素平面グリッド（正本 plot_complex_plane_N3_N40_stage123_v1.py と同一様式）"""
import csv
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'results')
CTRL = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'results')
NS = list(range(3, 17))


def series(path):
    out = {}
    for r in csv.DictReader(open(path)):
        if r['series'] == 'N':
            out.setdefault(int(r['N']), {})[int(r['step'])] = float(r['Hperp_frac'])
    return out


test = series(os.path.join(RES, 'timeseries_64bit_with124_N3_N40.csv'))
ctrl = series(os.path.join(CTRL, 'timeseries_64bit_with124_N3_N40.csv'))

fig, axs = plt.subplots(3, 5, figsize=(20, 11))
axs = axs.ravel()
for k, N in enumerate(NS):
    ax = axs[k]
    xs = sorted(test[N]); ax.semilogy(xs, [test[N][s] for s in xs], linewidth=1.1, label='R反射あり')
    xs = sorted(ctrl[N]); ax.semilogy(xs, [ctrl[N][s] for s in xs], linewidth=0.9, ls='--', color='gray', label='Control(理論床)')
    ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.grid(alpha=.25)
    ax.set_title(f'N={N}', fontsize=10)
    if k == 0: ax.legend(fontsize=7, loc='lower right')
    if k // 5 == 2: ax.set_xlabel('step')
    if k % 5 == 0: ax.set_ylabel('Hperp/H')
axs[14].axis('off'); axs[15].axis('off') if len(axs) > 15 else None
fig.suptitle('固定 R_{124,23} 反射散乱（対称分割合成）vs Control — 高対称理論床・den=N・500step・N=3..16', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_Hperp_R124_23_vs_control_N3_N16.png'), dpi=170)
plt.close(fig)


def grid(step, fname, title):
    fig, axs = plt.subplots(3, 5, figsize=(20, 12)); axs = axs.ravel()
    for k, N in enumerate(NS):
        ax = axs[k]
        d = np.load(os.path.join(RES, f'hm_N{N}_den_{N}_states_500.npz'))
        z = np.asarray(d['Z'][step], np.complex128)
        for w in z:
            ax.plot([0, w.real], [0, w.imag], color='tab:blue', linewidth=0.6, alpha=0.5)
        ax.plot(z.real, z.imag, 'o', ms=2.2, color='tab:red', alpha=0.85, linestyle='none')
        cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
        for (a, b), c in cnt.items():
            if c > 1:
                ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3), fontsize=5)
        r = float(np.abs(z).max()); lim = r * 1.15 if r > 0 else 1
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
        ax.grid(alpha=.25); ax.tick_params(labelsize=6); ax.ticklabel_format(style='sci', scilimits=(-2, 3))
        ax.set_title(f'N={N} (M={N*(N-1)//2})', fontsize=9)
    for k in range(14, 15):
        axs[k].axis('off')
    if len(axs) > 15: axs[15].axis('off')
    fig.suptitle(title, y=.995); fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=170); plt.close(fig)


grid(0, 'fig_complex_plane_step0_R124_23_N3_N16.png',
     'R124,23 反射系: complex plane at step 0（理論床親）; N=3..16, den=N')
grid(500, 'fig_complex_plane_step500_R124_23_N3_N16.png',
     'R124,23 反射系: complex plane at final step 500; N=3..16, den=N')
print('PLOTS DONE')
