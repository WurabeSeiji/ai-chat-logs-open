#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""修正1（den124）・修正2（allpairs）の標準図化（読み出しのみ）。
使い方: python3 plot_R124_23_variant_figs_v1.py den124|allpairs
出力: fig_Hperp_R124_23_{variant}_vs_control_N3_N16.png、
      fig_complex_plane_step0/step500_R124_23_{variant}_N3_N16.png
Control は理論床の den=124 系列（../対称親v2_500step走行_20260906/results、読み取り専用）。"""
import csv
import os
import sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

VARIANT = sys.argv[1]
assert VARIANT in ('den124', 'allpairs')
BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, f'results_{VARIANT}')
CTRL = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'results')
NS = list(range(3, 17))
LABEL = {'den124': 'R反射（頂点共有I）・den=124統一', 'allpairs': 'R反射（全対I・床免疫）・den=124統一'}[VARIANT]


def series124(path):
    out = {}
    for r in csv.DictReader(open(path)):
        if int(r['denominator']) == 124:
            out.setdefault(int(r['N']), {})[int(r['step'])] = float(r['Hperp_frac'])
    return out


test = series124(os.path.join(RES, 'timeseries_64bit_with124_N3_N40.csv'))
ctrl = series124(os.path.join(CTRL, 'timeseries_64bit_with124_N3_N40.csv'))

fig, axs = plt.subplots(3, 5, figsize=(20, 11)); axs = axs.ravel()
for k, N in enumerate(NS):
    ax = axs[k]
    xs = sorted(test[N]); ax.semilogy(xs, [test[N][s] for s in xs], linewidth=1.1, label=LABEL)
    xs = sorted(ctrl[N]); ax.semilogy(xs, [ctrl[N][s] for s in xs], linewidth=0.9, ls='--', color='gray', label='Control(理論床 den=124)')
    ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.grid(alpha=.25); ax.set_title(f'N={N}', fontsize=10)
    if k == 0: ax.legend(fontsize=6, loc='lower right')
    if k // 5 == 2: ax.set_xlabel('step')
    if k % 5 == 0: ax.set_ylabel('Hperp/H')
for k in (14, 15):
    if k < len(axs): axs[k].axis('off')
fig.suptitle(f'{LABEL} vs Control — 高対称理論床・500step・N=3..16', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, f'fig_Hperp_R124_23_{VARIANT}_vs_control_N3_N16.png'), dpi=170)
plt.close(fig)


def grid(step, fname, title):
    fig, axs = plt.subplots(3, 5, figsize=(20, 12)); axs = axs.ravel()
    for k, N in enumerate(NS):
        ax = axs[k]
        z = np.asarray(np.load(os.path.join(RES, f'hm_N{N}_den_124_states_500.npz'))['Z'][step], np.complex128)
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
    for k in (14, 15):
        if k < len(axs): axs[k].axis('off')
    fig.suptitle(title, y=.995); fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=170); plt.close(fig)


grid(0, f'fig_complex_plane_step0_R124_23_{VARIANT}_N3_N16.png', f'{LABEL}: complex plane at step 0; N=3..16')
grid(500, f'fig_complex_plane_step500_R124_23_{VARIANT}_N3_N16.png', f'{LABEL}: complex plane at final step 500; N=3..16')
print('PLOTS DONE', VARIANT)
