#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""§5 図B/C/D: 複素平面図（正本 plot_complex_plane_N3_N40_stage123_v1.py と同一様式:
原点からの線分・赤点・12桁丸め重複 x本数・等アスペクト）。データは full_N3_N40_sweep/states/ の den=N。
B: step0 全N 8x5 / C: step500 全N 8x5 / D: N=3..10 の step0|step500 左右対比（8行2列）。"""
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
ST = os.path.join(BASE, 'full_N3_N40_sweep', 'states')
OUT = os.path.join(BASE, 'full_N3_N40_sweep')


def load(N, step):
    d = np.load(os.path.join(ST, f'hm_N{N}_den_{N}_states_500.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    return np.asarray(d['Z'][step], dtype=np.complex128)


def panel(ax, z, title, ms=2.5):
    for w in z:
        ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    ax.plot(z.real, z.imag, 'o', ms=ms, color='tab:red', alpha=0.85, linestyle='none')
    cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    for (a, b), c in cnt.items():
        if c > 1:
            ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3), fontsize=5)
    r = float(np.abs(z).max()); lim = r * 1.15 if r > 0 else 1.0
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5); ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.grid(alpha=.25); ax.tick_params(labelsize=6); ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.set_title(title, fontsize=9)


def grid(step, fname, title):
    fig, axs = plt.subplots(8, 5, figsize=(20, 24)); axs = axs.ravel()
    for k, N in enumerate(range(3, 41)):
        z = load(N, step)
        panel(axs[k], z, f'N={N} (M={N*(N-1)//2})')
    for k in range(38, 40):
        axs[k].axis('off')
    fig.suptitle(title, y=.998); fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=180); plt.close(fig)


grid(0, 'fig_complex_plane_step0_makeparent_N3_N40.png',
     'make_parent parents sweep: complex plane at step 0 (den=N); N=3..40')
grid(500, 'fig_complex_plane_step500_makeparent_N3_N40.png',
     'make_parent parents sweep: complex plane at final step 500 (den=N); N=3..40')

fig, axs = plt.subplots(8, 2, figsize=(9, 34))
for i, N in enumerate(range(3, 11)):
    panel(axs[i, 0], load(N, 0), f'N={N} step0', ms=3)
    panel(axs[i, 1], load(N, 500), f'N={N} step500', ms=3)
fig.suptitle('make_parent parents: step0 vs step500 (den=N); N=3..10', y=.999)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_complex_plane_step0_step500_makeparent_N3_N10.png'), dpi=170)
plt.close(fig)
print('PLOTS DONE')
