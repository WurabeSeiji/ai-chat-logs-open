#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""合成親実験の複素平面図: step0（処理後の合成親そのもの）と step500。
色分け: ブロック内波（凝縮体1・2、青）／クロス波（赤）。ε=1e-8 のクロス波は step0 では原点近傍。"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']
BASE = os.path.dirname(os.path.abspath(__file__))
NS = [6, 8, 10, 12]

for step, tag in ((0, 'step0'), (500, 'step500')):
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.6))
    for ax, N in zip(axs, NS):
        d = np.load(os.path.join(BASE, 'results_composite', f'hm_N{N}_den_{N}_states_500.npz'))
        z = np.asarray(d['Z'], np.complex128)[step]
        k = N // 2
        ea, eb = np.triu_indices(N, k=1)
        cross = (ea < k) != (eb < k)
        ax.scatter(z[~cross].real, z[~cross].imag, s=26, c='tab:blue', label=f'intra 2x{k}-floor')
        ax.scatter(z[cross].real, z[cross].imag, s=26, c='tab:red', marker='x', label='cross (born)')
        r = np.max(np.abs(z)) * 1.15
        ax.set_xlim(-r, r); ax.set_ylim(-r, r); ax.set_aspect('equal')
        ax.axhline(0, color='gray', lw=.5); ax.axvline(0, color='gray', lw=.5)
        ax.set_title(f'N={N} = {k}+{k} composite, {tag}')
        ax.grid(alpha=.25)
        if N == NS[0]:
            ax.legend(fontsize=8, loc='upper right')
    fig.suptitle(f'Composite parents k+k (eps=1e-8, contact gauge 0, dtau=2pi/N): complex plane at {tag}', y=.98)
    fig.tight_layout()
    p = os.path.join(BASE, f'fig_complex_plane_{tag}_composite_kk.png')
    fig.savefig(p, dpi=180)
    plt.close(fig)
    print('wrote', os.path.basename(p))
print('ALL DONE')
