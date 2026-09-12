#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40 最終ステップ（step2000）の複素平面図・単図（読み出しのみ）。
N=64版 plot_complex_plane_final_N64_v1.py の忠実コピー（N/入出力名/タイトルのみ変更）。
入力: results/hm_N40_den_40_states_2000.npz の Z[2000]
出力: fig_complex_plane_step2000_N40.png"""
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(BASE, 'results', 'hm_N40_den_40_states_2000.npz'))
assert int(d['denominator']) == 40 and int(d['steps']) == 2000
z = np.asarray(d['Z'][2000], np.complex128)
N = 40; M = N * (N - 1) // 2
assert z.size == M
fig, ax = plt.subplots(figsize=(9, 9))
for w in z:
    ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.4, alpha=0.35)
ax.plot(z.real, z.imag, 'o', ms=1.8, color='tab:red', alpha=0.8, linestyle='none')
cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
for (a, b), c in cnt.items():
    if c > 1:
        ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3), fontsize=5)
r = float(np.abs(z).max()); lim = r * 1.15
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5); ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
ax.grid(alpha=.25); ax.tick_params(labelsize=8); ax.ticklabel_format(style='sci', scilimits=(-2, 3))
ax.set_xlabel('Re z'); ax.set_ylabel('Im z')
ax.set_title(f'N={N} (M={M}) make_parent parent: complex plane at final step 2000 (den=N)')
fig.tight_layout()
out = os.path.join(BASE, 'fig_complex_plane_step2000_N40.png')
fig.savefig(out, dpi=180)
print('saved:', out)
