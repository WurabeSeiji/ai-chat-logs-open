#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=64 親（step0 = Z0）の複素平面図・単図（読み出しのみ）。
様式は正本 plot_complex_plane_N3_N40_stage123_v1.py の draw_grid と同一
（原点からの線分・赤点・12桁丸め重複の x本数表記・等アスペクト）。
入力: parents_N64/parent_static_N00064_makeparent_20260905.npz の Z0
出力: fig_complex_plane_step0_N64.png"""
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
z = np.asarray(np.load(os.path.join(BASE, 'parents_N64',
                'parent_static_N00064_makeparent_20260905.npz'))['Z0'], np.complex128)
N = 64
M = N * (N - 1) // 2
assert z.size == M

fig, ax = plt.subplots(figsize=(9, 9))
for w in z:
    ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.4, alpha=0.35)
ax.plot(z.real, z.imag, 'o', ms=1.8, color='tab:red', alpha=0.8, linestyle='none')
cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
for (a, b), c in cnt.items():
    if c > 1:
        ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3), fontsize=5)
r = float(np.abs(z).max())
lim = r * 1.15
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
ax.grid(alpha=.25); ax.tick_params(labelsize=8)
ax.ticklabel_format(style='sci', scilimits=(-2, 3))
ax.set_xlabel('Re z'); ax.set_ylabel('Im z')
ax.set_title(f'N={N} (M={M}) make_parent static parent: complex plane at step 0')
fig.tight_layout()
out = os.path.join(BASE, 'fig_complex_plane_step0_N64.png')
fig.savefig(out, dpi=180)
print('saved:', out)
