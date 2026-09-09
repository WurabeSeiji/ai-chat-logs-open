#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=7+6 混合走行【90度解析床】の部分系別インフレーション図（2026-09-09、木原指示）。
plot_Hperp_mixed_parts_N7N6_v1.py の忠実コピー（コード無変更・BASE相対で本フォルダの90度解析床データに自動適用）。

union（全36成分の共通平面）でのインフレーション図は無意味のため、保存済み状態
results_mixed_simplex_N7N6/hm_mixed_N7N6_states_500.npz から、各部分系（N=7 の21成分・
N=6 の15成分）ごとに自分の step0 支配平面を基準に H⊥/H を計算し、2枚のインフレーション
図を作る。力学の再走行はしない（図化のみ）。plane/metrics は物理正本から SHA 照合で抽出。
図の様式は run_single_simplex_N9_v1.py 等のインフレーション図と同一。
"""
import hashlib
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical_funcs', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
plane, metrics = g['plane'], g['metrics']

d = np.load(os.path.join(BASE, 'results_mixed_simplex_N7N6', 'hm_mixed_N7N6_states_500.npz'))
Z = d['Z']; sims = d['sims']; offsets = d['offsets']; STEPS = int(d['steps'])

for r, N in enumerate(sims):
    s0, s1 = int(offsets[r]), int(offsets[r + 1])
    part = Z[:, s0:s1]
    M = s1 - s0
    p, q = plane(part[0])
    f = np.empty(STEPS + 1); cl = np.empty(STEPS + 1)
    for s in range(STEPS + 1):
        fr, h, c = metrics(part[s], p, q); f[s] = fr; cl[s] = c
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.semilogy(np.arange(STEPS + 1), np.maximum(f, 1e-34), lw=1.2)
    ax.set_xlim(0, STEPS); ax.set_ylim(1e-34, 3)
    ax.set_title(f'N={N} part of N7+N6 mixed (M={M}, H={M}/36): inflation Hperp/H')
    ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.grid(alpha=.25)
    fig.tight_layout()
    out = os.path.join(BASE, f'fig_Hperp_mixed_simplex_part_N{N}.png')
    fig.savefig(out, dpi=180); plt.close(fig)
    onset = int(np.flatnonzero(f > 0.05)[0]) if np.any(f > 0.05) else -1
    print(f'N={N} part: step0 f={f[0]:.3e} onset(>0.05)={onset} final f={f[-1]:.4f} '
          f'最大閉塞={cl.max():.2e} -> {os.path.basename(out)}', flush=True)
print('ALL DONE', flush=True)
