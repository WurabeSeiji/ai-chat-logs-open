#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ロックした刻み（q=8, 4）とロックなし（den=40）のインフレーション図と複素図（2026-10-20）。

インフレーション図：H⊥/H（元プログラムの metrics、床の平面基準）を τ（= step × Δτ）に対して片対数で描く。
複素図：各走行の step 0 と終端の 780 本の関係波 z_e を複素平面に描く（共回転は除かない。形だけを見る）。
入力：results/locked_N40_q{8,4}_states.npz、および元の N3_N40_stage123_sweep_20260905/results/hm_N40_den_40_states_500.npz。
"""
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'run_N40_locked_step_20261020.py')).read().split("if __name__")[0])
FIG = os.path.join(HERE, 'figures'); os.makedirs(FIG, exist_ok=True)
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

runs = []
for qord in (8, 4):
    d = np.load(os.path.join(HERE, 'results', f'locked_N40_q{qord}_states.npz'))
    runs.append((f'ロックあり q={qord}（Δτ=2π/({qord}σ)）', d['Z'].astype(np.complex128), float(d['dtau'])))
ref = np.load(os.path.join(SRC, 'results', 'hm_N40_den_40_states_500.npz'))['Z'].astype(np.complex128)
runs.append(('ロックなし den=40（Δτ=2π/40）', ref, 2 * math.pi / 40))
cols = ['tab:blue', 'tab:orange', 'tab:gray']

# --- インフレーション図 ---
fig, ax = plt.subplots(figsize=(8, 4.8))
for (name, S, dt), c in zip(runs, cols):
    p, q = plane(S[0]); h = np.array([metrics(z, p, q)[0] for z in S])
    ax.semilogy(np.arange(len(S)) * dt, np.maximum(h, 1e-33), color=c, lw=1.2, label=name)
ax.set_xlabel('τ（= step × Δτ）'); ax.set_ylabel('H⊥/H'); ax.set_ylim(1e-33, 1)
ax.set_title('N=40・M=780　インフレーション：ロックあり（q=8, 4）とロックなし（den=40）')
ax.grid(alpha=0.3); ax.legend(fontsize=9, loc='lower right')
fig.tight_layout()
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(FIG, f'fig1_inflation_locked_vs_unlocked_N40.{ext}'), dpi=150)
plt.close(fig)

# --- 複素図 ---
fig, axs = plt.subplots(2, 3, figsize=(13, 8.6))
lim = 1.15 * max(np.abs(S).max() for _, S, _ in runs)
for j, ((name, S, dt), c) in enumerate(zip(runs, cols)):
    for i, (t, lab) in enumerate(((0, 'step 0（床）'), (len(S) - 1, f'終端 step {len(S) - 1}'))):
        a = axs[i, j]; z = S[t]
        a.scatter(z.real, z.imag, s=6, color=c, alpha=0.7)
        a.axhline(0, color='k', lw=0.4); a.axvline(0, color='k', lw=0.4)
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect('equal')
        a.set_title(f'{name}\n{lab}', fontsize=9); a.set_xlabel('Re z'); a.set_ylabel('Im z')
fig.suptitle('N=40・M=780 の関係波の複素図（上：床、下：終端）', fontsize=12)
fig.tight_layout()
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(FIG, f'fig2_complex_plane_locked_vs_unlocked_N40.{ext}'), dpi=150)
plt.close(fig)
print('done')
