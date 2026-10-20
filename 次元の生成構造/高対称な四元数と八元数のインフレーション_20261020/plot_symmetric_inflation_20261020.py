#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""高対称な四元数と八元数のインフレーション実験（途中で停止した 20 組）の図（2026-10-20）。
図 1（インフレーション図）：results/series_*.npz の形の変化 g(t) を log10 で。g = 0 の組（そのままの配置）は描けないので注記だけ。
図 2（複素図）：T24 × C2（κ = +1、自己同型で回した配置）を symmetric_inflation_20261020.py と同じ条件（同じシード、同じ刻み、
  20000 歩）で流し直し、q = z₁ + z₂ j の z₁, z₂ の平面に初期（○）と終端（・）を描く。変位は 10⁻¹² の桁なので、
  拡大率を明記して矢印で描く。
"""
import importlib.util
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('si', os.path.join(HERE, 'symmetric_inflation_20261020.py'))
si = importlib.util.module_from_spec(spec); spec.loader.exec_module(si)

# 図 1
runs = ['T24_C2_k+1_rot', 'T24_C2_k-1_rot', 'O48_C1_k+1_rot', 'O48_C1_k-1_rot', 'T24_C1_k+1_rot', 'T24_C1_k-1_rot']
fig, ax = plt.subplots(figsize=(8, 5))
for n in runs:
    d = np.load(os.path.join(HERE, 'results', f'series_{n}.npz'))
    g = d['g']; m = g > 0
    ax.plot(d['t'][m], np.log10(g[m]), label=n.replace('_rot', '（回した配置）'))
ax.set_xlabel('時間 τ'); ax.set_ylabel('log₁₀ g(τ)（形の変化）')
ax.set_title('高対称な四元数のインフレーション図（種は丸め誤差だけ）\nそのままの配置（CROSS8、T24、O48 など）は g = 0 のまま（描けない）', fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figures', 'fig1_inflation_symmetric_quaternion.png'), dpi=150)
fig.savefig(os.path.join(HERE, 'figures', 'fig1_inflation_symmetric_quaternion.svg'))

# 図 2：T24 × C2、κ = +1、回した配置を同じ条件で流し直す
rng = np.random.default_rng(si.SEED)
qc = si.Q.configs(np.random.default_rng(si.SEED))
X = qc['CROSS8']; u = rng.normal(size=4); u /= np.linalg.norm(u)      # CROSS8 で消費した乱数を同じ順に進める
X = qc['T24']; N = len(X)
u = rng.normal(size=4); u /= np.linalg.norm(u)
X0 = np.array([si.Q.mul(si.Q.mul(u, x), u * si.Q.CONJ) for x in X])
Xg = np.random.default_rng(si.SEED + N).normal(size=(N, 4))
sref = float(max(np.linalg.norm(si.Q.F_C2(Xg, 1), axis=1) / np.linalg.norm(Xg, axis=1)))
dt = 2 * math.pi / (200 * sref); r0 = np.linalg.norm(X0, axis=1); Xt = X0.copy()
for _ in range(20000):
    Xt = si.rk4r(lambda Y: si.Q.F_C2(Y, 1), Xt, dt, r0)
G0 = X0 @ X0.T; gend = float(np.linalg.norm(Xt @ Xt.T - G0) / np.linalg.norm(G0))
d = Xt - X0; amp = 0.25 / max(np.abs(d).max(), 1e-300)
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
for k, (a, b, name) in enumerate(((0, 1, 'z₁ = q₀ + q₁ i'), (2, 3, 'z₂ = q₂ + q₃ i'))):
    ax = axs[k]
    ax.scatter(X0[:, a], X0[:, b], s=60, facecolors='none', edgecolors='C0', label='初期')
    ax.scatter(Xt[:, a], Xt[:, b], s=12, color='C3', label='終端')
    ax.quiver(X0[:, a], X0[:, b], d[:, a] * amp, d[:, b] * amp, angles='xy', scale_units='xy', scale=1, color='C2', width=0.004,
              label=f'変位 × {amp:.1e}')
    ax.set_aspect('equal'); ax.set_title(name); ax.grid(alpha=0.3); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
axs[0].legend(fontsize=8, loc='lower left')
fig.suptitle(f'複素図：T24 × C2（κ = +1、回した配置）、q = z₁ + z₂ j。20000 歩（τ = {20000 * dt:.0f}）、形の変化 g = {gend:.1e}', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'figures', 'fig2_complex_T24_C2.png'), dpi=150)
fig.savefig(os.path.join(HERE, 'figures', 'fig2_complex_T24_C2.svg'))
print('g_end', gend, 'max displacement', float(np.abs(d).max()))
