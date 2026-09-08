#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""無名化→復元 走行の複素平面図: step0（無名化→復元後の初期値）と step500。den=N 系列を表示。"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']
BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'results_anonymize_floor')
NS = [3, 4, 5, 6]

def plane(v):
    p = v.real.astype(np.float64).copy(); p /= np.linalg.norm(p)
    q = v.imag.astype(np.float64).copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    return p, q

# インフレーション図（den=N、Hperp/H）
fig, axs = plt.subplots(1, 4, figsize=(18, 4.5))
for ax, N in zip(axs, NS):
    S = np.asarray(np.load(os.path.join(RES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'], np.complex128)
    p, q = plane(S[0])
    f = np.empty(S.shape[0])
    for s in range(S.shape[0]):
        z = S[s]; zp = z - p * np.dot(p, z) - q * np.dot(q, z)
        f[s] = np.vdot(zp, zp).real / np.vdot(z, z).real
    ax.semilogy(np.arange(S.shape[0]), np.maximum(f, 1e-34), lw=1.0)
    ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.set_title(f'N={N} (M={N*(N-1)//2})')
    ax.grid(alpha=.25); ax.set_xlabel('step')
    if N == NS[0]:
        ax.set_ylabel('Hperp/H')
fig.suptitle('Anonymize->rebuild floor (theoretical floor, den=N): inflation Hperp/H', y=.99)
fig.tight_layout()
p_inf = os.path.join(BASE, 'fig_Hperp_anonymize_floor_N3_N6.png')
fig.savefig(p_inf, dpi=180); plt.close(fig)
print('wrote', os.path.basename(p_inf))

for step, tag in ((0, 'step0'), (500, 'step500')):
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.6))
    for ax, N in zip(axs, NS):
        d = np.load(os.path.join(RES, f'hm_N{N}_den_{N}_states_500.npz'))
        z = np.asarray(d['Z'], np.complex128)[step]
        ax.scatter(z.real, z.imag, s=30, c='tab:blue')
        r = max(np.max(np.abs(z)) * 1.15, 1e-3)
        ax.set_xlim(-r, r); ax.set_ylim(-r, r); ax.set_aspect('equal')
        ax.axhline(0, color='gray', lw=.5); ax.axvline(0, color='gray', lw=.5)
        ax.set_title(f'N={N} (M={N*(N-1)//2}), {tag}')
        ax.grid(alpha=.25)
    fig.suptitle(f'Anonymize->rebuild floor (theoretical floor, shuffle+reconstruct, den=N): complex plane at {tag}', y=.98)
    fig.tight_layout()
    p = os.path.join(BASE, f'fig_complex_plane_{tag}_anonymize_floor.png')
    fig.savefig(p, dpi=180); plt.close(fig)
    print('wrote', os.path.basename(p))
print('ALL DONE')
