#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M=15 × 2（2シンプレックス）走行（2026-09-08、木原指示）。N=6・M=15 を2つ、索引 r∈{0,1} で並走。

・状態を shape (R=2, M=15) に拡張（シンプレックス索引を1つ追加）。無名性は無視し明示索引で処理。
・各シンプレックスは同じ N=6 理論床から出発し、同じ力学（正本 one_step）で独立に回る。
・物理は物理正本 run_N3_N40_stage123_v1.py（SHA照合）の関数を「そのまま抽出」して使用（byte一致）。
  正本のモジュール本体（N=3..40 スイープ）は実行せず、関数定義部のみ exec で取り込む。
・図: インフレーションはシンプレックス無関係なはず → 全 30 成分を1系として H⊥/H を1枚。
      複素平面は step0 と終端 step を、両シンプレックス重ねて各1枚。
・出力データは results_two_simplex_N6/ に保存（対照 results_control_N6 は残す）。
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
PD = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')
OUT = os.path.join(BASE, 'results_two_simplex_N6')
os.makedirs(OUT, exist_ok=True)

N = 6
R = 2          # シンプレックス数
STEPS = 500

# --- 正本の関数だけを抽出（モジュール本体スイープは実行しない）---
src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
func_src = src.split('rows=[]; summaries=[]')[0]          # 関数定義部＋定数のみ
g = {'__name__': 'canonical_funcs', '__file__': ORIG}
exec(compile(func_src, ORIG, 'exec'), g)
one_step, adjacency, plane, metrics = g['one_step'], g['adjacency'], g['plane'], g['metrics']

# --- 初期状態: N=6 理論床を R=2 部にコピー（索引 r を追加）---
z0_single = np.array(np.load(os.path.join(PD, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],
                     dtype=np.complex128, copy=True)
M = z0_single.size
A = adjacency(N)
den = N
Z0 = np.stack([z0_single.copy() for _ in range(R)], axis=0)   # (R, M)

# --- 走行: 各シンプレックスを独立に1歩ずつ（同じ A・den・one_step）---
states = np.empty((STEPS + 1, R, M), dtype=np.complex128)
z = Z0.copy()
for t in range(STEPS + 1):
    states[t] = z
    if t < STEPS:
        for r in range(R):
            z[r] = one_step(z[r], A, den)
np.savez_compressed(os.path.join(OUT, f'hm_N{N}_R{R}_den_{den}_states_{STEPS}.npz'),
                    Z=states, N=np.int64(N), R=np.int64(R), M=np.int64(M),
                    denominator=np.int64(den), steps=np.int64(STEPS))
print(f'走行完了: 状態 shape={states.shape} 保存', flush=True)

# --- インフレーション図（全 30 成分を1系として、1枚）---
flat = states.reshape(STEPS + 1, R * M)          # (T, 30)
p, q = plane(flat[0])
f = np.empty(STEPS + 1); cl = np.empty(STEPS + 1)
for s in range(STEPS + 1):
    fr, h, c = metrics(flat[s], p, q); f[s] = fr; cl[s] = c
fig, ax = plt.subplots(figsize=(7.5, 4.6))
ax.semilogy(np.arange(STEPS + 1), np.maximum(f, 1e-34), lw=1.2)
ax.set_xlim(0, STEPS); ax.set_ylim(1e-34, 3)
ax.set_title(f'N={N} M={M} × {R} simplices (全体1系): inflation Hperp/H')
ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig(os.path.join(BASE, f'fig_Hperp_two_simplex_N{N}.png'), dpi=180); plt.close(fig)

# --- 複素平面図（step0 と終端、両シンプレックス重ね、各1枚）---
colors = ['tab:blue', 'tab:red']
markers = ['o', 'x']
for step, tag in ((0, 'step0'), (STEPS, f'step{STEPS}')):
    fig, ax = plt.subplots(figsize=(5.6, 5.6))
    for r in range(R):
        zz = states[step, r]
        ax.scatter(zz.real, zz.imag, s=48 if r == 0 else 60, c=colors[r], marker=markers[r],
                   label=f'simplex {r}', alpha=0.8)
    rr = np.max(np.abs(states[step])) * 1.15
    ax.set_xlim(-rr, rr); ax.set_ylim(-rr, rr); ax.set_aspect('equal')
    ax.axhline(0, color='gray', lw=.5); ax.axvline(0, color='gray', lw=.5); ax.grid(alpha=.25)
    ax.set_title(f'N={N} M={M} × {R} simplices: complex plane {tag}')
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, f'fig_complex_plane_{tag}_two_simplex_N{N}.png'), dpi=180); plt.close(fig)

print(f'図: fig_Hperp_two_simplex_N{N}.png, fig_complex_plane_step0/step{STEPS}_two_simplex_N{N}.png', flush=True)
print(f'onset(>0.05)={int(np.flatnonzero(f>0.05)[0]) if np.any(f>0.05) else -1} '
      f'final f={f[-1]:.4f} 最大閉塞={cl.max():.2e}', flush=True)
print('ALL DONE', flush=True)
