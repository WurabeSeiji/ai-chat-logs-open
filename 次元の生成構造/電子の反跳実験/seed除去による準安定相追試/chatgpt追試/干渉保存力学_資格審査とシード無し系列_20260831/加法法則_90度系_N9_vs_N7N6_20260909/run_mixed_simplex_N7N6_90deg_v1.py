#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=7 + N=6 混合2シンプレックス走行【90度系】（2026-09-09、木原指示・加法法則 N9=N7+N6 実験の合成側を0°/90°初期値で追試）。
ロールバック対照テスト_正本実装N3_N6_20260908/run_mixed_simplex_N7N6_v1.py の忠実コピーに最小変更のみ: 初期値を 90度系正規化親 parent_norm_*_90deg へ差替（＋出力先ラベル）。物理・構造・図・走行手順は無変更。

・3本同種 N=6 → 2本異種 (N=7, M=21) + (N=6, M=15)。各シンプレックスは自分の adjacency(N)・den=N で独立に回す。
・初期値 = make_normalized_parents_N7N6_90deg_v1.py が作った90度系正規化親（0°/90°、密度 H/M = 1/36 に統一、合計 H=1）。位相回転なし。
・状態は flat 連結 (STEPS+1, 36)（先頭21成分=N7、後15成分=N6）。図・保存の形式は基と同一。
・物理は物理正本 run_N3_N40_stage123_v1.py（SHA照合）。正本モジュール本体は実行せず関数定義部のみ取り込む。
・図（別名）: インフレーション1枚（全36成分1系）、step0 と終端 step の複素平面図 各1枚（2シンプレックス重ね）。
・データ・本ラッパーとも保存し再現性を確保。既存 results_* は残す。
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
PD = os.path.join(BASE, 'parents_mixed_N7N6')
OUT = os.path.join(BASE, 'results_mixed_simplex_N7N6')
os.makedirs(OUT, exist_ok=True)

SIMS = [7, 6]            # 各シンプレックスの N（M = N(N-1)/2 = 21, 15）
STEPS = 500

# --- 正本の関数だけを抽出（モジュール本体スイープは実行しない）---
src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical_funcs', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
one_step, adjacency, plane, metrics = g['one_step'], g['adjacency'], g['plane'], g['metrics']

# --- 初期状態: 正規化親（密度 H/M = 1/36）を flat 連結 ---
z0_list, A_list, den_list, Ms = [], [], [], []
for N in SIMS:
    z0 = np.array(np.load(os.path.join(PD, f'parent_norm_N{N:05d}_HtoM_over36_90deg_20260909.npz'))['Z0'],
                  dtype=np.complex128, copy=True)
    z0_list.append(z0)
    A_list.append(adjacency(N))
    den_list.append(N)
    Ms.append(z0.size)
M_total = sum(Ms)
offsets = np.cumsum([0] + Ms)          # [0, 21, 36]
Z0 = np.concatenate(z0_list)           # (36,)

# --- 走行（各シンプレックスは自分の A・den で独立）---
states = np.empty((STEPS + 1, M_total), dtype=np.complex128)
z = Z0.copy()
for t in range(STEPS + 1):
    states[t] = z
    if t < STEPS:
        for r, N in enumerate(SIMS):
            s0, s1 = offsets[r], offsets[r + 1]
            z[s0:s1] = one_step(z[s0:s1], A_list[r], den_list[r])
np.savez_compressed(os.path.join(OUT, f'hm_mixed_N7N6_states_{STEPS}.npz'),
                    Z=states, sims=np.array(SIMS, dtype=np.int64),
                    Ms=np.array(Ms, dtype=np.int64), offsets=np.array(offsets, dtype=np.int64),
                    denominators=np.array(den_list, dtype=np.int64), steps=np.int64(STEPS))
print(f'走行完了: 状態 shape={states.shape} sims={SIMS} M={Ms} den={den_list} 保存', flush=True)

# --- インフレーション図（全 36 成分を1系, 別名）---
flat = states
p, q = plane(flat[0])
f = np.empty(STEPS + 1); cl = np.empty(STEPS + 1)
for s in range(STEPS + 1):
    fr, h, c = metrics(flat[s], p, q); f[s] = fr; cl[s] = c
fig, ax = plt.subplots(figsize=(7.5, 4.6))
ax.semilogy(np.arange(STEPS + 1), np.maximum(f, 1e-34), lw=1.2)
ax.set_xlim(0, STEPS); ax.set_ylim(1e-34, 3)
ax.set_title(f'N=7+6 mixed 90deg0/90 (M=21+15=36, H∝M): inflation Hperp/H')
ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'fig_Hperp_mixed_simplex_N7N6.png'), dpi=180); plt.close(fig)

# --- 複素平面図（step0 と終端, 別名, 2シンプレックス重ね）---
colors = ['tab:blue', 'tab:red']; markers = ['o', 'x']; sizes = [48, 62]
for step, tag in ((0, 'step0'), (STEPS, f'step{STEPS}')):
    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    for r, N in enumerate(SIMS):
        s0, s1 = offsets[r], offsets[r + 1]
        zz = states[step, s0:s1]
        ax.scatter(zz.real, zz.imag, s=sizes[r], c=colors[r], marker=markers[r],
                   label=f'N={N} (M={Ms[r]}, H={Ms[r]}/36)', alpha=0.8)
    rr = np.max(np.abs(states[step])) * 1.15
    ax.set_xlim(-rr, rr); ax.set_ylim(-rr, rr); ax.set_aspect('equal')
    ax.axhline(0, color='gray', lw=.5); ax.axvline(0, color='gray', lw=.5); ax.grid(alpha=.25)
    ax.set_title(f'N=7+6 mixed 90deg0/90 (M=21+15=36, H∝M): complex plane {tag}')
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, f'fig_complex_plane_{tag}_mixed_simplex_N7N6.png'), dpi=180); plt.close(fig)

print(f'図(別名): fig_Hperp_mixed_simplex_N7N6.png, fig_complex_plane_step0/step{STEPS}_mixed_simplex_N7N6.png', flush=True)
print(f'onset(>0.05)={int(np.flatnonzero(f>0.05)[0]) if np.any(f>0.05) else -1} final f={f[-1]:.4f} 最大閉塞={cl.max():.2e}', flush=True)
print('ALL DONE', flush=True)
