#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=9 (M=36) 単一シンプレックス走行【90度系】（2026-09-09、木原指示・加法法則N9=N7+N6実験の単独側を0°/90°初期値で追試）。ロールバック対照テスト_正本実装N3_N6_20260908/run_single_simplex_N9_v1.py の忠実コピーに最小変更のみ: 初期値を staged → 90度系（0°/90°）親 parent_symmetric_N*_v2.npz へ差替（＋出力先ラベル）。物理・構造・図・走行手順は無変更。

・状態 shape (R=1, M=36)。simplex 0 = N=9 高対称初期値【90度系: 位相 0°/90° のみ】（最も対称性の高い初期値_20260906/parents_symmetric、H=1）、位相回転なし（PHASES=[0]）。
・シンプレックスは正本 one_step（関数を byte 一致で抽出）で回す。
・物理は物理正本 run_N3_N40_stage123_v1.py（SHA照合）。正本モジュール本体は実行せず関数定義部のみ取り込む。
・図（別名）: インフレーション1枚（全45成分1系）、step0 と終端 step の複素平面図 各1枚（N=6×3 系列と同形式）。
・データ・本ラッパーとも保存し再現性を確保。既存 results_* は残す。
"""
import hashlib
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
PD = os.path.join(BASE, '..', '最も対称性の高い初期値_20260906', 'parents_symmetric')
OUT = os.path.join(BASE, 'results_single_simplex_N9')
os.makedirs(OUT, exist_ok=True)

N = 9
R = 1
STEPS = 500
PHASES = [r * (2.0 * math.pi / N) for r in range(R)]     # 0, 2π/6, 2·2π/6 = 0°,60°,120°

# --- 正本の関数だけを抽出（モジュール本体スイープは実行しない）---
src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical_funcs', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
one_step, adjacency, plane, metrics = g['one_step'], g['adjacency'], g['plane'], g['metrics']

# --- 初期状態: simplex r = 理論床 × 大域位相 e^{i·r·2π/6} ---
z0_single = np.array(np.load(os.path.join(PD, f'parent_symmetric_N{N:05d}_v2.npz'))['Z0'],
                     dtype=np.complex128, copy=True)
M = z0_single.size
A = adjacency(N)
den = N
Z0 = np.stack([np.exp(1j * ph) * z0_single for ph in PHASES], axis=0)   # (R, M)

# --- 走行 ---
states = np.empty((STEPS + 1, R, M), dtype=np.complex128)
z = Z0.copy()
for t in range(STEPS + 1):
    states[t] = z
    if t < STEPS:
        for r in range(R):
            z[r] = one_step(z[r], A, den)
np.savez_compressed(os.path.join(OUT, f'hm_N{N}_R{R}_phase_den_{den}_states_{STEPS}.npz'),
                    Z=states, N=np.int64(N), R=np.int64(R), M=np.int64(M),
                    denominator=np.int64(den), steps=np.int64(STEPS),
                    phase_offsets=np.array(PHASES, dtype=np.float64))
print(f'走行完了: 状態 shape={states.shape} 大域位相offset={[round(math.degrees(p),1) for p in PHASES]}deg 保存', flush=True)

# --- インフレーション図（全 R·M=45 成分を1系, 別名）---
flat = states.reshape(STEPS + 1, R * M)
p, q = plane(flat[0])
f = np.empty(STEPS + 1); cl = np.empty(STEPS + 1)
for s in range(STEPS + 1):
    fr, h, c = metrics(flat[s], p, q); f[s] = fr; cl[s] = c
fig, ax = plt.subplots(figsize=(7.5, 4.6))
ax.semilogy(np.arange(STEPS + 1), np.maximum(f, 1e-34), lw=1.2)
ax.set_xlim(0, STEPS); ax.set_ylim(1e-34, 3)
ax.set_title(f'N={N} M={M} × {R} (単一simplex 90度系0/90): inflation Hperp/H')
ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig(os.path.join(BASE, f'fig_Hperp_single_simplex_N{N}.png'), dpi=180); plt.close(fig)

# --- 複素平面図（step0 と終端, 別名, 3シンプレックス重ね）---
colors = ['tab:blue', 'tab:red', 'tab:green']; markers = ['o', 'x', '+']; sizes = [48, 62, 80]
for step, tag in ((0, 'step0'), (STEPS, f'step{STEPS}')):
    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    for r in range(R):
        zz = states[step, r]
        ax.scatter(zz.real, zz.imag, s=sizes[r], c=colors[r], marker=markers[r],
                   label=f'simplex {r} (+{int(round(math.degrees(PHASES[r])))}°)', alpha=0.8)
    rr = np.max(np.abs(states[step])) * 1.15
    ax.set_xlim(-rr, rr); ax.set_ylim(-rr, rr); ax.set_aspect('equal')
    ax.axhline(0, color='gray', lw=.5); ax.axvline(0, color='gray', lw=.5); ax.grid(alpha=.25)
    ax.set_title(f'N={N} M={M} × {R} (単一simplex 90度系0/90): complex plane {tag}')
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, f'fig_complex_plane_{tag}_single_simplex_N{N}.png'), dpi=180); plt.close(fig)

print(f'図(別名): fig_Hperp_single_simplex_N{N}.png, fig_complex_plane_step0/step{STEPS}_single_simplex_N{N}.png', flush=True)
print(f'onset(>0.05)={int(np.flatnonzero(f>0.05)[0]) if np.any(f>0.05) else -1} final f={f[-1]:.4f} 最大閉塞={cl.max():.2e}', flush=True)
print('ALL DONE', flush=True)
