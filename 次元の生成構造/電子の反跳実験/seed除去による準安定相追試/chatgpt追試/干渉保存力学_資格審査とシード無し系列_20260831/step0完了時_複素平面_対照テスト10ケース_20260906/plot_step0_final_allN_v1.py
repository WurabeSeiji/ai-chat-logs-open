#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step0 と完了時（step10000）の複素平面読出し図・全ケース N=3..40（読み出しのみ・新規走行なし）。
plot_step0_final_10cases_v1.py の最小変更版: NS を代表10ケースから全38ケース N=3..40 に置換
しただけ（per-panel 描画コード・入力・SHAゲートは同一）。個別38図に加え、俯瞰用に
step0 と完了時の 8×5 グリッド（A の様式）も出す。
入力: ../N3_N40_long10000_20260905/results/hm_N{N}_den_{N}_states_10000.npz（SHA台帳照合）。
出力: fig_step0_final_N{N}.png × 38、fig_step0_grid_allN.png、fig_final_grid_allN.png。"""
import hashlib
import os
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
LONG = os.path.join(SERIES, 'N3_N40_long10000_20260905')
IN = os.path.join(LONG, 'results')

NS = list(range(3, 41))
FINAL = 10000

ledger = {}
with open(os.path.join(LONG, 'SHA256SUMS.txt')) as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]


def load(N, step):
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    path = os.path.join(LONG, rel)
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    d = np.load(path)
    assert int(d['denominator']) == N and int(d['steps']) == FINAL
    return np.asarray(d['Z'][step], dtype=np.complex128)


def draw_panel(ax, z, title, small=False):
    """10cases版と同一（A の per-panel 描画コード verbatim 流用）。small はグリッド用の縮小字。"""
    for w in z:
        ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
    cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    for (a, b), c in cnt.items():
        if c > 1:
            ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3),
                        fontsize=5 if small else 6, color='black')
    r = float(np.abs(z).max())
    lim = r * 1.15 if r > 0 else 1.0
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.grid(alpha=.25)
    ax.tick_params(labelsize=6 if small else 7)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.set_title(title, fontsize=9 if small else 10)
    if not small:
        ax.set_xlabel('Re z', fontsize=8); ax.set_ylabel('Im z', fontsize=8)


# (A) 各Nの個別図 38枚: [step0 ｜ 完了時]
for N in NS:
    M = N * (N - 1) // 2
    z0 = load(N, 0)
    zf = load(N, FINAL)
    assert z0.size == M and zf.size == M
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.6))
    draw_panel(axs[0], z0, f'N={N} (M={M}) — step 0')
    draw_panel(axs[1], zf, f'N={N} (M={M}) — final step {FINAL}')
    fig.suptitle(f'Stage1+2+3 (dt=2pi/N): complex-plane readout, N={N}  [step 0 vs final]', y=.99)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, f'fig_step0_final_N{N}.png'), dpi=170)
    plt.close(fig)
    print(f'N={N}: |z0|max={np.abs(z0).max():.3e} |zf|max={np.abs(zf).max():.3e}', flush=True)

# (B) 俯瞰グリッド 8×5（step0 / 完了時）— A の様式（38パネル＋余白2消し）
def draw_grid(step, fname, title):
    fig, axs = plt.subplots(8, 5, figsize=(20, 24))
    axs = axs.ravel()
    for k, N in enumerate(NS):
        draw_panel(axs[k], load(N, step), f'N={N} (M={N*(N-1)//2})', small=True)
    for k in range(len(NS), 40):
        axs[k].axis('off')
    fig.suptitle(title, y=.998)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=160)
    plt.close(fig)

draw_grid(0, 'fig_step0_grid_allN.png',
          'Stage1+2+3 (dt=2pi/N): complex-plane readout at step 0; N=3..40 (existing 10000-step data)')
draw_grid(FINAL, 'fig_final_grid_allN.png',
          f'Stage1+2+3 (dt=2pi/N): complex-plane readout at final step {FINAL}; N=3..40 (existing data)')
print('grids 保存: fig_step0_grid_allN.png / fig_final_grid_allN.png')
print('ALL DONE')
