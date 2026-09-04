#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""複素平面読出し図（読出し専用・データ生成なし）。
入力: results_2000steps/hm_N{N}_den_{N}_states_2000.npz（Δτ=2π/N 走行）の Z[step]。
N=3..33 を 7×5 グリッド（正本 PNG と同じ配置）で、M=N(N-1)/2 本の複素波 z_e を
複素平面（横=実部 a、縦=虚部 b）にプロットする。中心から各 (a,b) へ線分を引く。
各パネルは全体振幅（max|z|）でスケールを揃えるが、目盛は実際の値のまま。
同一複素数が複数本ある場合は (a,b) に x99 形式の小さい文字で本数を記す。
出力: step0 と step2000 の 2 ファイル。既存データは一切変更しない。"""
import os
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

ROOT = '/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/ChatGPT_denominator_controls_N3_N33_legacyparent_20260903'
IN = os.path.join(ROOT, 'results_2000steps')
OUT = os.path.join(ROOT, 'complex_plane_readout_step0_step2000_20260904')

def draw(step, fname, title):
    fig, axs = plt.subplots(7, 5, figsize=(20, 21))
    axs = axs.ravel()
    for k, N in enumerate(range(3, 34)):
        ax = axs[k]
        d = np.load(os.path.join(IN, f'hm_N{N}_den_{N}_states_2000.npz'))
        assert int(d['denominator']) == N and int(d['steps']) == 2000
        z = np.asarray(d['Z'][step], dtype=np.complex128)
        M = N * (N - 1) // 2
        assert z.size == M
        for w in z:
            ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
        ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
        # 同一複素数の本数表示（倍精度の同値判定は 12 桁丸めで行う）
        cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
        for (a, b), c in cnt.items():
            if c > 1:
                ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3),
                            fontsize=5, color='black')
        r = float(np.abs(z).max())
        lim = r * 1.15 if r > 0 else 1.0
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
        ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
        ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
        ax.grid(alpha=.25)
        ax.tick_params(labelsize=6)
        ax.ticklabel_format(style='sci', scilimits=(-2, 3))
        ax.set_title(f'N={N} (M={M})', fontsize=9)
        if k // 5 == 6: ax.set_xlabel('Re z', fontsize=7)
        if k % 5 == 0: ax.set_ylabel('Im z', fontsize=7)
    for k in range(31, 35):
        axs[k].axis('off')
    fig.suptitle(title, y=.995)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=180)
    plt.close(fig)

draw(0, 'fig_complex_plane_step0_N3_N33.png',
     'Complex-plane readout of M=N(N-1)/2 edge waves at step 0 (from dt=2pi/N runs, states_2000.npz); ticks in actual units; N=3..33')
draw(2000, 'fig_complex_plane_step2000_N3_N33.png',
     'Complex-plane readout of M=N(N-1)/2 edge waves at step 2000 (dt=2pi/N); ticks in actual units; N=3..33')
print('ALL DONE')
