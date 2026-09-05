#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""段1+2+3 スイープ N=3・Δτ=2π/3 のストロボ複素平面図（読み出しのみ・1枚）。
plot_complex_plane_N3_N40_stage123_v1.py（対照実験でグリッド図3枚の bit 一致を確認済み）の
忠実コピーからの最小変更版。グリッドのループ変数を「N=3..40（step固定）」から
「step=0,3,...,117（N=3固定・40コマ）」に置き換えただけで、各パネルの描画様式
（原点からの線分・全体振幅スケール 1.15r・実値目盛・12桁丸め重複の x本数表記）は同一。
入力: data/hm_N3_den_3_states_500.npz（本体スイープ results/ のコピー、SHA256 一致確認済み）。
Δτ=2π/3 なので 3 step ごとの標本化は時計1回転（位相2π）ごとのストロボに当たる。"""
import os
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(BASE, 'data')

N = 3

def load(step):
    d = np.load(os.path.join(IN, f'hm_N{N}_den_{N}_states_500.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    return np.asarray(d['Z'][step], dtype=np.complex128)

def draw_grid(steps, fname, title):
    fig, axs = plt.subplots(8, 5, figsize=(20, 24))
    axs = axs.ravel()
    for k, step in enumerate(steps):
        ax = axs[k]
        z = load(step)
        M = N * (N - 1) // 2
        assert z.size == M
        for w in z:
            ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
        ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
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
        ax.set_title(f'step {step}', fontsize=9)
        if k // 5 == 7: ax.set_xlabel('Re z', fontsize=7)
        if k % 5 == 0: ax.set_ylabel('Im z', fontsize=7)
    fig.suptitle(title, y=.998)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=180)
    plt.close(fig)

draw_grid(range(0, 120, 3), 'fig_complex_plane_N3_den3_steps0_117.png',
          'Stage1+2+3 sweep: N=3 (M=3), dt=2pi/3 — complex-plane readout every 3 steps (one clock turn per frame), step 0..117')
print('ALL DONE')
