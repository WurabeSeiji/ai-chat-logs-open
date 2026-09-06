#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""段1+2+3 スイープの複素平面読出し図（読み出しのみ・3枚）。
入力: results/hm_N{N}_den_{N}_states_500.npz（Δτ=2π/N 走行）の Z[0] と Z[500]。
(1) step0 グリッド、(2) 終了時（step500）グリッド:
    complex_plane_readout_step0_step2000_20260904/plot_complex_plane_step0_step2000.py と同一様式
    （原点からの線分・全体振幅スケール・実値目盛・12桁丸め重複の x本数表記）を N=3..40 の 8×5 に拡張。
(3) 終了時の拡大図グリッド: 各 N の最大角クラスター（粗解像度 1/100 群化、
    plot_complex_plane_N40_v1.py と同一アルゴリズム）を拡大し、本数・中心・割れ幅を題記。"""
import os
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(BASE, 'results')

def load(N, step):
    d = np.load(os.path.join(IN, f'hm_N{N}_den_{N}_states_500.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    return np.asarray(d['Z'][step], dtype=np.complex128)

def draw_grid(step, fname, title):
    fig, axs = plt.subplots(8, 5, figsize=(20, 24))
    axs = axs.ravel()
    for k, N in enumerate(range(3, 41)):
        ax = axs[k]
        z = load(N, step)
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
        ax.set_title(f'N={N} (M={M})', fontsize=9)
        if k // 5 == 7: ax.set_xlabel('Re z', fontsize=7)
        if k % 5 == 0: ax.set_ylabel('Im z', fontsize=7)
    for k in range(38, 40):
        axs[k].axis('off')
    fig.suptitle(title, y=.998)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=180)
    plt.close(fig)

draw_grid(0, 'fig_complex_plane_step0_N3_N40_stage123.png',
          'Stage1+2+3 sweep: complex-plane readout at step 0 (make_parent static parents, dt=2pi/N files); N=3..40')
draw_grid(500, 'fig_complex_plane_final_N3_N40_stage123.png',
          'Stage1+2+3 sweep: complex-plane readout at final step 500 (dt=2pi/N); N=3..40')

# (3) 終了時の最大角クラスター拡大グリッド
fig, axs = plt.subplots(8, 5, figsize=(20, 24))
axs = axs.ravel()
for k, N in enumerate(range(3, 41)):
    ax = axs[k]
    z = load(N, 500)
    amp = float(np.abs(z).max())
    coarse = {}
    for w in z:
        key = (round(float(w.real) / amp, 2), round(float(w.imag) / amp, 2))
        coarse.setdefault(key, []).append(w)
    mem = max(coarse.values(), key=len)
    zz = np.array(mem)
    c = zz.mean()
    dev = np.abs(zz - c)
    spread = float(dev.max())
    win = spread * 1.4 if spread > 0 else amp * 1e-12
    ax.plot(zz.real, zz.imag, 'o', ms=3, color='tab:red', alpha=0.8, linestyle='none')
    cnt = Counter((round(float(w.real), 15), round(float(w.imag), 15)) for w in zz)
    for (a, b), n in cnt.items():
        if n > 1:
            ax.annotate(f'x{n}', (a, b), textcoords='offset points', xytext=(3, 3),
                        fontsize=5, color='black')
    ax.set_xlim(c.real - win, c.real + win); ax.set_ylim(c.imag - win, c.imag + win)
    ax.set_aspect('equal')
    ax.grid(alpha=.25)
    ax.tick_params(labelsize=5)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3), useOffset=True)
    ax.set_title(f'N={N}: {len(zz)} waves, dev={spread:.2e} (|z|max={amp:.2e})', fontsize=7)
for k in range(38, 40):
    axs[k].axis('off')
fig.suptitle('Stage1+2+3 sweep: zoom into largest angle cluster at final step 500 (dt=2pi/N); N=3..40', y=.998)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_complex_plane_final_zoom_N3_N40_stage123.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
