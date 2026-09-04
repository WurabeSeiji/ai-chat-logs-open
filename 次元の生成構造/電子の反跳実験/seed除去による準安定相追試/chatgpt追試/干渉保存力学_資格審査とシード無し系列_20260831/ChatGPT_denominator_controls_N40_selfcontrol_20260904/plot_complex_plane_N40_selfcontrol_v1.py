#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40 自己対照走行（新相互作用・Δτ=2π/40）の複素平面読出し図（読み出しのみ）。
入力: results/hm_N40_den_40_states_500.npz の Z[0]（step0）と Z[500]（最大step）。
描画部は 自発的分裂予備実験_v1/N40_state_readout_20260904/plot_complex_plane_N40_v1.py と
同一様式（draw_plane・クラスター拡大とも同じアルゴリズム）。入力読出しのみ本npz形式に対応。
出力: step0 図・最大step図・凝縮部（角クラスター）拡大図の3枚。データは一切変更しない。"""
import os
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(BASE, 'results', 'hm_N40_den_40_states_500.npz')
d = np.load(IN)
Z = np.asarray(d['Z'], dtype=np.complex128)
STEPS = int(d['steps'])
assert Z.shape == (STEPS + 1, 780) and int(d['denominator']) == 40
Z0 = Z[0]
ZF = Z[STEPS]

def draw_plane(z, title, fname):
    fig, ax = plt.subplots(figsize=(8, 8))
    for w in z:
        ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.5, alpha=0.5)
    ax.plot(z.real, z.imag, 'o', ms=3, color='tab:red', alpha=0.85, linestyle='none')
    cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    for (a, b), c in cnt.items():
        if c > 1:
            ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(4, 4),
                        fontsize=6, color='black')
    r = float(np.abs(z).max())
    lim = r * 1.15 if r > 0 else 1.0
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.grid(alpha=.25)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.set_xlabel('Re z'); ax.set_ylabel('Im z')
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=180)
    plt.close(fig)

draw_plane(Z0, 'N=40 (M=780) selfcontrol step 0: hm_mp_free parent_v (new-interaction input)',
           'fig_N40_selfcontrol_complex_plane_step0.png')
draw_plane(ZF, f'N=40 (M=780) selfcontrol dtau=2pi/40: final step tau={STEPS}',
           'fig_N40_selfcontrol_complex_plane_final.png')

# 凝縮部（角クラスター）拡大: plot_complex_plane_N40_v1.py と同一アルゴリズム
# （粗解像度 1/100 で群化し、大きい順に4クラスターを各パネルで拡大）
amp = float(np.abs(ZF).max())
coarse = {}
for w in ZF:
    key = (round(float(w.real) / amp, 2), round(float(w.imag) / amp, 2))
    coarse.setdefault(key, []).append(w)
clusters = sorted(coarse.values(), key=len, reverse=True)[:4]
fig, axs = plt.subplots(2, 2, figsize=(13, 12))
axs = axs.ravel()
for k, mem in enumerate(clusters):
    ax = axs[k]
    zz = np.array(mem)
    c = zz.mean()
    dev = np.abs(zz - c)
    spread = float(dev.max())
    win = spread * 1.4 if spread > 0 else amp * 1e-12
    ax.plot(zz.real, zz.imag, 'o', ms=4, color='tab:red', alpha=0.8, linestyle='none')
    cnt = Counter((round(float(w.real), 15), round(float(w.imag), 15)) for w in zz)
    for (a, b), n in cnt.items():
        if n > 1:
            ax.annotate(f'x{n}', (a, b), textcoords='offset points', xytext=(4, 4),
                        fontsize=6, color='black')
    ax.set_xlim(c.real - win, c.real + win); ax.set_ylim(c.imag - win, c.imag + win)
    ax.set_aspect('equal')
    ax.grid(alpha=.25)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3), useOffset=True)
    ax.tick_params(labelsize=7)
    ax.set_title(f'cluster {k+1}: {len(zz)} waves, center=({c.real:+.6e},{c.imag:+.6e}),'
                 f' max dev={spread:.3e} (|z|max={amp:.3e})', fontsize=8)
    ax.set_xlabel('Re z', fontsize=8); ax.set_ylabel('Im z', fontsize=8)
for k in range(len(clusters), 4):
    axs[k].axis('off')
fig.suptitle(f'N=40 selfcontrol dtau=2pi/40 final step tau={STEPS}: zoom into angle clusters', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_N40_selfcontrol_complex_plane_final_cluster_zoom.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
