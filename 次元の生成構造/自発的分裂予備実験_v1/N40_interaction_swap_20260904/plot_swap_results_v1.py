#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""相互作用スワップ実験の図化（読み出しのみ）。
(1) f(τ) 比較図: 旧力学（位相のみ sin(Δθ)・σ正規化 Cayley）vs 新相互作用（振幅込み H・Δτ=2π/40）。
    両走行の step0 は bit 同一（検証済み）。
(2) 新相互作用の終状態 τ=8000 の複素平面図（N40_state_readout と同一様式）。
(3) 縮退角クラスター拡大図（同上）。
step0 図は N40_state_readout_20260904/fig_N40_complex_plane_step0.png と bit 同一の状態なので再掲しない。"""
import csv
import os
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
OLD_DIR = os.path.join(os.path.dirname(BASE), 'N40_state_readout_20260904', 'largeN_splitting_result_v1')
NEW_DIR = os.path.join(BASE, 'largeN_splitting_result_v1')

def load_fcurve(path):
    taus, fs = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            taus.append(float(row['tau'])); fs.append(float(row['f']))
    return np.array(taus), np.array(fs)

# (1) f(τ) 比較
t_old, f_old = load_fcurve(os.path.join(OLD_DIR, 'fcurve_N00040_delta1e-15_seed0.csv'))
t_new, f_new = load_fcurve(os.path.join(NEW_DIR, 'fcurve_N00040_delta1e-15_seed0_den40.csv'))
fig, ax = plt.subplots(figsize=(8.2, 5.2))
ax.semilogy(t_old, f_old, color='#1f77b4', lw=1.6,
            label='old interaction: phase-only sin(dtheta), sigma-normalized Cayley (tau_end=3511, f=0.166)')
ax.semilogy(t_new, f_new, color='#d62728', lw=1.6,
            label='swapped interaction: amplitude-weighted H, dtau=2pi/40 (tau_end=8000, f=0.0062)')
ax.set_xlabel(r'$\tau$')
ax.set_ylabel(r'dormant fraction $f(\tau)$ (complement projection)')
ax.set_title('N=40, identical step0 (bit-exact): interaction swap only')
ax.grid(alpha=0.3, which='both')
ax.legend(loc='lower right', fontsize=8)
ax.set_ylim(1e-31, 2.0)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_N40_fcurve_interaction_swap_compare.png'), dpi=160)
plt.close(fig)

# (2) 終状態の複素平面図（N40_state_readout の plot_complex_plane_N40_v1.py と同一様式）
d = np.load(os.path.join(NEW_DIR, 'states_N00040_delta1e-15_seed0_den40.npz'))
ZF = np.asarray(d['Zfinal'], dtype=np.complex128)
TAU = int(d['tau_final'])
assert ZF.size == 780

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

draw_plane(ZF, f'N=40 (M=780) swapped interaction (dtau=2pi/40): final step tau={TAU} (f=0.0062)',
           'fig_N40_swap_complex_plane_final.png')

# (3) 縮退角クラスター拡大図（同一アルゴリズム: 粗解像度 1/100 群化、上位4クラスター）
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
fig.suptitle(f'N=40 swapped interaction final step tau={TAU}: zoom into angle clusters', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_N40_swap_complex_plane_final_cluster_zoom.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
