#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 99 走行の初期状態（t=0）複素平面カタログ（2026-09-16）。

- 各走行の t=0 の全関係波 z_e(0) を 1 パネルの複素平面に散布。
- 振幅は正規化しない。全 99 走行の実測最大 |z_e(0)| に合わせた
  共通絶対スケール（±S）で全パネルを表示する。
- データは runs/*/states.npz を read-only 参照するのみ（再走行しない、書き込まない）。
- ページ 1: 主サーベイ 72 走行 = 行 L∈{8,10,12,16,20,24,32,40} × 列 9 組
  （奇奇 (1,3)(1,5)(3,5) / 偶偶 (2,4)(2,6)(4,6) / 奇偶 (1,2)(2,3)(3,6)）。
- ページ 2: gcd 系列 18 走行（行 L∈{12..40} × 列 (3,9)(4,8)(5,10)）＋
  den=40 対照 9 走行（行 L∈{12,24,32} × 列 (1,3)(2,4)(2,3)）。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, SURVEY)
from gen_manifest import MAIN_L, MAIN_PAIRS, GCD_L, GCD_PAIRS, CTRL_L, CTRL_PAIRS, CTRL_DEN  # noqa: E402


def load_z0(L, ma, mb, den):
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
        return np.array(d['Z'][0])


def collect():
    data = {}
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            data[(L, ma, mb, L)] = load_z0(L, ma, mb, L)
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            data[(L, ma, mb, L)] = load_z0(L, ma, mb, L)
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            data[(L, ma, mb, CTRL_DEN)] = load_z0(L, ma, mb, CTRL_DEN)
    return data


def panel(ax, z, S, label):
    ax.axhline(0, color='gray', lw=0.4, alpha=0.6)
    ax.axvline(0, color='gray', lw=0.4, alpha=0.6)
    ax.scatter(z.real, z.imag, s=3, c='tab:blue', alpha=0.7, linewidths=0)
    ax.set_xlim(-S, S)
    ax.set_ylim(-S, S)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(label, fontsize=7, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)


def main():
    data = collect()
    S = max(float(np.max(np.abs(z))) for z in data.values()) * 1.05
    Smax = S / 1.05
    print(f'共通絶対スケール: max|z_e(0)| = {Smax:.6f} → 表示範囲 ±{S:.4f}')

    # ページ 1: 主サーベイ 72（8 行 × 9 列）
    fig, axs = plt.subplots(len(MAIN_L), len(MAIN_PAIRS), figsize=(16.5, 15.5))
    for i, L in enumerate(MAIN_L):
        for j, (ma, mb) in enumerate(MAIN_PAIRS):
            panel(axs[i, j], data[(L, ma, mb, L)], S, f'L{L} ({ma},{mb})')
        axs[i, 0].set_ylabel(f'L={L}', fontsize=9)
    for j, (ma, mb) in enumerate(MAIN_PAIRS):
        cls = 'odd-odd' if j < 3 else ('even-even' if j < 6 else 'odd-even')
        axs[0, j].set_title(f'({ma},{mb})\n{cls}\nL{MAIN_L[0]}', fontsize=7, pad=2)
    fig.suptitle('初期状態 t=0（主サーベイ 72 走行、den=L）— 正規化なし・共通絶対スケール '
                 f'±{S:.3f}（全99走行の実測最大 |z|={Smax:.4f}）', fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    out1 = os.path.join(HERE, 'fig_initial_states_all99_page1_main72_20260916.png')
    fig.savefig(out1, dpi=170)
    plt.close(fig)
    print('page1 ->', os.path.basename(out1))

    # ページ 2: gcd 18（6 行 × 3 列）＋ den=40 対照 9（3 行 × 3 列）
    fig = plt.figure(figsize=(6.8, 19.0))
    gs = fig.add_gridspec(9, 3, hspace=0.32, wspace=0.12,
                          top=0.955, bottom=0.012, left=0.05, right=0.97)
    for i, L in enumerate(GCD_L):
        for j, (ma, mb) in enumerate(GCD_PAIRS):
            ax = fig.add_subplot(gs[i, j])
            panel(ax, data[(L, ma, mb, L)], S, f'L{L} ({ma},{mb}) den={L}')
    for i, L in enumerate(CTRL_L):
        for j, (ma, mb) in enumerate(CTRL_PAIRS):
            ax = fig.add_subplot(gs[6 + i, j])
            panel(ax, data[(L, ma, mb, CTRL_DEN)], S,
                  f'L{L} ({ma},{mb}) den={CTRL_DEN}')
    fig.suptitle('初期状態 t=0（gcd 系列 18 走行・上 6 行 ＋ den=40 対照 9 走行・下 3 行）\n'
                 f'正規化なし・共通絶対スケール ±{S:.3f}', fontsize=10, y=0.995)
    out2 = os.path.join(HERE, 'fig_initial_states_all99_page2_gcd18_ctrl9_20260916.png')
    fig.savefig(out2, dpi=170)
    plt.close(fig)
    print('page2 ->', os.path.basename(out2))
    print(f'panels total: {len(data)}')


if __name__ == '__main__':
    main()
