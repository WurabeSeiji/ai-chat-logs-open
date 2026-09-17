#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位相・振幅の時系列図・1 run = 1 図（2026-09-17 木原指示）。

木原指示: 「step10 までのステップを横軸に、縦軸は2つの値。左側が位相 −π..0..π、
右側が振幅（全体最大振幅で正規化）。波毎に、位相を実線、振幅を点線で表示。」

目的: 複素平面図では位相と振幅が混ざって回転が読めないため、各波の位相前進
（傾き）を直接測り、どの波群が同じ回転率を共有するかを見る。振幅点線は
r→0 通過（位相が跳ぶ／未定義になる瞬間）を同じ横軸で照合する。

- 位相: arg z_e(t) ∈ (−π, π]（unwrap しない生値）、実線、辺ごと1色（turbo）
- 振幅: |z_e(t)| / max_run|z_e(t)|（全4097ステップの実測最大で正規化）、点線、同色
- 区間: t=RANGE_T0..RANGE_T1（既定 0..10、コマンドライン引数で変更可）
- データは states.npz の read-only 参照のみ。1 run = 1 図（存在する run だけ）。

使い方: python3 plot_phase_amplitude_timeseries_20260917.py [t0] [t1]
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
from gen_manifest import MAIN_L, MAIN_PAIRS, GCD_L, GCD_PAIRS, CTRL_L, CTRL_PAIRS, CTRL_DEN, T  # noqa: E402

RANGE_T0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
RANGE_T1 = int(sys.argv[2]) if len(sys.argv) > 2 else 10


def have(L, ma, mb, den):
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    return os.path.exists(os.path.join(SURVEY, 'runs', rid, 'states.npz'))


def all_keys():
    keys = []
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            keys.append((L, ma, mb, L))
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            keys.append((L, ma, mb, L))
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            keys.append((L, ma, mb, CTRL_DEN))
    return keys


def main():
    found = 0
    for (L, ma, mb, den) in all_keys():
        if not have(L, ma, mb, den):
            continue
        found += 1
        rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
        with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
            Z = np.array(d['Z'])
            assert Z.shape[0] == T + 1
        M = Z.shape[1]
        ja, jb = np.triu_indices(L, k=1)
        amax_run = float(np.max(np.abs(Z)))          # 全体最大振幅（全4097ステップ）
        Zs = Z[RANGE_T0:RANGE_T1 + 1]
        tt = np.arange(RANGE_T0, RANGE_T1 + 1)
        cols = plt.cm.turbo(np.linspace(0.02, 0.98, M))

        fig, ax = plt.subplots(figsize=(11.0, 7.5))
        ax2 = ax.twinx()
        for e in range(M):
            ax.plot(tt, np.angle(Zs[:, e]), color=cols[e], lw=0.9, ls='-',
                    alpha=0.85, marker='o', markersize=2.2, zorder=3,
                    label=f'({ja[e]},{jb[e]})')
            ax2.plot(tt, np.abs(Zs[:, e]) / amax_run, color=cols[e], lw=0.9,
                     ls=':', alpha=0.85, marker='o', markersize=2.2, zorder=2)
        ax.set_xlim(RANGE_T0, RANGE_T1)
        ax.set_xticks(tt)
        ax.set_xlabel('step', fontsize=11)
        ax.set_ylim(-np.pi * 1.04, np.pi * 1.04)
        ax.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
        ax.set_yticklabels(['−π', '−π/2', '0', 'π/2', 'π'], fontsize=10)
        ax.set_ylabel('位相 arg z_e(t)（実線）', fontsize=11)
        ax.axhline(0, color='gray', lw=0.5, alpha=0.5)
        ax2.set_ylim(0, 1.04)
        ax2.set_ylabel(f'振幅 |z_e(t)| / {amax_run:.4f}（点線、全体最大で正規化）',
                       fontsize=11)
        ax.set_title(f'位相・振幅の時系列 t={RANGE_T0}..{RANGE_T1}（{rid}）— '
                     f'波毎に1色（実線=位相・点線=正規化振幅）', fontsize=12)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=7,
                  fontsize=7, frameon=False, title='辺 (j,k)', title_fontsize=8)
        fig.tight_layout(rect=[0, 0.10, 1, 1])
        outdir = os.path.join(SURVEY, 'plots', 'phase_amplitude')
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir,
                           f'fig_phase_amplitude_t{RANGE_T0}_{RANGE_T1}_{rid}_dense_rerun_20260917.png')
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f'{rid}: 全体最大振幅 {amax_run:.6f} -> {os.path.basename(out)}')
    if found == 0:
        print('states.npz が存在する run が無い — 図化せず終了')
    else:
        print(f'存在する run: {found}/99')


if __name__ == '__main__':
    main()
