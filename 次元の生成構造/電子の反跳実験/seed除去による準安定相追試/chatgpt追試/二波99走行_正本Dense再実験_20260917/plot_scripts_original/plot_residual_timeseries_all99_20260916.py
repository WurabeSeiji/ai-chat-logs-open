#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 99 走行の残差時系列カタログ（2026-09-16）。

木原指示 2026-09-16: 同じ配列で、横軸=ステップ t=0..4096、縦軸に 4 系列:
  ① Re Σ_e z_e(t)²（二乗閉塞残差の実部）  = 赤実線
  ② Im Σ_e z_e(t)²（同 虚部）            = 青実線
  ③ Re Σ_e z_e(t)（重心の実部）           = 赤点線
  ④ Im Σ_e z_e(t)（重心の虚部）           = 青点線

スケーリング（案B・木原決定）: 正規化・無次元化はしない。全 99 走行・全ステップ・
全 4 系列にわたる最大絶対値 1 つを基準にした全図共通の生値スケール ±Y とする。
小さい走行が見えないのは「実際に小さい」ことの正しい表現とする
（状態による差をそのまま見せる。パネル別正規化・上界正規化は読者を混乱させるため不採用）。

予告どおり Σz²=Q₂ は厳密保存量（K 実反対称 ⟹ exp(ΔτK) 実直交）なので①②は
水平線になり、動くのは重心③④のみ — 保存/非保存の対比がそのまま見える。

- ページ 1: 主サーベイ 72 走行（行 L × 列 9 組: 奇奇/偶偶/奇偶）
- ページ 2: gcd 系列 18 ＋ den=40 対照 9
- データは runs/*/states.npz の Z 全 4097 状態を read-only 参照のみ（再走行しない）。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, SURVEY)
from gen_manifest import MAIN_L, MAIN_PAIRS, GCD_L, GCD_PAIRS, CTRL_L, CTRL_PAIRS, CTRL_DEN, T  # noqa: E402


def load_series(L, ma, mb, den):
    """全ステップの (Σz(t), Σz(t)²) を返す（read-only）。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
        Z = d['Z']
        assert Z.shape[0] == T + 1
        s1 = Z.sum(axis=1)              # 重心 Σz(t)
        s2 = (Z ** 2).sum(axis=1)       # 二乗閉塞残差 Σz(t)² = Q₂(t)
    return np.asarray(s1), np.asarray(s2)


def panel(ax, s1, s2, Y, label, xlabels=False):
    tt = np.arange(len(s1))
    ax.axhline(0, color='gray', lw=0.4, alpha=0.5)
    ax.plot(tt, s2.real, color='tab:red', lw=0.8, ls='-', alpha=0.9, zorder=3)
    ax.plot(tt, s2.imag, color='tab:blue', lw=0.8, ls='-', alpha=0.9, zorder=3)
    ax.plot(tt, s1.real, color='tab:red', lw=0.7, ls=':', alpha=0.9, zorder=4)
    ax.plot(tt, s1.imag, color='tab:blue', lw=0.7, ls=':', alpha=0.9, zorder=4)
    # 実線が点線帯に埋もれても水準が読めるよう t=0 に ×（木原指示 2026-09-16）
    ax.scatter([0], [s2[0].real], marker='x', s=34, c='tab:red',
               linewidths=1.4, zorder=6, clip_on=False)
    ax.scatter([0], [s2[0].imag], marker='x', s=34, c='tab:blue',
               linewidths=1.4, zorder=6, clip_on=False)
    ax.set_xlim(0, T)
    ax.set_ylim(-Y, Y)
    ax.set_yticks([0])                    # 縦軸の原点位置の目盛（木原指示 2026-09-16）
    ax.set_yticklabels(['0'], fontsize=6)
    ax.set_xticks([0, T // 2, T])
    ax.set_xticklabels(['0', str(T // 2), str(T)] if xlabels else [], fontsize=6)
    ax.tick_params(axis='x', length=1.5)
    ax.tick_params(axis='y', length=3)
    ax.set_title(label, fontsize=8, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)


def legend_handles():
    return [
        Line2D([], [], color='tab:red', ls='-', lw=1.3, label='Re Σz²（二乗閉塞残差・実部）'),
        Line2D([], [], color='tab:blue', ls='-', lw=1.3, label='Im Σz²（二乗閉塞残差・虚部）'),
        Line2D([], [], color='tab:red', ls=':', lw=1.3, label='Re Σz（重心・実部）'),
        Line2D([], [], color='tab:blue', ls=':', lw=1.3, label='Im Σz（重心・虚部）'),
        Line2D([], [], color='tab:red', marker='x', ls='none', markersize=7,
               markeredgewidth=1.5, label='Re Σz²(0) の水準（t=0 マーク）'),
        Line2D([], [], color='tab:blue', marker='x', ls='none', markersize=7,
               markeredgewidth=1.5, label='Im Σz²(0) の水準（t=0 マーク）'),
    ]


def main():
    data = {}
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            data[(L, ma, mb, L)] = load_series(L, ma, mb, L)
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            data[(L, ma, mb, L)] = load_series(L, ma, mb, L)
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            data[(L, ma, mb, CTRL_DEN)] = load_series(L, ma, mb, CTRL_DEN)

    peak = {key: max(float(np.max(np.abs(s2.real))), float(np.max(np.abs(s2.imag))),
                     float(np.max(np.abs(s1.real))), float(np.max(np.abs(s1.imag))))
            for key, (s1, s2) in data.items()}
    kmax = max(peak, key=peak.get)
    Ymax = peak[kmax]
    Y = Ymax * 1.05
    kmin = min(peak, key=peak.get)
    print(f'全体最大（全 run・全ステップ・4 系列）: {Ymax:.4f} '
          f'(L{kmax[0]}_ma{kmax[1]}_mb{kmax[2]}_den{kmax[3]}) → 共通縦スケール ±{Y:.4f}')
    print(f'パネル内最大が最小の run: {peak[kmin]:.4f} '
          f'(L{kmin[0]}_ma{kmin[1]}_mb{kmin[2]}_den{kmin[3]}) — 全体スケールの '
          f'{100 * peak[kmin] / Ymax:.1f}%（見えないのは実際に小さいため）')

    # ページ 1: 主サーベイ 72（8 行 × 9 列）
    fig, axs = plt.subplots(len(MAIN_L), len(MAIN_PAIRS), figsize=(26.0, 13.5))
    for i, L in enumerate(MAIN_L):
        for j, (ma, mb) in enumerate(MAIN_PAIRS):
            s1, s2 = data[(L, ma, mb, L)]
            panel(axs[i, j], s1, s2, Y, f'L{L} ({ma},{mb})',
                  xlabels=(i == len(MAIN_L) - 1))
        axs[i, 0].set_ylabel(f'L={L}', fontsize=10)
    for j, (ma, mb) in enumerate(MAIN_PAIRS):
        cls = 'odd-odd' if j < 3 else ('even-even' if j < 6 else 'odd-even')
        axs[0, j].set_title(f'({ma},{mb})  {cls}\nL{MAIN_L[0]}', fontsize=8, pad=2)
    fig.suptitle('残差時系列 t=0..4096（主サーベイ 72 走行、den=L）— 生値・全図共通スケール '
                 f'±{Y:.1f}（全 99 走行・全ステップ・4 系列の最大 {Ymax:.1f}）。'
                 '実線=Σz²（赤=Re・青=Im、保存量 Q₂ ゆえ水平線）、点線=重心 Σz（赤=Re・青=Im）',
                 fontsize=12, y=0.995)
    fig.legend(handles=legend_handles(), loc='lower center', ncol=6, fontsize=9,
               frameon=False)
    fig.tight_layout(rect=[0, 0.028, 1, 0.972])
    out1 = os.path.join(HERE, 'fig_residual_timeseries_all99_page1_main72_20260916.png')
    fig.savefig(out1, dpi=130)
    plt.close(fig)
    print('page1 ->', os.path.basename(out1))

    # ページ 2: gcd 18（6 行 × 3 列）＋ den=40 対照 9（3 行 × 3 列）
    fig = plt.figure(figsize=(10.0, 16.5))
    gs = fig.add_gridspec(9, 3, hspace=0.55, wspace=0.10,
                          top=0.945, bottom=0.045, left=0.05, right=0.97)
    for i, L in enumerate(GCD_L):
        for j, (ma, mb) in enumerate(GCD_PAIRS):
            ax = fig.add_subplot(gs[i, j])
            s1, s2 = data[(L, ma, mb, L)]
            panel(ax, s1, s2, Y, f'L{L} ({ma},{mb}) den={L}')
    for i, L in enumerate(CTRL_L):
        for j, (ma, mb) in enumerate(CTRL_PAIRS):
            ax = fig.add_subplot(gs[6 + i, j])
            s1, s2 = data[(L, ma, mb, CTRL_DEN)]
            panel(ax, s1, s2, Y, f'L{L} ({ma},{mb}) den={CTRL_DEN}',
                  xlabels=(i == len(CTRL_L) - 1))
    fig.suptitle('残差時系列 t=0..4096（gcd 系列 18 走行・上 6 行 ＋ den=40 対照 9 走行・下 3 行）\n'
                 f'生値・全図共通スケール ±{Y:.1f}。実線=Σz²（赤=Re・青=Im、保存量）、'
                 '点線=重心 Σz（赤=Re・青=Im）', fontsize=10, y=0.99)
    fig.legend(handles=legend_handles(), loc='lower center', ncol=2, fontsize=8,
               frameon=False)
    out2 = os.path.join(HERE, 'fig_residual_timeseries_all99_page2_gcd18_ctrl9_20260916.png')
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print('page2 ->', os.path.basename(out2))
    print(f'panels total: {len(data)}')

    q2flat = max(float(np.max(np.abs(s2 - s2[0]))) for (_, s2) in data.values())
    g_move = max(float(np.max(np.abs(s1 - s1[0]))) for (s1, _) in data.values())
    print(f'Σz² の時間変動 max|Q2(t)−Q2(0)| = {q2flat:.3e}（保存の確認）')
    print(f'重心の時間変動 max|Σz(t)−Σz(0)| = {g_move:.4f}（非保存の確認）')


if __name__ == '__main__':
    main()
