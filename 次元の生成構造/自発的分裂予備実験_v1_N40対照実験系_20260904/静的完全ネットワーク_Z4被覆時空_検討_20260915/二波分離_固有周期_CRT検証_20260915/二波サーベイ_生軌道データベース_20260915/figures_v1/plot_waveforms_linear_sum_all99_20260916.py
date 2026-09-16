#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 99 走行の M 本辺波の完全線形合成波形カタログ（2026-09-16、図その2）。

定義（木原指示 2026-09-16）:
  各辺 e（差分 Δ_e、x_Δ=2πΔ/L）は生成波の平行移動コピーを 1 本の波として担う:
    w_e(x) = Re[A_e e^{i m_a x} + B_e e^{i m_b x}] = Re[z(x + x_Δ)]
    （A_e = c_a U^{m_aΔ}, B_e = c_b U^{m_bΔ} は初期状態カタログ v2 のフェーザ）
  全 M 本の完全線形干渉:
    W(x) = Σ_e w_e(x) = Re[(Σ_e A_e) e^{i m_a x} + (Σ_e B_e) e^{i m_b x}]
  Σ_e A_e, Σ_e B_e は M 辺の直接和で計算する（閉形式に依存しない）。

表示:
  - 横軸 x ∈ [−π, π]（差分座標 2πΔ/L の連続化、1 周期）
  - 細線 = 元の波形（振幅 1・真スケール）: ベース波 cos(m_a x−π/2L) 実線、
    倍音 cos(m_b x+π/2L) 点線。色は奇数倍音=青 / 偶数倍音=赤（v2 と同一）
  - 太黒線 = 合成波形 W(x)。振幅が O(L²)（数百）に達するため、この線のみ
    パネル枠最大に正規化し、真の振幅 A=max|W| をパネル右上に記載
  - 縦スケールは全 99 走行共通 ±S（S = 実測最大 |z(0)| × 1.05、v1/v2 と同一）
  - 配列は v1/v2 と同一（ページ1: 主 72、ページ2: gcd 18 ＋ den=40 対照 9）

データは runs/*/states.npz の Z[0] を read-only 参照のみ。成分分解 A_e+B_e が
データ z(0) と一致することを毎 run 検証する（再走行しない・書き込まない）。
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
from gen_manifest import MAIN_L, MAIN_PAIRS, GCD_L, GCD_PAIRS, CTRL_L, CTRL_PAIRS, CTRL_DEN  # noqa: E402

COL_ODD = 'tab:blue'    # 奇数倍音
COL_EVEN = 'tab:red'    # 偶数倍音
X = np.linspace(-np.pi, np.pi, 2001)


def pcol(m):
    return COL_ODD if m % 2 == 1 else COL_EVEN


def load_z0(L, ma, mb, den):
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
        return np.array(d['Z'][0])


def edge_phasor_sums(L, ma, mb, z0):
    """M 辺のフェーザ A_e, B_e を直接和し (Σ_e A_e, Σ_e B_e) を返す。
    A_e+B_e = データ z(0) の bit 一致を検証（read-only 整合性チェック）。"""
    j, k = np.triu_indices(L, 1)
    delta = (k - j).astype(np.int64)
    U = np.exp(2j * np.pi / L)
    A = np.exp(-1j * np.pi / (2 * L)) * U ** (ma * delta)
    B = np.exp(+1j * np.pi / (2 * L)) * U ** (mb * delta)
    resid = float(np.max(np.abs(A + B - z0)))
    return complex(np.sum(A)), complex(np.sum(B)), resid


def panel(ax, L, ma, mb, z0, S, label, xlabels=False):
    Sa, Sb, resid = edge_phasor_sums(L, ma, mb, z0)
    base = np.cos(ma * X - np.pi / (2 * L))       # 元のベース波（振幅 1）
    harm = np.cos(mb * X + np.pi / (2 * L))       # 元の倍音（振幅 1）
    W = (Sa * np.exp(1j * ma * X) + Sb * np.exp(1j * mb * X)).real
    Wmax = float(np.max(np.abs(W)))
    k = 0.97 * S / Wmax                            # 合成のみ枠内に正規化
    ax.axhline(0, color='gray', lw=0.4, alpha=0.5)
    ax.plot(X, base, color=pcol(ma), lw=0.7, ls='-', alpha=0.85, zorder=2)
    ax.plot(X, harm, color=pcol(mb), lw=0.7, ls=':', alpha=0.9, zorder=3)
    ax.plot(X, W * k, color='black', lw=1.5, alpha=0.9, zorder=4)
    ax.annotate(f'A={Wmax:.1f}', (0.98, 0.96), xycoords='axes fraction',
                ha='right', va='top', fontsize=5.5, color='black')
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-S, S)
    ax.set_yticks([])
    ax.set_xticks([-np.pi, 0, np.pi])
    ax.set_xticklabels(['−π', '0', 'π'] if xlabels else [], fontsize=6)
    ax.tick_params(length=1.5)
    ax.set_title(label, fontsize=8, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)
    return resid, Wmax


def legend_handles():
    return [
        Line2D([], [], color=COL_ODD, ls='-', lw=1.2, label='ベース波 m_a（奇数倍音・振幅1）'),
        Line2D([], [], color=COL_EVEN, ls='-', lw=1.2, label='ベース波 m_a（偶数倍音・振幅1）'),
        Line2D([], [], color=COL_ODD, ls=':', lw=1.2, label='倍音 m_b（奇数倍音・振幅1）'),
        Line2D([], [], color=COL_EVEN, ls=':', lw=1.2, label='倍音 m_b（偶数倍音・振幅1）'),
        Line2D([], [], color='black', lw=1.8,
               label='全 M 辺波の線形合成 W(x)（枠に正規化、A=真の振幅）'),
    ]


def main():
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
    S = max(float(np.max(np.abs(z))) for z in data.values()) * 1.05
    print(f'共通縦スケール: ±{S:.4f}（実測最大 |z(0)|={S/1.05:.6f}）')
    worst_resid = 0.0
    amps = []

    # ページ 1: 主サーベイ 72（8 行 × 9 列）
    fig, axs = plt.subplots(len(MAIN_L), len(MAIN_PAIRS), figsize=(26.0, 13.5))
    for i, L in enumerate(MAIN_L):
        for j, (ma, mb) in enumerate(MAIN_PAIRS):
            r, w = panel(axs[i, j], L, ma, mb, data[(L, ma, mb, L)], S,
                         f'L{L} ({ma},{mb})', xlabels=(i == len(MAIN_L) - 1))
            worst_resid = max(worst_resid, r)
            amps.append((f'L{L}_({ma},{mb})', w))
        axs[i, 0].set_ylabel(f'L={L}', fontsize=10)
    for j, (ma, mb) in enumerate(MAIN_PAIRS):
        cls = 'odd-odd' if j < 3 else ('even-even' if j < 6 else 'odd-even')
        axs[0, j].set_title(f'({ma},{mb})  {cls}\nL{MAIN_L[0]}', fontsize=8, pad=2)
    fig.suptitle('全 M 辺波の完全線形合成波形（主サーベイ 72 走行、den=L）— 各辺 e は生成波の'
                 '平行移動コピー Re[z(x+2πΔ/L)] を担う。太黒線=合成 W(x)（枠に正規化・A=真の振幅）、'
                 f'細線=元のベース波/倍音（振幅1・共通縦スケール ±{S:.3f}）', fontsize=12, y=0.995)
    fig.legend(handles=legend_handles(), loc='lower center', ncol=5, fontsize=9,
               frameon=False)
    fig.tight_layout(rect=[0, 0.028, 1, 0.972])
    out1 = os.path.join(HERE, 'fig_waveforms_linear_sum_all99_page1_main72_20260916.png')
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
            r, w = panel(ax, L, ma, mb, data[(L, ma, mb, L)], S,
                         f'L{L} ({ma},{mb}) den={L}')
            worst_resid = max(worst_resid, r)
            amps.append((f'L{L}_({ma},{mb})_gcd', w))
    for i, L in enumerate(CTRL_L):
        for j, (ma, mb) in enumerate(CTRL_PAIRS):
            ax = fig.add_subplot(gs[6 + i, j])
            r, w = panel(ax, L, ma, mb, data[(L, ma, mb, CTRL_DEN)], S,
                         f'L{L} ({ma},{mb}) den={CTRL_DEN}',
                         xlabels=(i == len(CTRL_L) - 1))
            worst_resid = max(worst_resid, r)
            amps.append((f'L{L}_({ma},{mb})_den40', w))
    fig.suptitle('全 M 辺波の完全線形合成波形（gcd 系列 18 走行・上 6 行 ＋ den=40 対照 9 走行・下 3 行）\n'
                 f'太黒線=合成 W(x)（枠に正規化・A=真の振幅）、細線=元の波形（振幅1・±{S:.3f}）',
                 fontsize=10, y=0.99)
    fig.legend(handles=legend_handles(), loc='lower center', ncol=2, fontsize=8,
               frameon=False)
    out2 = os.path.join(HERE, 'fig_waveforms_linear_sum_all99_page2_gcd18_ctrl9_20260916.png')
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print('page2 ->', os.path.basename(out2))

    print(f'panels total: {len(amps)}')
    print(f'成分分解とデータ z(0) の最大残差（全 99 run）: {worst_resid:.3e}')
    amps.sort(key=lambda t: t[1])
    print(f'合成振幅 A の範囲: 最小 {amps[0][1]:.2f} ({amps[0][0]}) / '
          f'最大 {amps[-1][1]:.2f} ({amps[-1][0]})')


if __name__ == '__main__':
    main()
