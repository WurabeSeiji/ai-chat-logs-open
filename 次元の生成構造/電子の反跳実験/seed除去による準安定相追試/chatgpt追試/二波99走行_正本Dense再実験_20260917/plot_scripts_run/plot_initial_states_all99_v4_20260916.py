#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 99 走行の初期状態（t=0）複素平面カタログ v4（2026-09-16、二乗閉塞残差追加版）。

v1 からの表示改良（木原指示 2026-09-16）:
 1. 同一値に縮退した辺の重複数 n で ◯（合成点 z の円）の面積を変える
 2. ◯ の横に重複辺数 n を整数で書く
 3. 各点を中心と細線で結ぶ。z(0)=c_a U^{m_aΔ}+c_b U^{m_bΔ} は 2 モードの和で
    あり単独の倍音に帰属しないため、「中心と結ぶ」はベクトル合成路
    原点 →(実線: ベース波成分 c_a U^{m_aΔ})→ 合成点 z（点線: 倍音成分）で実装
 4. 色: 奇数倍音=青 / 偶数倍音=赤（ベース波 m_a・倍音 m_b それぞれの偶奇で着色）
    線種: ベース波=実線 / 倍音=点線。ベース波成分の先端に小さな色点を打つ。
 5. v3 追加（木原指示 2026-09-16）: 各走行の M 波全ての重心 Σ_e z_e を緑×でマーク。
    計算式は高対称系列の verify() と同一（Z.sum()）。生の |Σz| は O(L²)（最大〜340）で
    枠外に飛ぶため、振幅スケールと同じ思想で「最も重心がズレている走行の max|Σz|」を
    基準に全 99 走行共通の線形圧縮で表示する（×位置 = Σz × 0.97S/max|Σz|）。
    真の |Σz| は各パネル左下に G=値 で記載。
 6. v4 追加（木原指示 2026-09-16）: この系は二乗閉塞 Σ_e z_e² が自動では閉じない
    （なお Σz² は保存量 Q₂ そのもの）ため、残差を外枠上に橙×で表示する:
    下枠の X 方向位置 = Re Σz²、左枠の Y 方向位置 = Im Σz²。二乗のままの値を使い、
    実部は全 99 走行の max|Re Σz²|、虚部は max|Im Σz²| をそれぞれ枠最大に対応させる
    全図共通の線形圧縮（重心×と同じ思想）。真値を各パネル右下に Q2=(Re,Im) で記載。

- データは runs/*/states.npz の Z[0] を read-only 参照のみ（再走行しない）。
  成分分解は初期条件式から計算し、和がデータの z(0) と一致することを毎 run 検証。
- 振幅正規化なし・全 99 走行共通の絶対スケール ±S（S=実測最大|z|×1.05）。
- ページ 1: 主サーベイ 72 走行（行 L × 列 9 組: 奇奇/偶偶/奇偶）
- ページ 2: gcd 系列 18 走行 ＋ den=40 対照 9 走行
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


def pcol(m):
    return COL_ODD if m % 2 == 1 else COL_EVEN


def load_z0(L, ma, mb, den):
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
        return np.array(d['Z'][0])


def decompose_groups(L, ma, mb, z0):
    """初期条件式による 2 成分分解と、同一値クラス（(m_aΔ, m_bΔ) mod L）ごとの
    重複辺数。和 = データ z(0) を検証して残差を返す。"""
    j, k = np.triu_indices(L, 1)
    delta = (k - j).astype(np.int64)
    U = np.exp(2j * np.pi / L)
    A = np.exp(-1j * np.pi / (2 * L)) * U ** (ma * delta)
    B = np.exp(+1j * np.pi / (2 * L)) * U ** (mb * delta)
    resid = float(np.max(np.abs(A + B - z0)))
    groups = {}
    for e in range(len(delta)):
        key = (int(ma * delta[e] % L), int(mb * delta[e] % L))
        groups.setdefault(key, []).append(e)
    assert sum(len(v) for v in groups.values()) == len(delta)
    return A, groups, resid


def panel(ax, L, ma, mb, z0, S, label, gdisp, gtrue, q2x, q2y, q2true):
    ax.axhline(0, color='gray', lw=0.4, alpha=0.5)
    ax.axvline(0, color='gray', lw=0.4, alpha=0.5)
    A, groups, resid = decompose_groups(L, ma, mb, z0)
    ca_col, cb_col = pcol(ma), pcol(mb)
    for key in sorted(groups):
        idx = groups[key]
        n = len(idx)
        a = A[idx[0]]
        z = z0[idx[0]]                       # 合成点はデータ値を使用
        ax.plot([0, a.real], [0, a.imag], color=ca_col, lw=0.5, ls='-',
                alpha=0.55, zorder=1)        # ベース波 m_a: 実線
        ax.plot([a.real, z.real], [a.imag, z.imag], color=cb_col, lw=0.5,
                ls=':', alpha=0.7, zorder=2)  # 倍音 m_b: 点線
        ax.scatter([a.real], [a.imag], s=4, c=ca_col, linewidths=0, zorder=3)
        s = float(np.clip(2.5 * n, 8, 300))  # ◯ 面積 ∝ 重複辺数（上限 300）
        ax.scatter([z.real], [z.imag], s=s, facecolors='white',
                   edgecolors='black', linewidths=0.7, zorder=4)
        ax.annotate(str(n), (z.real, z.imag), xytext=(3, 2),
                    textcoords='offset points', fontsize=4.5,
                    color='dimgray', zorder=5)
    ax.scatter([gdisp.real], [gdisp.imag], marker='x', s=42, c='tab:green',
               linewidths=1.4, zorder=6)
    ax.annotate(f'G={gtrue:.1f}', (0.03, 0.03), xycoords='axes fraction',
                ha='left', va='bottom', fontsize=5, color='tab:green', zorder=6)
    ax.scatter([q2x], [-0.94 * S], marker='x', s=30, c='tab:orange',
               linewidths=1.2, zorder=6)                 # 下枠: Re Σz²
    ax.scatter([-0.94 * S], [q2y], marker='x', s=30, c='tab:orange',
               linewidths=1.2, zorder=6)                 # 左枠: Im Σz²
    ax.annotate(f'Q2=({q2true.real:.1f},{q2true.imag:.1f})', (0.97, 0.03),
                xycoords='axes fraction', ha='right', va='bottom',
                fontsize=4.5, color='tab:orange', zorder=6)
    ax.set_xlim(-S, S)
    ax.set_ylim(-S, S)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(label, fontsize=8, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)
    return resid


def legend_handles():
    return [
        Line2D([], [], color=COL_ODD, ls='-', lw=1.2, label='ベース波 m_a（奇数倍音）'),
        Line2D([], [], color=COL_EVEN, ls='-', lw=1.2, label='ベース波 m_a（偶数倍音）'),
        Line2D([], [], color=COL_ODD, ls=':', lw=1.2, label='倍音 m_b（奇数倍音）'),
        Line2D([], [], color=COL_EVEN, ls=':', lw=1.2, label='倍音 m_b（偶数倍音）'),
        Line2D([], [], color='black', marker='o', ls='none', markersize=7,
               markerfacecolor='white', label='z(0)（面積∝重複辺数、数字=辺数）'),
        Line2D([], [], color='tab:green', marker='x', ls='none', markersize=8,
               markeredgewidth=1.6,
               label='重心 Σz（max|Σz| 基準の全図共通圧縮表示、G=真の|Σz|）'),
        Line2D([], [], color='tab:orange', marker='x', ls='none', markersize=7,
               markeredgewidth=1.4,
               label='二乗閉塞残差 Σz²（下枠=Re・左枠=Im、成分別 max 基準の共通圧縮、Q2=真値）'),
    ]


def have(L, ma, mb, den):
    """段階実行対応（木原承認 2026-09-17）: states.npz が存在する run だけ図化する。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    return os.path.exists(os.path.join(SURVEY, 'runs', rid, 'states.npz'))


def main():
    data = {}
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            if have(L, ma, mb, L):
                data[(L, ma, mb, L)] = load_z0(L, ma, mb, L)
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            if have(L, ma, mb, L):
                data[(L, ma, mb, L)] = load_z0(L, ma, mb, L)
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            if have(L, ma, mb, CTRL_DEN):
                data[(L, ma, mb, CTRL_DEN)] = load_z0(L, ma, mb, CTRL_DEN)
    if not data:
        print('states.npz が存在する run が無い — 図化せず終了')
        return
    print(f'存在する run: {len(data)}/99（欠損パネルは空欄表示）')
    S = max(float(np.max(np.abs(z))) for z in data.values()) * 1.05
    Smax = S / 1.05
    print(f'共通絶対スケール: max|z_e(0)| = {Smax:.6f} → 表示範囲 ±{S:.4f}')
    Gs = {key: complex(np.sum(z)) for key, z in data.items()}
    Gmax = max(abs(g) for g in Gs.values())
    kG = 0.97 * S / Gmax
    print(f'重心 max|Σz| = {Gmax:.4f} → ×共通圧縮係数 {kG:.6e}')
    Q2s = {key: complex(np.sum(z ** 2)) for key, z in data.items()}
    Q2re_max = max(abs(q.real) for q in Q2s.values())
    Q2im_max = max(abs(q.imag) for q in Q2s.values())
    kRe = 0.97 * S / Q2re_max
    kIm = 0.97 * S / Q2im_max
    print(f'二乗閉塞残差 max|Re Σz²| = {Q2re_max:.4f} / max|Im Σz²| = {Q2im_max:.4f} '
          f'→ 圧縮係数 Re {kRe:.6e} / Im {kIm:.6e}')
    worst_resid = 0.0

    # ページ 1: 主サーベイ 72（8 行 × 9 列）
    fig, axs = plt.subplots(len(MAIN_L), len(MAIN_PAIRS), figsize=(24.0, 21.5))
    for i, L in enumerate(MAIN_L):
        for j, (ma, mb) in enumerate(MAIN_PAIRS):
            key = (L, ma, mb, L)
            if key in data:
                r = panel(axs[i, j], L, ma, mb, data[key], S,
                          f'L{L} ({ma},{mb})', Gs[key] * kG, abs(Gs[key]),
                          Q2s[key].real * kRe, Q2s[key].imag * kIm, Q2s[key])
                worst_resid = max(worst_resid, r)
            else:
                axs[i, j].axis('off')
        axs[i, 0].set_ylabel(f'L={L}', fontsize=10)
    for j, (ma, mb) in enumerate(MAIN_PAIRS):
        cls = 'odd-odd' if j < 3 else ('even-even' if j < 6 else 'odd-even')
        axs[0, j].set_title(f'({ma},{mb})  {cls}\nL{MAIN_L[0]}', fontsize=8, pad=2)
    fig.suptitle('初期状態 t=0 v4（主サーベイ 72 走行、den=L）— 正規化なし・共通絶対スケール '
                 f'±{S:.3f}。実線=ベース波成分、点線=倍音成分（合成路）、'
                 '青=奇数倍音・赤=偶数倍音、◯面積∝重複辺数、'
                 f'緑×=重心 Σz（max|Σz|={Gmax:.1f} 共通圧縮）、'
                 f'橙×=Σz² 残差（下枠=Re/max {Q2re_max:.1f}・左枠=Im/max {Q2im_max:.1f}）',
                 fontsize=12, y=0.996)
    fig.legend(handles=legend_handles(), loc='lower center', ncol=4, fontsize=9,
               frameon=False)
    fig.tight_layout(rect=[0, 0.034, 1, 0.978])
    out1 = os.path.join(HERE, '..', 'plots', 'initial_complex', 'fig_initial_states_all99_v4_page1_main72_dense_rerun_20260917.png')
    fig.savefig(out1, dpi=130)
    plt.close(fig)
    print('page1 ->', os.path.basename(out1))

    # ページ 2: gcd 18（6 行 × 3 列）＋ den=40 対照 9（3 行 × 3 列）
    fig = plt.figure(figsize=(8.8, 25.0))
    gs = fig.add_gridspec(9, 3, hspace=0.30, wspace=0.14,
                          top=0.962, bottom=0.028, left=0.05, right=0.97)
    for i, L in enumerate(GCD_L):
        for j, (ma, mb) in enumerate(GCD_PAIRS):
            key = (L, ma, mb, L)
            if key in data:
                ax = fig.add_subplot(gs[i, j])
                r = panel(ax, L, ma, mb, data[key], S,
                          f'L{L} ({ma},{mb}) den={L}', Gs[key] * kG, abs(Gs[key]),
                          Q2s[key].real * kRe, Q2s[key].imag * kIm, Q2s[key])
                worst_resid = max(worst_resid, r)
    for i, L in enumerate(CTRL_L):
        for j, (ma, mb) in enumerate(CTRL_PAIRS):
            key = (L, ma, mb, CTRL_DEN)
            if key in data:
                ax = fig.add_subplot(gs[6 + i, j])
                r = panel(ax, L, ma, mb, data[key], S,
                          f'L{L} ({ma},{mb}) den={CTRL_DEN}', Gs[key] * kG, abs(Gs[key]),
                          Q2s[key].real * kRe, Q2s[key].imag * kIm, Q2s[key])
                worst_resid = max(worst_resid, r)
    fig.suptitle('初期状態 t=0 v4（gcd 系列 18 走行・上 6 行 ＋ den=40 対照 9 走行・下 3 行）\n'
                 f'正規化なし・共通絶対スケール ±{S:.3f}。実線=ベース波・点線=倍音、青=奇数・赤=偶数\n'
                 f'◯面積∝重複辺数、緑×=重心 Σz（max|Σz|={Gmax:.1f} 共通圧縮）\n'
                 f'橙×=Σz² 残差（下枠=Re/max {Q2re_max:.1f}・左枠=Im/max {Q2im_max:.1f}）',
                 fontsize=9, y=0.995)
    fig.legend(handles=legend_handles(), loc='lower center', ncol=2, fontsize=7,
               frameon=False)
    out2 = os.path.join(HERE, '..', 'plots', 'initial_complex', 'fig_initial_states_all99_v4_page2_gcd18_ctrl9_dense_rerun_20260917.png')
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print('page2 ->', os.path.basename(out2))
    print(f'panels total: {len(data)}')
    print(f'成分分解とデータ z(0) の最大残差（全 99 run）: {worst_resid:.3e}')
    gm = sorted((abs(g), f'L{k[0]}_ma{k[1]}_mb{k[2]}_den{k[3]}') for k, g in Gs.items())
    print(f'重心 |Σz| 範囲: 最小 {gm[0][0]:.4f} ({gm[0][1]}) / '
          f'最大 {gm[-1][0]:.4f} ({gm[-1][1]})')
    qm = sorted((abs(q), f'L{k[0]}_ma{k[1]}_mb{k[2]}_den{k[3]}') for k, q in Q2s.items())
    print(f'二乗閉塞残差 |Σz²| 範囲: 最小 {qm[0][0]:.4f} ({qm[0][1]}) / '
          f'最大 {qm[-1][0]:.4f} ({qm[-1][1]})')


if __name__ == '__main__':
    main()
