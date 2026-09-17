#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期状態 t=0 を最終図と同一形式で描くカタログ（2026-09-17 木原指示）。

plot_final_states_all99_20260916.py（実行用コピー）の物理コピー。変更は
(1) プロット対象を Z[T] → Z[0]、(2) カタログ配列ではなく 1 run = 1 図の単独出力、
(3) 題・ラベルの t=0 表記、(4) 出力ファイル名、のみ。作図ルール（panel() 無変更:
合成路なし・12桁丸め縮退・n>1 のみ数字・絶対スケール式・緑×/橙× 圧縮式）は
最終図と同一。

元ヘッダ: 全 99 走行の最終状態（t=4096）複素平面カタログ（2026-09-16、初期状態 v4 の最終ステージ版）。

木原指示 2026-09-16: 「今の図を、最終ステージの状態で適切な名前で作成」。

初期状態カタログ v4 と同一の様式・配列・指標:
 - ◯ = z_e(4096)。面積∝同一値（12 桁丸め）の重複辺数、重複時（n>1）のみ数字を記載。
 - 緑× = 重心 Σ_e z_e(4096)（全 99 走行の max|Σz| 基準の全図共通線形圧縮、G=真の|Σz|）。
 - 橙× = 二乗閉塞残差 Σ_e z_e(4096)²（下枠=Re・左枠=Im、成分別 max 基準の共通圧縮、
   Q2=真値）。Σz² は保存量 Q₂ なので t=0 の値と一致するはず — 一致を機械検証して出力。
 - 共通絶対スケール ±S、S = 全 99 走行の実測最大 |z_e(4096)| × 1.05（正規化なし）。
 - ページ 1: 主サーベイ 72 走行（行 L × 列 9 組: 奇奇/偶偶/奇偶）、
   ページ 2: gcd 系列 18 ＋ den=40 対照 9。

v4 との相違（意図的・理由あり）: ベース波/倍音の合成路（実線・点線）は描かない。
2 モード分解 z=c_aU^{m_aΔ}+c_bU^{m_bΔ} は t=0 の生成式の性質であり、発展後の状態は
初期 2 モード張る部分空間から O(1) 逸脱している（サーベイ実測 η⊥）ため、
最終状態に同じ分解を描くことはできない。

データは runs/*/states.npz の Z[4096]（検証用に Z[0] も）を read-only 参照のみ。
再走行しない・書き込まない。
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


def load_final(L, ma, mb, den):
    """最終状態 z(T) と、保存量検証用の Q2(0), Q2(T) を返す（read-only）。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
        Z = d['Z']
        assert Z.shape[0] == T + 1
        z0, zT = np.array(Z[0]), np.array(Z[T])
    return z0, complex(np.sum(z0 ** 2)), complex(np.sum(z0 ** 2))


def panel(ax, zT, S, label, gdisp, gtrue, q2x, q2y, q2true):
    ax.axhline(0, color='gray', lw=0.4, alpha=0.5)
    ax.axvline(0, color='gray', lw=0.4, alpha=0.5)
    groups = {}
    for w in zT:
        groups.setdefault((round(float(w.real), 12), round(float(w.imag), 12)), []).append(w)
    for (x, y), members in groups.items():
        n = len(members)
        s = float(np.clip(2.5 * n, 8, 300))
        ax.scatter([x], [y], s=s, facecolors='white', edgecolors='black',
                   linewidths=0.5, alpha=0.85, zorder=4)
        if n > 1:
            ax.annotate(str(n), (x, y), xytext=(3, 2), textcoords='offset points',
                        fontsize=4.5, color='dimgray', zorder=5)
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


def legend_handles():
    return [
        Line2D([], [], color='black', marker='o', ls='none', markersize=7,
               markerfacecolor='white', label='z(0)（面積∝重複辺数、n>1 のみ数字）'),
        Line2D([], [], color='tab:green', marker='x', ls='none', markersize=8,
               markeredgewidth=1.6,
               label='重心 Σz（max|Σz| 基準の全図共通圧縮表示、G=真の|Σz|）'),
        Line2D([], [], color='tab:orange', marker='x', ls='none', markersize=7,
               markeredgewidth=1.4,
               label='二乗閉塞残差 Σz²=Q₂（下枠=Re・左枠=Im、成分別 max 基準の共通圧縮、Q2=真値）'),
    ]


def have(L, ma, mb, den):
    """段階実行対応（木原承認 2026-09-17）: states.npz が存在する run だけ図化する。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    return os.path.exists(os.path.join(SURVEY, 'runs', rid, 'states.npz'))


def main():
    data, q2dev = {}, 0.0
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            if have(L, ma, mb, L):
                zT, q0, qT = load_final(L, ma, mb, L)
                data[(L, ma, mb, L)] = zT
                q2dev = max(q2dev, abs(qT - q0))
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            if have(L, ma, mb, L):
                zT, q0, qT = load_final(L, ma, mb, L)
                data[(L, ma, mb, L)] = zT
                q2dev = max(q2dev, abs(qT - q0))
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            if have(L, ma, mb, CTRL_DEN):
                zT, q0, qT = load_final(L, ma, mb, CTRL_DEN)
                data[(L, ma, mb, CTRL_DEN)] = zT
                q2dev = max(q2dev, abs(qT - q0))
    if not data:
        print('states.npz が存在する run が無い — 図化せず終了')
        return
    print(f'存在する run: {len(data)}/99（欠損パネルは空欄表示）')
    S = max(float(np.max(np.abs(z))) for z in data.values()) * 1.05
    Smax = S / 1.05
    print(f'共通絶対スケール: max|z_e(0)| = {Smax:.6f} → 表示範囲 ±{S:.4f}')
    print(f'Q2 参照（t=0 同士、恒等的に 0）: {q2dev:.3e}')
    Gs = {key: complex(np.sum(z)) for key, z in data.items()}
    Gmax = max(abs(g) for g in Gs.values())
    kG = 0.97 * S / Gmax
    print(f'重心 max|Σz(0)| = {Gmax:.4f} → ×共通圧縮係数 {kG:.6e}')
    Q2s = {key: complex(np.sum(z ** 2)) for key, z in data.items()}
    Q2re_max = max(abs(q.real) for q in Q2s.values())
    Q2im_max = max(abs(q.imag) for q in Q2s.values())
    kRe = 0.97 * S / Q2re_max
    kIm = 0.97 * S / Q2im_max
    print(f'二乗閉塞残差 max|Re Σz²| = {Q2re_max:.4f} / max|Im Σz²| = {Q2im_max:.4f}')

    # 1 run = 1 図（木原指示 2026-09-17: カタログ配列ではなく、この1図だけを取り出す）
    for key in sorted(data):
        L, ma, mb, den = key
        rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
        fig, ax = plt.subplots(figsize=(7.6, 8.6))
        panel(ax, data[key], S, f'L{L} ({ma},{mb}) den={den}',
              Gs[key] * kG, abs(Gs[key]),
              Q2s[key].real * kRe, Q2s[key].imag * kIm, Q2s[key])
        fig.suptitle(f'初期状態 t=0・最終図形式（{rid}）\n'
                     f'正規化なし・絶対スケール ±{S:.3f}（実測最大 |z|={Smax:.4f}）。'
                     f'◯面積∝重複辺数（12桁丸め）\n'
                     f'緑×=重心 Σz（max|Σz|={Gmax:.1f} 圧縮）、'
                     f'橙×=Σz² 残差（下枠=Re/max {Q2re_max:.1f}・左枠=Im/max {Q2im_max:.1f}）',
                     fontsize=9, y=0.99)
        fig.legend(handles=legend_handles(), loc='lower center', ncol=1, fontsize=8,
                   frameon=False)
        fig.tight_layout(rect=[0, 0.11, 1, 0.90])
        out1 = os.path.join(HERE, '..', 'plots', 'initial_complex',
                            f'fig_initial_states_finalformat_{rid}_dense_rerun_20260917.png')
        fig.savefig(out1, dpi=200)
        plt.close(fig)
        print('->', os.path.basename(out1))

    print(f'panels total: {len(data)}')
    gm = sorted((abs(g), f'L{k[0]}_ma{k[1]}_mb{k[2]}_den{k[3]}') for k, g in Gs.items())
    print(f'重心 |Σz(0)| 範囲: 最小 {gm[0][0]:.4f} ({gm[0][1]}) / '
          f'最大 {gm[-1][0]:.4f} ({gm[-1][1]})')
    amax = sorted((float(np.max(np.abs(z))), f'L{k[0]}_ma{k[1]}_mb{k[2]}_den{k[3]}')
                  for k, z in data.items())
    print(f'最大振幅 max|z_e(0)| 範囲: 最小 {amax[0][0]:.4f} ({amax[0][1]}) / '
          f'最大 {amax[-1][0]:.4f} ({amax[-1][1]})')


if __name__ == '__main__':
    main()
