#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 99 走行の形状指標 χ = |C²|/H カタログ 英語ラベル版（2026-09-16）。
plot_chi_all99_20260916.py の忠実コピー。ラベル文字列と出力ファイル名のみ変更。

定義（論文2 §6.4）: Z=x+iy に対し H=Z†Z、C²=Z^T Z は共に厳密保存され、
実部・虚部の Gram 行列が固定される。固有値 λ±=(H±|C²|)/2、形状指標

    χ = |C²|/H ∈ [0,1]

は各走行の保存ラベル。χ=0 は円形 2-frame、0<χ<1 は楕円、χ=1 は一直線縮退。
定理（全成分非零・線グラフ連結の下）: χ=1 ⟺ K=0（厳密固定点）。

- χ(0) は runs/*/states.npz の Z[0]（read-only）から計算。
- χ(t) の時間不変性は runs/*/invariants.csv（毎 step の H, Q2_re, Q2_im）で
  全 4097 step にわたり機械検証する（既知の repr 書式はパーサ側で剥がす。DB 不変更）。
- 表示: 主サーベイ 72 = 行 L × 列 9 組の行列、gcd 18 ＋ den=40 対照 9 = 右側 3 列行列。
  各セルに χ を 3 桁で記載、χ=1（凍結類）は赤枠で強調。色は共通 [0,1]。
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, SURVEY)
from gen_manifest import MAIN_L, MAIN_PAIRS, GCD_L, GCD_PAIRS, CTRL_L, CTRL_PAIRS, CTRL_DEN  # noqa: E402


def _f(s):
    s = s.strip()
    if s.startswith('np.float64('):
        s = s[11:-1]
    return float(s)


def chi_run(L, ma, mb, den):
    """χ(0) を states から、χ(t) の最大時間変動を invariants から返す（read-only）。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
        z0 = np.array(d['Z'][0])
    H0 = float(np.vdot(z0, z0).real)
    chi0 = abs(complex(np.sum(z0 ** 2))) / H0
    rows = []
    with open(os.path.join(SURVEY, 'runs', rid, 'invariants.csv')) as f:
        next(f)
        for line in f:
            p = line.split(',')
            rows.append((_f(p[1]), _f(p[2]), _f(p[3])))
    arr = np.array(rows)
    chi_t = np.hypot(arr[:, 1], arr[:, 2]) / arr[:, 0]
    return chi0, float(np.max(np.abs(chi_t - chi_t[0])))


def draw_matrix(ax, chis, row_labels, col_labels, frozen, title):
    m = np.array(chis)
    ax.imshow(m, cmap='viridis', vmin=0.0, vmax=1.0, aspect='auto')
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            c = 'black' if m[i, j] > 0.5 else 'white'
            ax.text(j, i, f'{m[i, j]:.3f}', ha='center', va='center',
                    fontsize=7.5, color=c)
            if frozen[i][j]:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor='red', lw=2.0))
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels(col_labels, fontsize=7)
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.tick_params(length=0)


def main():
    devmax = 0.0
    main_chi, main_frozen = [], []
    for L in MAIN_L:
        row, fro = [], []
        for ma, mb in MAIN_PAIRS:
            c, d = chi_run(L, ma, mb, L)
            row.append(c); fro.append((ma + mb) % L == 0); devmax = max(devmax, d)
        main_chi.append(row); main_frozen.append(fro)
    gcd_chi, gcd_frozen = [], []
    for L in GCD_L:
        row, fro = [], []
        for ma, mb in GCD_PAIRS:
            c, d = chi_run(L, ma, mb, L)
            row.append(c); fro.append((ma + mb) % L == 0); devmax = max(devmax, d)
        gcd_chi.append(row); gcd_frozen.append(fro)
    ctrl_chi, ctrl_frozen = [], []
    for L in CTRL_L:
        row, fro = [], []
        for ma, mb in CTRL_PAIRS:
            c, d = chi_run(L, ma, mb, CTRL_DEN)
            row.append(c); fro.append((ma + mb) % L == 0); devmax = max(devmax, d)
        ctrl_chi.append(row); ctrl_frozen.append(fro)

    allc = [c for rows_, _ in ((main_chi, 0), (gcd_chi, 0), (ctrl_chi, 0))
            for r in rows_ for c in r]
    n_one = sum(1 for c in allc if 1 - c < 1e-12)
    others = [c for c in allc if 1 - c >= 1e-12]
    print(f'χ=1（機械精度）: {n_one} / 99 走行（凍結類）')
    print(f'非凍結 94 走行の χ: 最小 {min(others):.6f} / 最大 {max(others):.6f}')
    print(f'χ(t) の最大時間変動（invariants 全 4097 step、全 99 run）: {devmax:.3e}')

    fig = plt.figure(figsize=(17.0, 7.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[3.0, 1.05], height_ratios=[2, 1],
                          wspace=0.16, hspace=0.35,
                          left=0.055, right=0.985, top=0.88, bottom=0.06)
    ax1 = fig.add_subplot(gs[:, 0])
    cls = ['odd-odd'] * 3 + ['even-even'] * 3 + ['odd-even'] * 3
    cols1 = [f'({a},{b})\n{c}' for (a, b), c in zip(MAIN_PAIRS, cls)]
    draw_matrix(ax1, main_chi, [f'L={L}' for L in MAIN_L], cols1, main_frozen,
                'Main survey, 72 runs (den=L)')
    ax2 = fig.add_subplot(gs[0, 1])
    draw_matrix(ax2, gcd_chi, [f'L={L}' for L in GCD_L],
                [f'({a},{b})' for a, b in GCD_PAIRS], gcd_frozen,
                'gcd series, 18 runs (den=L)')
    ax3 = fig.add_subplot(gs[1, 1])
    draw_matrix(ax3, ctrl_chi, [f'L={L}' for L in CTRL_L],
                [f'({a},{b})' for a, b in CTRL_PAIRS], ctrl_frozen,
                'den=40 controls, 9 runs')
    im = ax1.images[0]
    cb = fig.colorbar(im, ax=[ax1, ax2, ax3], fraction=0.025, pad=0.01)
    cb.set_label('χ', fontsize=11)
    fig.suptitle('Shape index χ = |C²|/H (all 99 runs; conserved) — fixed-Gram eigenvalues '
                 'λ±=(H±|C²|)/2; χ=1 ⟺ K=0 (red frame = 5 frozen runs, χ=1 exact)',
                 fontsize=12, y=0.965)
    out = os.path.join(HERE, 'fig_chi_all99_en_20260916.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print('fig ->', os.path.basename(out))


if __name__ == '__main__':
    main()
