#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# fine 系列（den=200L・runs_fine/・T=819200）版。正本カタログ版（存在runガード付き実行コピー）の
# 物理コピーに機械置換のみ適用（2026-09-17 木原指示。物理式・パネル描画ロジック無変更）。
"""全 99 走行のインフレーション図カタログ H⊥/H（2026-09-16）。

木原指示 2026-09-16: 横軸=ステップ、縦軸=インフレーション指標 H⊥/H。
配列は既存カタログと同一。作り込み禁止 — 正本の式・軸をそのまま使う。

指標の正本（無変更コピー）:
  干渉保存力学_資格審査とシード無し系列_20260831/make_parent型初期値_自己無撞着構造_20260907/
  plot3d_makeparent_selectedN_v1.py の series()（正本 metrics() と同一式・N16対照検証済み）:
    p = Re z0 / ‖Re z0‖
    q = ( Im z0 − (Im z0·p) p ) / ‖·‖
    f(t) = ⟨z⊥,z⊥⟩/⟨z,z⟩,  z⊥ = z − p(p·z) − q(q·z)
    fz = max(f, FLOOR),  FLOOR = 1e-31（表示床クリップ）
  閾値・場合分けは一切ない。凍結類 5 走行も Im z0 の機械残差(~1e-14)で
  そのまま数値的に通る（規約の追加をしない）。

軸の正本（N16_original_note_source_den16_step0_500_3D.html の layout から抽出）:
  縦軸: type=log, range=[1e-31, 1e0]
  目盛: 1e-31, 1e-25, 1e-20, 1e-15, 1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1e0
        （ticktext "10^-31" … "10^0"）、軸題 Hperp/H
  横軸: step, range [0, steps]
  f=H⊥/H≤1 が恒等的に成り立つため全データがこの範囲に必ず収まる。

- ページ 1: 主サーベイ 72 走行（行 L × 列 9 組）、ページ 2: gcd 18 ＋ den=40 対照 9。
- パネルが小さいため目盛ラベルは左端列・最下行のみ表示（目盛位置は全パネル正本どおり）。
- データは runs/*/states.npz の Z を read-only 参照のみ（再走行しない・書き込まない）。
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

DEN_MULT = 200
TF = T * DEN_MULT
RANGE_T0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
RANGE_T1 = int(sys.argv[2]) if len(sys.argv) > 2 else 1000  # 拡大区間（引数指定可）

FLOOR = 1e-31                                   # 正本と同値
YTICKS = [1e-31, 1e-25, 1e-20, 1e-15, 1e-10, 1e-05, 0.001, 0.01, 0.1, 1.0]
YLABELS = ['10^-31', '10^-25', '10^-20', '10^-15', '10^-10',
           '10^-5', '10^-3', '10^-2', '10^-1', '10^0']


def series_hperp(L, ma, mb, den):
    """正本 series() の H⊥/H 計算部を無変更コピー（入力の読み方のみ本 DB に適合）。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den * DEN_MULT}'
    d = np.load(os.path.join(SURVEY, 'runs_fine', rid, 'states.npz'))
    S = np.asarray(d['Z'], np.complex128)
    assert S.shape[0] == TF + 1
    z0 = S[0]
    p = z0.real.astype(np.float64).copy(); p /= np.linalg.norm(p)
    q = z0.imag.astype(np.float64).copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    Ts = S.shape[0]
    f = np.empty(Ts)
    for s in range(Ts):
        z = S[s]
        a = np.dot(p, z); b = np.dot(q, z)
        zp = z - p * a - q * b
        f[s] = np.vdot(zp, zp).real / np.vdot(z, z).real
    fz = np.maximum(f, FLOOR)
    return fz


def panel(ax, fz, label, xlabels=False, ylabels=False):
    ax.plot(np.arange(RANGE_T0, RANGE_T0 + len(fz)), fz, color='black', lw=0.6)
    ax.set_yscale('log')
    ax.set_xlim(RANGE_T0, RANGE_T1)
    ax.set_ylim(FLOOR, 1.0)
    ax.set_yticks(YTICKS)
    ax.set_yticklabels(YLABELS if ylabels else [], fontsize=4)
    ax.set_xticks([RANGE_T0, (RANGE_T0 + RANGE_T1) // 2, RANGE_T1])
    ax.set_xticklabels([str(RANGE_T0), str((RANGE_T0 + RANGE_T1) // 2), str(RANGE_T1)] if xlabels else [], fontsize=6)
    ax.tick_params(axis='y', length=1.5)
    ax.tick_params(axis='x', length=1.5)
    ax.set_title(label, fontsize=8, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)


def have(L, ma, mb, den):
    """段階実行対応（木原承認 2026-09-17）: states.npz が存在する run だけ図化する。"""
    rid = f'L{L}_ma{ma}_mb{mb}_den{den * DEN_MULT}'
    return os.path.exists(os.path.join(SURVEY, 'runs_fine', rid, 'states.npz'))


def main():
    data = {}
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            if have(L, ma, mb, L):
                data[(L, ma, mb, L)] = series_hperp(L, ma, mb, L)[RANGE_T0:RANGE_T1 + 1]
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            if have(L, ma, mb, L):
                data[(L, ma, mb, L)] = series_hperp(L, ma, mb, L)
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            if have(L, ma, mb, CTRL_DEN):
                data[(L, ma, mb, CTRL_DEN)] = series_hperp(L, ma, mb, CTRL_DEN)[RANGE_T0:RANGE_T1 + 1]
    if not data:
        print('states.npz が存在する run が無い — 図化せず終了')
        return
    print(f'存在する run: {len(data)}/99（欠損パネルは空欄表示）')
    f0 = sorted((fz[0], k) for k, fz in data.items())
    fe = sorted((fz[-1], k) for k, fz in data.items())
    print(f'H⊥/H step{RANGE_T0}: 最小 {f0[0][0]:.3e} {f0[0][1]} / 最大 {f0[-1][0]:.3e} {f0[-1][1]}')
    print(f'H⊥/H step{RANGE_T1}: 最小 {fe[0][0]:.3e} {fe[0][1]} / 最大 {fe[-1][0]:.3e} {fe[-1][1]}')

    # ページ 1: 主サーベイ 72（8 行 × 9 列）
    fig, axs = plt.subplots(len(MAIN_L), len(MAIN_PAIRS), figsize=(26.0, 13.5))
    for i, L in enumerate(MAIN_L):
        for j, (ma, mb) in enumerate(MAIN_PAIRS):
            if (L, ma, mb, L) in data:
                panel(axs[i, j], data[(L, ma, mb, L)], f'L{L} ({ma},{mb})',
                      xlabels=(i == len(MAIN_L) - 1), ylabels=(j == 0))
            else:
                axs[i, j].axis('off')
        axs[i, 0].set_ylabel(f'L={L}   Hperp/H', fontsize=8)
    for j, (ma, mb) in enumerate(MAIN_PAIRS):
        cls = 'odd-odd' if j < 3 else ('even-even' if j < 6 else 'odd-even')
        axs[0, j].set_title(f'({ma},{mb})  {cls}\nL{MAIN_L[0]}', fontsize=8, pad=2)
    fig.suptitle(f'インフレーション図 H⊥/H 拡大（fine 系列、den=200L、t={RANGE_T0}..{RANGE_T1}）— '
                 '正本 series() 式（p=Re z0、q=Im z0 直交化、床 10⁻³¹）・'
                 '正本軸 log [10⁻³¹, 10⁰]・正本目盛', fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0.012, 1, 0.972])
    out1 = os.path.join(HERE, '..', 'plots', 'inflation', f'fig_inflation_hperp_fine_zoom_t{RANGE_T0}_{RANGE_T1}_page1_20260917.png')
    fig.savefig(out1, dpi=130)
    plt.close(fig)
    print('page1 ->', os.path.basename(out1))

    # ページ 2: gcd 18（6 行 × 3 列）＋ den=40 対照 9（3 行 × 3 列）
    fig = plt.figure(figsize=(10.0, 16.5))
    gs = fig.add_gridspec(9, 3, hspace=0.55, wspace=0.16,
                          top=0.945, bottom=0.03, left=0.09, right=0.97)
    for i, L in enumerate(GCD_L):
        for j, (ma, mb) in enumerate(GCD_PAIRS):
            if (L, ma, mb, L) in data:
                ax = fig.add_subplot(gs[i, j])
                panel(ax, data[(L, ma, mb, L)], f'L{L} ({ma},{mb}) den={L}',
                      ylabels=(j == 0))
    for i, L in enumerate(CTRL_L):
        for j, (ma, mb) in enumerate(CTRL_PAIRS):
            if (L, ma, mb, CTRL_DEN) in data:
                ax = fig.add_subplot(gs[6 + i, j])
                panel(ax, data[(L, ma, mb, CTRL_DEN)], f'L{L} ({ma},{mb}) den={CTRL_DEN}',
                      xlabels=(i == len(CTRL_L) - 1), ylabels=(j == 0))
    fig.suptitle('インフレーション図 H⊥/H（gcd 系列 18 走行・上 6 行 ＋ den=40 対照 9 走行・下 3 行）\n'
                 '正本 series() 式・床 10⁻³¹・正本軸 log [10⁻³¹, 10⁰]・正本目盛',
                 fontsize=10, y=0.99)
    out2 = os.path.join(HERE, '..', 'plots', 'inflation', f'fig_inflation_hperp_fine_zoom_t{RANGE_T0}_{RANGE_T1}_page2_20260917.png')
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print('page2 ->', os.path.basename(out2))
    print(f'panels total: {len(data)}')


if __name__ == '__main__':
    main()
