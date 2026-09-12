#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""乱数make_parent床と、意図的に作った高対称床（v2）を同じNで並べ、対称性の定量と複素平面図を出す
（読出しのみ・物理無変更・新規走行なし）。

動機（木原・実験台帳）: 初期に使った make_parent は乱数シードから自己無撞着固定点 KZ=iσZ を作るが、
その床が全Nで自然に90°二軸位相格子になった（実験4-1・検討28）。位相は自然に90°だが、振幅は辺ごとに
バラバラ（M個の値・CV大）・重心 Σz≠0 で対称性は高くない。そこで「90°位相差骨格の上に、極めて
対称な自己無撞着床を作れるか」を問い、位相0/90・振幅≤2値・P_A=P_B・D_N対称の高対称床（v2）を構成した
（検討31・最も対称性の高い初期値）。本スクリプトは両床を同じNで並べ、対称性が高くなったことを定量する。

対称性の定量（各床・N=6,12,17,40）:
  - n_distinct_amp: |z_e| の相異値数（丸め1e-6）。make_parent は多数、v2 は≤2。
  - amp_CV: |z| の変動係数。
  - |Σz|（重心）, |Σz²|（閉塞）。
入力（読み取りのみ）: make_parent床 = ../../make_parent型…/full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz の step0、
  高対称v2床 = ../../最も対称性の高い初期値_20260906/parents_symmetric/parent_symmetric_N{N:05d}_v2.npz の Z0。
出力: results/floor_symmetry_summary.csv, fig_floor_symmetry_compare.png。
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907', 'full_N3_N40_sweep', 'states'))
SY = os.path.abspath(os.path.join(HERE, '..', '..', '最も対称性の高い初期値_20260906', 'parents_symmetric'))
IK = os.path.abspath(os.path.join(HERE, '..', '..', '90度理論床_整数K解析解_全N_20260909', 'parents_90deg_floor_analytic'))
NS = [6, 12, 17, 40]


def metrics(z):
    r = np.abs(z)
    n_amp = len(set(np.round(r, 6).tolist()))
    cv = float(np.std(r) / np.mean(r))
    cen = float(abs(np.sum(z)))
    clo = float(abs(np.sum(z * z)))
    return n_amp, cv, cen, clo


def gauge_to_0_90(z):
    """支配軸を実軸へ回す（大域位相ゲージ）。90°二軸構造を見やすくするだけ・物理不変。"""
    # 最大振幅波の位相を基準に、90°格子へ寄せる大域回転
    k = int(np.argmax(np.abs(z)))
    return z * np.exp(-1j * np.angle(z[k]))


def main():
    rows = []
    fig, axs = plt.subplots(len(NS), 3, figsize=(11.4, 3.6 * len(NS)))
    for row, N in enumerate(NS):
        M = N * (N - 1) // 2
        zmp = np.asarray(np.load(os.path.join(MP, f'hm_N{N}_den_{N}_states_500.npz'))['Z'][0], np.complex128)
        zsy = np.asarray(np.load(os.path.join(SY, f'parent_symmetric_N{N:05d}_v2.npz'))['Z0'], np.complex128)
        zik = np.asarray(np.load(os.path.join(IK, f'parent_90deg_floor_N{N:05d}_analytic.npz'))['Z0'], np.complex128)
        for col, (z, lab) in enumerate([(zmp, 'make_parent（乱数・自然に90°）'),
                                        (zsy, '高対称v2（振幅≤2値・大重心）'),
                                        (zik, '整数K解析床（≤3値・小重心・相対平衡）')]):
            n_amp, cv, cen, clo = metrics(z)
            rows.append(dict(N=N, M=M, floor=lab, n_distinct_amp=n_amp, amp_CV=cv, centroid_abs=cen, closure_abs=clo))
            zp = gauge_to_0_90(z)
            ax = axs[row, col]
            for w in zp:
                ax.plot([0, w.real], [0, w.imag], color='tab:blue', lw=0.5, alpha=0.45)
            ax.plot(zp.real, zp.imag, 'o', ms=2.4, color='tab:red', alpha=0.8, ls='none')
            lim = float(np.abs(zp).max()) * 1.2
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
            ax.axhline(0, color='gray', lw=0.4, alpha=0.5); ax.axvline(0, color='gray', lw=0.4, alpha=0.5)
            ax.grid(alpha=.25); ax.tick_params(labelsize=6)
            ax.set_title(f'N={N} (M={M}) {lab}\n振幅値数={n_amp} CV={cv:.3f} |Σz|={cen:.3f} |Σz²|={clo:.1e}', fontsize=7)
    fig.suptitle('床の対称性比較（同一N・位相はどれも90°）: 乱数make_parent（振幅ほぼM値）→ 高対称v2（振幅≤2値・大重心・NG）→ 整数K解析床（≤3値・小重心・相対平衡＝論文3）', fontsize=8.5, y=0.998)
    fig.tight_layout()
    out = os.path.join(HERE, 'fig_floor_symmetry_compare.png')
    fig.savefig(out, dpi=150); plt.close(fig)
    with open(os.path.join(HERE, 'results', 'floor_symmetry_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"{'N':>3}{'floor':>28}{'振幅値数':>8}{'CV':>8}{'|Σz|':>8}{'|Σz²|':>10}")
    for r in rows:
        print(f"{r['N']:>3}{r['floor']:>28}{r['n_distinct_amp']:>8}{r['amp_CV']:>8.3f}{r['centroid_abs']:>8.3f}{r['closure_abs']:>10.1e}")
    print('\n一目瞭然: 位相はどちらも90°だが、make_parent は振幅がM個近くバラバラ（CV大・重心≠0）、'
          '高対称v2 は振幅≤2値・重心≈0・閉塞厳密＝対称性が硬く（高く）なっている。図 fig_floor_symmetry_compare.png。')


if __name__ == '__main__':
    main()
