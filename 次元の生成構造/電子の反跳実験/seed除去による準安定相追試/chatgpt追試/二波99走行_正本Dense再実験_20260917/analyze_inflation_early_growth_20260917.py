#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""インフレーション指標 H⊥/H の早期成長解析（2026-09-17 木原指示・read-only 後処理）。

問い: 正本インフレーション図（横軸線形・縦軸 log）で立ち上がりが「一気に駆け上がる」
ように見えるのは拡大不足か、それとも床滞在段階が存在しないためか。

方法（丸め・間引きなし、fine 系列 den=200L の全 step を使用）:
 - 正本 series() と同式で H⊥/H を全 step 計算（p=Re z0 正規化、q=Im z0 直交化、
   f = <z⊥,z⊥>/<z,z>）。ベクトル化のみで式は無変更。
 - 早期値 f(1),f(2),...,f(2048) を生値で記録。
 - 初期成長則 f ∝ t^α の対数勾配 α を区間 [1,64]（事前固定）で最小二乗推定。
 - 図: 両対数（log t – log f）。t^2 成長は傾き2の直線、1step ジャンプは切片、
   飽和は折れとして分離される。正本の線形横軸図では潰れる構造を可視化する診断図。

出力: diagnostics/inflation_early_growth/
  early_growth_<rid>.csv（t, f_hperp 全 step）/ early_growth_summary.json /
  fig_inflation_early_growth_loglog_20260917.png
一次データ（runs_fine/*/states.npz）は read-only。正本 series() の式は変更しない。
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'diagnostics', 'inflation_early_growth')
FLOOR = 1e-31                                   # 正本と同値
FIT_LO, FIT_HI = 1, 64                          # 対数勾配の推定区間（事前固定）
PROBE = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]


def hperp_series(Z):
    """正本 series() の H⊥/H 計算部（式無変更・全 step ベクトル化）。"""
    z0 = Z[0]
    p = z0.real.astype(np.float64).copy(); p /= np.linalg.norm(p)
    q = z0.imag.astype(np.float64).copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    a = Z @ p                                    # <p,z>（p 実なので dot）
    b = Z @ q
    H = np.sum(np.abs(Z) ** 2, axis=1).real
    f = (H - np.abs(a) ** 2 - np.abs(b) ** 2) / H
    return f


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(HERE, 'runs_fine', '*', 'states.npz')))
    if not paths:
        print('runs_fine なし')
        sys.exit(1)
    summary = {}
    fig, ax = plt.subplots(figsize=(11.0, 7.5))
    cols = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, p in enumerate(paths):
        rid = os.path.basename(os.path.dirname(p))
        with np.load(p) as d:
            Z = np.array(d['Z'])
        f = hperp_series(Z)
        with open(os.path.join(OUTDIR, f'early_growth_{rid}.csv'), 'w') as fp:
            fp.write('t,f_hperp\n')
            fp.writelines(f'{t},{repr(x)}\n' for t, x in enumerate(f))
        tt = np.arange(FIT_LO, FIT_HI + 1)
        alpha = float(np.polyfit(np.log(tt), np.log(f[FIT_LO:FIT_HI + 1]), 1)[0])
        probe = {str(t): float(f[t]) for t in PROBE if t < len(f)}
        summary[rid] = {'f_step1': float(f[1]), 'alpha_fit_1_64': alpha,
                        'f_final_tail_mean': float(f[-8192:].mean()), 'probe': probe}
        tp = np.arange(1, min(4097, len(f)))
        ax.loglog(tp, np.maximum(f[1:len(tp) + 1], FLOOR), lw=0.8, color=cols[i],
                  label=f'{rid}（α={alpha:.2f}）')
        print(f'{rid}: f(1)={f[1]:.2e} α(1..64)={alpha:.3f} 飽和={f[-8192:].mean():.3f}')
    # 参照勾配 t^2
    tg = np.array([2, 200])
    ax.loglog(tg, 1e-5 * (tg / 1) ** 2, color='gray', ls='--', lw=1.0,
              label='参照勾配 f∝t² (弾道)')
    ax.set_xlabel('step t（対数）', fontsize=11)
    ax.set_ylabel('H⊥/H（対数）', fontsize=11)
    ax.set_title('インフレーション指標の早期成長（fine 系列 den=200L、両対数）\n'
                 '正本 series() 式・1step 目で床 10⁻³¹ から離脱、以後 f∝t^α（α≈2 弾道）で飽和',
                 fontsize=12)
    ax.grid(True, which='both', lw=0.3, alpha=0.4)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, 'fig_inflation_early_growth_loglog_20260917.png'), dpi=160)
    plt.close(fig)
    with open(os.path.join(OUTDIR, 'early_growth_summary.json'), 'w', encoding='utf-8') as f:
        json.dump({'fit_interval': [FIT_LO, FIT_HI], 'runs': summary}, f,
                  indent=1, ensure_ascii=False)
    print('-> diagnostics/inflation_early_growth/')


if __name__ == '__main__':
    main()
