#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""星凝縮の緩和時間解析（2026-09-17 木原指示・read-only 後処理、fine 系列）。

問い: 初期振幅（凝縮先の星の外の辺）が 0 に落ち着くまでの期間はどれだけか。

方法（丸め・間引きなし、全 819201 step）:
 1. 各 run の末尾窓（8192 step）で星エネルギー最大の頂点 v* を同定（凝縮先）。
 2. 星外エネルギー比 f_out(t) = Σ_{e∉star(v*)} |z_e(t)|² / Σ_e |z_e(t)|² を全 step で計算。
 3. 事前固定の水準 LEVELS について、
      first(t): f_out が初めて水準を下回った step
      settle(t): f_out が最後に水準以上だった step（以後二度と超えない＝落ち着き）
    を生値で記録する（水準は事前固定、結果を見て変更しない）。
 4. 図: (a) 3 run の f_out(t)（対数縦軸・全点）
        (b) run ごとの全辺振幅包絡 |z_e(t)|/amax（辺ごと1色・全点）

出力: diagnostics/star_condensation/
  f_out_<rid>.csv（t, f_out 全 step）/ settle_times.json /
  fig_star_condensation_fout_20260917.png / fig_edge_amplitudes_<rid>_20260917.png
一次データ（runs_fine/*/states.npz）は read-only。
"""
import glob
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'diagnostics', 'star_condensation')
LEVELS = [0.5, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8]   # 事前固定


def analyze(path):
    with np.load(path) as d:
        Z = np.array(d['Z'])
        L = int(d['L'])
    rid = os.path.basename(os.path.dirname(path))
    ja, jb = np.triu_indices(L, k=1)
    E = np.abs(Z) ** 2                              # (T+1, M)
    tail = E[-8192:]
    stars = {v: [e for e in range(len(ja)) if ja[e] == v or jb[e] == v]
             for v in range(L)}
    v_star = max(stars, key=lambda v: tail[:, stars[v]].sum())
    out = [e for e in range(len(ja)) if e not in stars[v_star]]
    f_out = E[:, out].sum(axis=1) / E.sum(axis=1)
    cross = {}
    for lv in LEVELS:
        below = f_out < lv
        first = int(np.argmax(below)) if below.any() else None
        above = np.nonzero(~below)[0]
        settle = (int(above[-1]) + 1) if len(above) and below.any() else (0 if below.all() else None)
        if settle is not None and settle > len(f_out) - 1:
            settle = None                            # 最終stepまで上回ったまま
        cross[repr(lv)] = {'first_below': first, 'settled_after': settle}
    return rid, L, v_star, f_out, cross, E, ja, jb


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(HERE, 'runs_fine', '*', 'states.npz')))
    if not paths:
        print('runs_fine なし')
        sys.exit(1)
    summary = {}
    figf, axf = plt.subplots(figsize=(13.0, 7.0))
    cols = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, p in enumerate(paths):
        t0 = time.time()
        rid, L, v_star, f_out, cross, E, ja, jb = analyze(p)
        with open(os.path.join(OUTDIR, f'f_out_{rid}.csv'), 'w') as f:
            f.write('t,f_out\n')
            f.writelines(f'{t},{repr(x)}\n' for t, x in enumerate(f_out))
        summary[rid] = {'star_vertex': v_star, 'f_out_final': float(f_out[-1]),
                        'f_out_min': float(f_out.min()), 'crossings': cross}
        axf.plot(np.arange(len(f_out)), np.maximum(f_out, 1e-31), lw=0.7,
                 color=cols[i], label=f'{rid}（星=頂点{v_star}）')
        # 辺別振幅包絡図
        amax = float(np.abs(np.sqrt(E)).max())
        fig, ax = plt.subplots(figsize=(13.0, 7.0))
        ecols = plt.cm.turbo(np.linspace(0.02, 0.98, E.shape[1]))
        for e in range(E.shape[1]):
            ax.plot(np.arange(E.shape[0]), np.sqrt(E[:, e]) / amax,
                    lw=0.4, alpha=0.7, color=ecols[e])
        ax.set_xlim(0, E.shape[0] - 1)
        ax.set_xlabel('step')
        ax.set_ylabel('|z_e(t)| / 全体最大')
        ax.set_title(f'{rid} — 全辺の振幅包絡（全 {E.shape[0]} step・全点）。凝縮先=頂点{v_star} の星')
        fig.tight_layout()
        fig.savefig(os.path.join(OUTDIR, f'fig_edge_amplitudes_{rid}_20260917.png'), dpi=150)
        plt.close(fig)
        print(f'{rid}: 星=頂点{v_star} f_out(末)={f_out[-1]:.3e} min={f_out.min():.3e} '
              f'({time.time() - t0:.0f}s)')
    axf.set_yscale('log')
    axf.set_xlabel('step')
    axf.set_ylabel('星外エネルギー比 f_out(t)')
    axf.set_title('星凝縮の緩和: 星外エネルギー比の全時系列（fine 系列 den=200L、対数縦軸・全点）\n'
                  f'水準線 = 事前固定 LEVELS {LEVELS}')
    for lv in LEVELS:
        axf.axhline(lv, color='gray', lw=0.4, ls=':', alpha=0.6)
    axf.legend(fontsize=9)
    figf.tight_layout()
    figf.savefig(os.path.join(OUTDIR, 'fig_star_condensation_fout_20260917.png'), dpi=160)
    plt.close(figf)
    with open(os.path.join(OUTDIR, 'settle_times.json'), 'w', encoding='utf-8') as f:
        json.dump({'levels': LEVELS, 'runs': summary}, f, indent=1, ensure_ascii=False)
    print('-> diagnostics/star_condensation/')


if __name__ == '__main__':
    main()
