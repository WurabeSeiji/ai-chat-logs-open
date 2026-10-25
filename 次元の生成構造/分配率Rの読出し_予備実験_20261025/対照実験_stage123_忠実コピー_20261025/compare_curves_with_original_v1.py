#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照実験（案B）の曲線照合（読み出しのみ）: 本フォルダの timeseries/summary CSV を原本（2026-09-05）の同じ (N, den) 行と比べる。

項目（走行ごと）: onset（H⊥/H>0.05 の最初の step）、f(0)・f(1)・final・max、log10 f の全 step での最大ずれ、
指数区間（−13 < log10 f < −3）の傾き（桁/step、最小二乗）、H_total の相対ドリフト最大、閉塞 |zᵀz|/H の最大。
出力: results/compare_curves_with_original_v1.json、results/compare_curves_table.md、results/fig_compare_curves_control_vs_original.png。解釈は含まない。"""
import csv
import json
import math
import os

import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = ('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/'
        'ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/'
        '干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905')
NS = (3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40)
TS = 'results/timeseries_64bit_with124_N3_N40.csv'


def load_ts(path):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            key = (int(r['N']), int(r['denominator']))
            d.setdefault(key, []).append((int(r['step']), float(r['Hperp_frac']), float(r['H_total']), float(r['global_closure'])))
    for k in d:
        d[k].sort()
    return d


def stats(rows):
    t = np.array([r[0] for r in rows]); f = np.array([r[1] for r in rows]); h = np.array([r[2] for r in rows]); c = np.array([r[3] for r in rows])
    lf = np.log10(np.maximum(f, 1e-300))
    ix = np.flatnonzero(f > 0.05)
    m = (lf > -13) & (lf < -3)
    slope = float(np.polyfit(t[m], lf[m], 1)[0]) if m.sum() >= 5 else None
    return dict(onset=int(ix[0]) if ix.size else -1, f0=float(f[0]), f1=float(f[1]), final=float(f[-1]), fmax=float(f.max()),
                slope_dec_per_step=slope, n_exp_points=int(m.sum()),
                Htot_rel_drift_max=float(np.max(np.abs(h - h[0]) / h[0])), closure_max=float(c.max())), lf


ctrl = load_ts(os.path.join(BASE, TS))
orig = load_ts(os.path.join(ORIG, TS))
out = {'runs': {}}
rows_md = ['| N | den | onset ctrl/orig | slope ctrl/orig (dec/step) | final ctrl/orig | max\\|Δlog10 f\\| | Htot drift ctrl/orig | closure max ctrl/orig | 親 |',
           '|---|---|---|---|---|---|---|---|---|']
fig, axs = plt.subplots(3, 4, figsize=(16, 11)); axs = axs.ravel()
for k, N in enumerate(NS):
    ax = axs[k]
    for den, col in ((N, 'tab:blue'), (124, 'tab:orange')):
        key = (N, den)
        if key not in ctrl or key not in orig:
            out['runs'][f'N{N}_den{den}'] = 'missing'; continue
        sc, lc = stats(ctrl[key]); so, lo = stats(orig[key])
        dl = float(np.max(np.abs(lc - lo)))
        rec = dict(control=sc, original=so, max_abs_dlog10f=dl,
                   parent_note=('別の解' if N in (3, 4, 5) else '末尾ビット差'))
        out['runs'][f'N{N}_den{den}'] = rec
        fs = lambda x: '—' if x is None else f'{x:.4f}'
        rows_md.append(f"| {N} | {den} | {sc['onset']} / {so['onset']} | {fs(sc['slope_dec_per_step'])} / {fs(so['slope_dec_per_step'])} | "
                       f"{sc['final']:.4f} / {so['final']:.4f} | {dl:.2f} | {sc['Htot_rel_drift_max']:.1e} / {so['Htot_rel_drift_max']:.1e} | "
                       f"{sc['closure_max']:.1e} / {so['closure_max']:.1e} | {rec['parent_note']} |")
        t = np.arange(len(lc))
        ax.semilogy(t, 10.0 ** lc, color=col, linewidth=1.2, label=f'ctrl den={den}')
        ax.semilogy(t, 10.0 ** lo, color=col, linewidth=1.0, linestyle='--', label=f'orig den={den}')
    ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.set_title(f'N={N}' + (' (parent differs)' if N in (3, 4, 5) else '')); ax.grid(alpha=.25)
    ax.legend(fontsize=6, loc='lower right')
    if k // 4 == 2: ax.set_xlabel('step')
    if k % 4 == 0: ax.set_ylabel('Hperp/H')
fig.suptitle('Control (solid, macOS 26.6.2) vs original 2026-09-05 (dashed): Hperp/H, den=N (blue) and 124 (orange)', y=.995)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'results', 'fig_compare_curves_control_vs_original.png'), dpi=160); plt.close(fig)
with open(os.path.join(BASE, 'results', 'compare_curves_with_original_v1.json'), 'w') as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
with open(os.path.join(BASE, 'results', 'compare_curves_table.md'), 'w') as fh:
    fh.write('\n'.join(rows_md) + '\n')
print('\n'.join(rows_md))
print('wrote results/compare_curves_with_original_v1.json, compare_curves_table.md, fig_compare_curves_control_vs_original.png')
