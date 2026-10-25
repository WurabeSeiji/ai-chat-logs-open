#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""二つの走行の H⊥/H 曲線照合（読み出しのみ）。対照実験_stage123_忠実コピー_20261025/compare_curves_with_original_v1.py の
入力を一般化したもの（A: 比較する走行、B: 基準の走行）。測る量は同じ: onset、f(0)・f(1)・final・max、log10 f の全 step 最大ずれ、
指数区間（−13 < log10 f < −3）の傾き、H_total 相対ドリフト最大、閉塞最大。
使い方: python compare_curves_two_runs_v1.py tsA.csv tsB.csv labelA labelB out_prefix N1 N2 ..."""
import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

tsA, tsB, labA, labB, out_prefix = sys.argv[1:6]
NS = [int(x) for x in sys.argv[6:]]


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


A = load_ts(tsA); B = load_ts(tsB)
out = {'A': tsA, 'B': tsB, 'runs': {}}
fs = lambda x: '—' if x is None else f'{x:.4f}'
md = [f'| N | den | onset {labA}/{labB} | slope {labA}/{labB} (dec/step) | f(0) {labA}/{labB} | f(1) {labA}/{labB} | final {labA}/{labB} | max\\|Δlog10 f\\| | Htot drift {labA}/{labB} | closure max {labA}/{labB} |',
      '|---|---|---|---|---|---|---|---|---|---|']
ncol = 4 if len(NS) > 6 else len(NS); nrow = int(np.ceil(len(NS) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.6 * nrow), squeeze=False); axs = axs.ravel()
for k, N in enumerate(NS):
    ax = axs[k]
    for den, col in ((N, 'tab:blue'), (124, 'tab:orange')):
        key = (N, den)
        if key not in A or key not in B:
            out['runs'][f'N{N}_den{den}'] = 'missing'; continue
        sa, la = stats(A[key]); sb, lb = stats(B[key])
        dl = float(np.max(np.abs(la - lb)))
        out['runs'][f'N{N}_den{den}'] = dict(A=sa, B=sb, max_abs_dlog10f=dl)
        md.append(f"| {N} | {den} | {sa['onset']} / {sb['onset']} | {fs(sa['slope_dec_per_step'])} / {fs(sb['slope_dec_per_step'])} | "
                  f"{sa['f0']:.1e} / {sb['f0']:.1e} | {sa['f1']:.1e} / {sb['f1']:.1e} | {sa['final']:.4f} / {sb['final']:.4f} | {dl:.2f} | "
                  f"{sa['Htot_rel_drift_max']:.1e} / {sb['Htot_rel_drift_max']:.1e} | {sa['closure_max']:.1e} / {sb['closure_max']:.1e} |")
        t = np.arange(len(la))
        ax.semilogy(t, 10.0 ** la, color=col, linewidth=1.2, label=f'{labA} den={den}')
        ax.semilogy(t, 10.0 ** lb, color=col, linewidth=1.0, linestyle='--', label=f'{labB} den={den}')
    ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.set_title(f'N={N}'); ax.grid(alpha=.25); ax.legend(fontsize=6, loc='lower right')
    if k % ncol == 0: ax.set_ylabel('Hperp/H')
    if k // ncol == nrow - 1: ax.set_xlabel('step')
for k in range(len(NS), nrow * ncol): axs[k].axis('off')
fig.suptitle(f'{labA} (solid) vs {labB} (dashed): Hperp/H, den=N (blue), den=124 (orange)', y=.995)
fig.tight_layout(); fig.savefig(out_prefix + '.png', dpi=160); plt.close(fig)
json.dump(out, open(out_prefix + '.json', 'w'), indent=1, ensure_ascii=False)
open(out_prefix + '.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md)); print('wrote', out_prefix + '.{png,json,md}')
