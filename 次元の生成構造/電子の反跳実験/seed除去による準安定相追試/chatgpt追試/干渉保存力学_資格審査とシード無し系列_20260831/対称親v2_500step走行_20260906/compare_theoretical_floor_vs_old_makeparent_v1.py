#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""理論床親（本フォルダ results/）と旧 make_parent 親（../N3_N40_stage123_sweep_20260905/results/）の
500 step 走行比較（指示書 §5・§7、読み出しのみ）。

比較項目: step0 / step1 / onset(>0.05) / 増幅期の指数成長率 / final / max / H_total / global_closure。
成長率: den=N 系列の timeseries から 1e-28 < H⊥/H < 1e-4 の窓で log10 f を最小二乗フィット（log10/step）。
出力: compare_theoretical_floor_vs_old_makeparent.csv、fig_compare_theoretical_floor_vs_old_makeparent.png"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

BASE = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.join(BASE, 'results')
OLD = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'results')


def load(resdir):
    S = {}
    for r in csv.DictReader(open(os.path.join(resdir, 'summary_64bit_with124_N3_N40.csv'))):
        S[(int(r['N']), r['series'])] = r
    T = {}
    CL = {}
    HT = {}
    for r in csv.DictReader(open(os.path.join(resdir, 'timeseries_64bit_with124_N3_N40.csv'))):
        key = (int(r['N']), r['series'])
        T.setdefault(key, {})[int(r['step'])] = float(r['Hperp_frac'])
        CL[key] = max(CL.get(key, 0.0), float(r['global_closure']))
        if r['step'] in ('0', '500'):
            HT.setdefault(key, {})[int(r['step'])] = float(r['H_total'])
    return S, T, CL, HT


def growth(ts):
    steps = sorted(ts)
    xs = [s for s in steps if 1e-28 < ts[s] < 1e-4]
    if len(xs) < 5:
        return float('nan'), len(xs)
    y = np.log10([ts[s] for s in xs])
    A = np.vstack([xs, np.ones(len(xs))]).T
    slope = float(np.linalg.lstsq(A, y, rcond=None)[0][0])
    return slope, len(xs)


Sn, Tn, CLn, HTn = load(NEW)
So, To, CLo, HTo = load(OLD)

rows = []
for N in range(3, 41):
    k = (N, 'N')
    gn, npt_n = growth(Tn[k]); go, npt_o = growth(To[k])
    rows.append(dict(
        N=N,
        new_step0=float(Sn[k]['initial']), old_step0=float(So[k]['initial']),
        new_step1=float(Sn[k]['step1']), old_step1=float(So[k]['step1']),
        new_onset=int(Sn[k]['onset_gt_0.05']), old_onset=int(So[k]['onset_gt_0.05']),
        new_growth_log10_per_step=gn, old_growth_log10_per_step=go,
        new_fit_points=npt_n, old_fit_points=npt_o,
        new_final=float(Sn[k]['final']), old_final=float(So[k]['final']),
        new_max=float(Sn[k]['max']), old_max=float(So[k]['max']),
        new_Htot0=HTn[k][0], new_Htot500=HTn[k][500],
        old_Htot0=HTo[k][0], old_Htot500=HTo[k][500],
        new_closure_max=CLn[k], old_closure_max=CLo[k],
    ))
out_csv = os.path.join(BASE, 'compare_theoretical_floor_vs_old_makeparent.csv')
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

Ns = [r['N'] for r in rows]
fig, axs = plt.subplots(2, 2, figsize=(13, 10))
ax = axs[0, 0]
ax.plot(Ns, [r['old_onset'] if r['old_onset'] >= 0 else np.nan for r in rows], 'o-', label='旧 make_parent 親')
ax.plot(Ns, [r['new_onset'] if r['new_onset'] >= 0 else np.nan for r in rows], 's-', label='理論床親')
floor_new = [r['N'] for r in rows if r['new_onset'] < 0]
ax.set_xlabel('N'); ax.set_ylabel('onset step (H⊥/H > 0.05)')
ax.set_title(f'点火時刻（den=N）。理論床の床維持: N={floor_new}')
ax.grid(alpha=.3); ax.legend()
ax = axs[0, 1]
ax.plot(Ns, [r['old_growth_log10_per_step'] for r in rows], 'o-', label='旧')
ax.plot(Ns, [r['new_growth_log10_per_step'] for r in rows], 's-', label='理論床')
ax.set_xlabel('N'); ax.set_ylabel('増幅率 [log10/step]'); ax.set_title('増幅期の指数成長率（窓 1e-28..1e-4）')
ax.grid(alpha=.3); ax.legend()
ax = axs[1, 0]
ax.plot(Ns, [r['old_final'] for r in rows], 'o-', label='旧 final')
ax.plot(Ns, [r['new_final'] for r in rows], 's-', label='理論床 final')
ax.axhline(1 / 6, color='gray', ls='--', lw=.8, label='1/6')
ax.set_xlabel('N'); ax.set_ylabel('H⊥/H final (step500)'); ax.set_title('飽和値'); ax.grid(alpha=.3); ax.legend()
ax = axs[1, 1]
ax.semilogy(Ns, [r['old_step1'] for r in rows], 'o-', label='旧 f(1)')
ax.semilogy(Ns, [r['new_step1'] for r in rows], 's-', label='理論床 f(1)')
ax.set_xlabel('N'); ax.set_ylabel('H⊥/H at step1'); ax.set_title('step1（床維持の確認、v3失敗系は 1e-1 級だった）')
ax.grid(alpha=.3); ax.legend()
fig.suptitle('理論床親 vs 旧 make_parent 親 — 同一物理正本・500 step・den=N', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_compare_theoretical_floor_vs_old_makeparent.png'), dpi=160)
print('wrote', out_csv)
print('N | onset new/old | growth new/old | final new/old')
for r in rows:
    print(f"{r['N']:2d} | {r['new_onset']:4d}/{r['old_onset']:4d} | {r['new_growth_log10_per_step']:.4f}/{r['old_growth_log10_per_step']:.4f} | {r['new_final']:.3f}/{r['old_final']:.3f}")
