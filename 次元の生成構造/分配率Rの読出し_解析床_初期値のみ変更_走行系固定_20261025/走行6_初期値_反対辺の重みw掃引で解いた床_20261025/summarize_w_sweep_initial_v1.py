#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期値の w 掃引の集計（読み出しのみ、走行系は全て走行1 と同一）。w=0（走行5′ results_w0 ＝ 走行1）、0.3（走行5′）、本走行 {0.5,0.7,0.9,0.95,0.99,1.0}。
N・den ごとに w に対する f(1)、台地の高さ（step 5〜50 の中央値）、onset、500 step で稼いだ桁 log10(f_500/f_1)、final、床の σ を並べる。図：log10 f(t) の曲線を w で重ねる（den=N、N=4,6,8,10）。解釈は含まない。"""
import csv, json, math, os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
R5 = os.path.join(BASE, '..', '走行5_初期値_反対辺に重みwで解いた床_w0_w0.3_20261025')
SRC = {0.0: (os.path.join(R5, 'results_w0'), os.path.join(R5, 'floors_init_w0')),
       0.3: (os.path.join(R5, 'results_w0.3'), os.path.join(R5, 'floors_init_w0.3'))}
for w in (0.5, 0.7, 0.9, 0.95, 0.99, 1.0):
    SRC[w] = (os.path.join(BASE, f'results_w{w}'), os.path.join(BASE, f'floors_init_w{w}'))
NS = (3, 4, 5, 6, 8, 10, 12, 15)


def load_ts(path):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d.setdefault((int(r['N']), int(r['denominator'])), {})[int(r['step'])] = float(r['Hperp_frac'])
    return d


out = {}; curves = {}
for w, (rdir, fdir) in SRC.items():
    ts = load_ts(os.path.join(rdir, 'timeseries_64bit_with124_N3_N40.csv'))
    man = json.load(open(os.path.join(fdir, 'manifest_floor_all.json')))
    for N in NS:
        sig = man[f'floor_90deg_all_N{N:05d}_analytic.npz']['sigma']
        for den in (N, 124):
            f = ts[(N, den)]; arr = np.array([f[s] for s in range(501)])
            ix = np.flatnonzero(arr > 0.05)
            out[f'w{w}_N{N}_den{den}'] = dict(f1=float(arr[1]), plateau_med_5_50=float(np.median(arr[5:51])), onset=int(ix[0]) if ix.size else -1,
                                             decades_1_to_500=float(math.log10(max(arr[-1], 1e-300)) - math.log10(max(arr[1], 1e-300))), final=float(arr[-1]), sigma=sig)
            curves[(w, N, den)] = arr

ws = sorted(SRC)
def table(key, fmt):
    rows = ['| N | den | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|---|' + '---|' * len(ws)]
    for N in NS:
        for den in (N, 124):
            rows.append(f'| {N} | {den} | ' + ' | '.join(fmt(out[f"w{w}_N{N}_den{den}"][key]) for w in ws) + ' |')
    return '\n'.join(rows)
md = ['## f(1)（第 1 歩の H⊥/H）', '', table('f1', lambda x: f'{x:.1e}'), '',
      '## 台地の高さ（step 5〜50 の中央値）', '', table('plateau_med_5_50', lambda x: f'{x:.1e}'), '',
      '## onset（H⊥/H > 0.05 の最初の step、−1 は未達）', '', table('onset', str), '',
      '## final（step 500）', '', table('final', lambda x: f'{x:.3f}'), '',
      '## 床の σ（den=N の行のみ）', '', '| N | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|' + '---|' * len(ws)]
for N in NS:
    md.append(f'| {N} | ' + ' | '.join(f'{out[f"w{w}_N{N}_den{N}"]["sigma"]:.2f}' for w in ws) + ' |')
open(os.path.join(BASE, 'results_summary', 'w_sweep_initial_summary.md'), 'w').write('\n'.join(md) + '\n')
json.dump(out, open(os.path.join(BASE, 'results_summary', 'w_sweep_initial_summary.json'), 'w'), indent=1)

fig, axs = plt.subplots(2, 4, figsize=(20, 8.5)); axs = axs.ravel()
for k, (N, den) in enumerate([(4, 4), (6, 6), (8, 8), (10, 10), (4, 124), (6, 124), (8, 124), (10, 124)]):
    ax = axs[k]
    for w in ws:
        ax.semilogy(np.arange(501), np.maximum(curves[(w, N, den)], 1e-300), linewidth=1.0, label=f'w={w}')
    ax.set_ylim(1e-34, 3); ax.set_xlim(0, 500); ax.set_title(f'N={N}, den={den}'); ax.grid(alpha=.3); ax.legend(fontsize=6, loc='lower right')
    if k % 4 == 0: ax.set_ylabel('Hperp/H')
    if k // 4 == 1: ax.set_xlabel('step')
fig.suptitle('RUN6 redo: INITIAL VALUE = floor solved with invisible-pair weight w; DYNAMICS = RUN1 unchanged (standard adjacency)', y=0.995)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'results_summary', 'fig_w_sweep_initial_curves.png'), dpi=150); plt.close(fig)
print('\n'.join(md)); print('wrote results_summary/')
