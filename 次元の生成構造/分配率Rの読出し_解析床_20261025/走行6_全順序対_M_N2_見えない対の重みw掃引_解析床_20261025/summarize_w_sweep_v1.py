#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""w 掃引の集計（読み出しのみ）。w=0（走行5 results_w0）、0.3（走行5 results_w0.3）、本走行 {0.5,0.7,0.9,0.95,0.99}、1.0（走行2 results_all_noadj）
の timeseries CSV と床の manifest から、N・den ごとに w に対する 指数区間の傾き（桁/step）、Δτσ で割った傾き（桁/rad）、
500 step で稼いだ桁 log10(f_500/f_1)、onset、final、床の σ・Δτσ を並べる。図：傾き vs w（den=N と 124 の 2 面）。解釈は含まない。"""
import csv, json, math, os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
UP = os.path.join(BASE, '..')
SRC = {
    0.0:  (os.path.join(UP, '走行5_全順序対_M_N2_見えない対に重みw_解析床_20261025', 'results_w0'),   os.path.join(UP, '走行5_全順序対_M_N2_見えない対に重みw_解析床_20261025', 'floors_w0')),
    0.3:  (os.path.join(UP, '走行5_全順序対_M_N2_見えない対に重みw_解析床_20261025', 'results_w0.3'), os.path.join(UP, '走行5_全順序対_M_N2_見えない対に重みw_解析床_20261025', 'floors_w0.3')),
    0.5:  (os.path.join(BASE, 'results_w0.5'),  os.path.join(BASE, 'floors_w0.5')),
    0.7:  (os.path.join(BASE, 'results_w0.7'),  os.path.join(BASE, 'floors_w0.7')),
    0.9:  (os.path.join(BASE, 'results_w0.9'),  os.path.join(BASE, 'floors_w0.9')),
    0.95: (os.path.join(BASE, 'results_w0.95'), os.path.join(BASE, 'floors_w0.95')),
    0.99: (os.path.join(BASE, 'results_w0.99'), os.path.join(BASE, 'floors_w0.99')),
    1.0:  (os.path.join(UP, '走行2_全順序対_M_N2_隣接なし_解析床_20261025', 'results_all_noadj'), os.path.join(UP, '走行2_全順序対_M_N2_隣接なし_解析床_20261025', 'floors_all_noadj')),
}
NS = (3, 4, 5, 6, 8, 10, 12, 15)


def load_ts(path):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d.setdefault((int(r['N']), int(r['denominator'])), []).append((int(r['step']), float(r['Hperp_frac'])))
    for k in d: d[k].sort()
    return d


def stats(rows):
    t = np.array([r[0] for r in rows]); f = np.array([r[1] for r in rows]); lf = np.log10(np.maximum(f, 1e-300))
    ix = np.flatnonzero(f > 0.05); m = (lf > -13) & (lf < -3)
    slope = float(np.polyfit(t[m], lf[m], 1)[0]) if m.sum() >= 5 else None
    return dict(onset=int(ix[0]) if ix.size else -1, f1=float(f[1]), final=float(f[-1]), fmax=float(f.max()),
                slope=slope, decades_1_to_500=float(lf[-1] - lf[1]))


out = {}
for w, (rdir, fdir) in SRC.items():
    ts = load_ts(os.path.join(rdir, 'timeseries_64bit_with124_N3_N40.csv'))
    man = json.load(open(os.path.join(fdir, 'manifest_floor_all.json')))
    for N in NS:
        sig = man[f'floor_90deg_all_N{N:05d}_analytic.npz']['sigma']
        for den in (N, 124):
            s = stats(ts[(N, den)]); dts = 2 * math.pi * sig / den
            s.update(sigma=sig, dtau_sigma=dts, slope_per_rad=(s['slope'] / dts if s['slope'] is not None else None))
            out[f'w{w}_N{N}_den{den}'] = s

ws = sorted(SRC)
fs = lambda x: '—' if x is None else f'{x:.3f}'
md = ['## 指数区間の傾き（桁/step）: 行＝(N, den)、列＝w', '', '| N | den | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|---|' + '---|' * len(ws)]
for N in NS:
    for den in (N, 124):
        md.append(f'| {N} | {den} | ' + ' | '.join(fs(out[f"w{w}_N{N}_den{den}"]["slope"]) for w in ws) + ' |')
md += ['', '## 500 step で稼いだ桁 log10(f_500/f_1)', '', '| N | den | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|---|' + '---|' * len(ws)]
for N in NS:
    for den in (N, 124):
        md.append(f'| {N} | {den} | ' + ' | '.join(f'{out[f"w{w}_N{N}_den{den}"]["decades_1_to_500"]:.1f}' for w in ws) + ' |')
md += ['', '## onset（H⊥/H > 0.05 の最初の step、−1 は未達）', '', '| N | den | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|---|' + '---|' * len(ws)]
for N in NS:
    for den in (N, 124):
        md.append(f'| {N} | {den} | ' + ' | '.join(str(out[f"w{w}_N{N}_den{den}"]["onset"]) for w in ws) + ' |')
md += ['', '## 床の σ と Δτσ（den=N）', '', '| N | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|' + '---|' * len(ws)]
for N in NS:
    md.append(f'| {N} | ' + ' | '.join(f'{out[f"w{w}_N{N}_den{N}"]["sigma"]:.2f} ({out[f"w{w}_N{N}_den{N}"]["dtau_sigma"]:.1f})' for w in ws) + ' |')
md += ['', '## 傾きを Δτσ で割ったもの（桁/rad）', '', '| N | den | ' + ' | '.join(f'w={w}' for w in ws) + ' |', '|---|---|' + '---|' * len(ws)]
for N in NS:
    for den in (N, 124):
        md.append(f'| {N} | {den} | ' + ' | '.join(fs(out[f"w{w}_N{N}_den{den}"]["slope_per_rad"]) for w in ws) + ' |')
open(os.path.join(BASE, 'results_summary', 'w_sweep_summary.md'), 'w').write('\n'.join(md) + '\n')
json.dump(out, open(os.path.join(BASE, 'results_summary', 'w_sweep_summary.json'), 'w'), indent=1)

fig, axs = plt.subplots(1, 3, figsize=(17, 4.8))
for k, (den_mode, title) in enumerate(((lambda N: N, 'den = N'), (lambda N: 124, 'den = 124'))):
    ax = axs[k]
    for N in NS:
        ys = [out[f'w{w}_N{N}_den{den_mode(N)}']['slope'] for w in ws]
        ax.plot([w for w, y in zip(ws, ys) if y is not None], [y for y in ys if y is not None], 'o-', label=f'N={N}')
    ax.set_xlabel('w (weight of invisible pairs)'); ax.set_ylabel('slope of exponential segment [dec/step]'); ax.set_title(title); ax.grid(alpha=.3); ax.legend(fontsize=7)
ax = axs[2]
for N in NS:
    ys = [out[f'w{w}_N{N}_den{N}']['decades_1_to_500'] for w in ws]
    ax.plot(ws, ys, 'o-', label=f'N={N}')
ax.set_xlabel('w'); ax.set_ylabel('decades gained, step 1 -> 500 (den=N)'); ax.set_title('total growth in 500 steps'); ax.grid(alpha=.3); ax.legend(fontsize=7)
fig.suptitle('RUN6: growth vs weight w of invisible (opposite-side) pairs; analytic floor, M=N^2, seedless', y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'results_summary', 'fig_w_sweep_growth.png'), dpi=150); plt.close(fig)
print('\n'.join(md)); print('wrote results_summary/')
