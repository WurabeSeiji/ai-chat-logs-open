#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""§6 比較: make_parent 型（本フォルダ）vs 高対称理論床（../対称親v2_500step走行_20260906/results、読取専用）。
den=N。比較量: onset / step1,10,100,500 / 初期 |Σz|,|Σz²| / 振幅CV / Z4位相残差 / 増幅率
（timeseries の 1e-28<f<1e-4 窓の log10 最小二乗勾配）。
出力: compare_makeparent_vs_theoretical_floor_N3_N40.csv、fig_compare_makeparent_vs_theoretical_floor_onset_growth.png"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

BASE = os.path.dirname(os.path.abspath(__file__))
MP_TS = os.path.join(BASE, 'full_N3_N40_sweep', 'timeseries_64bit_with124_N3_N40_makeparent.csv')
TF = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'results')
MP_PARENTS = os.path.join(BASE, 'parents_actual_N3_N40')
TF_PARENTS = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')
AUDIT = os.path.join(BASE, 'full_N3_N40_sweep', 'parent_structure_N3_N40.csv')


def series_f(path):
    out = {}
    for r in csv.DictReader(open(path)):
        if r['series'] != 'N':
            continue
        out.setdefault(int(r['N']), {})[int(r['step'])] = float(r['Hperp_frac'])
    return out


def growth(ts):
    xs = [s for s in sorted(ts) if 1e-28 < ts[s] < 1e-4]
    if len(xs) < 5:
        return float('nan')
    y = np.log10([ts[s] for s in xs])
    A = np.vstack([xs, np.ones(len(xs))]).T
    return float(np.linalg.lstsq(A, y, rcond=None)[0][0])


def parent_stats(pdir, N):
    z = np.asarray(np.load(os.path.join(pdir, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'], complex)
    th = np.angle(z); r = np.abs(z)
    phi = th[0]
    d = ((th - phi + np.pi / 4) % (np.pi / 2)) - np.pi / 4
    return abs(z.sum()), abs(z @ z), float(np.std(r) / np.mean(r)), float(np.max(np.abs(d)))


mp = series_f(MP_TS)
tf = series_f(os.path.join(TF, 'timeseries_64bit_with124_N3_N40.csv'))
audit = {int(r['N']): r for r in csv.DictReader(open(AUDIT))}

rows = []
for N in range(3, 41):
    m, t = mp[N], tf[N]
    om = next((s for s in sorted(m) if m[s] > 0.05), -1)
    ot = next((s for s in sorted(t) if t[s] > 0.05), -1)
    smp = parent_stats(MP_PARENTS, N)
    stf = parent_stats(TF_PARENTS, N)
    rows.append(dict(N=N,
                     mp_onset=om, tf_onset=ot,
                     mp_step1=m[1], tf_step1=t[1], mp_step10=m[10], tf_step10=t[10],
                     mp_step100=m[100], tf_step100=t[100], mp_step500=m[500], tf_step500=t[500],
                     mp_growth=growth(m), tf_growth=growth(t),
                     mp_abs_sum_z=smp[0], tf_abs_sum_z=stf[0],
                     mp_abs_sum_z2=smp[1], tf_abs_sum_z2=stf[1],
                     mp_amp_CV=smp[2], tf_amp_CV=stf[2],
                     mp_Z4_resid_rad=smp[3], tf_Z4_resid_rad=stf[3]))
out_csv = os.path.join(BASE, 'full_N3_N40_sweep', 'compare_makeparent_vs_theoretical_floor_N3_N40.csv')
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

Ns = [r['N'] for r in rows]
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]
ax.plot(Ns, [r['mp_onset'] if r['mp_onset'] >= 0 else np.nan for r in rows], 'o-', label='make_parent型')
ax.plot(Ns, [r['tf_onset'] if r['tf_onset'] >= 0 else np.nan for r in rows], 's-', label='高対称理論床')
ax.set_xlabel('N'); ax.set_ylabel('onset step (>0.05)'); ax.grid(alpha=.3); ax.legend()
ax.set_title('点火時刻（den=N。N=3: make_parent型=45歩で点火 / 理論床=床）')
ax = axs[1]
ax.plot(Ns, [r['mp_growth'] for r in rows], 'o-', label='make_parent型')
ax.plot(Ns, [r['tf_growth'] for r in rows], 's-', label='高対称理論床')
ax.set_xlabel('N'); ax.set_ylabel('増幅率 [log10/step]'); ax.grid(alpha=.3); ax.legend()
ax.set_title('増幅期の指数成長率（窓 1e-28..1e-4）')
fig.suptitle('make_parent型親 vs 高対称理論床親 — 同一物理正本・500step・den=N', y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'full_N3_N40_sweep', 'fig_compare_makeparent_vs_theoretical_floor_onset_growth.png'), dpi=160)
print('N | onset mp/tf | growth mp/tf | step500 mp/tf')
for r in rows[:8] + rows[-3:]:
    print(f"{r['N']:2d} | {r['mp_onset']:4d}/{r['tf_onset']:4d} | {r['mp_growth']:.4f}/{r['tf_growth']:.4f} | {r['mp_step500']:.3f}/{r['tf_step500']:.3f}")
print('COMPARE DONE')
