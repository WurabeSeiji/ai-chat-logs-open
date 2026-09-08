#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R反射系 vs Control の機械集計（読み出しのみ）。出力: compare_R124_23_vs_control_N3_N16.csv"""
import csv
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'results')
CTRL = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'results')


def load(path):
    ts, cl, ht = {}, {}, {}
    for r in csv.DictReader(open(path)):
        if r['series'] != 'N':
            continue
        N = int(r['N'])
        ts.setdefault(N, {})[int(r['step'])] = float(r['Hperp_frac'])
        cl[N] = max(cl.get(N, 0.0), float(r['global_closure']))
        ht.setdefault(N, {})[int(r['step'])] = float(r['H_total'])
    return ts, cl, ht


def growth(d):
    xs = [s for s in sorted(d) if 1e-28 < d[s] < 1e-4]
    if len(xs) < 5:
        return float('nan')
    y = np.log10([d[s] for s in xs])
    A = np.vstack([xs, np.ones(len(xs))]).T
    return float(np.linalg.lstsq(A, y, rcond=None)[0][0])


t_ts, t_cl, t_ht = load(os.path.join(RES, 'timeseries_64bit_with124_N3_N40.csv'))
c_ts, c_cl, c_ht = load(os.path.join(CTRL, 'timeseries_64bit_with124_N3_N40.csv'))
rows = []
for N in range(3, 17):
    d, c = t_ts[N], c_ts[N]
    rows.append(dict(N=N,
                     test_onset=next((s for s in sorted(d) if d[s] > 0.05), -1),
                     ctrl_onset=next((s for s in sorted(c) if c[s] > 0.05), -1),
                     test_growth=growth(d), ctrl_growth=growth(c),
                     test_f1=d[1], ctrl_f1=c[1], test_final=d[500], ctrl_final=c[500],
                     test_closure_max=t_cl[N], ctrl_closure_max=c_cl[N],
                     test_Htot_drift=abs(t_ht[N][500] - t_ht[N][0])))
out = os.path.join(BASE, 'compare_R124_23_vs_control_N3_N16.csv')
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('N | onset T/C | growth T/C | f(1) T | final T/C | closure_max T/C | Htot drift T')
for r in rows:
    print(f"{r['N']:2d} | {r['test_onset']:4d}/{r['ctrl_onset']:4d} | {r['test_growth']:.4f}/{r['ctrl_growth']:.4f} | "
          f"{r['test_f1']:.2e} | {r['test_final']:.3f}/{r['ctrl_final']:.3f} | {r['test_closure_max']:.2e}/{r['ctrl_closure_max']:.2e} | {r['test_Htot_drift']:.1e}")
