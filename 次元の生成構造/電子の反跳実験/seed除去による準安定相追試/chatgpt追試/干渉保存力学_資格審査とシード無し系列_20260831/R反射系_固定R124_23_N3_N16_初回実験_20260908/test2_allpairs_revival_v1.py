#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""テスト2（修正2=全対散乱の検証・読み出しのみ）: 床免疫・インフレーション復活・閉塞・124再帰。
出力: test2_allpairs_revival_results.csv"""
import csv
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CTRL = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'results')


def series(path, den):
    out = {}
    for r in csv.DictReader(open(path)):
        if int(r['denominator']) == den:
            out.setdefault(int(r['N']), {})[int(r['step'])] = float(r['Hperp_frac'])
    return out


def growth(d):
    xs = [s for s in sorted(d) if 1e-28 < d[s] < 1e-4]
    if len(xs) < 5:
        return float('nan')
    y = np.log10([d[s] for s in xs])
    A = np.vstack([xs, np.ones(len(xs))]).T
    return float(np.linalg.lstsq(A, y, rcond=None)[0][0])


test = series(os.path.join(BASE, 'results_allpairs', 'timeseries_64bit_with124_N3_N40.csv'), 124)
ctrl = series(os.path.join(CTRL, 'timeseries_64bit_with124_N3_N40.csv'), 124)
rows = []
for N in range(3, 17):
    d, c = test[N], ctrl[N]
    S = np.asarray(np.load(os.path.join(BASE, 'results_allpairs', f'hm_N{N}_den_124_states_500.npz'))['Z'], complex)
    h = np.einsum('ij,ij->i', S.conj(), S).real
    c2 = np.abs(np.einsum('ij,ij->i', S, S)) / h
    rec = abs(np.vdot(S[376], S[500])) / (np.linalg.norm(S[376]) * np.linalg.norm(S[500]))
    rows.append(dict(N=N, f1_test=d[1], f1_ctrl=c[1],
                     onset_test=next((s for s in sorted(d) if d[s] > 0.05), -1),
                     onset_ctrl=next((s for s in sorted(c) if c[s] > 0.05), -1),
                     growth_test=growth(d), growth_ctrl=growth(c),
                     f500_test=d[500], f500_ctrl=c[500],
                     c2_max=float(c2.max()), rec124_late=float(rec)))
with open(os.path.join(BASE, 'test2_allpairs_revival_results.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
