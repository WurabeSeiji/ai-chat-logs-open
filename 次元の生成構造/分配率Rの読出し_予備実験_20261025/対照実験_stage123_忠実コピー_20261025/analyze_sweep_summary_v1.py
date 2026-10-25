#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""スイープ結果の機械集計（読み出しのみ）。再現論文 第2章 §6 の数値の出所。
入力: results/summary_64bit_with124_N3_N40.csv, results/timeseries_64bit_with124_N3_N40.csv
出力: results/analysis_sweep_summary_v1.json
集計項目: 走行数・0.05交差数、f(0)/f(1)/final/max の範囲、閉塞列の初期値と全域最大、
系列別 onset 表（分母依存の記録）、未交差走行の一覧。解釈は含まない。"""
import csv
import json
import math
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'results')
PARENT_DIR = os.path.join(BASE, 'parents')

with open(os.path.join(RES, 'summary_64bit_with124_N3_N40.csv')) as fh:
    S = list(csv.DictReader(fh))

runs = len(S)
crossed = [r for r in S if int(r['onset_gt_0.05']) >= 0]
f0 = [float(r['initial']) for r in S]
f1 = [float(r['step1']) for r in S]
fin = [float(r['final']) for r in S]
mx = [float(r['max']) for r in S]

onset_by_series = {}
for r in S:
    onset_by_series.setdefault(r['series'], {})[int(r['N'])] = int(r['onset_gt_0.05'])
uncrossed = [[int(r['N']), r['series']] for r in S if int(r['onset_gt_0.05']) < 0]
uncrossed_final = [[int(r['N']), r['series'], float(r['final']), float(r['max'])]
                   for r in S if int(r['onset_gt_0.05']) < 0]

cl_max = 0.0
cl0_max = 0.0
with open(os.path.join(RES, 'timeseries_64bit_with124_N3_N40.csv')) as fh:
    for row in csv.DictReader(fh):
        c = float(row['global_closure'])
        if c > cl_max:
            cl_max = c
        if row['step'] == '0' and c > cl0_max:
            cl0_max = c

out = {
    'runs': runs,
    'crossed_0p05': len(crossed),
    'f0_min': min(f0), 'f0_max': max(f0),
    'f1_min': min(f1), 'f1_max': max(f1),
    'final_min': min(fin), 'final_max': max(fin),
    'max_hperp_min': min(mx), 'max_hperp_max': max(mx),
    'crossed_final_min': min(float(r['final']) for r in crossed),
    'crossed_final_max': max(float(r['final']) for r in crossed),
    'onset_min': min(int(r['onset_gt_0.05']) for r in crossed),
    'onset_max': max(int(r['onset_gt_0.05']) for r in crossed),
    'global_closure_step0_max': cl0_max,
    'global_closure_all_max': cl_max,
    'onset_by_series': onset_by_series,
    'uncrossed': uncrossed,
    'uncrossed_final': uncrossed_final,
}

# 親のσスペクトル最大値（第1章 npz の sigma）と回転角スケール ψ_max=2π·σ_max/den の実測
sigma_max_by_N = {}
psi_max_denN = {}
psi_max_den124 = {}
for N in range(3, 41):
    sig = np.load(os.path.join(PARENT_DIR, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['sigma']
    smax = float(np.max(sig))
    sigma_max_by_N[N] = smax
    psi_max_denN[N] = 2.0 * math.pi * smax / N
    psi_max_den124[N] = 2.0 * math.pi * smax / 124.0
out['sigma_max_by_N'] = sigma_max_by_N
out['psi_max_over_2pi_denN'] = {n: v / (2 * math.pi) for n, v in psi_max_denN.items()}
out['psi_max_over_2pi_den124'] = {n: v / (2 * math.pi) for n, v in psi_max_den124.items()}
path = os.path.join(RES, 'analysis_sweep_summary_v1.json')
with open(path, 'w') as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)
print(json.dumps({k: v for k, v in out.items() if k not in ('onset_by_series',)}, indent=2))
print('wrote', path)
