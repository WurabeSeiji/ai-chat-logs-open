#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""§4 集計: 正本出力（full_N3_N40_sweep/states/）から指示書指定の CSV 2 本を作る（読み出しのみ）。
- timeseries_64bit_with124_N3_N40_makeparent.csv: 正本 timeseries のコピー（名称のみ指定名）
- summary_64bit_with124_N3_N40_makeparent.csv: N, series, denominator, onset_gt_0.05,
  initial, step1, step10, step100, final, max, closure0, closure500（timeseries から機械集計）"""
import csv
import os
import shutil

BASE = os.path.dirname(os.path.abspath(__file__))
ST = os.path.join(BASE, 'full_N3_N40_sweep', 'states')
OUT = os.path.join(BASE, 'full_N3_N40_sweep')

src_ts = os.path.join(ST, 'timeseries_64bit_with124_N3_N40.csv')
shutil.copyfile(src_ts, os.path.join(OUT, 'timeseries_64bit_with124_N3_N40_makeparent.csv'))

data = {}
for r in csv.DictReader(open(src_ts)):
    key = (int(r['N']), r['series'], int(r['denominator']))
    d = data.setdefault(key, {})
    d[int(r['step'])] = (float(r['Hperp_frac']), float(r['global_closure']))
rows = []
for (N, series, den), d in sorted(data.items()):
    f = {s: v[0] for s, v in d.items()}
    onset = next((s for s in sorted(f) if f[s] > 0.05), -1)
    rows.append([N, series, den, onset, f[0], f[1], f[10], f[100], f[500], max(f.values()),
                 d[0][1], d[500][1]])
with open(os.path.join(OUT, 'summary_64bit_with124_N3_N40_makeparent.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['N', 'series', 'denominator', 'onset_gt_0.05', 'initial', 'step1', 'step10',
                'step100', 'final', 'max', 'closure0', 'closure500'])
    w.writerows(rows)
print('rows:', len(rows))
shutil.copyfile(os.path.join(ST, 'fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png'),
                os.path.join(OUT, 'fig_Hperp_makeparent_N3_N40_500.png'))
print('POSTPROCESS DONE')
