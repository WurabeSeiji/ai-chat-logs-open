#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""note記事の図説明の裏取り（読み出しのみ）。
終了時のリングで「大部分が同振幅・少数が極小振幅で散乱」する読みの事実確認:
(1) 終了時の振幅分布（中央値・最小・中央値の1/2未満と1/10未満の本数）
(2) 極小振幅の辺集合が step とともに入れ替わるか（一時的な打ち消し＝フェージングか、
    固定的な死んだ波か）を step 400/450/500 の集合重なりで判定。
入力: N3_N40_stage123_sweep_20260905/results/hm_N40_den_40_states_500.npz
出力: amplitude_fading_check_v1.json"""
import json
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, '..', '電子の反跳実験', 'seed除去による準安定相追試', 'chatgpt追試',
                   '干渉保存力学_資格審査とシード無し系列_20260831',
                   'N3_N40_stage123_sweep_20260905', 'results')

out = {}
for N in (3, 4, 5, 8, 12, 20, 30, 40):
    Z = np.load(os.path.join(RES, f'hm_N{N}_den_{N}_states_500.npz'))['Z']
    entry = {}
    for step in (400, 450, 500):
        a = np.abs(Z[step])
        med = float(np.median(a))
        entry[f'step{step}'] = {
            'median_amp': med,
            'min_amp': float(a.min()),
            'max_amp': float(a.max()),
            'n_below_half_median': int((a < 0.5 * med).sum()),
            'n_below_tenth_median': int((a < 0.1 * med).sum()),
            'n_total': int(a.size),
            'small_set_half': sorted(int(i) for i in np.flatnonzero(a < 0.5 * med)),
        }
    s400 = set(entry['step400']['small_set_half'])
    s450 = set(entry['step450']['small_set_half'])
    s500 = set(entry['step500']['small_set_half'])
    entry['overlap_400_500'] = len(s400 & s500)
    entry['overlap_450_500'] = len(s450 & s500)
    entry['union_size'] = len(s400 | s450 | s500)
    out[f'N{N}'] = entry
    e5 = entry['step500']
    print(f"N={N}: med={e5['median_amp']:.4f} min={e5['min_amp']:.4f} max={e5['max_amp']:.4f} "
          f"<half={e5['n_below_half_median']}/{e5['n_total']} <tenth={e5['n_below_tenth_median']} "
          f"| overlap400&500={entry['overlap_400_500']} union={entry['union_size']}")
path = os.path.join(BASE, 'amplitude_fading_check_v1.json')
with open(path, 'w') as fh:
    json.dump(out, fh, indent=2)
print('wrote', path)
