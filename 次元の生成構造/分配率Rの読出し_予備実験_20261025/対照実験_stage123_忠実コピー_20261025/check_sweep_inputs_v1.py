#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""スイープ入力ゲート（読み出しのみ）。
全 N=3..40 × 全分母の results npz について Z[0] が各 N の静的親 Z0 と bit 一致することを検証。
あわせて summary CSV から N ごとの交差状況を集計表示する。不一致が1つでもあれば exit 1。"""
import csv
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.join(BASE, 'parents')
RESULT_DIR = os.path.join(BASE, 'results')

ok = True
n_checked = 0
for N in range(3, 41):
    Z0 = np.load(os.path.join(PARENT_DIR, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0']
    dens = [N + o for o in (-2, -1, 0, 1, 2) if N + o > 0] + [124]
    for den in dens:
        p = os.path.join(RESULT_DIR, f'hm_N{N}_den_{den}_states_500.npz')
        same = bool(np.array_equal(np.load(p)['Z'][0], Z0))
        n_checked += 1
        if not same:
            print(f'MISMATCH: N={N} den={den}')
            ok = False
print(f'checked {n_checked} state files against 38 static parents')
print('INPUT GATE:', 'PASS' if ok else 'FAIL')

with open(os.path.join(RESULT_DIR, 'summary_64bit_with124_N3_N40.csv')) as fh:
    rows = list(csv.DictReader(fh))
crossed = sum(1 for r in rows if int(r['onset_gt_0.05']) >= 0)
print(f'summary: {len(rows)} runs, {crossed} crossed 0.05')
sys.exit(0 if ok else 1)
