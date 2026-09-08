#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ロールバック対照テストの bit 比較: 本フォルダ results/ の全 npz を
正本 N3_N40_stage123_sweep_20260905/results/ の同名ファイルと配列単位で厳密比較する。
判定: 全キー一致かつ全配列 np.array_equal（bit 同一。NaN 位置も equal_nan=True で一致要求）。
"""
import csv
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
MINE = os.path.join(BASE, 'results')
CANON = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'results')

rows = []
all_ok = True
files = sorted(f for f in os.listdir(MINE) if f.endswith('.npz'))
assert files, 'results/ に npz がない'
for fn in files:
    a = np.load(os.path.join(MINE, fn))
    b = np.load(os.path.join(CANON, fn))
    keys_ok = sorted(a.files) == sorted(b.files)
    arr_ok = keys_ok and all(np.array_equal(a[k], b[k], equal_nan=True) for k in a.files)
    maxdiff = ''
    if keys_ok and not arr_ok:
        ds = []
        for k in a.files:
            if a[k].shape != b[k].shape:
                ds.append(f'{k}:shape')
            elif not np.array_equal(a[k], b[k], equal_nan=True):
                d = np.nanmax(np.abs(np.asarray(a[k], dtype=np.complex128) - np.asarray(b[k], dtype=np.complex128)))
                ds.append(f'{k}:{d:.3e}')
        maxdiff = ';'.join(ds)
    ok = keys_ok and arr_ok
    all_ok &= ok
    rows.append({'file': fn, 'keys_match': keys_ok, 'bit_identical': ok, 'diff': maxdiff})
    print(f"{fn}: keys={'OK' if keys_ok else 'NG'} bit={'一致' if ok else ' 不一致 ' + maxdiff}")

with open(os.path.join(BASE, 'compare_bitwise_results.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['file', 'keys_match', 'bit_identical', 'diff'])
    w.writeheader()
    w.writerows(rows)
print(f"総合判定: {len(files)} ファイル中 bit一致 {sum(r['bit_identical'] for r in rows)} — "
      + ('全一致（正本再現）' if all_ok else '不一致あり'))
