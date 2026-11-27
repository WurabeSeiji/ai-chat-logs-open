#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""入力ゲート（指示書 §2）: parents_actual_N3_N40/ の38親を検査し、正本親と照合する。
- N=3..40 欠番なし・キー(N,v,g,Z0)・M=N(N-1)/2 整合
- parents_summary.csv との突合（N, M, parent_residual を親npzのresidualと比較）
- 正本 ../N3_N40_stage123_sweep_20260905/parents/ と SHA256 および全配列一致
不一致があれば最大絶対差を報告して exit 1（置換・丸め・再生成はしない）。"""
import csv
import hashlib
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(BASE, 'parents_actual_N3_N40')
CANON = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'parents')

summary = {int(r['N']): r for r in csv.DictReader(open(os.path.join(P, 'parents_summary.csv')))}
ok = True
sha_match = arr_match = 0
for N in range(3, 41):
    name = f'parent_static_N{N:05d}_makeparent_20260905.npz'
    path = os.path.join(P, name)
    if not os.path.exists(path):
        print(f'MISSING: {name}'); ok = False; continue
    d = np.load(path)
    M = N * (N - 1) // 2
    for k in ('N', 'v', 'g', 'Z0') if 'N' in d.files else ('n', 'v', 'g', 'Z0'):
        assert k in d.files, (N, k, d.files)
    nkey = int(d['N']) if 'N' in d.files else int(d['n'])
    assert nkey == N and d['Z0'].shape == (M,) and d['Z0'].dtype == np.complex128, (N, d['Z0'].shape, d['Z0'].dtype)
    s = summary[N]
    assert int(s['M']) == M, (N, s['M'], M)
    if 'residual' in d.files:
        dr = abs(float(d['residual']) - float(s['parent_residual']))
        if dr > 1e-18 and dr / max(float(s['parent_residual']), 1e-300) > 1e-6:
            print(f'N={N}: summary residual mismatch npz={float(d["residual"]):.6e} csv={s["parent_residual"]}'); ok = False
    cp = os.path.join(CANON, name)
    h1 = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    h2 = hashlib.sha256(open(cp, 'rb').read()).hexdigest()
    if h1 == h2:
        sha_match += 1; arr_match += 1
    else:
        c = np.load(cp)
        same = set(c.files) == set(d.files) and all(np.array_equal(c[k], d[k]) for k in c.files)
        if same:
            arr_match += 1
        else:
            md = max(float(np.max(np.abs(np.asarray(d[k], complex) - np.asarray(c[k], complex)))) if d[k].shape == c[k].shape else float('inf') for k in ('Z0',))
            print(f'N={N}: 正本と不一致 max|ΔZ0|={md:.3e}'); ok = False
print(f'欠番なし38本 / summary突合OK / 正本SHA一致 {sha_match}/38, 配列一致 {arr_match}/38')
print('INPUT GATE:', 'PASS' if ok and arr_match == 38 else 'FAIL')
sys.exit(0 if ok and arr_match == 38 else 1)
