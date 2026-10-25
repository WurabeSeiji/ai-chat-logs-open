#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照実験の照合（読み出しのみ）: 本フォルダの生成物を 2026-09-05 の原本パッケージの生成物と比較する。

比較対象
  (1) parents/parent_static_N*_makeparent_20260905.npz … v, g, Z0, sigma, residual の bit 一致（np.array_equal）
  (2) results/hm_N{N}_den_{den}_states_500.npz … Z（501×M）の bit 一致。不一致なら最初に食い違う step と |ΔZ| の最大
  (3) results/*.csv, results/analysis_sweep_summary_v1.json … SHA256 一致。CSV が不一致なら列ごとの数値差の最大
  (4) 図 png … SHA256（matplotlib のメタデータで異なり得るので参考値）
出力: results/compare_with_original_v1.json と標準出力の要約。解釈は含まない。"""
import csv
import hashlib
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = ('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/'
        'ai-chat-logs-open/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/'
        '干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def dens_of(N):
    return [N + o for o in (-2, -1, 0, 1, 2) if N + o > 0] + [124]


out = {'original_dir': ORIG, 'control_dir': BASE, 'parents': {}, 'states': {}, 'files': {}}

# (1) 親
pe = pm = 0
for N in range(3, 41):
    name = f'parent_static_N{N:05d}_makeparent_20260905.npz'
    a, b = os.path.join(BASE, 'parents', name), os.path.join(ORIG, 'parents', name)
    if not (os.path.exists(a) and os.path.exists(b)):
        out['parents'][N] = 'missing'; pm += 1; continue
    da, db = np.load(a), np.load(b)
    same = {k: bool(np.array_equal(da[k], db[k])) for k in ('v', 'g', 'Z0', 'sigma', 'residual')}
    same['max_abs_diff'] = {}
    for k in ('v', 'g', 'Z0', 'sigma', 'residual'):
        xa, xb = np.asarray(da[k]), np.asarray(db[k])
        same['max_abs_diff'][k] = float(np.max(np.abs(xa - xb))) if xa.shape == xb.shape else f'shape {list(xa.shape)} vs {list(xb.shape)}'
    out['parents'][N] = same
    if all(v for k, v in same.items() if k != 'max_abs_diff'):
        pe += 1
out['parents_bit_identical'] = pe
out['parents_missing'] = pm

# (2) 状態
se = sm = sd = 0
for N in range(3, 41):
    for den in dens_of(N):
        name = f'hm_N{N}_den_{den}_states_500.npz'
        a, b = os.path.join(BASE, 'results', name), os.path.join(ORIG, 'results', name)
        key = f'N{N}_den{den}'
        if not (os.path.exists(a) and os.path.exists(b)):
            out['states'][key] = 'missing'; sm += 1; continue
        Za, Zb = np.load(a)['Z'], np.load(b)['Z']
        if Za.shape != Zb.shape:
            out['states'][key] = {'identical': False, 'shape_mismatch': [list(Za.shape), list(Zb.shape)]}; sd += 1; continue
        if np.array_equal(Za, Zb):
            out['states'][key] = {'identical': True}; se += 1
        else:
            diff = np.abs(Za - Zb)
            rows = np.flatnonzero(diff.max(axis=1) > 0)
            out['states'][key] = {'identical': False, 'first_diff_step': int(rows[0]),
                                  'max_abs_diff_all': float(diff.max()), 'max_abs_diff_final': float(diff[-1].max()),
                                  'Z0_identical': bool(np.array_equal(Za[0], Zb[0]))}
            sd += 1
out['states_bit_identical'] = se
out['states_different'] = sd
out['states_missing'] = sm

# (3)(4) ファイル
for rel in ('results/timeseries_64bit_with124_N3_N40.csv', 'results/summary_64bit_with124_N3_N40.csv',
            'results/analysis_sweep_summary_v1.json', 'results/RUN_METADATA_N3_N40_stage123.json',
            'results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png',
            'fig_complex_plane_step0_N3_N40_stage123.png', 'fig_complex_plane_final_N3_N40_stage123.png',
            'fig_complex_plane_final_zoom_N3_N40_stage123.png', 'parents/parents_summary.csv'):
    a, b = os.path.join(BASE, rel), os.path.join(ORIG, rel)
    if not (os.path.exists(a) and os.path.exists(b)):
        out['files'][rel] = 'missing'; continue
    ha, hb = sha256(a), sha256(b)
    rec = {'sha256_identical': ha == hb, 'sha256_control': ha, 'sha256_original': hb}
    if ha != hb and rel.endswith('.csv'):
        with open(a) as fa, open(b) as fb:
            ra, rb = list(csv.reader(fa)), list(csv.reader(fb))
        rec['rows_control'], rec['rows_original'] = len(ra), len(rb)
        if len(ra) == len(rb) and ra[0] == rb[0]:
            colmax = {}
            for j, col in enumerate(ra[0]):
                m = 0.0
                for x, y in zip(ra[1:], rb[1:]):
                    try:
                        d = abs(float(x[j]) - float(y[j]))
                    except ValueError:
                        d = 0.0 if x[j] == y[j] else float('inf')
                    if d > m:
                        m = d
                colmax[col] = m
            rec['max_abs_diff_by_column'] = colmax
    out['files'][rel] = rec

os.makedirs(os.path.join(BASE, 'results'), exist_ok=True)
with open(os.path.join(BASE, 'results', 'compare_with_original_v1.json'), 'w') as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print(f"parents: bit-identical {pe}/38, missing {pm}")
print(f"states : bit-identical {se}, different {sd}, missing {sm} (of 228)")
for rel, rec in out['files'].items():
    print('file   :', rel, '->', rec if isinstance(rec, str) else ('IDENTICAL' if rec['sha256_identical'] else 'DIFFERENT'))
print('wrote results/compare_with_original_v1.json')
