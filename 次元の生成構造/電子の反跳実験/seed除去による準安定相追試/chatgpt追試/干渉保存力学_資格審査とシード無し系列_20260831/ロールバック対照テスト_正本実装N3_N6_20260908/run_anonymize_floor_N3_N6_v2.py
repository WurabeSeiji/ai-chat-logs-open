#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""無名化→復元 走行（2026-09-08、木原指示）: 高対称親（理論床）N=3..6 を無名化→復元してから走行。

初期値だけを anonymize_and_rebuild(z0) に通す（走行前にシャッフル→90度系列情報だけで完全リスト復元）。
N=3..6 は復元が元と bit 一致するので、走行結果も無名化しない場合と bit 一致するはず＝無名性ゲート。
力学は物理正本 run_N3_N40_stage123_v1.py（SHA照合・無変更）。親は 対称親v2/parents_symmetric_staged（理論床）。
"""
import hashlib
import os

import numpy as np

from anonymize_rebuild_lib_v2 import anonymize_and_rebuild

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
PD = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')
out = os.path.join(BASE, 'results_anonymize_floor')
os.makedirs(out, exist_ok=True)

# 事前ゲート: 無名化→復元が全 N で bit 一致することを確認（走行の正当性保証）
print('=== 事前ゲート: 無名化→復元 bit 一致確認 ===', flush=True)
for N in (3, 4, 5, 6):
    z0 = np.asarray(np.load(os.path.join(PD, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'], np.complex128)
    zc = anonymize_and_rebuild(z0, seed=20260908)
    assert np.array_equal(zc, z0), f'N={N}: 無名化→復元が bit 一致しない'
    print(f'  N={N}: bit 一致 OK', flush=True)

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

INJ = ("z0=anonymize_and_rebuild(np.array(np.load(os.path.join(PARENT_DIR,"
       "f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True),seed=20260908)")
REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')", f'PARENT_DIR={PD!r}'),
    ("for N in range(3,41):", "for N in range(3,7):"),
    ("    z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)",
     "    " + INJ),
    ("'N_range':[3,40]", "'N_range':[3,6]"),
    ("'input':'per-N static make_parent parents/parent_static_N*_makeparent_20260905.npz Z0 (N=40 bit-identical to canonical)'",
     "'input':'theoretical-floor parents (parents_symmetric_staged) passed through anonymize->rebuild (shuffle + reconstruct complete list from 90-degree phase + distance-class amplitude). N=3..6 bit-identical to un-anonymized.'"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（初期値を無名化→復元・力学は無変更）', flush=True)
g = {'__name__': 'stage123_anonymize_floor', '__file__': ORIG, 'anonymize_and_rebuild': anonymize_and_rebuild}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
