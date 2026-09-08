#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=6 (M=15) 対照実験（2026-09-08、木原指示）: 無名化処理をコメントアウトした素の走行。

新しい実験のベースライン。初期値は理論床（parents_symmetric_staged）をそのまま使用し、
無名化→復元は「行わない」（下の該当行をコメントアウトして残す）。
力学は物理正本 run_N3_N40_stage123_v1.py（SHA照合・無変更）、den=N・500 step、N=6 固定。
"""
import hashlib
import os

import numpy as np

from anonymize_rebuild_lib_v2 import anonymize_and_rebuild  # 無名化を戻す時のために import は残す

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
PD = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')
out = os.path.join(BASE, 'results_control_N6')
os.makedirs(out, exist_ok=True)

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

# 初期値: 無名化処理はコメントアウト（= 理論床をそのまま使用）
INJ = ("z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)\n"
       "        # 無名化処理はコメントアウト（対照実験）:\n"
       "        # z0=anonymize_and_rebuild(z0,seed=20260908)")

REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')", f'PARENT_DIR={PD!r}'),
    ("for N in range(3,41):", "for N in [6]:"),
    ("STEPS=500; OFFSETS=(-2,-1,0,1,2)", "STEPS=500; OFFSETS=(0,)"),
    ("pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]",
     "pairs=[(N,'N')]"),
    ("    z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)",
     "    " + INJ),
    ("'N_range':[3,40]", "'N_range':[6,6]"),
    ("'input':'per-N static make_parent parents/parent_static_N*_makeparent_20260905.npz Z0 (N=40 bit-identical to canonical)'",
     "'input':'theoretical-floor parent N=6 (parents_symmetric_staged), anonymize->rebuild DISABLED (control).'"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（N=6・無名化コメントアウト・力学無変更）', flush=True)
g = {'__name__': 'stage123_control_N6', '__file__': ORIG, 'anonymize_and_rebuild': anonymize_and_rebuild, 'np': np}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
