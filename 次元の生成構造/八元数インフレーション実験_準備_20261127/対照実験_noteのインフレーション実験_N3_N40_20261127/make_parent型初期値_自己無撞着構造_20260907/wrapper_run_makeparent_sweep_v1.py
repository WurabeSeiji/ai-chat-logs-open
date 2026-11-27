#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_parent 型親の全Nスイープ実行ラッパー（指示書 §4、物理正本無変更）。
- 正本 ../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py を SHA256 照合のうえ読み取り、
  PARENT_DIR / OUT の2行だけを置換（置換行数=2 をアサート）して exec。
- 親: parents_actual_N3_N40/（入力ゲート合格・正本と byte 一致）
- 出力: full_N3_N40_sweep/states/（正本の命名のまま。CSV改名・追加列・図は後段の集計プログラムで作る）"""
import hashlib
import os

import numpy as np  # noqa: F401（正本実行環境と同一であることの明示）

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
OLD_PD = "PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')"
OLD_OUT = "OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')"

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
pd = os.path.join(BASE, 'parents_actual_N3_N40')
out = os.path.join(BASE, 'full_N3_N40_sweep', 'states')
os.makedirs(out, exist_ok=True)
new = src.replace(OLD_PD, f'PARENT_DIR={pd!r}').replace(OLD_OUT, f'OUT={out!r}')
diff = sum(1 for a, b in zip(src.splitlines(), new.splitlines()) if a != b)
assert diff == 2 and len(src.splitlines()) == len(new.splitlines()), f'置換行数が2でない: {diff}'
g = {'__name__': 'stage123_makeparent_sweep', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
