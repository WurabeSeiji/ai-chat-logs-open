#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ロールバック対照テスト（2026-09-08、木原指示）: ラッパーを正本と同じ実装に戻し、N=3..6・500 step で
正本 N3_N40_stage123_sweep_20260905/results と同一結果が得られるかを確認する。

方式: 物理正本 run_N3_N40_stage123_v1.py（SHA照合）を読み取り、以下の宣言行だけを置換して exec。
 - OUT → 本フォルダ results/（正本 results の上書き防止）
 - N範囲 → 3..6（走行削減）、メタデータの N_range 表記も同値に
散乱注入なし・PARENT_DIR は正本のまま（正本 parents を使用）・STEPS=500 不変。物理は完全に正本と同一。
"""
import hashlib
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

out = os.path.join(BASE, 'results')
os.makedirs(out, exist_ok=True)

REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("for N in range(3,41):", "for N in range(3,7):"),
    ("'N_range':[3,40]", "'N_range':[3,6]"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（散乱注入なし・物理無変更）')
g = {'__name__': 'stage123_rollback_control', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
