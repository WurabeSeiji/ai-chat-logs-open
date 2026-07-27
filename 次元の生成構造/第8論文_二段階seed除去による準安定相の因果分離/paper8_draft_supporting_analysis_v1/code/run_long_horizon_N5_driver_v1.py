#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 草稿用補助実験: 完全無seed N=5 の t=110000 長時間対照

既存の長時間専用ラッパ run_seedless_natural_long_horizon_v1.py（N=40,300 で
検証済み）を一切編集せず import し、モジュール定数 NS / SAMPLE だけを
N=5 を含むよう実行時に拡張して run_long(5) を呼ぶ。

- 記録間隔は条件Aと同じ 25 step
- 出力は既存構造に従い
  paper7_seedless_natural_figures3_4_v1/outputs/long_horizon_110000/raw/N00005_...
- 既存出力があれば上書きせず停止する（原本ラッパの保護をそのまま利用）
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WRAPPER_DIR = os.path.normpath(
    os.path.join(HERE, "..", "..", "paper7_seedless_natural_figures3_4_v1")
)
sys.path.insert(0, WRAPPER_DIR)
os.chdir(WRAPPER_DIR)

import run_seedless_natural_long_horizon_v1 as lh  # noqa: E402

lh.NS = (5, 40, 300)
lh.SAMPLE = dict(lh.SAMPLE)
lh.SAMPLE[5] = 25

result = lh.run_long(5)
print("DONE run_long(5)")
print(result if isinstance(result, dict) else type(result))
