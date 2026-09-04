#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""静的親差し替え走行の入力ゲート（読み出しのみ）。
results_staticparent/ の全6分母 npz の Z[0] が、静的親ファイルの Z0 と bit 一致することを検証。
不一致が1つでもあれば exit 1。"""
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
STATIC_PARENT = '/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904/largeN_splitting_result_v1/parent_static_N40_makeparent_20260904.npz'

Z0 = np.load(STATIC_PARENT)['Z0']
ok = True
for den in (38, 39, 40, 41, 42, 124):
    p = os.path.join(BASE, 'results_staticparent', f'hm_N40_den_{den}_states_500.npz')
    same = bool(np.array_equal(np.load(p)['Z'][0], Z0))
    print(f'den={den}: Z[0] bit-identical to static Z0 = {same}')
    ok &= same
print('INPUT GATE:', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
