#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位相正規化版走行の入力ゲート（読み出しのみ）。
results_staticparent_sigmaclock/ の全6分母 npz の Z[0] が静的親 Z0 と bit 一致することを検証。
あわせて step1 の H⊥/H（ミスマッチ測定値）を分母ごとに表示する。不一致は exit 1。"""
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
STATIC_PARENT = '/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904/largeN_splitting_result_v1/parent_static_N40_makeparent_20260904.npz'

sp = np.load(STATIC_PARENT)
Z0 = sp['Z0']
p = Z0.real / np.linalg.norm(Z0.real)
q = Z0.imag - (Z0.imag @ p) * p
q = q / np.linalg.norm(q)
ok = True
for den in (38, 39, 40, 41, 42, 124):
    d = np.load(os.path.join(BASE, 'results_staticparent_sigmaclock', f'hm_N40_den_{den}_states_500.npz'))
    Z = d['Z']
    same = bool(np.array_equal(Z[0], Z0))
    z1 = Z[1]
    zp = z1 - p * np.dot(p, z1) - q * np.dot(q, z1)
    f1 = float(np.vdot(zp, zp).real / np.vdot(z1, z1).real)
    print(f'den={den}: Z[0] bit-identical = {same} | step1 Hperp/H = {f1:.6e}')
    ok &= same
print('INPUT GATE:', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
