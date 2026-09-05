#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ストロボ図の目視所見の数値検証（読み出しのみ）。
入力: data/hm_N3_den_3_states_500.npz（本体スイープのコピー、SHA256 一致確認済み）。
検証項目:
 (1) 交互性: 各 step s で ||Z[s+3] + Z[s]|| / ||Z[s]||（符号反転なら 0）
     と ||Z[s+6] − Z[s]|| / ||Z[s]||（周期6なら 0）を全区間で算出し、
     どの step から交互軌道に入るかを閾値 1e-6 で判定。
 (2) 等振幅性: 終盤（step 117）の |z| 3本の値と min/max 比。
出力: check_strobe_alternation_v1.json（判定はすべて機械適用）。"""
import json
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(BASE, 'data', 'hm_N3_den_3_states_500.npz'))
assert int(d['denominator']) == 3 and int(d['steps']) == 500
Z = np.asarray(d['Z'], dtype=np.complex128)

TOL = 1e-6
flip = {}   # s -> ||Z[s+3]+Z[s]|| / ||Z[s]||
per6 = {}   # s -> ||Z[s+6]-Z[s]|| / ||Z[s]||
for s in range(0, 495):
    nrm = float(np.linalg.norm(Z[s]))
    flip[s] = float(np.linalg.norm(Z[s + 3] + Z[s])) / nrm
    if s <= 494 - 6 + 6 and s + 6 <= 500:
        per6[s] = float(np.linalg.norm(Z[s + 6] - Z[s])) / nrm

def onset(res):
    """res[s] < TOL が s=最終まで連続して成立し始める最小 s"""
    keys = sorted(res)
    ok = [s for s in keys if res[s] < TOL]
    if not ok:
        return None
    last = keys[-1]
    s0 = last
    for s in reversed(keys):
        if res[s] < TOL:
            s0 = s
        else:
            break
    return s0

onset_flip = onset(flip)
onset_per6 = onset(per6)
amp117 = np.abs(Z[117])
out = {
    'tolerance': TOL,
    'onset_signflip_per_3steps': onset_flip,
    'onset_period6': onset_per6,
    'signflip_residual_samples': {str(s): flip[s] for s in (0, 30, 45, 48, 60, 90, 117, 300, 494) if s in flip},
    'period6_residual_samples': {str(s): per6[s] for s in (0, 30, 45, 48, 60, 90, 117, 300, 494) if s in per6},
    'amplitudes_step117': [float(a) for a in amp117],
    'amp_min_over_max_step117': float(amp117.min() / amp117.max()),
}
with open(os.path.join(BASE, 'check_strobe_alternation_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
