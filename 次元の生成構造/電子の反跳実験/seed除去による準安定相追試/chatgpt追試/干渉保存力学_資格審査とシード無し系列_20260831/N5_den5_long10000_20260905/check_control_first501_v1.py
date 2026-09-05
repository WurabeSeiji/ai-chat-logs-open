#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照ゲート: 10000歩走行の最初の501歩（Z[0..500]）が正本スイープの
hm_N5_den_5_states_500.npz と bit 一致することを検証（読み出しのみ）。
正本 npz は SHA256 で正本台帳と照合してから使う。timeseries も同区間で一致確認。
出力: check_control_first501_v1.json"""
import csv
import hashlib
import json
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
REL = 'results/hm_N5_den_5_states_500.npz'

ledger = {}
with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]
h = hashlib.sha256(open(os.path.join(PKG, REL), 'rb').read()).hexdigest()
assert ledger[REL] == h, 'INPUT GATE FAIL: canonical npz SHA mismatch'

can = np.load(os.path.join(PKG, REL))
new = np.load(os.path.join(BASE, 'results', 'hm_N5_den_5_states_10000.npz'))
Zc = np.asarray(can['Z'], dtype=np.complex128)
Zn = np.asarray(new['Z'], dtype=np.complex128)
assert Zc.shape == (501, 10) and Zn.shape == (10001, 10)
bit_identical_states = bool(np.array_equal(Zn[:501], Zc))

# timeseries の同区間一致（Hperp/H_total/closure、文字列レベルでなく float 一致）
def load_ts(path, filt):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            if filt(row):
                out[int(row['step'])] = (float(row['Hperp_frac']), float(row['H_total']), float(row['global_closure']))
    return out
ts_can = load_ts(os.path.join(PKG, 'results', 'timeseries_64bit_with124_N3_N40.csv'),
                 lambda r: r['N'] == '5' and r['denominator'] == '5')
ts_new = load_ts(os.path.join(BASE, 'results', 'timeseries_64bit_N5_den5_steps10000.csv'),
                 lambda r: True)
ts_identical = all(ts_can[s] == ts_new[s] for s in range(501))

out = {'canonical_sha256': h,
       'states_first501_bit_identical': bit_identical_states,
       'timeseries_first501_identical': ts_identical}
with open(os.path.join(BASE, 'check_control_first501_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
assert bit_identical_states and ts_identical, 'CONTROL FAIL'
