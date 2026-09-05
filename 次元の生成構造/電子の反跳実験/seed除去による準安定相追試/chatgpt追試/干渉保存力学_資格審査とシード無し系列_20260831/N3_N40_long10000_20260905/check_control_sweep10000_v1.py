#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照ゲート（全38系）: 10000歩スイープの最初の501歩が正本スイープと bit 一致すること、
および N=5 は既存10000歩走行（N5_den5_long10000_20260905）と全10001歩 bit 一致することを検証。
正本側 npz は SHA256 台帳と照合してから使う。
出力: check_control_sweep10000_v1.json（1系でも不一致なら assert で失敗）"""
import hashlib
import json
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
LONG5 = os.path.join(BASE, '..', 'N5_den5_long10000_20260905')

def ledger_of(pkg):
    led = {}
    with open(os.path.join(pkg, 'SHA256SUMS.txt')) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                led[parts[1]] = parts[0]
    return led

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

led = ledger_of(PKG)
out = {'first501_bit_identical': {}, 'N5_full_bit_identical': None}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    assert led[rel] == sha256(os.path.join(PKG, rel)), f'INPUT GATE FAIL: {rel}'
    Zc = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    Zn = np.asarray(np.load(os.path.join(BASE, 'results', f'hm_N{N}_den_{N}_states_10000.npz'))['Z'],
                    dtype=np.complex128)
    ok = bool(np.array_equal(Zn[:501], Zc))
    out['first501_bit_identical'][str(N)] = ok
    print(f'N={N}: first501 bit-identical = {ok}', flush=True)

led5 = ledger_of(LONG5)
rel5 = 'results/hm_N5_den_5_states_10000.npz'
assert led5[rel5] == sha256(os.path.join(LONG5, rel5)), 'INPUT GATE FAIL: N5 long npz'
Z5a = np.asarray(np.load(os.path.join(LONG5, rel5))['Z'], dtype=np.complex128)
Z5b = np.asarray(np.load(os.path.join(BASE, 'results', 'hm_N5_den_5_states_10000.npz'))['Z'],
                 dtype=np.complex128)
out['N5_full_bit_identical'] = bool(np.array_equal(Z5a, Z5b))
print('N=5 full 10001-step bit-identical to prior long run:', out['N5_full_bit_identical'])

with open(os.path.join(BASE, 'check_control_sweep10000_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
assert all(out['first501_bit_identical'].values()) and out['N5_full_bit_identical'], 'CONTROL FAIL'
print('ALL GATES PASSED')
