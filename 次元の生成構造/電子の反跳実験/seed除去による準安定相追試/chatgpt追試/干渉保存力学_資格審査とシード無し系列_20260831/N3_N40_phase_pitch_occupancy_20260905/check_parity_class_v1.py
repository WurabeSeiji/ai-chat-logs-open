#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""等角射影読みのパリティ構造検定（検討7-未検証(i)、読み出しのみ）。
問い: 終状態の二乗 z² の6方向占有は、単一パリティ類（占有方向が相互120°＝
純立方根枠＝直交3枠の影として整合）に落ちるか、両類混在（60°混じり）か。
方法: 各状態の二乗を6方向格子（オフセット最良フィット）に割り当て、
偶類 {0,2,4}・奇類 {1,3,5} の占有数から純度 purity = |n_even − n_odd|/M を算出。
格子割り当ての信頼度は E_total/ρ（証明書つき判定器と同じ量）で併記する。
対象: 全 N=3..40 の step0 / step500（正本スイープ、SHA ゲート付き）＋ N=5 の step10000。
出力: parity_class_table_v1.csv、check_parity_class_v1.json"""
import csv
import hashlib
import json
import math
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

def gate(pkg, rel):
    h = hashlib.sha256(open(os.path.join(pkg, rel), 'rb').read()).hexdigest()
    assert ledger_of(pkg)[rel] == h, f'INPUT GATE FAIL: {rel}'

def parity_of(z):
    s = z.astype(np.complex128) ** 2
    rho = float(np.median(np.abs(s)))
    u = np.angle(s) / (math.pi / 3.0)
    off = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    k = np.round(u - off).astype(int) % 6
    ideal = rho * np.exp(1j * (math.pi / 3.0) * (k + off))
    e_rel = float(np.sum(np.abs(s - ideal))) / rho
    a = [int(np.sum(k == j)) for j in range(6)]
    n_even = a[0] + a[2] + a[4]
    n_odd = a[1] + a[3] + a[5]
    Mn = len(z)
    return {'occ6': a, 'n_even': n_even, 'n_odd': n_odd,
            'purity': abs(n_even - n_odd) / Mn, 'E_over_rho': e_rel}

rows = []
detail = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(PKG, rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    r0 = parity_of(Z[0]); rf = parity_of(Z[500])
    rows.append([N, N * (N - 1) // 2, (N * (N - 1) // 2) % 3,
                 round(r0['purity'], 4), round(r0['E_over_rho'], 3),
                 round(rf['purity'], 4), round(rf['E_over_rho'], 3)])
    detail[N] = {'step0': r0, 'step500': rf}
gate(LONG5, 'results/hm_N5_den_5_states_10000.npz')
Z5 = np.asarray(np.load(os.path.join(LONG5, 'results/hm_N5_den_5_states_10000.npz'))['Z'],
                dtype=np.complex128)
detail['N5_step10000'] = parity_of(Z5[10000])

HEADER = ['N', 'M', 'M_mod_3', 'purity_step0', 'E_rho_step0', 'purity_final', 'E_rho_final']
with open(os.path.join(BASE, 'parity_class_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(HEADER); w.writerows(rows)
with open(os.path.join(BASE, 'check_parity_class_v1.json'), 'w') as f:
    json.dump({'table_header': HEADER, 'table': rows, 'detail': detail}, f, indent=2, default=str)

print(f"{'N':>3} {'M%3':>3} {'pur0':>6} {'E/ρ0':>7} {'purF':>6} {'E/ρF':>7}")
for r in rows:
    print(f'{r[0]:>3} {r[2]:>3} {r[3]:>6} {r[4]:>7} {r[5]:>6} {r[6]:>7}')
d5 = detail['N5_step10000']
print(f"N=5@10000: occ6={d5['occ6']} purity={d5['purity']:.4f} E/rho={d5['E_over_rho']:.3f}")
print('ALL DONE')
