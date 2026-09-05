#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""無名性（振幅エントロピー）の時間発展の全系列調査（読み出しのみ）。
背景（木原 2026-09-05）: 宇宙の初期値は振幅も含めて極めてコヒーレントで無名性が高く
情報量が極めて少ない波の集合だったのではないか、という予感。現在のシード（make_parent
静的親）は振幅が非等分（例 N=3: {1/2,1/2,1/√2}）で情報を持つため、初期値として
正しいのか疑わしい。
測定: 無名性の定量指標として正規化振幅エントロピー
  S(s) = −Σ_e p_e ln p_e / ln M、 p_e = |z_e(s)|²/Σ|z|²
（S=1 ⟺ 等振幅=振幅が個体識別情報を持たない最大無名状態）を全 N=3..40 の全501歩
＋N=5 の10000歩で算出。判定: (a) S(0) と S(final)、(b) 単調増加か（最大逆行幅）、
(c) 上昇が遷移窓に集中するか。
出力: anonymity_entropy_table_v1.csv、check_anonymity_entropy_v1.json"""
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

def entropy_series(Z):
    P = np.abs(Z) ** 2
    P = P / P.sum(axis=1, keepdims=True)
    Mn = Z.shape[1]
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.where(P > 0, np.log(P), 0.0)
    return (-(P * L).sum(axis=1)) / math.log(Mn)

rows = []
detail = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(PKG, rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    S = entropy_series(Z)
    dS = np.diff(S)
    worst_back = float(dS.min())
    rows.append([N, round(float(S[0]), 6), round(float(S[-1]), 6),
                 f'{1.0 - float(S[-1]):.2e}', f'{worst_back:.2e}',
                 int(np.argmax(dS)) if len(dS) else -1])
    detail[N] = {'S0': float(S[0]), 'S_final': float(S[-1]),
                 'one_minus_Sfinal': float(1.0 - S[-1]),
                 'max_single_step_rise_at': int(np.argmax(dS)),
                 'worst_backstep': worst_back,
                 'S_sampled': {str(s): float(S[s]) for s in range(0, 501, 50)}}
    print(f"N={N}: S0={S[0]:.4f} → S500={S[-1]:.6f} (1−S={1-S[-1]:.2e}) "
          f"最大上昇step={int(np.argmax(dS))} 最大逆行={worst_back:.2e}", flush=True)

gate(LONG5, 'results/hm_N5_den_5_states_10000.npz')
Z5 = np.asarray(np.load(os.path.join(LONG5, 'results/hm_N5_den_5_states_10000.npz'))['Z'],
                dtype=np.complex128)
S5 = entropy_series(Z5)
detail['N5_long10000'] = {'S0': float(S5[0]), 'S_10000': float(S5[-1]),
                          'one_minus_Sfinal': float(1.0 - S5[-1]),
                          'worst_backstep': float(np.diff(S5).min())}
print(f"N=5@10000: S0={S5[0]:.4f} → S={S5[-1]:.6f} (1−S={1-S5[-1]:.2e}) "
      f"最大逆行={float(np.diff(S5).min()):.2e}")

with open(os.path.join(BASE, 'anonymity_entropy_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'S0', 'S_final', 'one_minus_Sfinal', 'worst_backstep', 'max_rise_step'])
    w.writerows(rows)
with open(os.path.join(BASE, 'check_anonymity_entropy_v1.json'), 'w') as f:
    json.dump(detail, f, indent=2, default=str)
print('ALL DONE')
