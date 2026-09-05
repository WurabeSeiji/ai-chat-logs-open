#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""頂点スター閉塞（ガウス束縛候補）の全系列調査（読み出しのみ）。
問い（検討9のゲージ場読みから）: 波は K_N の辺に住むリンク変数であり、格子ゲージ理論なら
頂点に局所保存量（ガウス束縛）があるはず。候補は頂点スター閉塞
  G_v(s) = Σ_{e∋v} z_e(s)²  （頂点 v に接続する N−1 辺の二乗和）
である。測定: (a) 初期に G_v=0 か、(b) 時間発展で G_v は保存されるか
（値がゼロでなくても一定なら保存量）。
入力: 正本スイープの den=N npz（SHA ゲート付き）＋ N=5 の10000歩。
出力: vertex_star_table_v1.csv、check_vertex_star_closure_v1.json"""
import csv
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

def gate(pkg, rel):
    h = hashlib.sha256(open(os.path.join(pkg, rel), 'rb').read()).hexdigest()
    assert ledger_of(pkg)[rel] == h, f'INPUT GATE FAIL: {rel}'

def star_matrix(N):
    """頂点×辺の接続行列（K_N、辺は triu 順）"""
    a, b = np.triu_indices(N, k=1)
    S = np.zeros((N, len(a)))
    for e in range(len(a)):
        S[a[e], e] = 1.0
        S[b[e], e] = 1.0
    return S

def analyze(Z, N):
    """G_v(s) 全時系列。返り値: 統計 dict"""
    S = star_matrix(N)
    sq = Z ** 2                       # (steps+1, M)
    G = sq @ S.T                      # (steps+1, N) 複素
    absG = np.abs(G)
    scale = float(np.mean(np.abs(sq[0])))  # ρ 相当
    drift = np.abs(G - G[0])
    return {
        'scale': scale,
        'G0_max_over_v': float(absG[0].max()),
        'G0_max_rel': float(absG[0].max() / scale),
        'absG_max_over_all': float(absG.max()),
        'absG_max_rel': float(absG.max() / scale),
        'drift_max_rel': float(drift.max() / scale),
        'G0_per_vertex_rel': [float(x / scale) for x in absG[0]],
        'absG_final_rel': [float(x / scale) for x in absG[-1]],
    }

rows = []
detail = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(PKG, rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    r = analyze(Z, N)
    rows.append([N, f"{r['G0_max_rel']:.3e}", f"{r['absG_max_rel']:.3e}", f"{r['drift_max_rel']:.3e}"])
    detail[N] = r
    print(f"N={N}: |G_v(0)|/ρ max={r['G0_max_rel']:.3e}  全時間max={r['absG_max_rel']:.3e}  "
          f"ドリフトmax={r['drift_max_rel']:.3e}", flush=True)

gate(LONG5, 'results/hm_N5_den_5_states_10000.npz')
Z5 = np.asarray(np.load(os.path.join(LONG5, 'results/hm_N5_den_5_states_10000.npz'))['Z'],
                dtype=np.complex128)
r5 = analyze(Z5, 5)
detail['N5_long10000'] = r5
print(f"N=5@10000歩: |G_v(0)|/ρ max={r5['G0_max_rel']:.3e}  全時間max={r5['absG_max_rel']:.3e}  "
      f"ドリフトmax={r5['drift_max_rel']:.3e}")

with open(os.path.join(BASE, 'vertex_star_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'G0_max_rel', 'absG_max_rel', 'drift_max_rel'])
    w.writerows(rows)
with open(os.path.join(BASE, 'check_vertex_star_closure_v1.json'), 'w') as f:
    json.dump(detail, f, indent=2, default=str)
print('ALL DONE')
