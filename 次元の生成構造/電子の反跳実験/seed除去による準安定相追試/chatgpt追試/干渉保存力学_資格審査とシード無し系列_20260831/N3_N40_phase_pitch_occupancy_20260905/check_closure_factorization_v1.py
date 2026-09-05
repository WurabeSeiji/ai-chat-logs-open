#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""閉塞 Σz²=0 のブロック因数分解の検出（読み出しのみ）。
問い（木原 2026-09-05）: 全体閉塞は、どの程度の塊（部分集合）で成立しているか。
2ブロック（z_e²+z_f²=0）と3ブロック（3項零和）への完全分割を厳密探索で求める。
対象: N=3,4,5 の step0 / step500（N=5 は 10000歩終状態も）。M≤10 なので全探索可能。
方法: 二乗ベクトル s_e=z_e² について、未割当の最小添字を含む2元/3元部分集合で
|Σs| < tol·mean|s| のものを深さ優先で選び、完全分割を探す。tol は 1e-9（厳密級）から
1e-1 まで10倍刻みで緩め、成立した最小 tol と分割・最大残差を記録する。
出力: check_closure_factorization_v1.json"""
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

def find_partition(s, tol_abs, max_block=3):
    """s: 二乗値配列。サイズ 2..max_block のブロックへの完全分割を DFS で探す
    （小さいブロック優先）。返り値: ブロック(添字タプル)のリスト、または None。"""
    from itertools import combinations
    Mn = len(s)
    unassigned = list(range(Mn))

    def dfs(rem):
        if not rem:
            return []
        e0 = rem[0]
        rest = rem[1:]
        for size in range(2, max_block + 1):
            for combo in combinations(rest, size - 1):
                if abs(s[e0] + sum(s[c] for c in combo)) < tol_abs:
                    rr = [x for x in rest if x not in combo]
                    sub = dfs(rr)
                    if sub is not None:
                        return [(e0,) + tuple(int(c) for c in combo)] + sub
        return None

    return dfs(unassigned)

def analyze(z, label):
    s = z.astype(np.complex128) ** 2
    scale = float(np.mean(np.abs(s)))
    total = complex(np.sum(s))
    res = {'label': label, 'M': len(z), 'global_closure_abs': abs(total),
           'mean_abs_square': scale}
    for max_block in (3, 4, 5):
        for exp in range(-12, 0):
            tol = 10.0 ** exp * scale
            part = find_partition(list(s), tol, max_block)
            if part is not None:
                blocks = [tuple(int(i) for i in b) for b in part]
                resid = [float(abs(sum(s[list(b)]))) for b in blocks]
                res.update({'tol_rel': 10.0 ** exp, 'max_block_allowed': max_block,
                            'blocks': blocks,
                            'block_sizes': sorted((len(b) for b in blocks), reverse=True),
                            'max_block_residual_rel': max(r / scale for r in resid)})
                return res
    res.update({'tol_rel': None, 'blocks': None})
    return res

def zero_subsets(z, tol_rel=1e-9):
    """真部分集合で Σz²≈0 のものを全列挙（M≤12 用）。既約性の厳密判定。"""
    from itertools import combinations
    s = z.astype(np.complex128) ** 2
    scale = float(np.mean(np.abs(s)))
    Mn = len(z)
    hits = []
    for k in range(2, Mn):
        for combo in combinations(range(Mn), k):
            if abs(sum(s[c] for c in combo)) < tol_rel * scale:
                hits.append(list(combo))
    return hits

out = []
for N in (3, 4, 5):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(PKG, rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    out.append(analyze(Z[0], f'N={N} step0'))
    out.append(analyze(Z[500], f'N={N} step500'))
gate(LONG5, 'results/hm_N5_den_5_states_10000.npz')
Z5 = np.asarray(np.load(os.path.join(LONG5, 'results/hm_N5_den_5_states_10000.npz'))['Z'],
                dtype=np.complex128)
out.append(analyze(Z5[10000], 'N=5 step10000'))

# 既約性の厳密判定（tol 1e-9 の真部分集合零和の全列挙）
irr = {}
for N in (3, 4, 5):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    for step in (0, 500):
        hits = zero_subsets(Z[step])
        irr[f'N={N} step{step}'] = {'n_zero_proper_subsets': len(hits),
                                    'examples': hits[:6]}
irr['N=5 step10000'] = {'n_zero_proper_subsets': len(zero_subsets(Z5[10000])),
                        'examples': zero_subsets(Z5[10000])[:6]}
for k, v in irr.items():
    print(f"{k}: 零和真部分集合 {v['n_zero_proper_subsets']} 個 {v['examples']}")

with open(os.path.join(BASE, 'check_closure_factorization_v1.json'), 'w') as f:
    json.dump({'partitions': out, 'irreducibility_tol1e-9': irr}, f, indent=2)
for r in out:
    if r['blocks'] is not None:
        print(f"{r['label']}: sizes={r['block_sizes']} tol={r['tol_rel']:.0e} "
              f"max_resid={r['max_block_residual_rel']:.2e} blocks={r['blocks']}")
    else:
        print(f"{r['label']}: 分割見つからず（tol 1e-1 まで）")
print('ALL DONE')
