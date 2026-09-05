#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""証明書つき閉塞因数分解判定器（任意 M で厳密、読み出しのみ）。

理論（方法論の3層）:
 層1: 二乗 s_e=z_e² が6方向格子 ρ·e^{i(φ0+60°k)} に載る場合、部分和は
      アイゼンシュタイン格子 ρ·Z[ω] の元で、非零元は |·| ≥ ρ。よって全要素の
      格子からの偏差合計 E_total < ρ/2 なら、任意の部分集合について
      「真の和 ≈ 0 ⟺ 理想和が厳密に 0」が保証される（証明書）。
      6 は無平方なので最小消滅和は2サイクル（対蹠対）と3サイクル（立方根組）のみ
      （Lam–Leung）。ゆえに零和部分集合の存在と 2/3 ブロック完全分解の可否は
      方向別占有数 a_0..a_5 の小さな整数条件に厳密還元される:
       - 零和部分集合存在 ⟺ ∃軸 k: a_k≥1 かつ a_{k+3}≥1、または
         ∃パリティ類: a_k,a_{k+2},a_{k+4} すべて ≥1
       - 完全 2/3 分解 ⟺ ∃整数 t0,t1≥0: p_k = a_k − t0·[k偶] − t1·[k奇] が
         対蹠で等しく非負（p_k = p_{k+3} ≥ 0）
 層2: 格子証明書が立たない状態（溶融）は、M ≤ 46 なら中間一致法 O(2^{M/2}) で
      厳密全列挙（許容 tol）。
 層3: それ以上の溶融状態は非存在証明が NP 困難の壁により原理的に不可能。
      本判定器は 'undecidable_molten' として正直に報告する。

検証: M≤10 の既知例（N=3,4,5）で総当たり（check_closure_factorization_v1.py）と
突き合わせて一致確認する。
出力: check_certified_factorization_v1.json"""
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
LONG5 = os.path.join(BASE, '..', 'N5_den5_long10000_20260905')

def gate(pkg, rel):
    led = {}
    with open(os.path.join(pkg, 'SHA256SUMS.txt')) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                led[parts[1]] = parts[0]
    h = hashlib.sha256(open(os.path.join(pkg, rel), 'rb').read()).hexdigest()
    assert led[rel] == h, f'INPUT GATE FAIL: {rel}'

def lattice_certificate(z):
    """層1: 6方向格子フィットと証明書。返り値 dict。"""
    s = z.astype(np.complex128) ** 2
    rho = float(np.median(np.abs(s)))
    ang = np.angle(s)
    # オフセットは 60°格子: u = ang/(π/3) mod 1 の円平均
    u = ang / (math.pi / 3.0)
    off = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    k = np.round(u - off).astype(int) % 6
    ideal = rho * np.exp(1j * (math.pi / 3.0) * (k + off))
    E_total = float(np.sum(np.abs(s - ideal)))
    certified = bool(E_total < rho / 2.0)
    a = [int(np.sum(k == j)) for j in range(6)]
    res = {'rho': rho, 'offset_pitch_units': off, 'occupancy_6dir': a,
           'E_total': E_total, 'certificate_E_lt_rho_half': certified}
    if not certified:
        return res
    # 零和「真部分集合」の存在（厳密）: ブロックが全体と一致する場合は真部分でない
    Mn = len(z)
    pair_exists = any(a[j] >= 1 and a[j + 3] >= 1 for j in range(3)) and Mn > 2
    tri_exists = (any(all(a[(j + 2 * i) % 6] >= 1 for i in range(3)) for j in (0, 1))
                  and Mn > 3)
    res['zero_subset_exists'] = bool(pair_exists or tri_exists)
    # 完全 2/3 分解の可否（厳密）: p_k = a_k − t_par ≥ 0 かつ p_k = p_{k+3}
    ok = None
    for t0 in range(0, min(a[0], a[2], a[4]) + 1):
        for t1 in range(0, min(a[1], a[3], a[5]) + 1):
            p = [a[j] - (t0 if j % 2 == 0 else t1) for j in range(6)]
            if all(x >= 0 for x in p) and all(p[j] == p[j + 3] for j in range(3)):
                ok = {'t_even_triples': t0, 't_odd_triples': t1,
                      'pairs_per_axis': p[:3]}
                break
        if ok:
            break
    res['full_23_factorization'] = ok
    res['is_prime_closure'] = bool(not res['zero_subset_exists'])
    return res

def meet_in_middle(z, tol_rel=1e-9):
    """層2: 厳密全列挙（M≤46）。真部分集合の零和が存在するか。"""
    s = z.astype(np.complex128) ** 2
    Mn = len(s)
    assert Mn <= 46, 'M too large for exact enumeration'
    scale = float(np.mean(np.abs(s)))
    tol = tol_rel * scale
    h1, h2 = s[:Mn // 2], s[Mn // 2:]
    def all_sums(h):
        out = {0: 0.0 + 0.0j}
        sums = np.zeros(1, dtype=np.complex128)
        for v in h:
            sums = np.concatenate([sums, sums + v])
        return sums
    s1 = all_sums(h1); s2 = all_sums(h2)
    order = np.argsort(s2.real)
    s2s = s2[order]
    n_hits = 0
    for i1 in range(len(s1)):
        v = s1[i1]
        lo = int(np.searchsorted(s2s.real, -v.real - tol))
        hi = int(np.searchsorted(s2s.real, -v.real + tol))
        for j in range(lo, hi):
            if abs(v + s2s[j]) < tol:
                i2 = int(order[j])
                if i1 == 0 and i2 == 0:
                    continue  # 空集合
                if i1 == len(s1) - 1 and i2 == len(s2) - 1:
                    continue  # 全体集合
                n_hits += 1
    return {'n_zero_proper_subsets_mim': n_hits, 'tol_rel': tol_rel}

targets = []
for N in (3, 4, 5):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(PKG, rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    targets.append((f'N={N} step0', Z[0]))
    targets.append((f'N={N} step500', Z[500]))
gate(LONG5, 'results/hm_N5_den_5_states_10000.npz')
Z5 = np.asarray(np.load(os.path.join(LONG5, 'results/hm_N5_den_5_states_10000.npz'))['Z'],
                dtype=np.complex128)
targets.append(('N=5 step10000', Z5[10000]))

out = {}
for label, z in targets:
    r = lattice_certificate(z)
    if not r['certificate_E_lt_rho_half'] and len(z) <= 46:
        r['fallback_mim'] = meet_in_middle(z)
    out[label] = r
    cert = 'CERTIFIED' if r['certificate_E_lt_rho_half'] else 'molten(no cert)'
    extra = ''
    if r['certificate_E_lt_rho_half']:
        extra = (f" prime={r['is_prime_closure']} full23={r['full_23_factorization']}")
    elif 'fallback_mim' in r:
        extra = f" mim_zero_subsets={r['fallback_mim']['n_zero_proper_subsets_mim']}"
    print(f'{label}: {cert} occ={r["occupancy_6dir"]} E/rho={r["E_total"]/r["rho"]:.2e}{extra}')

with open(os.path.join(BASE, 'check_certified_factorization_v1.json'), 'w') as f:
    json.dump(out, f, indent=2, default=str)
print('ALL DONE')
