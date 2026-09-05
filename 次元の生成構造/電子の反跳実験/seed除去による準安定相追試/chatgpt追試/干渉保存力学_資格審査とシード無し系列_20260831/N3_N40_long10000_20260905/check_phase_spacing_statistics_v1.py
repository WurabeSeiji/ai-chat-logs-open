#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ガラス終状態の位相分布の規則性検定（読み出しのみ）。
問い（木原 2026-09-06）: 位相の分布にはルールがないのか。等間隔に並んでいないか。
方法（step 10000 の全 N）:
 (1) 位相クラスタ分割（円環ギャップ > 2° で分割）→ クラスタ中心と占有数
 (2) 等間隔検定: 中心が 360°/n_c の一様格子（オフセット最良フィット）に載るかの最大偏差
 (3) 間隔統計: 正規化ギャップ s_j = g_j·n_c/360° の CV
     （等間隔=0／独立一様（ポアソン）≈1／反発型（Wigner系）≈0.5）
 (4) フーリエ構造因子 S_k = |Σ_e e^{ikθ_e}|²/M（k=1..60）の卓越モード
出力: phase_spacing_table_v1.csv、check_phase_spacing_statistics_v1.json"""
import csv
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))

def clusters(ph, gap_thresh):
    order = np.argsort(ph)
    ps = ph[order]
    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    cut = np.flatnonzero(gaps > gap_thresh)
    if cut.size == 0:
        return [(float(np.angle(np.mean(np.exp(1j * ph)))), len(ph))]
    groups = []
    start = (cut[-1] + 1) % len(ps)
    idx = list(order[start:]) + list(order[:start])
    bounds = sorted(((c - start) % len(ps)) for c in cut)
    prev = 0
    for b in bounds:
        mem = idx[prev:b + 1]
        c = float(np.angle(np.sum(np.exp(1j * ph[mem]))))
        groups.append((c, len(mem)))
        prev = b + 1
    return groups

def grid_dev_deg(centers, pitch):
    u = np.array([c / pitch for c in centers])
    mean = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    return math.degrees(float(max(abs((x - mean + 0.5) % 1.0 - 0.5) for x in u) * pitch))

rows = []
detail = {}
for N in range(3, 41):
    d = np.load(os.path.join(BASE, 'results', f'hm_N{N}_den_{N}_states_10000.npz'))
    z = np.asarray(d['Z'][10000], dtype=np.complex128)
    M = len(z)
    ph = np.angle(z)
    g = clusters(ph, math.radians(2.0))
    cs = np.sort(np.array([c for c, _ in g]))
    occ = sorted((n for _, n in g), reverse=True)
    n_c = len(cs)
    gaps = np.diff(np.concatenate([cs, [cs[0] + 2 * math.pi]]))
    s = gaps * n_c / (2 * math.pi)
    cv = float(s.std() / s.mean()) if n_c > 1 else 0.0
    dev_uniform = grid_dev_deg(list(cs), 2 * math.pi / n_c) if n_c > 1 else 0.0
    # フーリエ構造因子（波全体、占有重みそのまま）
    ks = np.arange(1, 61)
    S = np.array([abs(np.sum(np.exp(1j * k * ph))) ** 2 / M for k in ks])
    top = ks[np.argsort(S)[::-1][:5]]
    rows.append([N, M, n_c, round(cv, 3), round(dev_uniform, 2),
                 '|'.join(map(str, occ[:6])), '|'.join(map(str, top)),
                 round(float(S[1]), 1)])
    detail[N] = {'n_clusters': n_c, 'gap_cv': cv, 'dev_from_uniform_deg': dev_uniform,
                 'occupancy': occ, 'centers_deg': [round(math.degrees(c) % 360, 2) for c in cs],
                 'normalized_gaps': [round(float(x), 3) for x in np.sort(s)],
                 'top_fourier_k': [int(k) for k in top],
                 'S_k2': float(S[1]), 'S_first12': [round(float(x), 1) for x in S[:12]]}
    print(f"N={N}: n_c={n_c} gapCV={cv:.3f} 等間隔偏差={dev_uniform:.1f}° "
          f"occ={occ[:5]} topk={list(top)} S(k=2)={S[1]:.0f}", flush=True)

HEADER = ['N', 'M', 'n_clusters', 'gap_cv', 'dev_uniform_deg', 'occ_top6', 'top_fourier_k', 'S_k2']
with open(os.path.join(BASE, 'phase_spacing_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(HEADER); w.writerows(rows)
with open(os.path.join(BASE, 'check_phase_spacing_statistics_v1.json'), 'w') as f:
    json.dump(detail, f, indent=2, default=str)
print('ALL DONE')
