#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=5 長時間走行（10000歩）終状態の位相幾何の読み出し。
入力: ../N5_den5_long10000_20260905/results/hm_N5_den_5_states_10000.npz
（同フォルダ SHA256SUMS.txt と照合）。クラスタ・ギャップ・占有・振幅と、
60° 格子および 36°(=180/5) 格子への最良フィット偏差を機械算出。
出力: check_N5_long_geometry_v1.json"""
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
LONG = os.path.join(BASE, '..', 'N5_den5_long10000_20260905')
REL = 'results/hm_N5_den_5_states_10000.npz'

ledger = {}
with open(os.path.join(LONG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]
h = hashlib.sha256(open(os.path.join(LONG, REL), 'rb').read()).hexdigest()
assert ledger[REL] == h, 'INPUT GATE FAIL'

d = np.load(os.path.join(LONG, REL))
Z = np.asarray(d['Z'], dtype=np.complex128)

def clusters(z, gap_thresh):
    ph = np.angle(z)
    order = np.argsort(ph)
    ps = ph[order]
    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    cut = np.flatnonzero(gaps > gap_thresh)
    if cut.size == 0:
        return [(float(np.angle(np.mean(np.exp(1j * ph)))), list(range(len(z))))]
    groups = []
    start = (cut[-1] + 1) % len(ps)
    idx = list(order[start:]) + list(order[:start])
    bounds = sorted(((c - start) % len(ps)) for c in cut)
    prev = 0
    for b in bounds:
        mem = idx[prev:b + 1]
        c = float(np.angle(np.sum(np.exp(1j * ph[mem]))))
        groups.append((c, mem))
        prev = b + 1
    return groups

def grid_fit(centers, pitch):
    u = np.array([c / pitch for c in centers])
    mean = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    dev = np.array([abs((x - mean + 0.5) % 1.0 - 0.5) for x in u]) * pitch
    return mean % 1.0, float(dev.max())

out = {}
for step in (500, 2000, 5000, 10000):
    z = Z[step]
    g = clusters(z, math.radians(9.0))
    cs = sorted(math.degrees(c) % 360.0 for c, _ in g)
    gaps = [round((cs[(i + 1) % len(cs)] - cs[i]) % 360.0, 3) for i in range(len(cs))]
    amps = np.abs(z)
    _, dev60 = grid_fit([math.radians(c) for c in cs], math.pi / 3)
    _, dev36 = grid_fit([math.radians(c) for c in cs], math.pi / 5)
    out[str(step)] = {
        'n_clusters': len(g),
        'centers_deg': [round(c, 3) for c in cs],
        'gaps_deg': gaps,
        'occupancy': sorted((len(m) for _, m in g), reverse=True),
        'amp_min': float(amps.min()), 'amp_max': float(amps.max()),
        'dev_from_60deg_grid_deg': math.degrees(dev60),
        'dev_from_36deg_grid_deg': math.degrees(dev36),
    }
with open(os.path.join(BASE, 'check_N5_long_geometry_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
