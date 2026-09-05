#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""軸構造（方向 = 位相 mod 180°）の検定（読み出しのみ）。
仮説（木原 2026-09-05）: 初期 90° ピッチ = 2軸（2次元）構造、終状態 60° ピッチ =
3軸（3方向）構造への移行ではないか。
方法: 各波の軸角 = arg(z_e²)/2（mod 180°）。z² の位相をクラスタ分割すれば軸が数えられる。
測定: 全 N の step0 / step500（＋N=5 は 10000歩）で、軸の本数・軸間隔の 90°/60° 格子
フィット偏差・軸ごとの波数。3軸等分配には 3|M が必要（M=N(N−1)/2、N≡2 mod 3 で不成立）
という算術も表に併記。
出力: check_axis_structure_v1.json"""
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
LONG = os.path.join(BASE, '..', 'N5_den5_long10000_20260905')

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

def clusters(ph, gap_thresh):
    """位相列 ph (rad, mod 2π想定) をクラスタ分割 → [(中心rad, メンバー添字), ...]"""
    order = np.argsort(ph)
    ps = ph[order]
    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    cut = np.flatnonzero(gaps > gap_thresh)
    if cut.size == 0:
        return [(float(np.angle(np.mean(np.exp(1j * ph)))), list(range(len(ph))))]
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
    return float(dev.max())

def axes_of(z, thresh_doubled_rad):
    """軸クラスタ: 2θ 空間でクラスタ → 軸角(deg, mod 180)・軸ごとの波数"""
    ph2 = np.angle(z ** 2)
    g = clusters(ph2, thresh_doubled_rad)
    axes = sorted((math.degrees(c) / 2.0 % 180.0, len(m)) for c, m in g)
    return axes

out = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(PKG, rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    ax0 = axes_of(Z[0], 0.4)
    axf = axes_of(Z[500], math.pi / (2 * N))
    # 軸間隔フィット: 2θ 空間ピッチ 180°(=軸90°) / 120°(=軸60°)
    dev0_90 = grid_fit([math.radians(2 * a) for a, _ in ax0], math.pi)
    devf_60 = grid_fit([math.radians(2 * a) for a, _ in axf], 2 * math.pi / 3)
    out[str(N)] = {
        'M': M, 'M_mod_3': M % 3,
        'step0_axes': [(round(a, 2), c) for a, c in ax0],
        'step0_n_axes': len(ax0),
        'step0_dev_from_90deg_axes_deg': math.degrees(dev0_90) / 2.0,
        'final_axes': [(round(a, 2), c) for a, c in axf],
        'final_n_axes': len(axf),
        'final_dev_from_60deg_axes_deg': math.degrees(devf_60) / 2.0,
    }
# N=5 長時間走行の終状態
gate(LONG, 'results/hm_N5_den_5_states_10000.npz')
Z5 = np.asarray(np.load(os.path.join(LONG, 'results/hm_N5_den_5_states_10000.npz'))['Z'],
                dtype=np.complex128)
ax5 = axes_of(Z5[10000], 0.15)
out['N5_step10000'] = {'axes': [(round(a, 2), c) for a, c in ax5], 'n_axes': len(ax5)}

with open(os.path.join(BASE, 'check_axis_structure_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
for N in (3, 4, 5, 6, 7, 8):
    v = out[str(N)]
    print(f"N={N} M={v['M']} (M mod 3 = {v['M_mod_3']}): step0 axes {v['step0_axes']} "
          f"dev90={v['step0_dev_from_90deg_axes_deg']:.2e}° | final axes {v['final_axes']} "
          f"dev60={v['final_dev_from_60deg_axes_deg']:.2f}°")
print('N=5@10000 axes:', out['N5_step10000']['axes'])
print('step0 n_axes all N:', sorted({out[str(N)]['step0_n_axes'] for N in range(3, 41)}))
print('ALL DONE')
