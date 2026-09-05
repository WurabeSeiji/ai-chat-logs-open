#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""軸結晶化の全系列・長時間分析（N=3..40、den=N、10000歩、読み出しのみ）。
検定する予言（木原 2026-09-05）: 終状態は「3軸・60°・軸あたり M/3 本の等分配」に
結晶化する。3|M は N≡0,1 (mod 3) でのみ成立し、N≡2 (mod 3)（N=5,8,11,…）は
3軸等分配に入れないため別構造（N=5 の実測では2軸型）へ向かうはず。
測定（snapshot: step 0,500,1000,...,10000）:
 (1) 軸数（軸角 = arg(z²)/2 のクラスタ）、60°軸格子への最大偏差、軸別波数
 (2) 3軸等分配フラグ: 軸数=3 かつ 全軸の波数=M/3 かつ 60°偏差<0.5°
 (3) 等振幅相対偏差、回転数 x(s)（末尾フレーム）と 1/1 距離
出力: axis_crystallization_table_v1.csv（Nごと1行・最終snapshot）、
      axis_crystallization_v1.json（全snapshot軌跡）、fig_axis_crystallization_N3_N40.png"""
import csv
import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SNAPS = list(range(0, 10001, 500))

def clusters(ph, gap_thresh):
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

def grid_dev(centers, pitch):
    u = np.array([c / pitch for c in centers])
    mean = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    return float(max(abs((x - mean + 0.5) % 1.0 - 0.5) for x in u) * pitch)

rows = []
traj = {}
for N in range(3, 41):
    d = np.load(os.path.join(BASE, 'results', f'hm_N{N}_den_{N}_states_10000.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 10000
    Z = np.asarray(d['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    tr = []
    for s in SNAPS:
        z = Z[s]
        g = clusters(np.angle(z ** 2), math.pi / (3 * N))
        counts = sorted((len(m) for _, m in g), reverse=True)
        dev60 = math.degrees(grid_dev([c for c, _ in g], 2 * math.pi / 3)) / 2.0
        amps = np.abs(z); target = float(np.linalg.norm(z)) / math.sqrt(M)
        eq = float(np.max(np.abs(amps - target)) / target)
        three_eq = bool(len(g) == 3 and M % 3 == 0 and all(c == M // 3 for c in counts) and dev60 < 0.5)
        tr.append({'step': s, 'n_axes': len(g), 'counts': counts[:8], 'dev60_deg': dev60,
                   'eqdev': eq, 'three_axis_equal': three_eq})
    # 末尾の回転数
    s = 10000 - N
    ip = np.vdot(Z[s], Z[s + N])
    x = (float(np.angle(ip)) / (2 * math.pi)) % 1.0
    dist1 = abs((x + 0.5) % 1.0 - 0.5)
    last = tr[-1]
    rows.append([N, M, M % 3, last['n_axes'], '|'.join(map(str, last['counts'])),
                 round(last['dev60_deg'], 3), f"{last['eqdev']:.3e}",
                 last['three_axis_equal'], f'{dist1:.3e}'])
    traj[str(N)] = tr
    print('done N', N, flush=True)

HEADER = ['N', 'M', 'M_mod_3', 'n_axes_10000', 'axis_counts_10000', 'dev60_deg_10000',
          'eqdev_10000', 'three_axis_equal_10000', 'dist_to_1_last']
with open(os.path.join(BASE, 'axis_crystallization_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(HEADER); w.writerows(rows)
with open(os.path.join(BASE, 'axis_crystallization_v1.json'), 'w') as f:
    json.dump({'snapshots': SNAPS, 'table_header': HEADER, 'table': rows, 'trajectories': traj},
              f, indent=2, default=str)

fig, axs = plt.subplots(3, 1, figsize=(12, 14), sharex=False)
colors = {0: 'tab:green', 1: 'tab:red'}
ax = axs[0]
for N in range(3, 41):
    M3 = (N * (N - 1) // 2) % 3
    ax.plot(SNAPS, [t['n_axes'] for t in traj[str(N)]], linewidth=0.8,
            color=colors[M3], alpha=0.6)
ax.axhline(3, color='gray', linestyle='--', linewidth=0.8)
ax.set_ylabel('n_axes'); ax.set_title('Axis count vs step (green: 3|M, red: M=1 mod 3)')
ax.grid(alpha=.3)
ax = axs[1]
for N in range(3, 41):
    M3 = (N * (N - 1) // 2) % 3
    ax.semilogy(SNAPS, [max(t['dev60_deg'], 1e-13) for t in traj[str(N)]], linewidth=0.8,
                color=colors[M3], alpha=0.6)
ax.set_ylabel('dev from 60-deg axis grid [deg]'); ax.grid(alpha=.3)
ax = axs[2]
Ns = list(range(3, 41))
ax.bar(Ns, [1 if r[7] else 0 for r in rows],
       color=[colors[(n * (n - 1) // 2) % 3] for n in Ns])
ax.set_ylabel('3-axis equal @10000'); ax.set_xlabel('N'); ax.grid(alpha=.3)
fig.suptitle('Axis crystallization, N=3..40, dt=2pi/N, 10000 steps (readout only)', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_axis_crystallization_N3_N40.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
