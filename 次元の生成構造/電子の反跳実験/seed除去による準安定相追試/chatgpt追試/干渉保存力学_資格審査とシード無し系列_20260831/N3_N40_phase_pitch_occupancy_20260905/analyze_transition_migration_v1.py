#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""遷移期の振幅保護の検証（N=3, 4、読み出しのみ）。
問い（木原）: 90°格子から終状態格子への移行で、(i) 波1本ごとの振幅は保護されるか、
(ii) 初期の位相位置グループごとの合計エネルギー Σ|z|² は保護されるか。
入力: ../N3_N40_stage123_sweep_20260905/results/hm_N{N}_den_{N}_states_500.npz
（SHA256 入力ゲート付き）。波の同一性は成分添字で厳密に追跡できる。
出力: transition_migration_v1.json、fig_transition_migration_N3_N4.png"""
import hashlib
import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')

ledger = {}
with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

def clusters_idx(z, gap_thresh):
    ph = np.angle(z)
    order = np.argsort(ph)
    ps = ph[order]
    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    cut = np.flatnonzero(gaps > gap_thresh)
    groups = []
    if cut.size == 0:
        return [list(range(len(z)))]
    start = (cut[-1] + 1) % len(ps)
    idx = list(order[start:]) + list(order[:start])
    bounds = sorted(((c - start) % len(ps)) for c in cut)
    prev = 0
    for b in bounds:
        groups.append(idx[prev:b + 1])
        prev = b + 1
    return groups

WINDOW = {3: (0, 90), 4: (0, 110)}
out = {}
fig, axs = plt.subplots(2, 2, figsize=(14, 9))
for col, N in enumerate((3, 4)):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    h = hashlib.sha256(open(os.path.join(PKG, rel), 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    d = np.load(os.path.join(PKG, rel))
    Z = np.asarray(d['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    s0, s1 = WINDOW[N]
    groups0 = clusters_idx(Z[0], 0.2)
    steps = np.arange(s0, s1 + 1)
    amps = np.abs(Z[s0:s1 + 1])  # (steps, M)
    gnorm = np.array([[float(np.sum(np.abs(Z[s][g]) ** 2)) for g in groups0] for s in steps])

    ax = axs[0][col]
    for e in range(M):
        ax.plot(steps, amps[:, e], linewidth=1.0)
    ax.set_title(f'N={N}: per-wave |z_e| through transition (M={M})')
    ax.grid(alpha=.3); ax.set_ylabel('|z_e|')
    ax = axs[1][col]
    for gi, g in enumerate(groups0):
        ax.plot(steps, gnorm[:, gi], linewidth=1.2,
                label=f'group{gi} (n={len(g)}, init={gnorm[0][gi]:.4f})')
    ax.set_title(f'N={N}: initial-cluster group energies sum|z|^2')
    ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_xlabel('step'); ax.set_ylabel('sum |z|^2')

    per_wave_delta = float(np.max(np.abs(amps[-1] - amps[0])))
    group_delta = [float(gnorm[-1][gi] - gnorm[0][gi]) for gi in range(len(groups0))]
    out[str(N)] = {
        'groups0_sizes': [len(g) for g in groups0],
        'group_energy_initial': [float(x) for x in gnorm[0]],
        'group_energy_final_window': [float(x) for x in gnorm[-1]],
        'group_energy_delta': group_delta,
        'max_group_energy_delta_abs': float(max(abs(x) for x in group_delta)),
        'per_wave_amp_initial': [float(a) for a in amps[0]],
        'per_wave_amp_final_window': [float(a) for a in amps[-1]],
        'max_per_wave_amp_change': per_wave_delta,
    }
with open(os.path.join(BASE, 'transition_migration_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
fig.suptitle('Amplitude protection through the transition (readout only)', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_transition_migration_N3_N4.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
