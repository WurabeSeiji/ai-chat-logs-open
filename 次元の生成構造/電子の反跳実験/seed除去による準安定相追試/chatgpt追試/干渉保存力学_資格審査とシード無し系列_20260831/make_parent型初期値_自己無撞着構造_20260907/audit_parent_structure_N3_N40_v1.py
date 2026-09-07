#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""親構造監査（指示書 §3、全N・分類のみ・親は変更しない）。
式は analyze_makeparent_parent_structure_20260907.py（ChatGPT 版 N3..7）と同一で、
N=3..40 へ拡張し振幅CVを追加。出力: full_N3_N40_sweep/parent_structure_N3_N40.csv
振幅種数の許容誤差: 12桁丸め（np.round(r,12) の一意数）。"""
import csv
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(BASE, 'parents_actual_N3_N40')
OUT = os.path.join(BASE, 'full_N3_N40_sweep')
os.makedirs(OUT, exist_ok=True)


def adjacency(N):
    ea, eb = np.triu_indices(N, k=1)
    M = len(ea)
    A = np.zeros((M, M))
    for e in range(M):
        sh = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        sh[e] = False
        A[e, sh] = 1
    return A


rows = []
for N in range(3, 41):
    z = np.asarray(np.load(os.path.join(P, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'], complex)
    th = np.angle(z); r = np.abs(z)
    A = adjacency(N)
    W = A * np.sin(th[None, :] - th[:, None]) ** 2
    lam = float(np.vdot(r, W @ r).real / np.vdot(r, r).real)
    res = float(np.linalg.norm(W @ r - lam * r) / np.linalg.norm(r))
    phi = th[0]
    d = ((th - phi + np.pi / 4) % (np.pi / 2)) - np.pi / 4
    q = np.rint((th - phi - d) / (np.pi / 2)).astype(int) % 4
    p = [float(np.sum(r[q == k] ** 2)) for k in range(4)]
    rows.append([N, len(z), abs(z @ z), abs(np.sum(z)),
                 float(np.max(np.abs(d))),
                 len(np.unique(np.round(r, 12))), float(np.std(r) / np.mean(r)),
                 lam, res, *p, p[0] + p[2], p[1] + p[3]])
with open(os.path.join(OUT, 'parent_structure_N3_N40.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'M', 'abs_sum_z2', 'abs_sum_z', 'max_Z4_phase_residual_rad',
                'amplitude_species_1e-12', 'amplitude_CV',
                'W_rayleigh_sigma', 'W_rayleigh_rel_residual',
                'P0', 'P1', 'P2', 'P3', 'P_even', 'P_odd'])
    w.writerows(rows)
for r_ in rows[:6] + rows[-2:]:
    print(f"N={r_[0]:2d} |Σz²|={r_[2]:.1e} |Σz|={r_[3]:.3f} Z4残差={r_[4]:.2e}rad 振幅種={r_[5]:3d} CV={r_[6]:.3f} σ={r_[7]:.4f} 残差={r_[8]:.1e} P_even={r_[13]:.4f} P_odd={r_[14]:.4f}")
print('AUDIT DONE')
