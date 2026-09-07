#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=100 親の構造監査（audit_parent_structure_N3_N40_v1.py と同一式・読み出しのみ）。
出力: parents_N100/audit_parent_N100_v1.json"""
import json
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(BASE, 'parents_N100', 'parent_static_N00100_makeparent_20260905.npz'))
z = np.asarray(d['Z0'], complex); N = 100; M = len(z)
th = np.angle(z); r = np.abs(z)
ea, eb = np.triu_indices(N, k=1)
A = np.zeros((M, M))
for e in range(M):
    sh = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); sh[e] = False; A[e, sh] = 1
W = A * np.sin(th[None, :] - th[:, None]) ** 2
lam = float(np.vdot(r, W @ r).real / np.vdot(r, r).real)
res = float(np.linalg.norm(W @ r - lam * r) / np.linalg.norm(r))
phi = th[0]; dv = ((th - phi + np.pi / 4) % (np.pi / 2)) - np.pi / 4
q = np.rint((th - phi - dv) / (np.pi / 2)).astype(int) % 4
p = [float(np.sum(r[q == k] ** 2)) for k in range(4)]
out = dict(N=N, M=M, parent_residual=float(d['residual']),
           norm2=float(np.vdot(z, z).real), abs_sum_z2=float(abs(z @ z)), abs_sum_z=float(abs(z.sum())),
           max_Z4_phase_residual_rad=float(np.max(np.abs(dv))),
           four_position_occupancy=[int(np.sum(q == k)) for k in range(4)],
           amplitude_species_1e12=int(len(set(np.round(r, 12)))), amplitude_CV=float(r.std() / r.mean()),
           amp_min=float(r.min()), amp_max=float(r.max()), one_over_sqrtM=float(1 / np.sqrt(M)),
           W_rayleigh_sigma=lam, W_rayleigh_rel_residual=res,
           P0123=p, P_even=p[0] + p[2], P_odd=p[1] + p[3])
with open(os.path.join(BASE, 'parents_N100', 'audit_parent_N100_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
