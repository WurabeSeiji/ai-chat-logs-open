#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=5 異常系の収束レート測定（読み出しのみ・延長実験の規模設計用）。
入力: ../N3_N40_stage123_sweep_20260905/results/hm_N5_den_5_states_500.npz（SHA ゲート付き）
測定:
 (1) 各フレーム s=0,5,...,495 の回転数 x(s) と、1/1 および 1/2 からの円距離
 (2) 剛体残差 r(s)、等振幅相対偏差 eqdev(s)
 (3) 後半窓での log10 減衰勾配（steps/decade）と、倍精度床 1e-15 到達に要する歩数の外挿
出力: check_N5_convergence_v1.json"""
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
REL = 'results/hm_N5_den_5_states_500.npz'

ledger = {}
with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]
h = hashlib.sha256(open(os.path.join(PKG, REL), 'rb').read()).hexdigest()
assert ledger[REL] == h, 'INPUT GATE FAIL'

d = np.load(os.path.join(PKG, REL))
assert int(d['denominator']) == 5 and int(d['steps']) == 500
Z = np.asarray(d['Z'], dtype=np.complex128)
N = 5; M = 10

def circdist(a, b):
    return abs((a - b + 0.5) % 1.0 - 0.5)

rows = []
for s in range(0, 496, 5):
    ip = np.vdot(Z[s], Z[s + N]); ph = float(np.angle(ip))
    nrm = float(np.linalg.norm(Z[s]))
    r = float(np.linalg.norm(Z[s + N] - np.exp(1j * ph) * Z[s])) / nrm
    x = (ph / (2 * math.pi)) % 1.0
    amps = np.abs(Z[s]); target = nrm / math.sqrt(M)
    eq = float(np.max(np.abs(amps - target)) / target)
    rows.append((s, x, circdist(x, 0.0), circdist(x, 0.5), r, eq))

def slope_steps_per_decade(vals, steps):
    """log10(vals) の最小二乗勾配から steps/decade（正=減衰）"""
    y = np.log10(np.maximum(vals, 1e-300)); x = np.asarray(steps, float)
    a = np.polyfit(x, y, 1)[0]
    return float(-1.0 / a) if a < 0 else float('inf')

late = [row for row in rows if row[0] >= 250]
steps = [row[0] for row in late]
d1_vals = [row[2] for row in late]
r_vals = [row[4] for row in late]
eq_vals = [row[5] for row in late]
out = {
    'frames_sampled': [{'step': s, 'x_turns': x, 'dist_to_1': d1, 'dist_to_half': dh,
                        'rigid_residual': r, 'eqmod_rel_dev': eq}
                       for (s, x, d1, dh, r, eq) in rows if s % 25 == 0 or s >= 450],
    'late_window_start': 250,
    'steps_per_decade': {
        'dist_to_1': slope_steps_per_decade(d1_vals, steps),
        'rigid_residual': slope_steps_per_decade(r_vals, steps),
        'eqmod_rel_dev': slope_steps_per_decade(eq_vals, steps),
    },
    'current_at_495': {'dist_to_1': rows[-1][2], 'rigid_residual': rows[-1][4], 'eqmod_rel_dev': rows[-1][5]},
}
# 床 1e-15 到達の外挿（各量、レートが有限の場合）
proj = {}
for key, cur in (('dist_to_1', rows[-1][2]), ('rigid_residual', rows[-1][4]), ('eqmod_rel_dev', rows[-1][5])):
    rate = out['steps_per_decade'][key]
    if math.isfinite(rate) and cur > 1e-15:
        proj[key] = int(495 + rate * (math.log10(cur) - (-15.0)))
    else:
        proj[key] = None
out['projected_step_to_1e-15'] = proj
with open(os.path.join(BASE, 'check_N5_convergence_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
