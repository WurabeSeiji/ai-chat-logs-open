#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回転数・振幅・飽和値の代数値ロックの数値検証（読み出しのみ）。
入力: data/hm_N3_den_3_states_500.npz（正本コピー）、
      ../N3_N40_stage123_sweep_20260905/results/timeseries_64bit_with124_N3_N40.csv（H⊥/H 引用）
検証項目（すべて機械比較、tol は倍精度床の視認用）:
 (1) 増幅期（frame 0..30）の時計1回転あたり剛体回転角 = (√2−1)·360° か
 (2) 飽和後（frame 60..114）の回転角 = ±180° か
 (3) step0 振幅 = {1/2, 1/2, 1/√2} か、飽和後振幅 = 1/√3 か
 (4) 飽和後 H⊥/H = 1/6 か（step 51..120）
出力: check_rotation_lock_v1.json"""
import csv
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(BASE, 'data', 'hm_N3_den_3_states_500.npz'))
assert int(d['denominator']) == 3 and int(d['steps']) == 500
Z = np.asarray(d['Z'], dtype=np.complex128)

def mean_phase(z):
    return np.angle(np.sum(np.exp(1j * np.angle(z))))

def rot3(s):
    dphi = math.degrees(mean_phase(Z[s + 3]) - mean_phase(Z[s]))
    return (dphi + 180.0) % 360.0 - 180.0

SQRT2_LOCK = (math.sqrt(2.0) - 1.0) * 360.0
early = [rot3(s) for s in range(0, 31, 3)]
late = [abs(rot3(s)) for s in range(60, 115, 3)]

amps0 = np.sort(np.abs(Z[0]))
amps117 = np.abs(Z[117])
target0 = np.array([0.5, 0.5, 1.0 / math.sqrt(2.0)])

hp = []
csv_path = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'results',
                        'timeseries_64bit_with124_N3_N40.csv')
with open(csv_path) as f:
    for row in csv.DictReader(f):
        if row['N'] == '3' and row['denominator'] == '3' and 51 <= int(row['step']) <= 120:
            hp.append(float(row['Hperp_frac']))

out = {
    'lock_sqrt2_deg': SQRT2_LOCK,
    'early_rot3_mean_deg': float(np.mean(early)),
    'early_rot3_max_absdev_from_lock_deg': float(max(abs(r - SQRT2_LOCK) for r in early)),
    'late_rot3_max_absdev_from_180_deg': float(max(abs(r - 180.0) for r in late)),
    'amps_step0': [float(a) for a in amps0],
    'amps_step0_max_absdev_from_half_half_invsqrt2': float(np.max(np.abs(amps0 - target0))),
    'amps_step117_max_absdev_from_invsqrt3': float(np.max(np.abs(amps117 - 1.0 / math.sqrt(3.0)))),
    'Hperp_sat_max_absdev_from_1over6': float(max(abs(h - 1.0 / 6.0) for h in hp)),
    'n_sat_samples': len(hp),
}
with open(os.path.join(BASE, 'check_rotation_lock_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
