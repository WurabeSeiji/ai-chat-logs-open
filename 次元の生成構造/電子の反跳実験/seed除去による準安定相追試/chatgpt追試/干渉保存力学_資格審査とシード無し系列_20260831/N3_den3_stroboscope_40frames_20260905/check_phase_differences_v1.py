#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ストロボ図の位相差読み取り（90°→60°移行・定常入り時期）の数値検証（読み出しのみ）。
入力:
  data/hm_N3_den_3_states_500.npz — 状態（本体スイープのコピー、SHA256 一致確認済み）
  ../N3_N40_stage123_sweep_20260905/results/timeseries_64bit_with124_N3_N40.csv
    — H⊥/H の正本読み出し値（引用のみ。再実装しない）
検証項目:
 (1) 各 step の3波の円環上の隣接位相差（ソート後の3ギャップ、和=360°）と振幅
 (2) ギャップ組が {60°,60°,240°} に収束する step（閾値 0.5° と 1e-9° の2段判定）
 (3) 3 step（時計1回転）あたりの全体回転角（平均位相の前進）の時間発展
 (4) H⊥/H の対応値（step 0..120、正本 CSV から引用）
出力: check_phase_differences_v1.json（判定はすべて機械適用）。"""
import csv
import json
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(BASE, 'data', 'hm_N3_den_3_states_500.npz'))
assert int(d['denominator']) == 3 and int(d['steps']) == 500
Z = np.asarray(d['Z'], dtype=np.complex128)

def gaps_deg(z):
    ph = np.sort(np.degrees(np.angle(z)))
    g = np.diff(np.concatenate([ph, [ph[0] + 360.0]]))
    return np.sort(g)  # 昇順3ギャップ（和=360）

TARGET = np.array([60.0, 60.0, 240.0])

gap_by_step = {}
amp_by_step = {}
for s in range(0, 121):
    gap_by_step[s] = [float(x) for x in gaps_deg(Z[s])]
    amp_by_step[s] = [float(a) for a in np.sort(np.abs(Z[s]))]

def onset(tol):
    """s=120 まで max|gaps−TARGET|<tol が連続成立し始める最小 s（0..120内）"""
    s0 = None
    for s in reversed(range(0, 121)):
        if max(abs(np.array(gap_by_step[s]) - TARGET)) < tol:
            s0 = s
        else:
            break
    return s0

# 3 step ごとの全体回転角: 平均位相（振幅重み無し円平均）の前進を deg で
def mean_phase(z):
    return np.angle(np.sum(np.exp(1j * np.angle(z))))

rot3_by_frame = {}
for s in range(0, 118, 3):
    dphi = np.degrees(mean_phase(Z[s + 3]) - mean_phase(Z[s]))
    rot3_by_frame[s] = float((dphi + 180.0) % 360.0 - 180.0)

# H⊥/H（正本 CSV から N=3, den=3, step 0..120 を引用）
hperp = {}
csv_path = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'results',
                        'timeseries_64bit_with124_N3_N40.csv')
with open(csv_path) as f:
    for row in csv.DictReader(f):
        if row['N'] == '3' and row['denominator'] == '3':
            s = int(row['step'])
            if s <= 120:
                hperp[s] = float(row['Hperp_frac'])

KEY_STEPS = [0, 1, 2, 3, 6, 9, 12, 24, 30, 36, 42, 45, 46, 47, 48, 49, 51, 60, 90, 117]
out = {
    'target_gaps_deg': [60.0, 60.0, 240.0],
    'onset_gaps60_tol_0p5deg': onset(0.5),
    'onset_gaps60_tol_1e-9deg': onset(1e-9),
    'gaps_deg_key_steps': {str(s): gap_by_step[s] for s in KEY_STEPS},
    'amps_key_steps': {str(s): amp_by_step[s] for s in KEY_STEPS},
    'rotation_per_3steps_deg_frames': {str(s): rot3_by_frame[s] for s in sorted(rot3_by_frame)},
    'Hperp_key_steps': {str(s): hperp[s] for s in KEY_STEPS if s in hperp},
}
with open(os.path.join(BASE, 'check_phase_differences_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
