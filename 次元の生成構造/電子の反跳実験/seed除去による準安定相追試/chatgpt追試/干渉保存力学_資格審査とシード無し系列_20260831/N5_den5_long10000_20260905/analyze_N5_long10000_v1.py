#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=5・Δτ=2π/5・10000歩走行の読み出し分析（判定は走行後のここでのみ行う）。
入力: results/hm_N5_den_5_states_10000.npz（本フォルダの走行、対照ゲート済み）
測定（フレーム s=0,5,...,9995）:
 (1) 回転数 x(s)=arg⟨Z[s],Z[s+5]⟩/2π、1/1・1/2 からの円距離、剛体残差 r(s)
 (2) 等振幅相対偏差 eqdev(s)、H⊥/H（timeseries から）
 (3) 最終フレームの x の最良有理近似（q≤48）、dist_to_1 の全区間最小とその step
出力: n5_long10000_frames_v1.csv（全フレーム）、analyze_N5_long10000_v1.json（要約）、
      fig_N5_long10000_convergence.png（収束診断 2面図）"""
import csv
import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(BASE, 'results', 'hm_N5_den_5_states_10000.npz'))
assert int(d['denominator']) == 5 and int(d['steps']) == 10000
Z = np.asarray(d['Z'], dtype=np.complex128)
N = 5; M = 10

hperp = np.full(10001, np.nan)
with open(os.path.join(BASE, 'results', 'timeseries_64bit_N5_den5_steps10000.csv')) as f:
    for row in csv.DictReader(f):
        hperp[int(row['step'])] = float(row['Hperp_frac'])

def circdist(a, b):
    return abs((a - b + 0.5) % 1.0 - 0.5)

def best_rational(x, qmax):
    x = x % 1.0
    best = (None, None, 2.0)
    for q in range(1, qmax + 1):
        p = round(x * q)
        res = abs(x - p / q)
        if res < best[2] - 1e-18:
            best = (p % q if q > 1 else p, q, res)
    return best

frames = list(range(0, 10001 - N, N))
rows = []
for s in frames:
    ip = np.vdot(Z[s], Z[s + N]); ph = float(np.angle(ip))
    nrm = float(np.linalg.norm(Z[s]))
    r = float(np.linalg.norm(Z[s + N] - np.exp(1j * ph) * Z[s])) / nrm
    x = (ph / (2 * math.pi)) % 1.0
    amps = np.abs(Z[s]); target = nrm / math.sqrt(M)
    eq = float(np.max(np.abs(amps - target)) / target)
    rows.append((s, x, circdist(x, 0.0), circdist(x, 0.5), r, eq, float(hperp[s])))

with open(os.path.join(BASE, 'n5_long10000_frames_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['step', 'x_turns', 'dist_to_1', 'dist_to_half', 'rigid_residual', 'eqmod_rel_dev', 'Hperp_frac'])
    w.writerows(rows)

arr = np.array([(r[0], r[2], r[4], r[5]) for r in rows])
late = arr[arr[:, 0] >= 200]  # 飽和後
imin = int(np.argmin(late[:, 1]))
xf = rows[-1][1]
bp, bq, bres = best_rational(xf, 48)
out = {
    'final_frame': {'step': rows[-1][0], 'x_turns': xf, 'dist_to_1': rows[-1][2],
                    'rigid_residual': rows[-1][4], 'eqmod_rel_dev': rows[-1][5], 'Hperp': rows[-1][6]},
    'best_rational_final': {'p': bp, 'q': bq, 'residual_turns': bres},
    'min_dist_to_1_after_200': {'step': int(late[imin, 0]), 'value': float(late[imin, 1])},
    'dist_to_1_range_last_5000': [float(np.min(late[late[:, 0] >= 5000][:, 1])),
                                  float(np.max(late[late[:, 0] >= 5000][:, 1]))],
    'eqdev_range_last_5000': [float(np.min(late[late[:, 0] >= 5000][:, 3])),
                              float(np.max(late[late[:, 0] >= 5000][:, 3]))],
    'Hperp_mean_last_5000': float(np.nanmean(hperp[5000:])),
    'Hperp_std_last_5000': float(np.nanstd(hperp[5000:])),
}

# 減衰形の同定: steps 2000..9995 の log-log 勾配（冪則 t^α の α。指数減衰なら
# log-log では下に凸で直線に乗らない——残差も併記）
fitwin = arr[arr[:, 0] >= 2000]
for name, col in (('dist_to_1', 1), ('rigid_residual', 2), ('eqmod_rel_dev', 3)):
    lx = np.log10(fitwin[:, 0]); ly = np.log10(np.maximum(fitwin[:, col], 1e-300))
    coef = np.polyfit(lx, ly, 1)
    resid = float(np.sqrt(np.mean((np.polyval(coef, lx) - ly) ** 2)))
    out[f'powerlaw_{name}'] = {'exponent': float(coef[0]), 'rms_residual_log10': resid}
with open(os.path.join(BASE, 'analyze_N5_long10000_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))

fig, axs = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
ax = axs[0]
ax.semilogy(arr[:, 0], np.maximum(arr[:, 1], 1e-18), label='dist(x, 1/1) [turns]', linewidth=0.8)
ax.semilogy(arr[:, 0], np.maximum(arr[:, 2], 1e-18), label='rigid residual', linewidth=0.8)
ax.semilogy(arr[:, 0], np.maximum(arr[:, 3], 1e-18), label='equimodularity dev', linewidth=0.8)
ax.grid(alpha=.3); ax.legend(); ax.set_title('N=5, dt=2pi/5, 10000 steps: convergence diagnostics (readout only)')
ax = axs[1]
ax.plot(np.arange(10001), hperp, linewidth=0.6)
ax.grid(alpha=.3); ax.set_xlabel('step'); ax.set_ylabel('Hperp/H')
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_N5_long10000_convergence.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
