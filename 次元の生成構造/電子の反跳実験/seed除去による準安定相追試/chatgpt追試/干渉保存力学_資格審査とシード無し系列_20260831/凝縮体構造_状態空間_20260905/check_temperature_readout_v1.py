#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""温度読出し曲線 T_read の実測（読み出しのみ）。
定義は熱力学読出し論文（DOI 10.5281/zenodo.22240034、note「絶対零度の幾何学」）:
  1/T_read = ΔS/ΔE。
本検定での具体化: S = 振幅エントロピー −Σp ln p（nats、p_e=|z_e|²/Σ|z|²）、
E = H⊥/H（休眠比、正本 timeseries から引用）。軌道に沿って中心差分で
1/T(s) = dS/dE を算出し、増幅期（T=∞ 予想: S凍結・E急増）→ 遷移（有限T）→
終状態（T→0 = 絶対零度の結晶）の熱史を数値化する。N=5 は例外挙動を記録。
出力: check_temperature_readout_v1.json"""
import csv
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')

def gate(rel):
    led = {}
    with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                led[parts[1]] = parts[0]
    h = hashlib.sha256(open(os.path.join(PKG, rel), 'rb').read()).hexdigest()
    assert led[rel] == h, f'INPUT GATE FAIL: {rel}'

hperp = {}
with open(os.path.join(PKG, 'results', 'timeseries_64bit_with124_N3_N40.csv')) as f:
    for row in csv.DictReader(f):
        if row['series'] == 'N':
            hperp.setdefault(int(row['N']), np.full(501, np.nan))[int(row['step'])] = float(row['Hperp_frac'])

out = {}
for N in (3, 4, 5, 6, 7, 10, 20, 40):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    P = np.abs(Z) ** 2
    P = P / P.sum(axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.where(P > 0, np.log(P), 0.0)
    S = -(P * L).sum(axis=1)
    E = hperp[N]
    # 中心差分（stride 2）で 1/T = dS/dE
    invT = np.full(501, np.nan)
    for s in range(2, 499):
        dE = E[s + 2] - E[s - 2]
        dS = S[s + 2] - S[s - 2]
        if abs(dE) > 0:
            invT[s] = dS / dE
    # 特徴点: 増幅期中盤（onset/2）、遷移（S上昇最大step）、終端
    dS1 = np.diff(S)
    st_tr = int(np.argmax(dS1))
    ix = np.flatnonzero(E > 0.05)
    onset = int(ix[0]) if ix.size else -1
    ramp_s = max(4, onset // 2)
    out[str(N)] = {
        'onset': onset, 'transition_step': st_tr,
        'invT_ramp_mid': float(invT[ramp_s]) if not math.isnan(invT[ramp_s]) else None,
        'invT_transition': float(invT[st_tr]) if not math.isnan(invT[st_tr]) else None,
        'invT_final_mean_last50': float(np.nanmean(invT[445:499])),
        'S_total_rise': float(S[-1] - S[0]),
        'E_final': float(E[500]),
        'sample_invT': {str(s): (float(invT[s]) if not math.isnan(invT[s]) else None)
                        for s in (10, 20, 30, 40, st_tr - 3, st_tr, st_tr + 3, 100, 200, 400, 480)
                        if 2 <= s <= 498},
    }
    print(f"N={N}: onset={onset} 遷移step={st_tr} | 1/T: 増幅期中盤={out[str(N)]['invT_ramp_mid']:.3e} "
          f"遷移={out[str(N)]['invT_transition']:.3e} 終端平均={out[str(N)]['invT_final_mean_last50']:.3e}",
          flush=True)

with open(os.path.join(BASE, 'check_temperature_readout_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print('ALL DONE')
