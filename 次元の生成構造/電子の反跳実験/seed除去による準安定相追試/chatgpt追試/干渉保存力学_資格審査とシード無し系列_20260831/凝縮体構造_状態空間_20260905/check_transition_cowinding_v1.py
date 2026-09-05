#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""遷移期の共巻き構造の検定（読み出しのみ）。
問い: 凝縮体=零閉塞＋共巻きの定義（実験6）の下で、遷移（旧凝縮体の死→新凝縮体の形成）
の間に、系は (a) 全体が一体のまま巻数レートを連続変化させるのか、
(b) 旧レートの部分と新レートの部分が共存する（相共存型）のか。
方法: N=3, 4 の遷移窓で、波1本ごとの瞬時位相速度
  v_e(s) = arg(z_e(s+1)/z_e(s)) / 2π  [回転/step]
を全ステップで算出し、ステップごとの波間の広がり（max−min）と、
旧レート（固有値予測 w₀·Δτ/2π）・新レート帯への分岐の有無を機械判定する。
出力: check_transition_cowinding_v1.json"""
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

WINDOW = {3: (20, 70), 4: (25, 80)}
RAMP_RATE = {3: math.sqrt(2) / 3.0, 4: 2.0 * math.sqrt(2) / 4.0}  # w0·Δτ/2π = |w0|/den

out = {}
for N in (3, 4):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    gate(rel)
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    s0, s1 = WINDOW[N]
    rows = []
    for s in range(s0, s1):
        v = np.angle(Z[s + 1] / Z[s]) / (2 * math.pi)  # 波ごと、mod 1 の円環上の量
        # 円環スプレッド: mod 1 でソートし最大ギャップの補（±0.5 巻き付きの誤検出を防ぐ）
        u = np.sort(v % 1.0)
        gaps = np.diff(np.concatenate([u, [u[0] + 1.0]]))
        spread = float(1.0 - gaps.max())
        rows.append({'step': s,
                     'v_spread_circular': spread,
                     'v_per_wave': [round(float(x), 6) for x in v]})
    spread = np.array([r['v_spread_circular'] for r in rows])
    # 分岐判定: 遷移窓内で spread の最大と、その step
    imax = int(np.argmax(spread))
    out[str(N)] = {
        'window': [s0, s1], 'M': M,
        'ramp_rate_predicted': RAMP_RATE[N],
        'max_spread': float(spread.max()),
        'max_spread_step': rows[imax]['step'],
        'spread_before_transition_mean': float(spread[:5].mean()),
        'spread_after_transition_mean': float(spread[-5:].mean()),
        'per_step': rows,
    }
    print(f"N={N}: ramp予測 {RAMP_RATE[N]:.4f} 回転/step | spread: 遷移前 "
          f"{out[str(N)]['spread_before_transition_mean']:.2e} → 最大 {spread.max():.4f} "
          f"(step {rows[imax]['step']}) → 遷移後 {out[str(N)]['spread_after_transition_mean']:.4f}")
    # 遷移最大時の波別レートの分布（分岐の直接読み）
    print('  step', rows[imax]['step'], 'の波別レート:', rows[imax]['v_per_wave'])

with open(os.path.join(BASE, 'check_transition_cowinding_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print('ALL DONE')
