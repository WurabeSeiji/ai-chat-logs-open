#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""検討10-5 の検定: 結晶化域での2時計レート比の系普遍性（読み出しのみ）。
仮説: t が単一の導出時間なら、質量時計（ロック接近の減衰レート）と電荷時計
（巻数レート）の比 = 「結晶化1桁あたりの巻数」は N によらぬ定数のはず。
方法:
 (1) 各 N（den=N、結晶化した指数ロック系）で、飽和後フレーム s の回転数
     x(s)=arg⟨Z[s],Z[s+N]⟩/2π の 1/1 からの円距離 d(s) を計算し、
     log10 d の直線区間（1e-12 < d < 0.3、3点以上・1.5桁以上）をフィット
     → steps/decade_lock（質量時計レート）。
 (2) 電荷時計レート: q = w_dom/N（素）と q_eff = |q − round(q)|（エイリアス、
     ストロボ領域で機能する読み・実験6-4）。w_dom は回転数ロック調査の
     保存出力（SHA台帳済みフォルダ）から引用。
 (3) 比 r = q_eff × steps/decade_lock（結晶化1桁あたりのエイリアス巻数）を
     全成立系で並べ、CV を steps/decade_lock 単独の CV と比較。
     さらに増幅域の定数（実験6-4: ストロボ領域 0.399 回転/桁）とも比較する。
出力: check_clock_ratio_v1.json"""
import csv
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
RL = os.path.join(BASE, '..', 'N3_N40_rotation_lock_analysis_20260905')

ledger = {}
with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

rl = json.load(open(os.path.join(RL, 'rotation_lock_analysis_v1.json')))
wdom = {r[0]: abs(float(r[8])) for r in rl['table']}
onset = {r[0]: int(r[2]) for r in rl['table']}

def circdist(a, b):
    return abs((a - b + 0.5) % 1.0 - 0.5)

rows = []
for N in range(4, 41):
    if N == 5:
        continue  # 冪則（周辺的）系は指数フィットの対象外
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    assert ledger[rel] == sha256(os.path.join(PKG, rel)), f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    sc = onset[N]
    if sc < 0:
        continue
    pts = []
    for s in range(0, 501 - N, N):
        if s < sc + 10:
            continue
        x = (float(np.angle(np.vdot(Z[s], Z[s + N]))) / (2 * math.pi)) % 1.0
        d = circdist(x, 0.0)
        if 1e-12 < d < 0.3:
            pts.append((s, d))
    if len(pts) < 3:
        continue
    xs = np.array([p[0] for p in pts], float)
    ys = np.log10([p[1] for p in pts])
    if ys.max() - ys.min() < 1.5:
        continue
    slope = float(np.polyfit(xs, ys, 1)[0])
    if slope >= 0:
        continue
    spd_lock = -1.0 / slope
    q = wdom[N] / N
    q_eff = abs(q - round(q))
    rows.append({'N': N, 'onset': sc, 'n_pts': len(pts),
                 'decades_span': round(float(ys.max() - ys.min()), 2),
                 'lock_steps_per_decade': round(spd_lock, 3),
                 'q': round(q, 6), 'q_eff': round(q_eff, 6),
                 'turns_alias_per_decade': round(q_eff * spd_lock, 4),
                 'turns_raw_per_decade': round(q * spd_lock, 4)})

def stats(key):
    a = np.array([r[key] for r in rows])
    return {'n': len(a), 'mean': float(a.mean()), 'std': float(a.std()),
            'cv': float(a.std() / a.mean()), 'min': float(a.min()), 'max': float(a.max())}

out = {'rows': rows,
       'stats_lock_steps_per_decade': stats('lock_steps_per_decade'),
       'stats_turns_alias_per_decade': stats('turns_alias_per_decade'),
       'stats_turns_raw_per_decade': stats('turns_raw_per_decade'),
       'ramp_reference_alias_turns_per_decade': 0.399}
with open(os.path.join(BASE, 'check_clock_ratio_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(f"{'N':>3} {'onset':>5} {'pts':>4} {'span':>5} {'spd_lock':>9} {'q_eff':>8} {'alias_t/dec':>11} {'raw_t/dec':>10}")
for r in rows:
    print(f"{r['N']:>3} {r['onset']:>5} {r['n_pts']:>4} {r['decades_span']:>5} "
          f"{r['lock_steps_per_decade']:>9} {r['q_eff']:>8} {r['turns_alias_per_decade']:>11} {r['turns_raw_per_decade']:>10}")
for k in ('stats_lock_steps_per_decade', 'stats_turns_alias_per_decade', 'stats_turns_raw_per_decade'):
    s = out[k]
    print(f"{k}: n={s['n']} mean={s['mean']:.4f} CV={s['cv']:.3f} [{s['min']:.3f},{s['max']:.3f}]")
print('ALL DONE')
