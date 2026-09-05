#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""「時間=巻数」仮説の検定（読み出しのみ・新規走行なし）。

仮説（木原 2026-09-05）: 処理ステップ τ は時間に比例しない。時間の実体は
全体位相の 2π に対する巻数である。
検定可能な帰結: もし時間=巻数なら、増幅レートを「ステップあたり」で測ると
N・den に依存してばらつくが、「巻数あたり」で測ると普遍定数に collapse するはず。

方法:
 (1) 各走行（N=3..40 × den∈{N-2,N-1,N,N+1,N+2,124}、正本 timeseries から）で
     log10 H⊥/H の直線区間（f∈[max(1e-28,100·f(1)), 1e-6] かつ3桁以上）を最小二乗
     フィットし steps/decade を得る。
 (2) 増幅期は親が固有状態（占有1.0、rotation_lock 調査で確立）なので、
     巻数レート = Δτ·|w_dom|/2π = |w_dom|/den [回転/step]。w_dom は初期生成子
     H=i·K の固有値（エンジンは正本の逐語コピー、入力は SHA ゲート付き npz）。
 (3) turns/decade = steps/decade × |w_dom|/den を全走行で算出し、
     steps/decade と turns/decade の変動係数（CV=std/mean）を比較する。
     時間=巻数なら CV(turns/decade) ≪ CV(steps/decade)。
 (4) 検証: N=3 den=3 の巻数レートを状態から直接積算（±180°補正の逐次巻数）し、
     増幅期 √2/3 回転/step と一致するか、ロック後のレートも実測する。
出力: winding_rates_table_v1.csv、analyze_winding_rates_v1.json、
      fig_winding_collapse.png"""
import csv
import hashlib
import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
RES = os.path.join(PKG, 'results')

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

# ---- エンジン（run_N3_N40_stage123_v1.py 14-22行の逐語コピー） ----
def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A
def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
# ---- コピーここまで ----

# w_dom（占有最大モードの固有値）を全 N で
wdom = {}
Z0_cache = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    assert ledger[rel] == sha256(os.path.join(PKG, rel)), f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(os.path.join(PKG, rel))['Z'], dtype=np.complex128)
    Z0_cache[N] = Z
    A = adjacency(N)
    H = H_of(np.exp(1j * np.angle(Z[0])), A); H = (1j * np.imag(H)).astype(np.complex128, copy=False)
    w, V = np.linalg.eigh(H)
    occ = np.abs(V.conj().T @ Z[0]) ** 2
    wdom[N] = abs(float(w[int(np.argmax(occ))]))

# timeseries 全読み
series = {}
with open(os.path.join(RES, 'timeseries_64bit_with124_N3_N40.csv')) as f:
    for row in csv.DictReader(f):
        key = (int(row['N']), int(row['denominator']))
        series.setdefault(key, np.full(501, np.nan))[int(row['step'])] = float(row['Hperp_frac'])

rows = []
for (N, den), fvals in sorted(series.items()):
    f1 = fvals[1]
    lo = max(1e-28, 100.0 * f1)
    hi = 1e-6
    idx = [t for t in range(1, 501) if lo < fvals[t] < hi]
    if len(idx) < 10:
        continue
    span = math.log10(fvals[idx[-1]] / fvals[idx[0]])
    if span < 3.0:
        continue
    coef = np.polyfit(idx, np.log10([fvals[t] for t in idx]), 1)
    if coef[0] <= 0:
        continue
    steps_per_decade = 1.0 / float(coef[0])
    winding_per_step = wdom[N] / den
    turns_per_decade = steps_per_decade * winding_per_step
    rows.append([N, den, round(steps_per_decade, 4), round(wdom[N], 10),
                 round(winding_per_step, 8), round(turns_per_decade, 4),
                 len(idx), round(span, 2)])

arr_s = np.array([r[2] for r in rows]); arr_t = np.array([r[5] for r in rows])
cv = {'n_runs_fit': len(rows),
      'steps_per_decade': {'mean': float(arr_s.mean()), 'std': float(arr_s.std()),
                           'cv': float(arr_s.std() / arr_s.mean()),
                           'min': float(arr_s.min()), 'max': float(arr_s.max())},
      'turns_per_decade': {'mean': float(arr_t.mean()), 'std': float(arr_t.std()),
                           'cv': float(arr_t.std() / arr_t.mean()),
                           'min': float(arr_t.min()), 'max': float(arr_t.max())}}

# (4) N=3 den=3 の直接巻数積算（検証）
Z = Z0_cache[3]
adv = []
for s in range(500):
    d = np.angle(np.vdot(Z[s], Z[s + 1]))
    adv.append(float(d))
adv = np.array(adv)  # rad/step、|adv|<π なので巻数は逐次和で一意
Wcum = np.cumsum(adv) / (2 * math.pi)
ramp_rate = float(np.mean(adv[:30])) / (2 * math.pi)
lock_rate = float(np.mean(adv[60:500])) / (2 * math.pi)
direct = {'ramp_turns_per_step_measured': ramp_rate,
          'ramp_turns_per_step_predicted_sqrt2_over_3': math.sqrt(2) / 3.0,
          'ramp_dev': abs(ramp_rate - math.sqrt(2) / 3.0),
          'lock_turns_per_step_measured': lock_rate,
          'total_winding_turns_500steps': float(Wcum[-1])}

out = {'cv_comparison': cv, 'direct_winding_N3_den3': direct,
       'wdom': {str(k): v for k, v in wdom.items()}}
with open(os.path.join(BASE, 'analyze_winding_rates_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
with open(os.path.join(BASE, 'winding_rates_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'den', 'steps_per_decade', 'w_dom', 'winding_per_step',
                'turns_per_decade', 'n_fit_points', 'decades_span'])
    w.writerows(rows)
print(json.dumps(cv, indent=2))
print(json.dumps(direct, indent=2))

fig, axs = plt.subplots(2, 1, figsize=(12, 10))
Ns = [r[0] for r in rows]; dens = [r[1] for r in rows]
c124 = ['tab:orange' if d == 124 else 'tab:blue' for d in dens]
axs[0].scatter(Ns, [r[2] for r in rows], c=c124, s=18)
axs[0].set_ylabel('steps / decade'); axs[0].set_yscale('log')
axs[0].set_title('Growth rate per STEP (blue: den near N, orange: den=124)')
axs[0].grid(alpha=.3)
axs[1].scatter(Ns, [r[5] for r in rows], c=c124, s=18)
axs[1].set_ylabel('winding turns / decade'); axs[1].set_yscale('log')
axs[1].set_title('Growth rate per WINDING TURN (hypothesis: collapse if time = winding)')
axs[1].set_xlabel('N'); axs[1].grid(alpha=.3)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_winding_collapse.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
