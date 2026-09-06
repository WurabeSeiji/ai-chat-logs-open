#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A対照再現（シード有り δ=0.01 / 真空 δ=0）と現行シード無し系列のインフレーション図
（読み出しのみ・新規走行なし）。
plot_inflation_N40_stage2init_v1.py（本フォルダに SHA 一致で忠実コピー）の描画様式を verbatim
流用（semilogy・縦軸 H⊥/H・横軸 step・lw=1.4・ylim 1e-34..2・grid・legend lower right）。
変更は入力の差し替えのみ:
 (a) 論文A対照再現 npz の m_f2（物質＝シード有り）/ v_f2（真空＝δ=0）時系列。
     f₂ は「新しい方向の割合」で、現行系列の H⊥/H と同じ定義量。
 (b) 現行シード無し系列 N=12 den=12（10,000歩）の H⊥/H = 1 − (|p·Z|²+|q·Z|²)/‖Z‖²
     （実験11 check_inplane_direction_conservation_v1.py と同一の p,q 構成、SHA台帳照合）。
木原の問い「4パネル図ではインフレーションが起きていないように見える」の検証図。
出力: fig_inflation_control_vs_seedless_zoom.png（0..2500歩）、同 _full.png（全域）、
      inflation_control_vs_seedless_v1.json（開始値・到達step・増幅期傾き）。"""
import hashlib
import json
import os

import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(SERIES))))
REPRO = os.path.join(ROOT, '対照実験_N掃引1to20_三系_v2_対照再現_20260906')
LONG = os.path.join(SERIES, 'N3_N40_long10000_20260905')

N = 12


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


# (a) 論文A対照再現の f₂ 時系列
srcs = {
    'paperA mixed ctl (seed d=0.01) [m_f2]':
        ('nsweep_mixed_T42000_d0.01_rep-controlrep1_N12_v2.npz', 'm_f2'),
    'paperA neutral ctl (seed d=0.01) [m_f2]':
        ('nsweep_neutral_T42000_d0.01_rep-controlrep1_N12_v2.npz', 'm_f2'),
    'paperA vacuum (d=0, no seed) [v_f2]':
        ('nsweep_vacuum_N12_v2.npz', 'v_f2'),
}
curves = {}
meta = {}
for label, (fn, key) in srcs.items():
    path = os.path.join(REPRO, fn)
    f = np.asarray(np.load(path)[key], dtype=float)
    curves[label] = (np.arange(f.size, dtype=float), f)
    meta[label] = {'file': fn, 'key': key, 'sha256': sha(path)}

# (b) 現行シード無し系列 N=12 の H⊥/H（SHA台帳ゲート）
ledger = {}
with open(os.path.join(LONG, 'SHA256SUMS.txt')) as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]
rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
path = os.path.join(LONG, rel)
assert ledger[rel] == sha(path), f'INPUT GATE FAIL: {rel}'
Z = np.asarray(np.load(path)['Z'], dtype=np.complex128)
rp, ip = Z[0].real, Z[0].imag
p = rp / np.linalg.norm(rp)
q = ip - (ip @ p) * p
q = q / np.linalg.norm(q)
h1 = (np.abs(Z @ p) ** 2 + np.abs(Z @ q) ** 2) / np.sum(np.abs(Z) ** 2, axis=1)
f_seedless = np.abs(1.0 - h1)
label_sl = f'this series N={N} den={N} seedless (d~1e-15) [Hperp/H]'
curves[label_sl] = (np.arange(f_seedless.size, dtype=float), f_seedless)
meta[label_sl] = {'file': rel, 'key': '1-h1', 'sha256': ledger[rel]}


# 数値要約（開始値・到達step・増幅期の傾き）
def summarize(f):
    f = np.asarray(f, float)
    def first(th):
        i = int(np.argmax(f >= th))
        return i if f[i] >= th else None
    m = (f > 1e-14) & (f < 1e-4)
    idx = np.where(m)[0]
    slope = float(np.polyfit(idx, np.log10(f[idx]), 1)[0]) if len(idx) >= 5 else None
    return {'f_step0': float(f[0]), 'reach_1e-10': first(1e-10), 'reach_1e-6': first(1e-6),
            'reach_1e-3': first(1e-3), 'reach_1e-1': first(1e-1), 'final': float(f[-1]),
            'ramp_slope_log10_per_step_(1e-14..1e-4)': slope,
            'steps_per_decade': (1.0 / slope) if slope else None,
            'n_points_in_ramp_window': int(len(idx))}

summary = {lab: summarize(f) for lab, (t, f) in curves.items()}
for lab, s in summary.items():
    print(f"{lab}: f0={s['f_step0']:.2e} | 1e-10→{s['reach_1e-10']} 1e-3→{s['reach_1e-3']} "
          f"1e-1→{s['reach_1e-1']} | ramp点数={s['n_points_in_ramp_window']} "
          f"slope={s['ramp_slope_log10_per_step_(1e-14..1e-4)']} step/桁={s['steps_per_decade']}")


def draw(fname, xlim, title):
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for lab, (t, f) in curves.items():
        ax.semilogy(t, f, lw=1.4, label=lab)
    ax.set_xlabel('step')
    ax.set_ylabel('Hperp/H (complement projection)  [= paperA f2]')
    ax.set_title(title)
    ax.grid(alpha=0.3, which='both')
    ax.legend(loc='lower right', fontsize=8)
    ax.set_ylim(1e-34, 2.0)
    if xlim:
        ax.set_xlim(*xlim)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, fname), dpi=160)
    plt.close(fig)
    print('wrote', fname)

draw('fig_inflation_control_vs_seedless_zoom.png', (0, 2500),
     f'N={N}: paperA control (seed d=0.01 vs vacuum d=0) vs this series seedless — H⊥/H rise, steps 0..2500')
draw('fig_inflation_control_vs_seedless_full.png', None,
     f'N={N}: paperA control (seed d=0.01 vs vacuum d=0) vs this series seedless — H⊥/H rise, full range')

with open(os.path.join(BASE, 'inflation_control_vs_seedless_v1.json'), 'w') as fh:
    json.dump({'N': N, 'inputs': meta, 'summary': summary}, fh, indent=2, ensure_ascii=False)
print('ALL DONE')
