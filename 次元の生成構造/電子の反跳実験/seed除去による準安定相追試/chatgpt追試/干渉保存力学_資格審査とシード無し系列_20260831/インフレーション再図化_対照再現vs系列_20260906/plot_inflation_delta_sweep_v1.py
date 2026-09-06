#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A正本のシード強度 δ 掃引（5系 × 8段・N=12・T42000）のインフレーション図
（読み出しのみ・新規走行なし・正本は読み取りのみで不可侵）。
plot_inflation_control_vs_seedless_v1.py（= N40 現行様式の verbatim）から入力だけ差し替え:
正本 `../../../../../対照実験_N掃引1to20_三系_v2/nsweep_{mode}_T42000_d{δ}_N12_v2.npz` の m_f2。
木原の問い「より少ないシードからの系列は？」→ 正本に完備。δ 依存性を1枚に重ね描きし、
増幅率（増幅窓 10⁻¹⁴..10⁻⁴ の log₁₀ 傾き）と到達時刻の δ 依存を要約図で示す。
出力: fig_inflation_delta_sweep_5modes.png（2×3 パネル・各系に全δ重ね描き）、
      fig_inflation_delta_sweep_summary.png（傾き・ランプ点数・10⁻¹到達 vs δ）、
      inflation_delta_sweep_v1.json（全 (mode, δ) の統計と入力SHA）。"""
import glob
import hashlib
import json
import os
import re

import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(SERIES))))
CANON = os.path.join(ROOT, '対照実験_N掃引1to20_三系_v2')   # 正本（読み取りのみ）

MODES = ['neutral', 'electron', 'mixed', 'fermion_family', 'boson_family']
N = 12


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


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


data = {}     # mode -> {delta: f2}
meta = {}
stats = {}
for mode in MODES:
    files = sorted(glob.glob(os.path.join(CANON, f'nsweep_{mode}_T42000_d*_N{N}_v2.npz')))
    files = [x for x in files if 'rep-' not in os.path.basename(x)]
    data[mode] = {}
    meta[mode] = {}
    stats[mode] = {}
    for fn in files:
        d = float(re.search(r'_d([0-9.e-]+)_N12', os.path.basename(fn)).group(1))
        f = np.asarray(np.load(fn)['m_f2'], dtype=float)
        data[mode][d] = f
        meta[mode][repr(d)] = {'file': os.path.basename(fn), 'key': 'm_f2', 'sha256': sha(fn)}
        stats[mode][repr(d)] = summarize(f)
    print(f'{mode}: δ = ' + ' '.join(f'{d:g}' for d in sorted(data[mode])), flush=True)

# 図1: 2×3 パネル、各系に全 δ を重ね描き（描画様式は N40 現行様式 verbatim）
fig, axs = plt.subplots(2, 3, figsize=(16.4, 10.4))
axs = axs.ravel()
for k, mode in enumerate(MODES):
    ax = axs[k]
    for d in sorted(data[mode]):
        t = np.arange(data[mode][d].size, dtype=float)
        ax.semilogy(t, data[mode][d], lw=1.4, label=f'd={d:g}')
    ax.set_xlabel('step')
    ax.set_ylabel('Hperp/H (= paperA f2)')
    ax.set_title(f'{mode}  N={N} T42000: seed delta sweep')
    ax.grid(alpha=0.3, which='both')
    ax.legend(loc='lower right', fontsize=8)
    ax.set_ylim(1e-34, 2.0)
    ax.set_xlim(0, 2500)
axs[5].axis('off')
fig.suptitle('paperA canonical: inflation (Hperp/H rise) vs seed strength delta — 5 seed modes, N=12, steps 0..2500', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_inflation_delta_sweep_5modes.png'), dpi=150)
plt.close(fig)
print('wrote fig_inflation_delta_sweep_5modes.png')

# 図2: 要約（傾き・ランプ点数・10⁻¹到達 vs δ）
fig, axs = plt.subplots(1, 3, figsize=(16.4, 4.8))
for mode in MODES:
    ds = sorted(data[mode])
    sl = [stats[mode][repr(d)]['ramp_slope_log10_per_step_(1e-14..1e-4)'] for d in ds]
    npts = [stats[mode][repr(d)]['n_points_in_ramp_window'] for d in ds]
    t1 = [stats[mode][repr(d)]['reach_1e-1'] for d in ds]
    axs[0].semilogx([d for d, s in zip(ds, sl) if s is not None],
                    [s for s in sl if s is not None], 'o-', lw=1.4, label=mode)
    axs[1].loglog(ds, [max(n, 0.5) for n in npts], 'o-', lw=1.4, label=mode)
    axs[2].semilogx(ds, t1, 'o-', lw=1.4, label=mode)
axs[0].axhline(0.018, color='gray', ls=':', lw=1)
axs[0].set_xlabel('seed delta'); axs[0].set_ylabel('ramp slope  log10(f2)/step'); axs[0].set_title('growth rate (window 1e-14..1e-4)')
axs[1].set_xlabel('seed delta'); axs[1].set_ylabel('# points in ramp window'); axs[1].set_title('ramp length (0 = no inflation)')
axs[2].set_xlabel('seed delta'); axs[2].set_ylabel('step at f2 = 1e-1'); axs[2].set_title('arrival time')
for ax in axs:
    ax.grid(alpha=0.3, which='both'); ax.legend(fontsize=8)
fig.suptitle('paperA canonical delta sweep: rate is seed-independent for delta <= 1e-4; ramp collapses at delta ~ 1e-3..1e-2', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_inflation_delta_sweep_summary.png'), dpi=150)
plt.close(fig)
print('wrote fig_inflation_delta_sweep_summary.png')

with open(os.path.join(BASE, 'inflation_delta_sweep_v1.json'), 'w') as fh:
    json.dump({'N': N, 'canonical_dir': CANON, 'inputs': meta, 'stats': stats}, fh, indent=2, ensure_ascii=False)
print('ALL DONE')
