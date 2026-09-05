#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成子スペクトルの流れとランクの解析（読み出しのみ・新規走行なし）。

背景（木原 2026-09-06）: データの固有値・固有ベクトル・ランクの基本解析が未実施。
実施済みは初期生成子のスペクトルのみ（回転数ロック調査・倍音走査）。本解析は:
 (1) スペクトル流: K(θ(t)) の σ スペクトル（グラム縮約 [S3] 経由、2N×2N で厳密）を
     軌道に沿って標本化。遷移前後で固有値がどう動くか、σ_max・縮退対の運命。
 (2) グラム階数 r_G(t): gram_reduce の保持ランク（τ_G=1e-12）。
 (3) 軌道の実効ランク: 窓ごとの状態行列（標本×M）の特異値から
     参加比 R_eff = (Σs²)²/Σs⁴ と数値ランク（s > 1e-8 s_max の個数）。
     窓 = 増幅期 [0, onset−10] / 遷移 [onset−10, onset+40] / 終盤 [9000, 10000]。
 (4) 垂直成長の実効ランク: 親平面 Π=span(p,q) の補直交成分軌道の同量
     （「成長は新1平面か」の直接判定）。
入力: 10,000歩スイープ npz（SHA台帳照合）＋回転数ロック調査の onset。
出力: spectral_flow_rank_table_v1.csv、analyze_spectral_flow_rank_v1.json、
      fig_spectral_flow_smallN.png"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[5]
ENGINE = ROOT / "時間軸Q軸とフェルミオンの生成構造" / "検証_対照実験" / "第5論文原本_自発的分裂予備実験_v1"
sys.path.insert(0, str(ENGINE))
sys.path.insert(0, str(ENGINE / "exact_lowN_eigenspectrum_v2" / "code"))
from run_n_scaling_lowrank_v1 import LowRankSystem
from run_n300_dimension_saturation_v2 import gram_reduce

SWEEP = BASE.parent / "N3_N40_long10000_20260905"
RL = BASE.parent / "N3_N40_rotation_lock_analysis_20260905"
NS = [3, 4, 5, 6, 7, 10, 40]

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

rl = json.load(open(RL / 'rotation_lock_analysis_v1.json'))
onset = {r[0]: int(r[2]) for r in rl['table']}

def eff_rank(Mrows):
    """複素行列（標本×M）の特異値から (数値ランク, 参加比 R_eff, 上位特異値)"""
    s = np.linalg.svd(Mrows, compute_uv=False)
    if s.size == 0 or s[0] == 0:
        return 0, 0.0, []
    nrank = int(np.sum(s > 1e-8 * s[0]))
    p = s ** 2
    reff = float(p.sum() ** 2 / np.sum(p ** 2))
    return nrank, reff, [float(x) for x in s[:8]]

rows = []
flow = {}
for N in NS:
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'], dtype=np.complex128)
    sys_lr = LowRankSystem(N)
    sc = onset[N]
    # (1)(2) スペクトル流: 遷移周辺は細かく、全域は粗く
    samples = sorted(set(list(range(0, min(2 * sc + 50, 1001), 10))
                         + list(range(0, 10001, 250)) + [10000]))
    spec = []
    for s in samples:
        gr = gram_reduce(sys_lr, Z[s])
        mu = np.sort(gr['mu'])[::-1]
        pos = mu[mu > 1e-9]
        spec.append({'step': s, 'r_G': int(gr['diag']['r_G']) if 'diag' in gr else int(len(gr['lam_r'])),
                     'sigma_top8': [float(x) for x in pos[:8]],
                     'n_pos': int(pos.size)})
    flow[str(N)] = spec
    # スペクトルの前後比較
    sig0 = spec[0]['sigma_top8']
    sigF = spec[-1]['sigma_top8']
    # (3) 軌道の実効ランク（複素行ベクトルの行列）
    w_ramp = Z[0:max(sc - 10, 2)]
    w_tr = Z[max(sc - 10, 0):min(sc + 40, 10001)]
    w_late = Z[9000:10001:10]
    nr_r, re_r, _ = eff_rank(w_ramp[::max(1, len(w_ramp) // 200)])
    nr_t, re_t, _ = eff_rank(w_tr)
    nr_l, re_l, _ = eff_rank(w_late)
    # (4) 垂直成長の実効ランク
    z0 = Z[0]
    p = z0.real / np.linalg.norm(z0.real)
    q = z0.imag - (z0.imag @ p) * p; q = q / np.linalg.norm(q)
    def perp(rows_):
        return rows_ - np.outer(rows_ @ p, p) - np.outer(rows_ @ q, q)
    nrp_r, rep_r, sp_r = eff_rank(perp(w_ramp[::max(1, len(w_ramp) // 200)]))
    nrp_t, rep_t, sp_t = eff_rank(perp(w_tr))
    nrp_l, rep_l, sp_l = eff_rank(perp(w_late))
    rows.append([N, sc,
                 round(re_r, 2), round(re_t, 2), round(re_l, 2),
                 nr_r, nr_t, nr_l,
                 round(rep_r, 2), round(rep_t, 2), round(rep_l, 2),
                 nrp_r, nrp_t, nrp_l,
                 spec[0]['n_pos'], spec[-1]['n_pos'],
                 round(sig0[0], 6) if sig0 else 0, round(sigF[0], 6) if sigF else 0])
    print(f"N={N}: onset={sc} | 軌道R_eff 増幅期={re_r:.2f} 遷移={re_t:.2f} 終盤={re_l:.2f} | "
          f"垂直R_eff {rep_r:.2f}/{rep_t:.2f}/{rep_l:.2f}（数値ランク {nrp_r}/{nrp_t}/{nrp_l}） | "
          f"σ_max {sig0[0] if sig0 else 0:.4f}→{sigF[0] if sigF else 0:.4f} 正枝本数 {spec[0]['n_pos']}→{spec[-1]['n_pos']}",
          flush=True)

HEADER = ['N', 'onset',
          'Reff_orbit_ramp', 'Reff_orbit_trans', 'Reff_orbit_late',
          'nrank_orbit_ramp', 'nrank_orbit_trans', 'nrank_orbit_late',
          'Reff_perp_ramp', 'Reff_perp_trans', 'Reff_perp_late',
          'nrank_perp_ramp', 'nrank_perp_trans', 'nrank_perp_late',
          'n_pos_sigma_step0', 'n_pos_sigma_final', 'sigma_max_step0', 'sigma_max_final']
with open(BASE / 'spectral_flow_rank_table_v1.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(HEADER); w.writerows(rows)
with open(BASE / 'analyze_spectral_flow_rank_v1.json', 'w') as fh:
    json.dump({'table_header': HEADER, 'table': rows, 'spectral_flow': flow},
              fh, indent=2, default=str)

fig, axs = plt.subplots(2, 3, figsize=(16, 9))
for ax, N in zip(axs.ravel(), [3, 4, 5, 6, 7, 10]):
    sp = flow[str(N)]
    ts = [r['step'] for r in sp]
    kmax = max(len(r['sigma_top8']) for r in sp)
    for k in range(kmax):
        ys = [r['sigma_top8'][k] if k < len(r['sigma_top8']) else np.nan for r in sp]
        ax.plot(ts, ys, lw=0.9)
    ax.axvline(onset[N], color='k', ls=':', lw=0.8)
    ax.set_xlim(0, min(2 * onset[N] + 200, 2000))
    ax.set_title(f'N={N}: sigma spectrum flow (top branches)')
    ax.grid(alpha=.3)
fig.suptitle('Generator spectrum flow along trajectory (gram-reduced, exact); dotted = onset', y=.995)
fig.tight_layout()
fig.savefig(BASE / 'fig_spectral_flow_smallN.png', dpi=140)
plt.close(fig)
print('ALL DONE')
