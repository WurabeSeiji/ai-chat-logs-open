#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""インフレーション図の「対数軸で等間隔の準振動」と「一周期で円周に張り付く」を、既存データで測る（2026-10-19）。

入力（読み取りのみ）：N3_N40_stage123_sweep_20260905/results/hm_N40_den_{den}_states_500.npz（den = 38〜42, 124）。
  図（K4 論文 図4）の元データ make_parent 型スイープの npz は削除済み。stage123 の den=40 は図の H⊥/H 系列と一致
  （step500 で 0.1005、対数で最大 0.002 桁差）するので、同じ走行として使う。
読み出し：readout_copied_20261019.py（図の生成コードの忠実コピー）。新しい走行はしない。

測る量（事前に固定）：
  growth      : log10(H⊥/H) が (−28, −6) にある区間の一次回帰の傾き（桁/step）
  omega_app   : 床（onset−5 まで）での一歩あたりの見かけの回転 arg Σ z_{t+1} conj(z_t) の中央値（rad/step）
  P_rot       : 2π/|omega_app|（見かけの一回転のステップ数）
  P_wiggle    : 共回転フレーム W の位相の全波平均（床区間）の FFT ピーク周期（ハニング窓）
  dur         : log10(H⊥/H) が −3 を超えてから最終値（末尾 50 step 平均）−0.1 に達するまでのステップ数
  sigma_est   : omega_app から逆算した σ（σΔτ ≡ omega_app mod 2π、σ≈N 付近の枝）
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from readout_copied_20261019 import series_from_states, corotating_W  # noqa: E402

R = os.path.join(HERE, '..', '電子の反跳実験', 'seed除去による準安定相追試', 'chatgpt追試',
                 '干渉保存力学_資格審査とシード無し系列_20260831', 'N3_N40_stage123_sweep_20260905', 'results')
N = 40
DENS = [38, 39, 40, 41, 42, 124]
OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
FIG = os.path.join(HERE, 'figures'); os.makedirs(FIG, exist_ok=True)

rows = []
for den in DENS:
    S = np.load(os.path.join(R, f'hm_N{N}_den_{den}_states_500.npz'))['Z'].astype(np.complex128)
    X, Y, fz = series_from_states(S)
    W, th0, sign, fv, other, onset, ok = corotating_W(S, X, Y, fz, strict=False)
    lg = np.log10(fz)
    idx = np.where((lg > -28) & (lg < -6))[0]
    growth = float(np.polyfit(idx, lg[idx], 1)[0])
    pre = S[:max(onset - 5, 3)]
    omega_app = float(np.median(np.angle(np.sum(pre[1:] * np.conj(pre[:-1]), axis=1))))
    P_rot = 2 * np.pi / abs(omega_app)
    ph = np.unwrap(np.angle(W[:max(onset - 5, 3)]), axis=0); ph -= ph.mean(0)
    x = ph.mean(1); x -= x.mean()
    F = np.abs(np.fft.rfft(x * np.hanning(len(x)))); fr = np.fft.rfftfreq(len(x))
    k = int(np.argmax(F[1:]) + 1); P_wiggle = float(1 / fr[k])
    wiggle_amp = float(np.std(x) * np.sqrt(2))
    fin = float(lg[-50:].mean())
    t3 = int(np.argmax(lg > -3)); tf = int(np.argmax(lg > fin - 0.1))
    dtau = 2 * np.pi / den
    sigma_est = float((omega_app + 2 * np.pi * np.round(N * dtau / (2 * np.pi))) / dtau)
    rows.append(dict(den=den, floor_static_ok=bool(ok), onset=onset, growth_dec_per_step=growth,
                     growth_dec_per_tau=growth / dtau, omega_app=omega_app, P_rot=P_rot,
                     P_wiggle=P_wiggle, wiggle_amp_rad=wiggle_amp, t_lg_gt_m3=t3, t_final=tf,
                     dur=tf - t3, dur_over_P_rot=(tf - t3) / P_rot, sigma_est=sigma_est,
                     decades_per_wiggle=growth * P_wiggle, final_log10_Hperp=fin))

json.dump(rows, open(os.path.join(OUT, 'wiggle_lock_N40.json'), 'w'), ensure_ascii=False, indent=1)
with open(os.path.join(OUT, 'wiggle_lock_N40.md'), 'w', encoding='utf-8') as f:
    f.write('| den | 床静止判定 | ω_app (rad/step) | σ 推定 | 2π/|ω_app| | 波打ち周期 | 振幅 (rad) | 成長 (桁/step) | 成長 (桁/τ) | 遷移 dur (step) | dur/P |\n')
    f.write('|---|---|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        f.write(f"| {r['den']} | {'可' if r['floor_static_ok'] else '不可'} | {r['omega_app']:+.4f} | {r['sigma_est']:.3f} | "
                f"{r['P_rot']:.1f} | {r['P_wiggle']:.1f} | {r['wiggle_amp_rad']:.4f} | {r['growth_dec_per_step']:.4f} | "
                f"{r['growth_dec_per_tau']:.3f} | {r['dur']} | {r['dur_over_P_rot']:.2f} |\n")

ok = [r for r in rows if r['floor_static_ok']]
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
pr = [r['P_rot'] for r in ok]; pw = [r['P_wiggle'] for r in ok]
ax[0].loglog(pr, pw, 'o'); m = [min(pr) * 0.8, max(pr) * 1.2]; ax[0].loglog(m, m, 'k--', lw=0.8)
for r in ok:
    ax[0].annotate(f"den={r['den']}", (r['P_rot'], r['P_wiggle']), fontsize=8)
ax[0].set_xlabel('2π/|ω_app| (steps)'); ax[0].set_ylabel('wiggle period (steps)')
ax[0].set_title('wiggle period = apparent rotation period')
ax[1].plot([r['P_rot'] for r in ok], [r['dur'] for r in ok], 's')
for r in ok:
    ax[1].annotate(f"den={r['den']}", (r['P_rot'], r['dur']), fontsize=8)
ax[1].plot([0, 45], [0, 45], 'k--', lw=0.8, label='dur = 1 period')
ax[1].set_xlabel('2π/|ω_app| (steps)'); ax[1].set_ylabel('lock-in duration (steps)')
ax[1].set_title('lock-in duration vs period'); ax[1].legend(fontsize=8)
fig.tight_layout()
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(FIG, f'fig1_wiggle_period_and_lock_N40.{ext}'), dpi=140)
print(open(os.path.join(OUT, 'wiggle_lock_N40.md'), encoding='utf-8').read())
