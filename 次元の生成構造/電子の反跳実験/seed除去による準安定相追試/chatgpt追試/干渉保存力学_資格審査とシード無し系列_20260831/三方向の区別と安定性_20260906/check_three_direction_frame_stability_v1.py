#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三方向の区別と安定性 (1): 状態の複素平面（偏光面）の向きと法線の安定性（読み出しのみ）。

三方向の定義（木原 2026-09-06 指示: 初期条件に錨を置く）:
  X̂ = p = Re Z0 / ‖Re Z0‖（初期実軸）
  Ŷ = q = Im Z0 を p に直交化・正規化（初期虚軸）。90°位相差で立つ複素平面 Π0 = span(p,q)。
  Ẑ = Π0 の法線（N=3 では R³ の p×q と生成子回転軸 n が一致するはず——これも検定）。
検定:
 (a) 定義の厳密性: step0 の |cos∠(Re Z0, Im Z0)|（90°からのずれ）と ‖Im Z0‖/‖Re Z0‖。
 (b) 偏光面の向きの安定性: Πs(t) = span(Re Z(t), Im Z(t)) と Π0 の主角（最大主角、度）。
     グローバル位相 Z→e^{iα}Z で span は不変（ゲージ不変量）。全 10001 歩。
     縮退警戒: 正規化後の [Re, Im] 行列の最小特異値も記録（面が潰れていないか）。
 (c) N=3 のみ: 生成子軸 n(t) = (K_32, K_13, K_21)/‖·‖（3次元では回転軸＝面法線）と
     n(0) の角、および n(0) と p×q の角（法線同定の確認）。
入力: N3_N40_long10000_20260905 の npz（SHA台帳照合）。新規走行なし。
出力: check_three_direction_frame_stability_v1.{json,csv,png}
"""
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
from run_n_scaling_lowrank_v1 import LowRankSystem

SWEEP = BASE.parent / "N3_N40_long10000_20260905"
RL = BASE.parent / "N3_N40_rotation_lock_analysis_20260905"
NS = [3, 4, 5, 6, 10, 40]

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

rl = json.load(open(RL / 'rotation_lock_analysis_v1.json'))
onset = {r[0]: int(r[2]) for r in rl['table']}

out = {}
rows = []
fig, axs = plt.subplots(2, 3, figsize=(16, 8))
for ax, N in zip(axs.ravel(), NS):
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'], dtype=np.complex128)
    T, M = Z.shape
    rp, ip = Z[0].real, Z[0].imag
    cos90 = float(abs(rp @ ip) / (np.linalg.norm(rp) * np.linalg.norm(ip)))
    ratio = float(np.linalg.norm(ip) / np.linalg.norm(rp))
    p = rp / np.linalg.norm(rp)
    q = ip - (ip @ p) * p
    q = q / np.linalg.norm(q)
    B0 = np.stack([p, q], axis=1)  # M×2 正規直交
    # (b) 偏光面の主角（全歩）
    maxang = np.empty(T)
    smin = np.empty(T)
    for t in range(T):
        A = np.stack([Z[t].real, Z[t].imag], axis=1)
        nrm = np.linalg.norm(A, axis=0)
        smin[t] = float(np.linalg.svd(A / nrm, compute_uv=False)[-1]) if nrm.min() > 0 else 0.0
        Bt, _ = np.linalg.qr(A)
        sv = np.clip(np.linalg.svd(B0.T @ Bt, compute_uv=False), 0.0, 1.0)
        maxang[t] = np.degrees(np.arccos(sv[-1]))
    sc = onset[N]
    ramp_end = max(sc - 10, 1)
    res = {'M': M, 'onset': sc,
           'def_cos_90deg_offset': cos90, 'def_amp_ratio_im_over_re': ratio,
           'plane_maxangle_deg': {
               'ramp_max_(0..onset-10)': float(maxang[:ramp_end].max()),
               'transition_max': float(maxang[:min(sc + 200, T)].max()),
               'global_max': float(maxang.max()),
               'argmax_step': int(np.argmax(maxang)),
               'final_mean_last1000': float(maxang[-1000:].mean()),
               'final_std_last1000': float(maxang[-1000:].std())},
           'plane_smin_min': float(smin.min())}
    # (c) N=3: 生成子軸
    if N == 3:
        sys_lr = LowRankSystem(3)
        n_axis = np.empty((T, 3))
        for t in range(T):
            sys_lr.set_theta(np.angle(Z[t]))
            K = np.column_stack([sys_lr.kmatvec(np.eye(3)[:, j]) for j in range(3)])
            v = np.array([K[2, 1], K[0, 2], K[1, 0]])
            nv = np.linalg.norm(v)
            n_axis[t] = v / nv if nv > 0 else 0.0
        ang_n = np.degrees(np.arccos(np.clip(np.abs(n_axis @ n_axis[0]), 0, 1)))
        pxq = np.cross(p, q)
        ang_n0_pxq = float(np.degrees(np.arccos(np.clip(abs(pxq @ n_axis[0]), 0, 1))))
        res['axis_N3'] = {'angle_n0_vs_pxq_deg': ang_n0_pxq,
                          'ramp_max_deg': float(ang_n[:ramp_end].max()),
                          'global_max_deg': float(ang_n.max()),
                          'argmax_step': int(np.argmax(ang_n)),
                          'final_mean_last1000_deg': float(ang_n[-1000:].mean())}
        ax.plot(ang_n, lw=0.7, color='tab:red', label='generator axis n(t) vs n(0)')
    out[str(N)] = res
    rows.append([N, M, sc, f'{cos90:.2e}', f'{ratio:.6f}',
                 f"{res['plane_maxangle_deg']['ramp_max_(0..onset-10)']:.3e}",
                 f"{res['plane_maxangle_deg']['global_max']:.3f}",
                 res['plane_maxangle_deg']['argmax_step'],
                 f"{res['plane_maxangle_deg']['final_mean_last1000']:.3f}",
                 f"{res['plane_maxangle_deg']['final_std_last1000']:.2e}"])
    ax.plot(maxang, lw=0.7, color='tab:blue', label='polarization plane vs Π0')
    ax.axvline(sc, color='k', ls=':', lw=0.8)
    ax.set_title(f'N={N}: max principal angle (deg)')
    ax.set_xscale('log')
    ax.grid(alpha=.3)
    ax.legend(fontsize=7)
    print(f"N={N}: 90°ずれ={cos90:.1e} 振幅比Im/Re={ratio:.6f} | 面の主角: 増幅期max="
          f"{res['plane_maxangle_deg']['ramp_max_(0..onset-10)']:.2e}° 全域max="
          f"{res['plane_maxangle_deg']['global_max']:.3f}°(step {res['plane_maxangle_deg']['argmax_step']}) "
          f"終盤平均={res['plane_maxangle_deg']['final_mean_last1000']:.3f}°±"
          f"{res['plane_maxangle_deg']['final_std_last1000']:.1e}", flush=True)
    if N == 3:
        a3 = res['axis_N3']
        print(f"  N=3軸: n(0)とp×qの角={a3['angle_n0_vs_pxq_deg']:.2e}° 増幅期max={a3['ramp_max_deg']:.2e}° "
              f"全域max={a3['global_max_deg']:.3f}°(step {a3['argmax_step']}) 終盤={a3['final_mean_last1000_deg']:.3f}°",
              flush=True)

HEADER = ['N', 'M', 'onset', 'cos_90deg_offset', 'amp_ratio_im_re',
          'plane_ang_ramp_max_deg', 'plane_ang_global_max_deg', 'argmax_step',
          'plane_ang_final_mean_deg', 'plane_ang_final_std_deg']
with open(BASE / 'check_three_direction_frame_stability_v1.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(HEADER); w.writerows(rows)
with open(BASE / 'check_three_direction_frame_stability_v1.json', 'w') as fh:
    json.dump(out, fh, indent=2)
fig.suptitle('Polarization-plane orientation drift vs Π0 = span(Re Z0, Im Z0); dotted = onset', y=.995)
fig.tight_layout()
fig.savefig(BASE / 'check_three_direction_frame_stability_v1.png', dpi=140)
print('ALL DONE')
