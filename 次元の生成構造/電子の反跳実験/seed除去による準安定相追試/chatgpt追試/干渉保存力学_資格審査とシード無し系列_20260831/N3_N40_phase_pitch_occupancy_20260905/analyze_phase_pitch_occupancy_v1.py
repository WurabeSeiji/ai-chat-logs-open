#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位相ピッチ・占有数・格子オフセットの全系列調査（N=3..40、den=N、読み出しのみ）。

検証する仮説（木原 2026-09-05）:
 (H-a) 終状態の位相は pitch=180°/N の格子にロックする（N=3 の 60° の一般形）。
 (H-b) 占有数は1位置あたり (N−1)/2 本。奇数 N でのみ整数なので、偶数 N は等分配が
       不可能——偶数系では格子オフセットが 1/2 ピッチずつずれる（旋回する）のではないか。
 (H-c) 初期状態は 90° ピッチ格子。増幅中は占有数・振幅分布とも凍結か。

入力: ../N3_N40_stage123_sweep_20260905/results/hm_N{N}_den_{N}_states_500.npz
      （SHA256 入力ゲートで正本台帳と照合）、summary CSV（onset 引用）。
測定（すべて機械判定）:
 (1) step0: 位相クラスタ（円環ギャップ分割）→ 位置・占有数・振幅、90°格子への
     最良フィット偏差。
 (2) 中間ステップ（onset の 60%）: 占有数と振幅多重集合が step0 と一致するか（凍結検査）。
 (3) 終状態（step500）: クラスタ → 180°/N 格子への最良フィット偏差、占有数分布、
     振幅の等振幅偏差。
 (4) 格子オフセットのドリフト: 飽和後の各フレーム s（N 歩=時計1回転ごと）の格子
     オフセット u(s)（ピッチ単位、mod 1）と隣接フレーム差 Δu。1/2 ピッチ交互なら
     Δu ≈ 0.5 が出る。
出力: phase_pitch_table_v1.csv、phase_pitch_analysis_v1.json、fig_phase_pitch_N3_N40.png"""
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
SAT_MARGIN = 50

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

onset = {}
with open(os.path.join(RES, 'summary_64bit_with124_N3_N40.csv')) as f:
    for row in csv.DictReader(f):
        if row['series'] == 'N':
            onset[int(row['N'])] = int(row['onset_gt_0.05'])

def clusters(z, gap_thresh):
    """円環上の位相クラスタ分割。返り値: [(中心位相rad, メンバー添字list), ...]"""
    ph = np.angle(z)
    order = np.argsort(ph)
    ps = ph[order]
    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    cut = np.flatnonzero(gaps > gap_thresh)
    if cut.size == 0:
        return [(float(np.angle(np.mean(np.exp(1j * ph)))), list(range(len(z))))]
    groups = []
    start = (cut[-1] + 1) % len(ps)
    idx = list(order[start:]) + list(order[:start])
    bounds = sorted(((c - start) % len(ps)) for c in cut)
    prev = 0
    for b in bounds:
        mem = idx[prev:b + 1]
        c = float(np.angle(np.sum(np.exp(1j * ph[mem]))))
        groups.append((c, mem))
        prev = b + 1
    return groups

def grid_fit(centers, pitch):
    """centers を offset+k*pitch 格子にフィット。返り値: (offset_in_pitch mod1, 最大偏差 rad)"""
    u = np.array([c / pitch for c in centers])
    mean = float(np.angle(np.sum(np.exp(2j * math.pi * u)))) / (2 * math.pi)
    dev = np.array([abs((x - mean + 0.5) % 1.0 - 0.5) for x in u]) * pitch
    return mean % 1.0, float(dev.max())

rows = []
detail = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    assert ledger[rel] == sha256(os.path.join(PKG, rel)), f'INPUT GATE FAIL: {rel}'
    d = np.load(os.path.join(PKG, rel))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    Z = np.asarray(d['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    pitch = math.pi / N

    # (1) step0
    g0 = clusters(Z[0], 0.2)
    occ0 = sorted((len(m) for _, m in g0), reverse=True)
    off90, dev90 = grid_fit([c for c, _ in g0], math.pi / 2)
    amps0 = np.abs(Z[0])
    amp_groups0 = [{'center_deg': math.degrees(c) % 360.0, 'count': len(m),
                    'amp_min': float(amps0[m].min()), 'amp_max': float(amps0[m].max())}
                   for c, m in sorted(g0, key=lambda t: t[0])]

    # (2) 中間（凍結検査）
    smid = max(1, int(0.6 * onset[N]))
    gm = clusters(Z[smid], 0.2)
    occm = sorted((len(m) for _, m in gm), reverse=True)
    ampm = np.abs(Z[smid])
    freeze_occ = bool(occ0 == occm)
    freeze_amp = float(np.max(np.abs(np.sort(ampm) - np.sort(amps0))))

    # (3) 終状態
    gf = clusters(Z[500], pitch / 4)
    occf = sorted((len(m) for _, m in gf), reverse=True)
    offf, devf = grid_fit([c for c, _ in gf], pitch)
    ampsf = np.abs(Z[500])
    target = float(np.linalg.norm(Z[500])) / math.sqrt(M)
    eqdev = float(np.max(np.abs(ampsf - target)) / target)

    # (4) 格子オフセットのドリフト（飽和後フレーム）
    sc = onset[N]
    frames = [s for s in range(0, 501 - N, N) if sc >= 0 and s >= sc + SAT_MARGIN]
    du_list = []
    u_prev = None
    for s in frames:
        gs = clusters(Z[s], pitch / 4)
        u, _ = grid_fit([c for c, _ in gs], pitch)
        if u_prev is not None:
            du_list.append((u - u_prev) % 1.0)
        u_prev = u
    du_arr = np.array(du_list) if du_list else np.array([np.nan])
    du_c = np.angle(np.nanmean(np.exp(2j * math.pi * du_arr))) / (2 * math.pi) % 1.0 if du_list else float('nan')
    du_absmean = float(np.nanmean(np.minimum(du_arr % 1.0, 1.0 - du_arr % 1.0))) if du_list else float('nan')

    per_pos = (N - 1) / 2
    rows.append([N, M, N % 2, onset[N],
                 len(g0), '|'.join(map(str, occ0)), math.degrees(dev90),
                 freeze_occ, freeze_amp,
                 len(gf), '|'.join(map(str, occf)), math.degrees(devf),
                 per_pos, bool(len(gf) == N and all(c == (N - 1) // 2 for c in occf)),
                 eqdev, du_c, du_absmean, len(frames)])
    final_groups = [{'center_deg': math.degrees(c) % 360.0, 'count': len(m),
                     'amp_min': float(ampsf[m].min()), 'amp_max': float(ampsf[m].max())}
                    for c, m in sorted(gf, key=lambda t: t[0])]
    detail[N] = {'step0_groups': amp_groups0,
                 'final_groups': final_groups,
                 'final_occupancy': occf, 'final_n_clusters': len(gf),
                 'final_grid_offset_pitch_units': offf,
                 'du_per_frame_pitch_units': [float(x) for x in du_arr[:20]] if du_list else []}
    print('done N', N, flush=True)

HEADER = ['N', 'M', 'parity', 'onset',
          'n_clusters_0', 'occupancy_0', 'dev90_deg_0',
          'freeze_occupancy', 'freeze_amp_maxdev',
          'n_clusters_f', 'occupancy_f', 'devPitch_deg_f',
          'target_per_pos', 'exact_N_by_halfNm1',
          'eqmod_rel_dev_f', 'du_mean_pitch', 'du_absmean_pitch', 'n_late_frames']
with open(os.path.join(BASE, 'phase_pitch_table_v1.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(HEADER); w.writerows(rows)
with open(os.path.join(BASE, 'phase_pitch_analysis_v1.json'), 'w') as f:
    json.dump({'table_header': HEADER, 'table': rows, 'detail': detail}, f, indent=2, default=str)

tab = {r[0]: r for r in rows}
Ns = list(range(3, 41))
odd = [n for n in Ns if n % 2 == 1]; even = [n for n in Ns if n % 2 == 0]
fig, axs = plt.subplots(4, 1, figsize=(12, 16))
ax = axs[0]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    ax.semilogy(grp, [max(tab[n][11], 1e-12) for n in grp], 'o', color=c, label=lab)
ax.set_ylabel('max dev from pi/N grid [deg]'); ax.set_title('Final-state fit to 180/N-pitch lattice')
ax.grid(alpha=.3); ax.legend()
ax = axs[1]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    ax.plot(grp, [tab[n][9] for n in grp], 'o', color=c, label=lab)
ax.plot(Ns, Ns, '--', color='gray', linewidth=0.8, label='N')
ax.plot(Ns, [2 * n for n in Ns], ':', color='gray', linewidth=0.8, label='2N')
ax.set_ylabel('final n_clusters'); ax.grid(alpha=.3); ax.legend()
ax = axs[2]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    xs = [n for n in grp if tab[n][17] >= 2]
    ax.plot(xs, [tab[n][16] for n in xs], 'o', color=c, label=lab)
ax.axhline(0.5, color='gray', linewidth=0.8, linestyle='--')
ax.set_ylabel('|du| per clock turn [pitch]'); ax.set_title('Lattice-offset drift per clock turn (0.5 = half-pitch alternation)')
ax.grid(alpha=.3); ax.legend()
ax = axs[3]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    ax.plot(grp, [tab[n][4] for n in grp], 'o', color=c, label=lab, alpha=0.6)
ax.set_ylabel('step0 n_clusters'); ax.set_xlabel('N'); ax.grid(alpha=.3); ax.legend()
fig.suptitle('Phase pitch / occupancy / lattice drift across N=3..40 (den=N, readout only)', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_phase_pitch_N3_N40.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
