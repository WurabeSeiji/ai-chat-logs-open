#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""シード型（ボゾン帯／フェルミオン帯）が効く物理量の横並び比較（読み出しのみ・新規走行なし）。
論文A正本 `対照実験_N掃引1to20_三系_v2/` の δ=1e-2・T42000・N=12 の5系:
  neutral(1セル k=1,η=0) / electron(1セル k=1,η=3) / mixed(8セル) /
  fermion_family(5セル k=1) / boson_family(3セル k=6)、＋ 真空腕（δ=0, v_ 配列）。
木原の問い「ボゾン型/フェルミオン型はインフレーションに効かないが読出しには効く——
他にどの物理に差が出るか」に、記録済み全観測量（99キー）を窓 [21000,42000] 中央値で
横並びにして答える。加えて奇数帯パワーの**全步厳密零**（偶奇選択則）を検定する。
出力: mode_physics_table_v1.csv / check_mode_physics_diff_v1.json /
      fig_mode_bands.png（帯別パワー）/ fig_mode_ledger.png（128セル帳簿）/
      fig_mode_timeseries.png（r・Q̂・担い手・奇数帯の時系列）"""
import csv
import hashlib
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(SERIES))))
CANON = os.path.join(ROOT, '対照実験_N掃引1to20_三系_v2')   # 正本（読み取りのみ）

FILES = {
    'neutral':        'nsweep_neutral_T42000_N12_v2.npz',
    'electron':       'nsweep_electron_T42000_d0.01_N12_v2.npz',
    'mixed':          'nsweep_mixed_T42000_N12_v2.npz',
    'fermion_family': 'nsweep_fermion_family_T42000_d0.01_N12_v2.npz',
    'boson_family':   'nsweep_boson_family_T42000_d0.01_N12_v2.npz',
}
WIN = (21000, 42000)   # 宣言窓 [T//2, T]

# 比較する量（窓中央値）。(キー, 表示名, 分類)
SCALARS = [
    ('m_f2', 'f2 = H⊥/H（空間・インフレ後）', '空間'),
    ('m_n_eff', 'n_eff（実効次元）', '空間'),
    ('m_align', 'align', '空間'),
    ('rec_m_odd_power', 'P_odd（奇数帯パワー）', '帯偶奇'),
    ('rec_m_even_power', 'P_even（偶数帯パワー）', '帯偶奇'),
    ('rec_m_k3_power', 'P_k3（相棒帯3）', '帯偶奇'),
    ('rec_m_pump_power', 'P_pump（k=2）', '帯偶奇'),
    ('rec_m_total_power', 'P_total', '帯偶奇'),
    ('rec_m_r_mean', 'r = sin²θ（混合率・node平均）', '電荷/α'),
    ('rec_m_r_nopump', 'r_nopump（物質側）', '電荷/α'),
    ('rec_m_absdist_alpha', '|r − R_α|', '電荷/α'),
    ('rec_m_absdist_alpha_nopump', '|r_nopump − R_α|', '電荷/α'),
    ('rec_m_q_hat', 'Q̂（電荷読出し）', '電荷/α'),
    ('rec_m_dom_m', 'm̂（優勢巻き）', '電荷/α'),
    ('rec_m_readable', '可読率', '電荷/α'),
    ('rec_m_conc_k3', '相棒帯3 巻き集中度', '電荷/α'),
    ('rec_m_conc_partner', '相棒 巻き集中度', '電荷/α'),
    ('rec_m_partner_q_hat', '相棒 Q̂', '電荷/α'),
    ('rec_m_partner_readable', '相棒 可読率', '電荷/α'),
    ('rec_m_partner_power', '相棒パワー', '電荷/α'),
    ('rec_m_excl_ratio', '排他比（狙い/非狙い）', '帳簿'),
    ('rec_m_target_power', '狙いセル総パワー', '帳簿'),
    ('rec_m_nontarget_power', '非狙いセル総パワー', '帳簿'),
    ('rec_m_seed_power', 'シードセルパワー', '帳簿'),
    ('m_carrier_power', '時間の担い手パワー', '時間'),
    ('m_omega_hat', 'ω̂', '時間'),
    ('m_phi_weight', 'φ重み', '時間'),
    ('m_coherence', 'コヒーレンス', '時間'),
    ('m_closure', '閉塞残差（全体）', '閉塞'),
    ('m_cond_closure', '閉塞残差（凝縮体セル）', '閉塞'),
    ('m_seed_closure', '閉塞残差（シードセル）', '閉塞'),
]
EXACT_ZERO_KEYS = ['rec_m_odd_power', 'rec_m_k3_power', 'rec_m_odd_amp_max', 'rec_m_partner_power']


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def wmed(x):
    w = np.asarray(x, float)[WIN[0]:WIN[1]]
    w = w[np.isfinite(w)]
    return float(np.median(w)) if w.size else float('nan')


data, meta = {}, {}
for mode, fn in FILES.items():
    path = os.path.join(CANON, fn)
    data[mode] = np.load(path, allow_pickle=False)
    F = data[mode].files
    meta[mode] = {'file': fn, 'sha256': sha(path),
                  'seed_cells': data[mode]['seed_cells_index'].tolist() if 'seed_cells_index' in F else '(旧形式: 未記録)',
                  'seed_delta': data[mode]['seed_cells_delta'].tolist() if 'seed_cells_delta' in F else None,
                  'n_keys': len(F)}
    print(f"{mode}: cells={meta[mode]['seed_cells']}", flush=True)

modes = list(FILES)
cols = modes + ['vacuum(δ=0)']
table = []
for key, name, cat in SCALARS:
    row = [cat, name, key]
    for mode in modes:
        row.append(wmed(data[mode][key]) if key in data[mode].files else float('nan'))
    vkey = key.replace('rec_m_', 'rec_v_').replace('m_', 'v_', 1) if key.startswith('m_') else key.replace('rec_m_', 'rec_v_')
    row.append(wmed(data['neutral'][vkey]) if vkey in data['neutral'] else float('nan'))
    table.append(row)

# 厳密零検定（全42000步の最大値）
exact = {}
for mode in modes:
    exact[mode] = {}
    for key in EXACT_ZERO_KEYS:
        if key not in data[mode].files:
            exact[mode][key] = {'max_all_steps': None, 'exact_zero': None, 'note': 'key absent (旧形式)'}
            continue
        a = np.asarray(data[mode][key], float)
        a = a[np.isfinite(a)]
        exact[mode][key] = {'max_all_steps': float(np.max(np.abs(a))) if a.size else None,
                            'exact_zero': bool(a.size and np.max(np.abs(a)) == 0.0)}
# 真空腕
exact['vacuum(δ=0)'] = {}
for key in EXACT_ZERO_KEYS:
    vk = key.replace('rec_m_', 'rec_v_')
    if vk not in data['neutral'].files:
        exact['vacuum(δ=0)'][key] = {'max_all_steps': None, 'exact_zero': None, 'note': 'key absent'}
        continue
    a = np.asarray(data['neutral'][vk], float); a = a[np.isfinite(a)]
    exact['vacuum(δ=0)'][key] = {'max_all_steps': float(np.max(np.abs(a))) if a.size else None,
                                 'exact_zero': bool(a.size and np.max(np.abs(a)) == 0.0)}

# 帯別パワー（窓中央値）と帳簿最終スナップショット
bands = {mode: np.median(np.asarray(data[mode]['rec_m_bands'], float)[WIN[0]:WIN[1]], axis=0) for mode in modes if 'rec_m_bands' in data[mode].files}
if 'rec_v_bands' in data['neutral'].files:
    bands['vacuum(δ=0)'] = np.median(np.asarray(data['neutral']['rec_v_bands'], float)[WIN[0]:WIN[1]], axis=0)
ledger = {mode: np.asarray(data[mode]['rec_m_ledger'], float)[-1] for mode in modes if 'rec_m_ledger' in data[mode].files}
if 'rec_v_ledger' in data['neutral'].files:
    ledger['vacuum(δ=0)'] = np.asarray(data['neutral']['rec_v_ledger'], float)[-1]
occ = {}
for k, L in ledger.items():
    tot = L.sum()
    occ[k] = {'n_cells_gt_1e-12_total': int(np.sum(L > 1e-12 * tot)),
              'n_cells_exact_zero': int(np.sum(L == 0.0)),
              'odd_band_rows_all_zero': bool(np.all(L[1::2] == 0.0))}

# 表示
hdr = ['分類', '量', 'key'] + cols
print('\n' + ' | '.join(f'{c:>14s}' for c in hdr[1:2] + cols))
for row in table:
    vals = ' | '.join(f'{v:14.4e}' if isinstance(v, float) else f'{v:>14s}' for v in row[3:])
    print(f'[{row[0]}] {row[1]:28s} {vals}')
print('\n厳密零検定（全步 max）:')
for k in exact:
    print(' ', k, {kk: ('0（厳密）' if v['exact_zero'] else ('(欠)' if v['max_all_steps'] is None else f"{v['max_all_steps']:.2e}")) for kk, v in exact[k].items()})
print('\n帳簿占有:', json.dumps(occ, ensure_ascii=False))

with open(os.path.join(BASE, 'mode_physics_table_v1.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(hdr); w.writerows(table)
with open(os.path.join(BASE, 'check_mode_physics_diff_v1.json'), 'w') as fh:
    json.dump({'window': WIN, 'inputs': meta, 'table_header': hdr, 'table': table,
               'exact_zero_tests': exact, 'bands_window_median': {k: v.tolist() for k, v in bands.items()},
               'ledger_final_occupancy': occ}, fh, indent=2, ensure_ascii=False)

# 図1: 帯別パワー
fig, ax = plt.subplots(figsize=(11, 5))
xs = np.arange(16); wdt = 0.13
for i, k in enumerate([c for c in cols if c in bands]):
    ax.bar(xs + (i - 2.5) * wdt, np.maximum(bands[k], 1e-40), wdt, label=k)
ax.set_yscale('log'); ax.set_xticks(xs); ax.set_xlabel('band k (odd = fermion, even = boson)')
ax.set_ylabel('band power (window median)'); ax.grid(alpha=.3, axis='y'); ax.legend(fontsize=8)
ax.set_title('Band powers by seed mode (N=12, T42000, delta=1e-2, window [21000,42000])')
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'fig_mode_bands.png'), dpi=150); plt.close(fig)

# 図2: 帳簿最終スナップショット（log10）
fig, axs = plt.subplots(2, 3, figsize=(15, 8)); axs = axs.ravel()
for i, k in enumerate([c for c in cols if c in ledger]):
    L = ledger[k]; im = axs[i].imshow(np.log10(np.maximum(L, 1e-40)), aspect='auto', origin='lower',
                                       vmin=-40, vmax=0, cmap='viridis')
    axs[i].set_title(f'{k}: ledger P[k,η] final (log10)', fontsize=9)
    axs[i].set_xlabel('η (winding)'); axs[i].set_ylabel('k (band)')
fig.colorbar(im, ax=axs.tolist(), shrink=.6)
fig.savefig(os.path.join(BASE, 'fig_mode_ledger.png'), dpi=140); plt.close(fig)

# 図3: 時系列（r, Q̂, 担い手, P_odd）
fig, axs = plt.subplots(2, 2, figsize=(14, 8)); axs = axs.ravel()
ts = np.arange(1, 42001)
for mode in modes:
    axs[0].plot(ts, data[mode]['rec_m_r_mean'], lw=0.9, label=mode)
    axs[1].plot(ts, data[mode]['rec_m_q_hat'], lw=0.9, label=mode)
    axs[2].semilogy(ts, np.maximum(data[mode]['m_carrier_power'], 1e-40), lw=0.9, label=mode)
    if 'rec_m_odd_power' in data[mode].files:
        axs[3].semilogy(ts, np.maximum(data[mode]['rec_m_odd_power'], 1e-40), lw=0.9, label=mode)
axs[0].axhline(0.697177928, color='k', ls=':', lw=1, label='alpha root R_alpha')
axs[0].set_title('r = sin^2 theta (mixing rate)'); axs[1].set_title('Q_hat (charge readout)')
axs[2].set_title('carrier power (time carrier)'); axs[3].set_title('P_odd (odd-band power)')
for a in axs:
    a.set_xlabel('step'); a.grid(alpha=.3); a.legend(fontsize=7)
fig.suptitle('Seed-mode dependence of readouts (N=12, T42000, delta=1e-2)', y=.995)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'fig_mode_timeseries.png'), dpi=140); plt.close(fig)
print('ALL DONE')
