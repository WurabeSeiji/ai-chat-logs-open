#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SVD 解析（仕様書 v1.1 §3–§8。判定条件は凍結済み・変更禁止）。

- 軌道行列 X=[Z(0),...,Z(T)] ∈ C^{M×(T+1)}、前処理なし（U(1) 位相不変性により
  位相 align 不要、H 保存により列ノルム一定）。
- 全区間と時間窓 W1=[0,50], W2=[50,200], W3=[200,500]（境界含む・凍結）。
- 観測量: 特異値スペクトル、E_r、r90/r99/r999、R_eff、principal angles
  （個数は min(2, r)。r=2 主解析と r=r99 副解析）。
- 事前登録判定（第一対象 = 主刻み AB のみに機械適用）:
    (a) 低ランク再凝縮: r99(W3) < r99(W1) かつ r99(W3) <= 6
    (b) 多モード拡散:  r99(W3) >= r99(W1) かつ R_eff(W3) >= min(M,|W3|)/2
    (c) 中間: いずれにも該当しない
- 正値対照: L6_den6_B で r99=1 を期待（代数 rank>1 は不合格としない）。
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PHASEA = os.path.abspath(os.path.join(HERE, '..', 'PhaseA_実装一式_v1_20260915'))
sys.path.insert(0, PHASEA)
from generator_phaseA import analysis_basis  # noqa: E402

STATES_DIR = os.path.join(HERE, 'results', 'states')
SV_DIR = os.path.join(HERE, 'results', 'singular_values')
os.makedirs(SV_DIR, exist_ok=True)

WINDOWS = {'full': (0, 500), 'W1': (0, 50), 'W2': (50, 200), 'W3': (200, 500)}
PRIMARY = ['L6_den6_AB', 'L12_den12_AB']
CONTROL = ['L6_den40_AB', 'L12_den40_AB']
POSITIVE = ['L6_den6_B']


def rank_indices(s):
    e = s ** 2
    tot = float(np.sum(e))
    cum = np.cumsum(e) / tot
    def r_at(th):
        return int(np.searchsorted(cum, th) + 1)
    reff = float(tot ** 2 / np.sum(e ** 2))
    return {'r90': r_at(0.90), 'r99': r_at(0.99), 'r999': r_at(0.999),
            'R_eff': reff, 'E_2': float(cum[1]) if len(cum) > 1 else 1.0,
            'E_4': float(cum[3]) if len(cum) > 3 else 1.0}


def principal_angles(Qb, Ur):
    s = np.linalg.svd(Qb.conj().T @ Ur, compute_uv=False)
    s = np.clip(s, 0.0, 1.0)
    n = min(Qb.shape[1], Ur.shape[1])
    return [float(np.arccos(x)) for x in s[:n]]


def analyze(run_id):
    d = np.load(os.path.join(STATES_DIR, f'{run_id}_states.npz'))
    Z = d['Z']            # (T+1, M)
    L = int(d['L'])
    M = Z.shape[1]
    B = analysis_basis(L)
    Qb, _ = np.linalg.qr(B)
    out = {'run_id': run_id, 'L': L, 'M': M}
    for wname, (t0, t1) in WINDOWS.items():
        X = Z[t0:t1 + 1].T                        # M × |W|
        U, s, _ = np.linalg.svd(X, full_matrices=False)
        idx = rank_indices(s)
        r99 = idx['r99']
        pa2 = principal_angles(Qb, U[:, :2])       # r=2 主解析（角 2 個）
        pa_r = principal_angles(Qb, U[:, :r99])    # r=r99 副解析（角 min(2,r99) 個）
        np.savetxt(os.path.join(SV_DIR, f'{run_id}_{wname}_sv.csv'),
                   np.column_stack([np.arange(1, len(s) + 1), s]),
                   delimiter=',', header='k,sigma_k', comments='')
        out[wname] = {**idx, 'ncols': X.shape[1],
                      'sigma_head': [float(x) for x in s[:8]],
                      'pa_r2_rad': pa2, 'pa_r99_rad': pa_r}
    return out


def main():
    results = {}
    for run_id in PRIMARY + CONTROL + POSITIVE:
        results[run_id] = analyze(run_id)
        r = results[run_id]
        print(f"{run_id:>14}: r99 full={r['full']['r99']} W1={r['W1']['r99']} "
              f"W2={r['W2']['r99']} W3={r['W3']['r99']}  "
              f"R_eff(W3)={r['W3']['R_eff']:.2f}  "
              f"pa_r2(full)=[{r['full']['pa_r2_rad'][0]:.3f},{r['full']['pa_r2_rad'][1]:.3f}]rad", flush=True)

    # 事前登録判定の機械適用（第一対象のみ）
    verdicts = {}
    for run_id in PRIMARY:
        r = results[run_id]
        M = r['M']
        w3len = WINDOWS['W3'][1] - WINDOWS['W3'][0] + 1
        a = (r['W3']['r99'] < r['W1']['r99']) and (r['W3']['r99'] <= 6)
        b = (r['W3']['r99'] >= r['W1']['r99']) and (r['W3']['R_eff'] >= 0.5 * min(M, w3len))
        verdicts[run_id] = 'a_low_rank_recondensation' if a else (
            'b_multimode_diffusion' if b else 'c_intermediate')
        print(f'判定 {run_id}: {verdicts[run_id]}')

    # 正値対照
    pos = results['L6_den6_B']
    pos_ok = pos['full']['r99'] == 1
    print(f"正値対照 L6_den6_B: r99(full)={pos['full']['r99']} "
          f"（期待 1）→ {'PASS' if pos_ok else 'FAIL'}")

    meta = {'windows': WINDOWS,
            'criteria_frozen': '(a) r99(W3)<r99(W1) and r99(W3)<=6; '
                               '(b) r99(W3)>=r99(W1) and R_eff(W3)>=min(M,|W3|)/2; (c) otherwise',
            'principal_angle_count': 'min(2, r)',
            'positive_control': {'run': 'L6_den6_B', 'expect': 'r99=1', 'pass': pos_ok}}
    with open(os.path.join(HERE, 'results', 'svd_summary.json'), 'w', encoding='utf-8') as f:
        json.dump({'meta': meta, 'results': results, 'verdicts': verdicts}, f, indent=1)
    print('SVD ANALYSIS DONE')


if __name__ == '__main__':
    main()
