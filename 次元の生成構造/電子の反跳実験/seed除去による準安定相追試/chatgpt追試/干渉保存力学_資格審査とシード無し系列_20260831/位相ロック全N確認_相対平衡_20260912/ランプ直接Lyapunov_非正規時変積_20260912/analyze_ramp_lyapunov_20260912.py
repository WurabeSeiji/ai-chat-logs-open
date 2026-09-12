#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""追試A: インフレーション実ランプの床値超過（~18-31%）が、時変・非可換ヤコビアンの
時間順序積の増幅で説明できるかの直接検証（読出しのみ・物理無変更・新規走行なし）。

木原指示（2026-09-12）: 論文6は「実ランプが床固有値予測 2log10|μ_max| を~18-31%上回る分は
時変非正規Jの時間順序積が担う」を物理的解釈として書いた。直接Lyapunov照合が未実施だったので
既存軌道データで確定させる。

方法（N=6,8,10、既存 den=N 500step 軌道）:
  J(z) = one_step の 2M×2M 実ヤコビアン（有限差分、analyze_inflation_linear_instability と同一構成）。
  同一の微小 δ を床平面直交成分 H⊥ で追跡する:
   (a) 床凍結: J0=J(z_0) を固定して δ を発展 → 傾き slope_frozen（床固有値予測 2log10|μ_max| に一致するはず）。
   (b) 時変積: 実軌道 φ_t 上の J(z_t) を順に掛けた時間順序積 ∏J(z_t) で δ を発展 → 傾き slope_tv。
   (c) 実測: 実軌道そのものの H⊥/H ランプ傾き slope_meas。
  判定: slope_tv ≈ slope_meas かつ slope_tv > slope_frozen なら「超過分＝時変非正規積の増幅」を確認。
        slope_tv ≈ slope_frozen なら時変効果は無く別要因。
正本 run_N3_N40_stage123_v1.py（SHA256 照合）から adjacency/one_step を import（忠実）。
"""
import csv
import hashlib
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.abspath(os.path.join(HERE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py'))
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
STATES = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907', 'full_N3_N40_sweep', 'states'))

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
adjacency, one_step = g['adjacency'], g['one_step']

rv = lambda z: np.concatenate([z.real, z.imag])
cp = lambda v, M: v[:M] + 1j * v[M:]


def jac(z, A, den, M, eps=1e-8):
    base = rv(one_step(z, A, den)); n = 2 * M; J = np.zeros((n, n))
    for k in range(n):
        v = rv(z).copy(); v[k] += eps
        J[:, k] = (rv(one_step(cp(v, M), A, den)) - base) / eps
    return J


def fit(x, y):
    return float(np.polyfit(x, y, 1)[0])


def main():
    rows = []
    for N in (6, 8, 10):
        A = adjacency(N); M = N * (N - 1) // 2; den = N; n = 2 * M
        S = np.asarray(np.load(os.path.join(STATES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'], np.complex128)
        z0 = S[0]
        p = z0.real / np.linalg.norm(z0.real); q = z0.imag - np.dot(z0.imag, p) * p; q /= np.linalg.norm(q)
        P = np.eye(M) - np.outer(p, p) - np.outer(q, q)
        perp2 = lambda dz: (lambda dr: np.vdot(P @ dr, P @ dr).real)(dz[:M] + 1j * dz[M:])
        # (c) 実測ランプ
        hp = np.array([(lambda zp: np.vdot(zp, zp).real / np.vdot(z, z).real)(z - p * np.dot(p, z) - q * np.dot(q, z)) for z in S])
        lo = int(np.argmax(hp > 1e-18)); hi = int(np.argmax(hp > 1e-8))
        w0, w1 = max(lo, 12), hi
        sl_meas = fit(np.arange(w0, w1), np.log10(hp[w0:w1]))
        # 床凍結 J0
        J0 = jac(z0, A, den, M); amax = float(np.max(np.abs(np.linalg.eigvals(J0))))
        slope_floor_pred = 2 * math.log10(amax)
        rng = np.random.default_rng(0); d0 = 1e-6 * rng.standard_normal(n); d0 /= np.linalg.norm(d0)
        T = hi + 10
        # (a) 床凍結 J0 を固定して δ を発展（N=6,8,10 の T では overflow しない）
        ha = np.empty(T); dd = d0.copy()
        for t in range(T):
            ha[t] = perp2(dd); dd = J0 @ dd
        sl_frozen = fit(np.arange(w0, w1), np.log10(ha[w0:w1] / ha[w0]))
        # (b) 時変積: 実軌道 φ_t 上の J(z_t) を順に掛ける
        hb = np.empty(T); dd = d0.copy()
        for t in range(T):
            hb[t] = perp2(dd)
            Jt = jac(S[min(t, len(S) - 1)], A, den, M)
            dd = Jt @ dd
        sl_tv = fit(np.arange(w0, w1), np.log10(hb[w0:w1] / hb[w0]))
        ratio_meas = sl_meas / slope_floor_pred
        ratio_tv = sl_tv / slope_floor_pred
        rows.append(dict(N=N, window=f'[{w0},{w1}]', slope_floor_pred=slope_floor_pred,
                         slope_frozen_dev=sl_frozen, slope_tv_prod=sl_tv, slope_meas=sl_meas,
                         tv_over_floor=ratio_tv, meas_over_floor=ratio_meas,
                         tv_vs_meas=sl_tv / sl_meas))
        print(f'N={N} [{w0},{w1}]: 床予測2log10|μ|={slope_floor_pred:.4f}  床凍結δ発展={sl_frozen:.4f}  '
              f'時変積={sl_tv:.4f}  実測={sl_meas:.4f}  | 時変/床={ratio_tv:.3f} 実測/床={ratio_meas:.3f} '
              f'時変/実測={sl_tv/sl_meas:.3f}')
    with open(os.path.join(HERE, 'results', 'ramp_lyapunov_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    print('\n判定: 時変積の傾きが実測ランプに一致し床凍結を上回れば、超過分は時変非正規Jの時間順序積の増幅。')


if __name__ == '__main__':
    main()
