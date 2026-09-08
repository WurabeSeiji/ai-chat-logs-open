#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""テスト1（修正1の検証・読み出しのみ）: den=124 統一時計での 4 検査。
(1) 蹴りの den 非依存性: f(1) を den=N 版（results/）と比較（予測: ほぼ同値）
(2) 124歩再帰: A124(t)=|<z_t,z_{t+124}>|/(‖z_t‖‖z_{t+124}‖) の max/mean（t≤376）。奇数Nは248歩再帰も（t≤252）
(3) 閉塞 C2(t) の型: 全期間max と「後半平均/前半max」で 開きっぱなし vs 往復 を判別。周期124のFFTピークも記録
(4) Control（理論床 den=124、既存 ../対称親v2_500step走行_20260906/results）との終値比較
出力: test1_den124_sync_results.csv"""
import csv
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
R124 = os.path.join(BASE, 'results_den124')
RN = os.path.join(BASE, 'results')
CTRL = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'results')


def f1_of(path, den):
    for r in csv.DictReader(open(os.path.join(path, 'timeseries_64bit_with124_N3_N40.csv'))):
        pass
    out = {}
    for r in csv.DictReader(open(os.path.join(path, 'timeseries_64bit_with124_N3_N40.csv'))):
        if int(r['denominator']) == den if den else r['series'] == 'N':
            pass
    return out


rows = []
for N in range(3, 17):
    d = np.load(os.path.join(R124, f'hm_N{N}_den_124_states_500.npz'))
    S = np.asarray(d['Z'], complex)
    dn = np.load(os.path.join(RN, f'hm_N{N}_den_{N}_states_500.npz'))
    Sn = np.asarray(dn['Z'], complex)
    # (1) f(1): 初期平面への垂直比を直接計算（正本 metrics と同一式）
    def f_of(Sarr, t):
        z0 = Sarr[0]
        p = z0.real.copy(); p /= np.linalg.norm(p)
        q = z0.imag.copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
        z = Sarr[t]
        zp = z - p * np.dot(p, z) - q * np.dot(q, z)
        return float(np.vdot(zp, zp).real / np.vdot(z, z).real)
    f1_124 = f_of(S, 1); f1_N = f_of(Sn, 1)
    # (2) 再帰
    def recur(Sarr, P):
        vals = [abs(np.vdot(Sarr[t], Sarr[t + P])) / (np.linalg.norm(Sarr[t]) * np.linalg.norm(Sarr[t + P]))
                for t in range(0, Sarr.shape[0] - P)]
        return float(np.max(vals)), float(np.mean(vals)), float(vals[-1])
    r124max, r124mean, r124last = recur(S, 124)
    r248max, r248mean, r248last = recur(S, 248)
    # (3) 閉塞
    h = np.einsum('ij,ij->i', S.conj(), S).real
    c2 = np.abs(np.einsum('ij,ij->i', S, S)) / h
    half = c2[250:]
    # 周期124成分
    x = c2 - c2.mean()
    F = np.abs(np.fft.rfft(x))
    freqs = np.fft.rfftfreq(len(x))
    k124 = np.argmin(np.abs(freqs - 1 / 124))
    p124 = float(F[k124] / (F.sum() + 1e-300))
    rows.append(dict(N=N, f1_den124=f1_124, f1_denN=f1_N, f1_ratio=f1_124 / f1_N if f1_N > 0 else float('inf'),
                     rec124_max=r124max, rec124_mean=r124mean, rec124_last=r124last,
                     rec248_max=r248max, rec248_last=r248last,
                     c2_max=float(c2.max()), c2_late_mean=float(half.mean()), c2_final=float(c2[-1]),
                     c2_period124_weight=p124,
                     f_final=f_of(S, 500)))
with open(os.path.join(BASE, 'test1_den124_sync_results.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('N | f1 124/N(比) | 再帰124 max/last | 再帰248 last | C2 max/late/final | 124周期重み | f(500)')
for r in rows:
    print(f"{r['N']:2d} | {r['f1_den124']:.2e}/{r['f1_denN']:.2e}({r['f1_ratio']:.2f}) | "
          f"{r['rec124_max']:.4f}/{r['rec124_last']:.4f} | {r['rec248_last']:.4f} | "
          f"{r['c2_max']:.2e}/{r['c2_late_mean']:.2e}/{r['c2_final']:.2e} | {r['c2_period124_weight']:.3f} | {r['f_final']:.3f}")
