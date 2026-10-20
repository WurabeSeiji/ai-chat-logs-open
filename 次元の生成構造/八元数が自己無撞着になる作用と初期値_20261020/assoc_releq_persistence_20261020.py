#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ASSOC の g₂ 相対平衡の持続性の検査（2026-10-20）。入力：results/assoc_stage2.json（assoc_stage2_20261020.py の出力）。

解析：ASSOC が G₂ 共変（g ∈ G₂ で F(gX) = g F(X)）なら、F(X) = DX（D ∈ g₂）の解 X は X(t) = exp(tD)X で厳密に持続する。
(a) 共変性：乱数の D ∈ g₂ から g = exp(D) を作り、‖F(gX) − gF(X)‖/‖F(X)‖。g が自己同型であることも ‖g(xy) − (gx)(gy)‖ で確認。
(b) 持続：見つかった解の全部を、刻み dt と dt/4 で同じ時間 T = 4000·dt 積分し、形の変化 g(t) を記録。
    厳密解 exp(tD)X との差も記録。形の崩れが刻みによらず同じ時刻に起き、g(t) が指数的に育つなら、解のまわりの不安定さ
    （丸め誤差がシード）による離脱。刻みを細かくして崩れが消えるなら積分誤差。
    指数の当てはめ：log10 g(t) が 1e-12 から 1e-2 の区間の傾き（桁/時間）。
"""
import importlib.util
import json
import math
import os

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('st2', os.path.join(HERE, 'assoc_stage2_20261020.py'))
st2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(st2)
F_assoc, mul, rk4, G2, releq_g2 = st2.F_assoc, st2.mul, st2.rk4, st2.G2, st2.releq_g2


def run(X, dt, steps):
    G0 = X @ X.T; Y = X.copy(); gs = []
    for _ in range(steps):
        Y = rk4(lambda Z: F_assoc(Z, 1), Y, dt)
        gs.append(float(np.linalg.norm(Y @ Y.T - G0) / np.linalg.norm(G0)))
    return Y, np.array(gs)


def slope(t, g):
    lg = np.log10(np.maximum(g, 1e-300)); m = (lg > -12) & (lg < -2)
    if m.sum() < 10:
        return None
    return float(np.polyfit(t[m], lg[m], 1)[0])


if __name__ == '__main__':
    rng = np.random.default_rng(20261020)
    st = json.load(open(os.path.join(HERE, 'results', 'assoc_stage2.json')))
    out = dict(equivariance=[], persistence=[])
    sols = [r for r in st['g2_relative_equilibria'] if 'X' in r]
    for r in sols[:5]:
        X = np.array(r['X'])
        D = sum(c * M for c, M in zip(rng.normal(size=14), G2)); g = expm(D)
        x, y = rng.normal(size=(2, 8))
        out['equivariance'].append(dict(N=r['N'], trial=r['trial'],
                                        automorphism_err=float(np.linalg.norm(g @ mul(x, y) - mul(g @ x, g @ y))),
                                        equivariance_rel=float(np.linalg.norm(F_assoc(X @ g.T, 1) - F_assoc(X, 1) @ g.T) / np.linalg.norm(F_assoc(X, 1)))))
        print('eq', json.dumps(out['equivariance'][-1]), flush=True)
    for r in sols:
        X = np.array(r['X']); F = F_assoc(X, 1); rg, D = releq_g2(X, F)
        dt = r['check_dt']; steps = 4000; T = dt * steps
        row = dict(N=r['N'], trial=r['trial'], releq_g2=rg, freqs=r['freqs'])
        for k, (h, n) in enumerate(((dt, steps), (dt / 4, steps * 4))):
            Y, gs = run(X, h, n); t = np.arange(1, n + 1) * h
            ex = X @ expm(T * D).T
            row[f'dt{"" if k == 0 else "_quarter"}'] = dict(shape_max=float(gs.max()), shape_end=float(gs[-1]),
                                                           diff_from_exact_end=float(np.linalg.norm(Y - ex) / np.linalg.norm(X)),
                                                           t_shape_gt_1e_6=float(t[np.argmax(gs > 1e-6)]) if (gs > 1e-6).any() else None,
                                                           growth_dec_per_time=slope(t, gs))
        out['persistence'].append(row)
        print('pe', json.dumps(row), flush=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'assoc_releq_persistence.json'), 'w'), ensure_ascii=False, indent=1)
