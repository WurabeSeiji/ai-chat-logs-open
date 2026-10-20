#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""高対称な四元数と八元数のインフレーション実験を τ = 5000 まで走らせる版（2026-10-20）。

著者の指示：テスト回数が少なすぎる。最低 τ = 5000 程度までは走らせて状態を確認する。
symmetric_inflation_20261020.py（20000 歩、τ ≈ 130〜290 で停止）の配置・作用・種（そのまま／自己同型で回す）・刻み・射影つき RK4 を
そのまま読み込み、変えるのは (1) 歩数を τ_end / dt にする、(2) 記録の間隔を全体で約 5000 点にする、(3) 組ごとに並列に走らせる、だけ。
組：四元数 CROSS8, T24, O48, I120, TWO_C × C1, C2 × κ = ±1 × (そのまま、回した) の 40 組（引数 quat）。
八元数は引数 oct で別に走らせる（ASSOC の高速版を使う場合はその検証を済ませてから）。
出力：results_long/series_<組>.npz（t, g, v）、results_long/summary_<組>.json、results_long/summary_all.json。
"""
import importlib.util
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('si', os.path.join(HERE, 'symmetric_inflation_20261020.py'))
si = importlib.util.module_from_spec(spec); spec.loader.exec_module(si)
OUT = os.path.join(HERE, 'results_long')
TAU_END = 5000.0
NREC = 5000


def quat_jobs():
    rng = np.random.default_rng(si.SEED)
    qc = si.Q.configs(np.random.default_rng(si.SEED))
    jobs = []
    for name in ('CROSS8', 'T24', 'O48', 'I120', 'TWO_C'):
        X = qc[name]; N = len(X)
        u = rng.normal(size=4); u /= np.linalg.norm(u)                 # symmetric_inflation_20261020.py と同じ乱数の順
        Xrot = np.array([si.Q.mul(si.Q.mul(u, x), u * si.Q.CONJ) for x in X])
        for aname in ('C1', 'C2'):
            for kappa in (1, -1):
                for tag, X0 in (('asis', X), ('rot', Xrot)):
                    jobs.append(('quat', name, aname, kappa, tag, X0))
    return jobs


def run_job(job):
    kind, name, aname, kappa, tag, X0 = job
    N = len(X0)
    if kind == 'quat':
        Fa = si.Q.F_C1 if aname == 'C1' else si.Q.F_C2
        Xg = np.random.default_rng(si.SEED + N).normal(size=(N, 4))
        extra = si.q_extra
    sref = float(max(np.linalg.norm(Fa(Xg, 1), axis=1) / np.linalg.norm(Xg, axis=1)))
    F = lambda Y: Fa(Y, kappa)
    dt = 2 * math.pi / (200 * sref); steps = int(math.ceil(TAU_END / dt)); every = max(1, steps // NREC)
    r0 = np.linalg.norm(X0, axis=1); sx = float(r0.sum()); G0 = X0 @ X0.T; nG = np.linalg.norm(G0)
    X = X0.copy(); T, G, V = [], [], []
    for t in range(1, steps + 1):
        X = si.rk4r(F, X, dt, r0)
        if t % every == 0 or t == steps:
            T.append(t * dt); G.append(float(np.linalg.norm(X @ X.T - G0) / nG)); V.append(float(np.linalg.norm(F(X)) / sx))
    label = f'{name}_{aname}_k{kappa:+d}_{tag}'
    tt = np.array(T); g = np.maximum(np.array(G), 1e-300); lg = np.log10(g)
    m = (lg > -13) & (lg < -3)
    slope = float(np.polyfit(tt[m], lg[m], 1)[0]) if m.sum() >= 10 else None
    i3 = np.where(lg > -3)[0]
    # 後半の育ち方：log g と log t の傾き（べき）
    h = len(tt) // 2; mm = g[h:] > 1e-290
    power = float(np.polyfit(np.log10(tt[h:][mm]), lg[h:][mm], 1)[0]) if mm.sum() >= 10 else None
    out = dict(label=label, N=N, dt=dt, steps=steps, T=steps * dt, F0=float(np.linalg.norm(F(X0)) / sx),
               g_at_tau={str(tau): float(g[np.argmin(np.abs(tt - tau))]) for tau in (100, 300, 1000, 3000, 5000)},
               g_max=float(g.max()), growth_dec_per_time=slope, power_late=power,
               t_reach_1e_3=float(tt[i3[0]]) if len(i3) else None, speed_end=float(V[-1]), **extra(X))
    np.savez_compressed(os.path.join(OUT, f'series_{label}.npz'), t=tt, g=np.array(G), v=np.array(V))
    json.dump(out, open(os.path.join(OUT, f'summary_{label}.json'), 'w'), ensure_ascii=False, indent=1)
    print(json.dumps({k: (float('%.4g' % v) if isinstance(v, float) else v) for k, v in out.items()}, ensure_ascii=False), flush=True)
    return out


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else 'quat'
    jobs = quat_jobs() if which == 'quat' else []
    jobs.sort(key=lambda j: -len(j[5]))                                   # 重い組から
    with Pool(8) as p:
        res = p.map(run_job, jobs, chunksize=1)
    json.dump(res, open(os.path.join(OUT, f'summary_all_{which}.json'), 'w'), ensure_ascii=False, indent=1)
