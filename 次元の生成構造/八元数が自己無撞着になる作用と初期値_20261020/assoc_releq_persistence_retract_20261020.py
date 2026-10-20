#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""g₂ 相対平衡の持続性の再確認：各歩のあとで半径を戻す射影（retraction）を入れた版（2026-10-20）。
assoc_releq_persistence_20261020.py をそのまま読み込み、積分の一歩を「RK4 → 各 |x_e| を初期値に戻す」に置き換える一点だけを変える
（claude.ai `定式化の突き合わせと走行の直し方_20261020.md` §3(2) の指摘。厳密解は |x_e| を保つので射影は正当）。
出力：results/assoc_releq_persistence_retract.json（同じ形式）と、射影なしの結果との離れる速さの比較。
"""
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('per', os.path.join(HERE, 'assoc_releq_persistence_20261020.py'))
per = importlib.util.module_from_spec(spec); spec.loader.exec_module(per)


def run_retract(X, dt, steps):
    r0 = np.linalg.norm(X, axis=1); G0 = X @ X.T; Y = X.copy(); gs = []
    for _ in range(steps):
        Y = per.rk4(lambda Z: per.F_assoc(Z, 1), Y, dt)
        Y = Y * (r0 / np.linalg.norm(Y, axis=1))[:, None]
        gs.append(float(np.linalg.norm(Y @ Y.T - G0) / np.linalg.norm(G0)))
    return Y, np.array(gs)


if __name__ == '__main__':
    per.run = run_retract
    from scipy.linalg import expm
    st = json.load(open(os.path.join(HERE, 'results', 'assoc_stage2.json')))
    old = {(r['N'], r['trial']): r for r in json.load(open(os.path.join(HERE, 'results', 'assoc_releq_persistence.json')))['persistence']}
    out = []
    for r in [r for r in st['g2_relative_equilibria'] if 'X' in r]:
        X = np.array(r['X']); F = per.F_assoc(X, 1); rg, D = per.releq_g2(X, F)
        dt = r['check_dt']; steps = 4000; T = dt * steps
        row = dict(N=r['N'], trial=r['trial'], releq_g2=rg, freqs=r['freqs'])
        for k, (h, n) in enumerate(((dt, steps), (dt / 4, steps * 4))):
            Y, gs = run_retract(X, h, n); t = np.arange(1, n + 1) * h
            key = 'dt' if k == 0 else 'dt_quarter'
            row[key] = dict(shape_max=float(gs.max()), norm_drift=float(np.max(np.abs(np.linalg.norm(Y, axis=1) - 1.0))),
                            diff_from_exact_end=float(np.linalg.norm(Y - X @ expm(T * D).T) / np.linalg.norm(X)),
                            t_shape_gt_1e_6=float(t[np.argmax(gs > 1e-6)]) if (gs > 1e-6).any() else None,
                            growth_dec_per_time=per.slope(t, gs),
                            growth_without_retract=old[(r['N'], r['trial'])][key]['growth_dec_per_time'],
                            shape_max_without_retract=old[(r['N'], r['trial'])][key]['shape_max'])
        out.append(row)
        print(json.dumps(row), flush=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'assoc_releq_persistence_retract.json'), 'w'), ensure_ascii=False, indent=1)
