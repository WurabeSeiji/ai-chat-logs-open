# -*- coding: utf-8 -*-
"""段 1 の事後診断（判定ではない。走行後に追加）：
事前固定の判定 (a) は |ΔΣz²| を「初期の」ノルムで割っていたため、ブーストの固有ベクトルが毎ステップ e^{±η} 倍に伸縮する
η≠0 の場合に意味を失った（オーバーフローを含む）。ここでは各ステップの Σz² を「そのときの」ノルム v†v で割った値と、
固有ベクトルの伸縮率を記録する。"""
import json, itertools
import numpy as np
from scipy.linalg import expm
exec(open('stage1_single_point_algebra.py').read().split("CASES_THETA")[0])
res = []
for uR, uL, eta in itertools.product('ijk', 'ijk', [0.3, 0.8, 1.5]):
    theta = 2*np.pi*23/124; G = theta*Rm[uR] + 1j*eta*Lm[uL]; U = expm(G)
    mu, V = np.linalg.eig(G)
    for e in range(4):
        v = V[:, e]/np.linalg.norm(V[:, e]); worst = 0.0; steps = 0
        for step in range(1, 249):
            v = U @ v; nv = np.vdot(v, v).real
            if not np.isfinite(nv) or nv > 1e250: break
            worst = max(worst, abs(v @ v)/nv); steps = step
        res.append(dict(uR=uR, uL=uL, eta=eta, eig_real=float(mu[e].real), steps_before_overflow=steps, max_abs_sumz2_over_current_norm=float(worst)))
json.dump(res, open('stage1_posthoc_relative_closure.json', 'w'), indent=1)
print('max |Σz²|/(current norm) over all:', max(r['max_abs_sumz2_over_current_norm'] for r in res))
print('growth per step = exp(Re eigenvalue); Re eigenvalues seen:', sorted({round(r['eig_real'],6) for r in res}))
print('min steps before overflow cut:', min(r['steps_before_overflow'] for r in res))
