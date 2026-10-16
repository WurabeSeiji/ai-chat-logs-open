# -*- coding: utf-8 -*-
"""段 2 の事後診断（判定ではない。走行後に追加）：ブーストの数値の条件。
ブーストは次数 n の空間でノルムを最大 e^{nη} 程度に伸ばすので、倍精度の丸め誤差は（伸びたノルム）に比例する。
B（本当の二乗）の相対誤差、可換性の差を伸びたノルムで割った値、右の回転 124 回の差の内訳を記録する。"""
import json, numpy as np
from s3_harmonics import *
rng = np.random.default_rng(20261016); N = 18
C0 = [rng.normal(size=(n+1, n+1)) + 1j*rng.normal(size=(n+1, n+1)) for n in range(N+1)]
rows = []
for eta in (0.1, 0.3, 0.8):
    for Nmax in (6, 12, 18):
        Cs = C0[:Nmax+1]; b0 = true_square(Cs); h0 = norm_h(Cs)
        Cb = translate_left(Cs, qexp(1j*eta, 'j')); growth = norm_h(Cb)/h0
        bR = qexp(2*np.pi*23/124, 'i')
        A1 = rotate_right(translate_left(Cs, qexp(1j*eta, 'i')), bR); A2 = translate_left(rotate_right(Cs, bR), qexp(1j*eta, 'i'))
        scale = max(np.max(np.abs(a)) for a in A1)
        rows.append(dict(eta=eta, N=Nmax, B_rel_err=float(abs(true_square(Cb)-b0)/abs(b0)), norm_growth=float(growth),
                         predicted_eps_scale=float(2.2e-16*growth*h0/abs(b0)),
                         commute_abs=float(max(np.max(np.abs(a-b)) for a, b in zip(A1, A2))), commute_rel_to_max_coef=float(max(np.max(np.abs(a-b)) for a, b in zip(A1, A2))/scale)))
C = C0; b = qexp(2*np.pi*23/124, 'i'); per_deg = []
for k in range(124): C = rotate_right(C, b)
per_deg = [float(np.max(np.abs(C[n]-C0[n]))/np.max(np.abs(C0[n]))) for n in range(N+1)]
out = dict(boost_rows=rows, right124_rel_err_per_degree=per_deg)
json.dump(out, open('stage2_posthoc_boost_conditioning.json', 'w'), indent=1)
for r in rows: print(r)
print('right rotation 124 steps, relative error per degree (max):', max(per_deg))
