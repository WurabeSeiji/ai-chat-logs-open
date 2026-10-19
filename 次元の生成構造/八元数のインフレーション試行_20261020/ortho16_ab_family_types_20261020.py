#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一つの作用（振幅込み L_amp、K = Q Lᵀ Qᵀ）で、相対平衡・周期的・準周期的の自己無撞着状態を一つの族で作る（2026-10-20）。

配置：正軸体 16 の ±1 を半径 a、±ℓ（ℓ = e₄）を半径 b にし、残り ±e_j は半径 1。S = diag(2a², …, 2b², …, 2)。
解析：ℓ は (1, ℓ) を結ぶので、この平面の振動数は 2ab、残りの三つの平面は 2。
  a = b → [L,S] = 0 → 相対平衡。ab が有理数 → 周期的（例 a = 2, b = 1 で振動数 4 と 2、周期 π）。
  ab が無理数 → 準周期的（例 a = √2, b = 1 で 2√2 と 2）。
確認：厳密解 Q(t) = Q₀ exp(−tLS)（=刻みの一歩 exp(ΔτK) と同じ）で、形（グラム行列）と状態の戻りを測る。
"""
import json
import math
import os

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'octonion_inflation_trials_20261020.py')).read().split("if __name__")[0]
     .replace("HERE = os.path.dirname(os.path.abspath(__file__))", f"HERE = {HERE!r}"))
L = LL


def family(a, b):
    X = np.array([s * E[k] for k in range(8) for s in (1, -1)])
    X[0:2] *= a; X[8:10] *= b
    return X


def run(a, b, T, n=4000):
    Q0 = family(a, b); S = Q0.T @ Q0; G0 = Q0 @ Q0.T
    K = kernel(Q0, 'L_amp'); dt = T / n
    Q = Q0.copy(); U = expm(dt * K)          # K は保存されない一般の場合に備え毎回作り直す
    shape, state = [], []
    for _ in range(n):
        Q = expm(dt * kernel(Q, 'L_amp')) @ Q
        shape.append(float(np.linalg.norm(Q @ Q.T - G0) / np.linalg.norm(G0)))
        state.append(float(np.linalg.norm(Q - Q0) / np.linalg.norm(Q0)))
    w, V = np.linalg.eigh(S); Sh = V @ np.diag(np.sqrt(w)) @ V.T
    fr = sorted(set(np.round(np.abs(np.linalg.eigvals(Sh @ L @ Sh)), 9).tolist()))
    exact = Q0 @ expm(-T * L @ S)
    return dict(a=a, b=b, T=T, comm_L_S=float(np.linalg.norm(L @ S - S @ L)), freqs=fr,
                shape_max=max(shape), shape_at_T=shape[-1], state_at_T=state[-1], state_min_after_half=min(state[n // 2:]),
                stepped_vs_exact=float(np.linalg.norm(Q - exact) / np.linalg.norm(Q0)))


if __name__ == '__main__':
    res = [run(1.0, 1.0, math.pi), run(2.0, 1.0, math.pi), run(math.sqrt(2), 1.0, 20 * math.pi, n=20000)]
    for r in res:
        print(json.dumps(r, ensure_ascii=False), flush=True)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'results', 'ortho16_ab_family_types.json'), 'w'), ensure_ascii=False, indent=1)
