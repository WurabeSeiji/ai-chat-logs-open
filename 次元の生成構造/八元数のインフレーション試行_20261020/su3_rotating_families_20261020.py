#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 3 段：ℓ を選んだ（SU(3) 共変の）過去の型の作用で、正軸体 16 と E₈ の根 240 の回転する自己無撞着状態を調べる（2026-10-20）。

作用：K_ef = A_ef⟨ℓ x_e, x_f⟩（ℓ = e₄、全員対全員、Q ← exp(ΔτK)Q）。振幅込み（amp）と位相のみ（ph）。
  一つの複素平面 span(1, ℓ) では過去の複素数の力学に一致する。
測る量：relEq = ‖Ω + Ωᵀ‖/‖Ω‖（Ω = Q⁺KQ、0 なら回転する相対平衡）、Ω の回転の振動数（複素平面ごと）、
  丸め誤差から走らせたときの形の変化 ‖QQᵀ − Q₀Q₀ᵀ‖/‖Q₀Q₀ᵀ‖（指数的に育てばインフレーションの資格あり）。
正軸体：{±1}・{±ℓ} の半径 a, b、残り 12 個（ℓ で組になる三つの複素平面）の半径 c を動かす。
E₈：span(1, ℓ) への射影の大きさで分けた類（4・48・128・60 個）の半径を動かす。
"""
import itertools
import json
import math
import os

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')).read().split("results = {")[0])
ELL = E[4]
LL = np.array([omul(ELL, E[j]) for j in range(8)]).T


def K_of(Q, ph=False):
    U = Q / np.linalg.norm(Q, axis=1, keepdims=True) if ph else Q
    M = len(Q); A = np.ones((M, M)) - np.eye(M)
    return A * ((U @ LL.T) @ U.T)


def releq(Q, ph=False):
    KQ = K_of(Q, ph) @ Q; Om = np.linalg.pinv(Q) @ KQ
    w = np.sort(np.abs(np.linalg.eigvals(Om).imag))[::2]
    return float(np.linalg.norm(Om + Om.T) / np.linalg.norm(Om)), [float(x) for x in w]


def evolve_shape(Q, ph, steps, every):
    K0 = K_of(Q, ph); sig = float(np.max(np.abs(np.linalg.eigvals(K0)))); dt = 2 * math.pi / (200 * sig)
    G0 = Q @ Q.T; out = []
    for t in range(steps):
        Q = expm(dt * K_of(Q, ph)) @ Q
        if t % every == every - 1:
            out.append(float(np.linalg.norm(Q @ Q.T - G0) / np.linalg.norm(G0)))
    return out


def ortho(a, b, c):
    return np.array([s * E[0] * a for s in (1, -1)] + [s * ELL * b for s in (1, -1)]
                    + [s * E[i] * c for i in (1, 2, 3, 5, 6, 7) for s in (1, -1)])


def e8():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si, sj in itertools.product((1, -1), repeat=2):
            v = np.zeros(8); v[i] = si; v[j] = sj; R.append(v)
    for signs in itertools.product((1, -1), repeat=8):
        if np.prod(signs) > 0:
            R.append(np.array(signs) / 2)
    return np.array(R) / math.sqrt(2)


if __name__ == '__main__':
    res = {'ortho_family': [], 'ortho_evolution': {}, 'e8_family': [], 'e8_evolution': {}}
    for a, b, c in [(1, 1, 1), (1, 2, 1), (2, 1, 1), (1, 1, 2), (1, 1, 0.5), (1.3, 0.7, 1.1)]:
        r, w = releq(ortho(a, b, c)); res['ortho_family'].append(dict(a=a, b=b, c=c, relEq=r, freqs=w))
    for ph in (False, True):
        for c in (0.5, 1.0, math.sqrt(2), 2.0):
            res['ortho_evolution'][f"{'ph' if ph else 'amp'}_c{c:.3f}"] = evolve_shape(ortho(1, 1, c), ph, 8000, 500)
    R = e8(); par = np.round(np.sqrt(R[:, 0] ** 2 + R[:, 4] ** 2), 4); cls = np.unique(par)
    res['e8_classes'] = {str(float(c)): int((par == c).sum()) for c in cls}
    for k in [-1] + list(range(len(cls))):
        Q = R.copy()
        if k >= 0:
            Q[par == cls[k]] *= 1.3
        r, w = releq(Q); res['e8_family'].append(dict(scaled_class=(None if k < 0 else float(cls[k])), relEq=r, freqs=w))
    for ph in (False, True):
        res['e8_evolution']['ph' if ph else 'amp'] = evolve_shape(R.copy(), ph, 3000, 300)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'results', 'su3_rotating_families.json'), 'w'), ensure_ascii=False, indent=1)
    print(json.dumps({k: res[k] for k in ('ortho_family', 'e8_classes', 'e8_family')}, ensure_ascii=False))
