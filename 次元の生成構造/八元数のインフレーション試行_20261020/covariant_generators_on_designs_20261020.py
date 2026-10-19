#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G₂ で共変・大きさを保つ・3 次の作用の族を、高対称な球面デザイン（正軸体 16、E₈ の根 240）に当てる（2026-10-20）。

作用（どれも八元数の積・左右の掛け算と内積だけで作る。参照の向き ℓ を使わない）：
  F_A : 要素ごとの SO(8) 回転 dx_e = Σ_f (x_f|x_e|² − x_e⟨x_e, x_f⟩)（各 |x_e| を保つ。SO(8) で共変）
  F_B : 要素ごとの左の掛け算 dx_e = u_e x_e、u_e = Σ_f Im(x_f x̄_e)（各 |x_e| を保つ。G₂ で共変）
  F_C : 要素どうしの振幅の移し替え dx_e = Σ_f s_ef x_f、s_ef = Re(x_e)Im(x_f) + Re(x_f)Im(x_e)
        （s は (e,f) について対称な純虚数なので、ブロック L(s_ef) の全体は反対称で、全ノルムを保つ。G₂ で共変）
  事実（表現論）：SO(n)（n ≥ 3）にも G₂ にも、二つのベクトルから作る反対称で不変なスカラーは無い。
  したがって過去の型（反対称スカラー K_ef を番号の側に掛ける）は、参照の向き ℓ なしでは作れない。
配置：ORTHO16 = {±1, ±e₁…±e₇}、E8_240 = E₈ の根（実軸を第 0 成分に取った標準座標、/√2 で単位長）。
出力：各作用の |F|、接線成分、全ノルムの変化、|Re| の類ごとの動径の変化率。
"""
import itertools
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')).read().split("results = {")[0])


def Im(x):
    y = x.copy(); y[0] = 0.0; return y


def F_A(X):
    F = np.zeros_like(X)
    for e in range(len(X)):
        for f in range(len(X)):
            if f != e:
                F[e] += X[f] * (X[e] @ X[e]) - X[e] * (X[e] @ X[f])
    return F


def F_B(X):
    F = np.zeros_like(X)
    for e in range(len(X)):
        u = np.zeros(8)
        for f in range(len(X)):
            if f != e:
                u += Im(omul(X[f], oconj(X[e])))
        F[e] = omul(u, X[e])
    return F


def F_C(X):
    F = np.zeros_like(X)
    for e in range(len(X)):
        for f in range(len(X)):
            if f != e:
                F[e] += omul(X[e, 0] * Im(X[f]) + X[f, 0] * Im(X[e]), X[f])
    return F


def ortho16():
    return np.array([s * E[a] for a in range(8) for s in (1, -1)])


def e8_240():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si, sj in itertools.product((1, -1), repeat=2):
            v = np.zeros(8); v[i] = si; v[j] = sj; R.append(v)
    for signs in itertools.product((1, -1), repeat=8):
        if np.prod(signs) > 0:
            R.append(np.array(signs) / 2)
    return np.array(R) / math.sqrt(2)


def report(name, X, scale_by_class=None):
    out = {}
    cls = np.round(np.abs(X[:, 0]) / np.max(np.abs(X[:, 0])), 3)
    for gn, Fn in (('A', F_A), ('B', F_B), ('C', F_C)):
        F = Fn(X)
        radial = np.array([X[e] @ F[e] / (X[e] @ X[e]) for e in range(len(X))])
        out[gn] = dict(F_norm=float(np.linalg.norm(F)), tangential=float(np.linalg.norm(F - X * radial[:, None])),
                       d_total_norm=float(np.sum(X * F)),
                       radial_rate_by_absRe_class={str(float(c)): float(radial[cls == c].mean()) for c in np.unique(cls)},
                       class_counts={str(float(c)): int(np.sum(cls == c)) for c in np.unique(cls)})
    return out


if __name__ == '__main__':
    res = {}
    for name, X in (('ORTHO16', ortho16()), ('E8_240', e8_240())):
        res[name] = report(name, X)
        # 半径を類ごとに変えても F_A, F_B が 0 のままか（連続な族か）
        Xs = X.copy(); cls = np.round(np.abs(X[:, 0]) / np.max(np.abs(X[:, 0])), 3)
        for k, c in enumerate(np.unique(cls)):
            Xs[cls == c] *= (1.0 + 0.37 * (k + 1))
        res[name]['rescaled_by_class'] = {g: float(np.linalg.norm(Fn(Xs))) for g, Fn in (('A', F_A), ('B', F_B))}
        print(name, json.dumps(res[name], ensure_ascii=False))
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'results', 'covariant_generators_on_designs.json'), 'w'), ensure_ascii=False, indent=1)
