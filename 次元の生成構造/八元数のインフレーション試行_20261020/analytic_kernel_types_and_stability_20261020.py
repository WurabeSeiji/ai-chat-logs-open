#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""番号の側の作用 K_ef = ⟨ℓ u_e, u_f⟩ の解析の数値確認（2026-10-20）。

(1) 振幅込み（L_amp）：K = Q Lᵀ Qᵀ（対角は ⟨ℓu,u⟩ = 0 で自分以外の制限は効かない）、S = QᵀQ は保存、
    解は厳密に Q(t) = Q₀ exp(−tLS)。型は [L,S] と S^{1/2} L S^{1/2} の固有振動数で決まる。
    確認する量：‖K − QLᵀQᵀ‖、‖[L,S]‖、振動数の集合、S の対角。
(2) 位相のみ（L_ph）：dq_e/dt = T L q_e/|q_e|、T = Σ_f |q_f| u_f u_fᵀ。大きさがそろい [L,T₀] = 0 の配置は相対平衡
    q_e(t) = exp(T₀Lt) x_e。回転系での線形化は厳密に dε_e = δT L x_e − a_e T₀ L x_e
    （a_e = x_e·ε_e、δT = Σ_f (ε_f x_fᵀ + x_f ε_fᵀ − a_f x_f x_fᵀ)）。この線形写像を基底に当てて行列を作り（差分なし）、固有値を出す。
作用と配置は octonion_inflation_trials_20261020.py と covariant_generators_on_designs_20261020.py の定義をそのまま読み込む。
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, 'octonion_inflation_trials_20261020.py')).read().split("if __name__")[0]
     .replace("HERE = os.path.dirname(os.path.abspath(__file__))", f"HERE = {HERE!r}"))
from covariant_generators_on_designs_20261020 import e8_240   # noqa: E402

L = LL


def amp_check(Q):
    S = Q.T @ Q; K = kernel(Q, 'L_amp')
    w, V = np.linalg.eigh(S); Sh = V @ np.diag(np.sqrt(np.maximum(w, 0))) @ V.T
    fr = sorted(set(np.round(np.abs(np.linalg.eigvals(Sh @ L @ Sh)), 9).tolist()))
    return dict(K_minus_QLtQt=float(np.linalg.norm(K - Q @ L.T @ Q.T)), comm_L_S=float(np.linalg.norm(L @ S - S @ L)),
                freqs=fr, S_diag=np.round(np.diag(S), 9).tolist())


def ph_linop(X):
    M = len(X); T0 = X.T @ X

    def apply(Ep):
        a = np.sum(Ep * X, 1)
        dT = Ep.T @ X + X.T @ Ep - (X * a[:, None]).T @ X
        return (dT @ L @ X.T).T - a[:, None] * (T0 @ L @ X.T).T
    n = M * 8; J = np.zeros((n, n))
    for i in range(n):
        e = np.zeros(n); e[i] = 1.0; J[:, i] = apply(e.reshape(M, 8)).ravel()
    return J


def poly(N):
    return np.array([math.cos(2 * math.pi * k / N) * E[0] + math.sin(2 * math.pi * k / N) * ELL for k in range(N)])


if __name__ == '__main__':
    cs = configs()
    amp_cfg = dict(ORTHO16=cs['ORTHO16'], NULL28=cs['NULL28'], CROSS4=cs['CROSS4'], E8_240=e8_240())
    ph_cfg = dict(ORTHO16=cs['ORTHO16'], CROSS4=cs['CROSS4'], E8_240=e8_240(), **{f'poly{N}': poly(N) for N in (3, 4, 5, 6, 8)})
    res = dict(amp={}, ph={})
    for n, Q in amp_cfg.items():
        res['amp'][n] = amp_check(Q); print('amp', n, json.dumps(res['amp'][n]), flush=True)
    for n, X in ph_cfg.items():
        S = X.T @ X; J = ph_linop(X); ev = np.linalg.eigvals(J)
        im = sorted(set(np.round(np.abs(ev.imag), 6).tolist()))
        # 冪零（ジョルダン）部分の確認：J の固有値 0 の代数的重複度と幾何的重複度
        rank0 = np.linalg.matrix_rank(J, tol=1e-9); rank0sq = np.linalg.matrix_rank(J @ J, tol=1e-9)
        res['ph'][n] = dict(M=len(X), comm_L_S=float(np.linalg.norm(L @ S - S @ L)), max_Re=float(ev.real.max()),
                            max_absIm=float(np.abs(ev.imag).max()), distinct_absIm=im[:12], rank_J=int(rank0), rank_J2=int(rank0sq))
        print('ph ', n, json.dumps(res['ph'][n]), flush=True)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'results', 'analytic_kernel_types_and_stability.json'), 'w'), ensure_ascii=False, indent=1)
