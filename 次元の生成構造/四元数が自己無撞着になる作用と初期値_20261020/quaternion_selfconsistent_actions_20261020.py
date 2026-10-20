#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""四元数の初期値と、交換子で作る作用の組で、自己無撞着になるかとその型を調べる（2026-10-20）。

枠組み：四元数の塊は力学から自然に出ないので初期値として与える。自己無撞着は作用と初期値の組で決まる。
代数：四元数の積は ../八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py の qmul をそのまま読み込む（構造定数に変換）。
材料：交換子 [a,b] = ab − ba（ℂ では 0、ℍ で初めて残る。八元数の結合子に当たる一段下の材料）。w = x/|x|、κ = ±1。
  C1：dx_e = κ|x_e| [w_e, P_e]、P_e = Σ_{f≠e} w_f（相手について 1 次）
  C2：dx_e = κ|x_e| Σ_{f≠e} ⟨w_e, w_f⟩ [w_e, w_f] = κ|x_e| [w_e, S_w w_e]、S_w = Σ_f w_f w_fᵀ（相手について 2 次）
  どちらも ⟨[w_e, ·], w_e⟩ = 0 で各 |x_e| を保つ。共変性は Aut(ℍ) = SO(3)（共役 q ↦ g q ḡ）。
解析の予想（C1）：[w, P] = 2 Im w × Im P。Σ_e |x_e|·(Im w_e × Im P)/|x_e| … w の虚部は Im P のまわりに歳差し、
  dP/dt = Σ_e d w_e/dt = 2κ Σ Im w_e × Im P = 2κ Im P × Im P = 0 で P は保存。よって全ての配置が相対平衡
  x_e(t) = exp(−tκ ad_P) x_e、振動数 2|Im P|。全要素が一つの複素平面上、または Im P = 0 なら静止。
初期値：CROSS8（±1, ±i, ±j, ±k）、T24（2T、24 胞体）、O48（2O）、I120（2I、600 胞体）、
  TWO_C（平面 (1,i) の 1 の 3 乗根 ∪ 平面 (1,j) の 1 の 3 乗根、N = 6）、ONE_C（平面 (1,i) の 1 の 6 乗根、対照）、
  RAND4, RAND8（SEED からの乱数）。
判定（各組）：F0（初期の ‖F‖/Σ|x|、< 1e-12 で静止）、relEq（F_e = D x_e、D ∈ so(3)=Aut(ℍ) の微分、最小二乗の残差、< 1e-9 で相対平衡）、
  so(4) でも同じ残差を記録。走行：RK4 ＋各歩のあと |x_e| を戻す射影、dt = 2π/(200σ)、σ = max_e |F_e|/|x_e|、4000 歩。
  g(t) = ‖G(t) − G₀‖/‖G₀‖（G = XXᵀ）、後半 9 割の g の最小（戻り）、1e-8 ずらした二軌道の形の差の育ち方（不規則の判定）、
  終端の relEq と速さ。型：静止 / 相対平衡 / 周期的（g が 1e-8 以下に戻る） / 漸近的（relEq か速さが 1e-3 倍以下に減る）/
  不規則（二軌道の差が 4 桁以上育つ） / 準周期的（以上のどれでもない）。
(2) C2 の相対平衡を初期値と D（so(3)）を同時に解く：N = 3..8、各 8 試行、各 |x_e| = 1。解の部分代数（ℂ か ℍ）と持続性。
乱数シード 20261020。
"""
import json
import math
import os

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')
exec(open(SRC).read().split("results = {")[0])          # qmul, qconj（SEED・rng も入るので後で上書き）
SEED = 20261020
STEPS = 4000
I4 = np.eye(4)
C = np.array([[qmul(I4[a], I4[b]) for b in range(4)] for a in range(4)])
CONJ = np.array([1.0, -1.0, -1.0, -1.0])


def mul(x, y):
    return np.einsum('...a,...b,abc->...c', x, y, C)


def comm(x, y):
    return mul(x, y) - mul(y, x)


def control():
    rng = np.random.default_rng(1)
    x, y = rng.normal(size=(2, 50, 4))
    return float(max(np.linalg.norm(mul(x[i], y[i]) - qmul(x[i], y[i])) for i in range(50)))


def F_C1(X, kappa):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; P = W.sum(0)
    return kappa * r[:, None] * comm(W, (P[None, :] - W))


def F_C2(X, kappa):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; S = W.T @ W
    return kappa * r[:, None] * comm(W, W @ S.T)


ACTIONS = {'C1': F_C1, 'C2': F_C2}

# so(3) = Aut(ℍ) の微分：ad_u(x) = ux − xu（u 純虚）。so(4)：6 生成子
SO3 = []
for k in range(1, 4):
    SO3.append(np.array([comm(I4[k], I4[j]) for j in range(4)]).T)
SO4 = []
for i in range(4):
    for j in range(i + 1, 4):
        M = np.zeros((4, 4)); M[i, j] = -1; M[j, i] = 1; SO4.append(M)


def releq(X, F, basis):
    nF = np.linalg.norm(F)
    if nF < 1e-300:
        return 0.0, None
    A = np.stack([(X @ B.T).ravel() for B in basis], axis=1)
    c, *_ = np.linalg.lstsq(A, F.ravel(), rcond=None)
    return float(np.linalg.norm(F.ravel() - A @ c) / nF), sum(ci * B for ci, B in zip(c, basis))


def subalgebra_dim(X, tol=1e-8):
    V = np.array(X, float)
    for _ in range(3):
        prods = [mul(a, b) for a in V for b in V]
        V = np.vstack([V, prods, I4[0][None, :]])
        u, s, vt = np.linalg.svd(V, full_matrices=False)
        V = vt[:int(np.sum(s > tol * s[0]))]
    return len(V)


def commutator_max(X):
    W = X / np.linalg.norm(X, axis=1, keepdims=True)
    return float(max(np.linalg.norm(comm(a, b)) for a in W for b in W))


def configs(rng):
    c = {}
    c['CROSS8'] = np.array([s * I4[a] for a in range(4) for s in (1, -1)])
    half = [np.array(v) / 2 for v in np.array(np.meshgrid(*[[1, -1]] * 4)).T.reshape(-1, 4)]
    c['T24'] = np.vstack([c['CROSS8'], half])
    pairs = [(s1 * I4[a] + s2 * I4[b]) / math.sqrt(2) for a in range(4) for b in range(a + 1, 4) for s1 in (1, -1) for s2 in (1, -1)]
    c['O48'] = np.vstack([c['T24'], pairs])
    phi = (1 + math.sqrt(5)) / 2
    import itertools
    ev = []
    base = [0, 1, phi, 1 / phi]
    for perm in itertools.permutations(range(4)):
        sign = 1
        p = list(perm)
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]:
                    sign = -sign
        if sign < 0:
            continue
        for s in itertools.product((1, -1), repeat=3):
            v = np.zeros(4); vals = [base[0], s[0] * base[1], s[1] * base[2], s[2] * base[3]]
            for k in range(4):
                v[p[k]] = vals[k]
            ev.append(v / 2)
    c['I120'] = np.unique(np.round(np.vstack([c['T24'], ev]), 12), axis=0)
    roots3 = [math.cos(2 * math.pi * k / 3) * I4[0] + math.sin(2 * math.pi * k / 3) * I4[1] for k in range(3)]
    roots3j = [math.cos(2 * math.pi * k / 3) * I4[0] + math.sin(2 * math.pi * k / 3) * I4[2] for k in range(3)]
    c['TWO_C'] = np.array(roots3 + roots3j)
    c['ONE_C'] = np.array([math.cos(2 * math.pi * k / 6) * I4[0] + math.sin(2 * math.pi * k / 6) * I4[1] for k in range(6)])
    c['RAND4'] = rng.normal(size=(4, 4))
    c['RAND8'] = rng.normal(size=(8, 4))
    return c


def rk4r(F, X, dt, r0):
    k1 = F(X); k2 = F(X + 0.5 * dt * k1); k3 = F(X + 0.5 * dt * k2); k4 = F(X + dt * k3)
    Y = X + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return Y * (r0 / np.linalg.norm(Y, axis=1))[:, None]


def classify(r):
    if r['F0'] < 1e-12:
        return '静止（固定点）'
    if r['relEq_so3'] < 1e-9:
        return '相対平衡'
    if r['g_min_late'] < 1e-8:
        return '周期的'
    if r['relEq_end'] / r['relEq_so3'] < 1e-3 or r['speed_end'] / r['F0'] < 1e-3:
        return '漸近的'
    if r['sep_growth_decades'] > 4:
        return '不規則'
    return '準周期的'


def trial(name, X0, aname, kappa):
    F = lambda X: ACTIONS[aname](X, kappa)
    F0 = F(X0); sx = float(np.linalg.norm(X0, axis=1).sum())
    r3, D3 = releq(X0, F0, SO3); r4, _ = releq(X0, F0, SO4)
    r = dict(config=name, action=aname, kappa=kappa, N=len(X0), F0=float(np.linalg.norm(F0) / sx),
             relEq_so3=r3, relEq_so4=r4, subalgebra_dim=subalgebra_dim(X0), commutator_max=commutator_max(X0))
    if r['F0'] < 1e-12 or r3 < 1e-9:
        if D3 is not None:
            w = np.abs(np.linalg.eigvals(D3).imag); r['freqs'] = sorted(set(np.round(w[w > 1e-9], 9).tolist()))
        r['type'] = classify(r); return r
    r0 = np.linalg.norm(X0, axis=1)
    sigma = float(max(np.linalg.norm(F0, axis=1) / r0)); dt = 2 * math.pi / (200 * sigma)
    G0 = X0 @ X0.T; X = X0.copy(); Xp = X0.copy(); Xp[0, 1] += 1e-8; Xp[0] *= r0[0] / np.linalg.norm(Xp[0])
    gs = []; sep0 = float(np.linalg.norm(Xp @ Xp.T - G0))
    for _ in range(STEPS):
        X = rk4r(F, X, dt, r0); Xp = rk4r(F, Xp, dt, r0)
        gs.append(float(np.linalg.norm(X @ X.T - G0) / np.linalg.norm(G0)))
    sepE = float(np.linalg.norm(Xp @ Xp.T - X @ X.T))
    Fe = F(X)
    r.update(sigma=sigma, dt=dt, steps=STEPS, g_max=float(max(gs)), g_min_late=float(min(gs[STEPS // 10:])),
             relEq_end=releq(X, Fe, SO3)[0], speed_end=float(np.linalg.norm(Fe) / sx),
             sep_growth_decades=float(math.log10(max(sepE, 1e-300) / sep0)),
             subalgebra_dim_end=subalgebra_dim(X))
    r['type'] = classify(r)
    return r


def solve_C2_releq(rng):
    out = []
    for N in range(3, 9):
        for t in range(8):
            p0 = np.concatenate([rng.normal(size=4 * N), rng.normal(size=3) * 0.1])

            def res(p):
                X = p[:4 * N].reshape(N, 4); X = X / np.linalg.norm(X, axis=1, keepdims=True)
                D = sum(ci * B for ci, B in zip(p[4 * N:], SO3))
                F = F_C2(X, 1); nF = max(np.linalg.norm(F), 1e-12)
                return np.concatenate([((F - X @ D.T) / nF).ravel(), [max(0.0, 0.05 - nF / N) * 10]])
            sol = least_squares(res, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
            X = sol.x[:4 * N].reshape(N, 4); X = X / np.linalg.norm(X, axis=1, keepdims=True)
            F = F_C2(X, 1); rg, D = releq(X, F, SO3)
            row = dict(N=N, trial=t, releq_so3=rg, speed=float(np.linalg.norm(F) / N), subalgebra_dim=subalgebra_dim(X),
                       commutator_max=commutator_max(X))
            row['solved'] = bool(rg < 1e-8 and row['speed'] > 1e-6)
            if row['solved']:
                w = np.abs(np.linalg.eigvals(D).imag); row['freq'] = float(w.max())
                r0 = np.ones(N); sig = float(max(np.linalg.norm(F, axis=1))); dt = 2 * math.pi / (400 * sig)
                G0 = X @ X.T; Y = X.copy(); gmax = 0.0
                for _ in range(4000):
                    Y = rk4r(lambda Z: F_C2(Z, 1), Y, dt, r0); gmax = max(gmax, float(np.linalg.norm(Y @ Y.T - G0) / np.linalg.norm(G0)))
                row.update(check_dt=dt, check_shape_max=gmax, X=X.tolist())
            out.append(row)
            print('C2rel', json.dumps({k: v for k, v in row.items() if k != 'X'}), flush=True)
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(SEED)
    res = dict(seed=SEED, control_mul_vs_qmul=control(), trials=[])
    for n, X0 in configs(rng).items():
        for a in ('C1', 'C2'):
            for kappa in (1, -1):
                r = trial(n, X0, a, kappa); res['trials'].append(r)
                print(json.dumps({k: (float('%.4g' % v) if isinstance(v, float) else v) for k, v in r.items()}, ensure_ascii=False), flush=True)
    res['C2_releq_solve'] = solve_C2_releq(rng)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'results', 'quaternion_selfconsistent_actions.json'), 'w'), ensure_ascii=False, indent=1)
