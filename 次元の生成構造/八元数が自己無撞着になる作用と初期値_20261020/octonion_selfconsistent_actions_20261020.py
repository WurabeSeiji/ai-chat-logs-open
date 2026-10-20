#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""八元数の集合が自己無撞着になる作用と初期値の組を探す（2026-10-20）。

公理：基本量は八元数（代数は公理）。作用は要素どうしの積（左右の掛け算）だけで作り、参照の向き ℓ は置かない（無名性）。
代数：八元数の積は ../八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py の omul をそのまま読み込み、
      構造定数 C[a,b,:] = omul(E_a, E_b) に変換して一括計算に使う（対照：omul との一致を確認して結果に記録）。

作用（どれも各 |x_e| を厳密に保つ。w = x/|x|、κ = ±1 は結合の符号）：
  PAIR_L : dx_e = κ |x_e| Σ_{f≠e} Im(w_f w̄_e) w_e        （要素の対だけ。Artin の定理で結合的）
  PAIR_R : dx_e = κ |x_e| Σ_{f≠e} w_e Im(w̄_e w_f)        （右版）
           解析：どちらも |x_e| (P_e − ⟨P_e, w_e⟩ w_e)、P_e = Σ_{f≠e} w_f（球面上の蔵本型＝Lohe 型。八元数の構造を使わない対照）
  ASSOC  : dx_e = κ |x_e| Σ_{f,g} ⟨w_e w_f, w_g⟩ [w_e, w_f, w_g]
           結合子 [a,b,c] = (ab)c − a(bc)。三つの要素が一つの積に入る、八元数でしか 0 にならない項。
           ⟨[x,f,g], x⟩ = 0（ケーリーの 4 形式の反対称性）なので大きさを保つ。重みなしの和は [w_e, P, P] = 0 に潰れるので、
           f, g について反対称な重み ⟨w_e w_f, w_g⟩ を付ける。
  解析：各相手の次数が 2 以下の作用は、相手を二次のモーメント S = Σ w wᵀ を通してしか見ない。S ∝ I（2-デザイン）なら
        ASSOC は各項が二つの要素で生成される結合的な部分代数に入り 0 になる（E₈・正軸体で確認）。

初期値（固定。乱数は SEED から）：
  N4      : {√3, e₁, e₂, e₄}（本当の二乗の閉塞 Σo² = 0 だけ。e₁, e₂, e₄ は結合的でない三つ組）
  N5      : Σo = 0 かつ Σo² = 0 の解（SEED からの初期点を最小二乗で閉塞へ射影）
  NULL28  : ±(1 ± e_a)/√2
  ORTHO_AB: 正軸体 16 の ±1 を半径 2、±e₄ を半径 1（二次のモーメントが異方的）
  RAND8, RAND16 : SEED からの一般の配置（対称性なし）
  E8, ORTHO16   : 2-デザイン（ASSOC が 0 になることの確認）
判定（各組、走行の生データから機械的に）：
  F0       : 初期の |F|/Σ|x|（< 1e-12 なら静止＝固定点）
  relEq    : F_e = D x_e（共通の D ∈ so(8)）への最小二乗の残差 ‖F − D X‖/‖F‖（< 1e-9 なら相対平衡）
  走行     : RK4、刻み dt = 2π/(200σ)、σ = max_e |F_e|/|x_e|、STEPS 歩。
             g(t) = ‖G(t) − G₀‖/‖G₀‖（G = XXᵀ、形の変化）、relEq(t)、speed(t) = ‖F‖/Σ|x|、assoc(t) = 結合子の大きさの最大
  型       : 静止 / 相対平衡 / 周期的（後半 9 割で g が 1e-8 以下に戻る） / 漸近的（relEq(終)/relEq(始) < 1e-3 か speed(終)/speed(始) < 1e-3）
             / 準周期的（以上のどれでもない、有界で減衰しない）。いずれも暫定の自動判定で、生データを results に残す。
  八元数のまま：assoc（終）> 1e-6 なら結合的な部分代数（ℂ・ℍ）に落ちていない。
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import least_squares

STEPS = 6000
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')
exec(open(SRC).read().split("results = {")[0])          # omul, oconj, E（読み込み元の SEED・rng も入るので、この後で上書きする）
SEED = 20261020
sys.path.insert(0, os.path.join(HERE, '..', '八元数のインフレーション試行_20261020'))
from covariant_generators_on_designs_20261020 import e8_240   # noqa: E402

C = np.array([[omul(E[a], E[b]) for b in range(8)] for a in range(8)])   # C[a,b,c]
CONJ = np.array([1.0] + [-1.0] * 7)


def mul(x, y):                       # x: (..., 8), y: (..., 8)，同じ形か放送できる形
    return np.einsum('...a,...b,abc->...c', x, y, C)


def control_check(rng):
    x, y = rng.normal(size=(2, 200, 8))
    return float(max(np.linalg.norm(mul(x[i], y[i]) - omul(x[i], y[i])) for i in range(200)))


def im(x):
    y = x.copy(); y[..., 0] = 0.0; return y


def F_pair_L(X, kappa):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; N = len(X)
    Wc = W * CONJ
    out = np.zeros_like(X)
    for e in range(N):
        prods = im(mul(W, Wc[e][None, :]))                  # Im(w_f w̄_e)
        prods[e] = 0.0
        out[e] = mul(prods.sum(0), W[e]) * r[e]
    return kappa * out


def F_pair_R(X, kappa):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; N = len(X)
    Wc = W * CONJ
    out = np.zeros_like(X)
    for e in range(N):
        prods = im(mul(Wc[e][None, :], W))                  # Im(w̄_e w_f)
        prods[e] = 0.0
        out[e] = mul(W[e], prods.sum(0)) * r[e]
    return kappa * out


def F_assoc(X, kappa):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; N = len(X)
    WW = mul(W[:, None, :], W[None, :, :])                  # WW[e,f] = w_e w_f
    out = np.zeros_like(X)
    for e in range(N):
        left = mul(WW[e][:, None, :], W[None, :, :])         # (w_e w_f) w_g   [f,g]
        right = mul(W[e][None, None, :], WW)                 # w_e (w_f w_g)   [f,g]
        A = left - right
        wgt = WW[e] @ W.T                                     # ⟨w_e w_f, w_g⟩  [f,g]
        out[e] = np.einsum('fg,fgc->c', wgt, A) * r[e]
    return kappa * out


ACTIONS = {'PAIR_L': F_pair_L, 'PAIR_R': F_pair_R, 'ASSOC': F_assoc}


def so8_basis():
    B = []
    for i in range(8):
        for j in range(i + 1, 8):
            M = np.zeros((8, 8)); M[i, j] = -1; M[j, i] = 1; B.append(M)
    return B


SO8 = so8_basis()


def releq_resid(X, F):
    nF = np.linalg.norm(F)
    if nF < 1e-300:
        return 0.0
    A = np.stack([(X @ Bk.T).ravel() for Bk in SO8], axis=1)
    c, *_ = np.linalg.lstsq(A, F.ravel(), rcond=None)
    return float(np.linalg.norm(F.ravel() - A @ c) / nF)


def assoc_max(X):
    W = X / np.linalg.norm(X, axis=1, keepdims=True)
    WW = mul(W[:, None, :], W[None, :, :])
    m = 0.0
    for e in range(len(W)):
        A = mul(WW[e][:, None, :], W[None, :, :]) - mul(W[e][None, None, :], WW)
        m = max(m, float(np.max(np.linalg.norm(A, axis=2))))
    return m


def configs(rng):
    c = {}
    c['N4'] = np.array([math.sqrt(3) * E[0], E[1], E[2], E[4]])

    def res(v):
        O = v.reshape(5, 8); s1 = O.sum(0); s2 = sum(mul(o, o) for o in O)
        return np.concatenate([s1, s2])
    sol = least_squares(res, rng.normal(size=40), xtol=1e-15, ftol=1e-15, gtol=1e-15)
    c['N5'] = sol.x.reshape(5, 8)
    c['NULL28'] = np.array([s * (E[0] + t * E[a]) / math.sqrt(2) for a in range(1, 8) for t in (1, -1) for s in (1, -1)])
    X = np.array([s * E[a] for a in range(8) for s in (1, -1)]); X[0:2] *= 2.0
    c['ORTHO_AB'] = X
    c['RAND8'] = rng.normal(size=(8, 8))
    c['RAND16'] = rng.normal(size=(16, 8))
    return c


def rk4(F, X, dt):
    k1 = F(X); k2 = F(X + 0.5 * dt * k1); k3 = F(X + 0.5 * dt * k2); k4 = F(X + dt * k3)
    return X + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def classify(r):
    if r['F0'] < 1e-12:
        return '静止（固定点）'
    if r['relEq0'] < 1e-9:
        return '相対平衡'
    if r['g_min_late'] < 1e-8:
        return '周期的'
    if r['relEq_end'] / r['relEq0'] < 1e-3 or r['speed_end'] / r['F0'] < 1e-3:
        return '漸近的'
    return '準周期的（または不規則）'


def trial(name, X0, aname, kappa):
    F = lambda X: ACTIONS[aname](X, kappa)
    F0 = F(X0); sx = float(np.linalg.norm(X0, axis=1).sum())
    r = dict(config=name, action=aname, kappa=kappa, N=len(X0), F0=float(np.linalg.norm(F0) / sx),
             relEq0=releq_resid(X0, F0), assoc0=assoc_max(X0))
    if r['F0'] < 1e-12 or r['relEq0'] < 1e-9:
        r['type'] = classify(r); return r, None
    sigma = float(max(np.linalg.norm(F0, axis=1) / np.linalg.norm(X0, axis=1)))
    dt = 2 * math.pi / (200 * sigma)
    G0 = X0 @ X0.T; X = X0.copy(); rows = []
    rec = max(1, STEPS // 300)
    for t in range(1, STEPS + 1):
        X = rk4(F, X, dt)
        g = float(np.linalg.norm(X @ X.T - G0) / np.linalg.norm(G0))
        if t % rec == 0 or t == STEPS:
            Fx = F(X)
            rows.append([t * dt, g, releq_resid(X, Fx), float(np.linalg.norm(Fx) / sx)])
        else:
            rows.append([t * dt, g, np.nan, np.nan])
    rows = np.array(rows)
    late = rows[STEPS // 10:]
    r.update(sigma=sigma, dt=dt, steps=STEPS, T=STEPS * dt, g_max=float(np.nanmax(rows[:, 1])),
             g_min_late=float(np.nanmin(late[:, 1])), relEq_end=float(rows[-1, 2]), speed_end=float(rows[-1, 3]),
             assoc_end=assoc_max(X),
             norm_drift=float(np.max(np.abs(np.linalg.norm(X, axis=1) - np.linalg.norm(X0, axis=1)))))
    r['type'] = classify(r)
    r['octonion_kept'] = bool(r['assoc_end'] > 1e-6)
    return r, rows


if __name__ == '__main__':
    rng = np.random.default_rng(SEED)
    out = dict(seed=SEED, steps=STEPS, control_mul_vs_omul=control_check(rng))
    cfg = configs(rng)
    # 初期値の性質
    out['configs'] = {}
    for n, X in cfg.items():
        S = X.T @ X
        out['configs'][n] = dict(N=len(X), sum_o=float(np.linalg.norm(X.sum(0))),
                                 sum_o2=float(np.linalg.norm(sum(mul(o, o) for o in X))),
                                 S_aniso=float(np.linalg.norm(S - np.trace(S) / 8 * np.eye(8)) / np.linalg.norm(S)),
                                 assoc_max=assoc_max(X))
    # 解析の確認：PAIR は Lohe 型と一致、2-デザインで ASSOC = 0
    X = cfg['RAND16']; r = np.linalg.norm(X, axis=1); W = X / r[:, None]
    P = W.sum(0)[None, :] - W
    lohe = r[:, None] * (P - np.sum(P * W, 1)[:, None] * W)
    out['check_pair_equals_lohe'] = dict(L=float(np.linalg.norm(F_pair_L(X, 1) - lohe)), R=float(np.linalg.norm(F_pair_R(X, 1) - lohe)))
    ortho16 = np.array([s * E[a] for a in range(8) for s in (1, -1)])
    E8 = e8_240()
    out['check_assoc_zero_on_2designs'] = dict(ORTHO16=float(np.linalg.norm(F_assoc(ortho16, 1))),
                                               E8_240=float(np.linalg.norm(F_assoc(E8, 1))))
    print(json.dumps({k: v for k, v in out.items()}, ensure_ascii=False), flush=True)
    out['trials'] = []
    series = {}
    for n, X0 in cfg.items():
        for aname in ('PAIR_L', 'ASSOC'):                    # PAIR_R は PAIR_L と恒等（上の確認）
            for kappa in (1, -1):
                r, rows = trial(n, X0, aname, kappa)
                out['trials'].append(r)
                if rows is not None:
                    series[f'{n}_{aname}_{"p" if kappa > 0 else "m"}'] = rows
                print(json.dumps({k: (float('%.3e' % v) if isinstance(v, float) else v) for k, v in r.items()}, ensure_ascii=False), flush=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'selfconsistent_actions.json'), 'w'), ensure_ascii=False, indent=1)
    np.savez_compressed(os.path.join(HERE, 'results', 'selfconsistent_actions_series.npz'), **series)
