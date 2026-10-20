#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 3 段（2026-10-20）：(A) g₂ 相対平衡の D の中心化部分環、(B) 閉塞を課した g₂ 相対平衡と k・安定性の照合、
(C) 枠ホロノミーと「干渉の順序」の恒等式（claude.ai `../八元数の干渉と二重スリット_20261020.md` §7 の 3）。
作用・g₂・積分は assoc_stage2_20261020.py（と、それが読み込む本体）をそのまま使う。乱数シード 20261020。

(A) 第 2 段で見つかった 59 解の D ∈ g₂ について、中心化部分環 {A ∈ g₂ : [A, D] = 0} の次元。g₂ の一般の元なら 2（階数）。
    振動数 (ω, ω, 0) の 3 例が特別か（次元が 2 より大きいか）を見る。
(B) 閉塞（Σx = 0、Σ x∘x = 0、八元数の積の二乗）と大きさ Σ|x_e|² = N を課し、F(X) = DX（D ∈ g₂）を X と D について最小二乗で解く。
    N = 4..9、各 8 試行。k = 虚部の行列 X[:,1:] の階数（../八元数閉塞と関係の分解_20261019 と同じ定義）。
    解は射影つき RK4（各歩のあと |x_e| を戻す）で dt と dt/4 の二通りに積分し、形の変化と離れる速さ、終端の閉塞の残差を記録する。
(C) 単位八元数 o_e, o_f, o_g、一歩 a = s_ef = o_f ō_e、b = s_fg、c = s_ge、枠 w。
    恒等式：c(b(aw)) − ((cb)a)w = −[c,b,a]w − [c, ba, w] − c[b,a,w]（[x,y,z] = (xy)z − x(yz)）
    を確かめ、小さな三角形（o_e = 1, o_f = exp(εA), o_g = exp(εB)）で
      |H − I|（枠ホロノミー）、|[c,b,a]|（板どうしの括り方の差）、|[c,ba,w]| + |c[b,a,w]|（枠を含む項）
    の ε 依存の傾き（両対数）を測る。
"""
import importlib.util
import json
import math
import os

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('st2', os.path.join(HERE, 'assoc_stage2_20261020.py'))
st2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(st2)
F_assoc, mul, rk4, G2, releq_g2, assoc_max, E = st2.F_assoc, st2.mul, st2.rk4, st2.G2, st2.releq_g2, st2.assoc_max, st2.E
SEED = 20261020
CONJ = np.array([1.0] + [-1.0] * 7)


def centralizer_dim(D, tol=1e-8):
    M = np.stack([(A @ D - D @ A).ravel() for A in G2], axis=1)
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s < tol * max(s[0], 1e-300)))


def part_A():
    st = json.load(open(os.path.join(HERE, 'results', 'assoc_stage2.json')))
    out = []
    for r in st['g2_relative_equilibria']:
        if 'X' not in r:
            continue
        X = np.array(r['X']); _, D = releq_g2(X, F_assoc(X, 1))
        out.append(dict(N=r['N'], trial=r['trial'], freqs=r['freqs'], centralizer_dim=centralizer_dim(D)))
        print('A', json.dumps(out[-1]), flush=True)
    return out


def k_of(X):
    s = np.linalg.svd(X[:, 1:], compute_uv=False)
    return 0 if s[0] < 1e-300 else int(np.sum(s > 1e-9 * s[0]))


def run_retract(X, dt, steps):
    r0 = np.linalg.norm(X, axis=1); G0 = X @ X.T; Y = X.copy(); gs = []
    for _ in range(steps):
        Y = rk4(lambda Z: F_assoc(Z, 1), Y, dt); Y = Y * (r0 / np.linalg.norm(Y, axis=1))[:, None]
        gs.append(float(np.linalg.norm(Y @ Y.T - G0) / np.linalg.norm(G0)))
    return Y, np.array(gs)


def slope(t, g):
    lg = np.log10(np.maximum(g, 1e-300)); m = (lg > -12) & (lg < -2)
    return float(np.polyfit(t[m], lg[m], 1)[0]) if m.sum() >= 10 else None


def closure_resid(X):
    return float(np.linalg.norm(X.sum(0)) + np.linalg.norm(sum(mul(x, x) for x in X)))


def part_B(rng):
    out = []
    for N in range(4, 10):
        for trial in range(8):
            p0 = np.concatenate([rng.normal(size=8 * N), rng.normal(size=14) * 0.1])

            def res(p):
                X = p[:8 * N].reshape(N, 8)
                D = sum(ci * M for ci, M in zip(p[8 * N:], G2))
                F = F_assoc(X, 1); nF = max(np.linalg.norm(F), 1e-12)
                return np.concatenate([((F - X @ D.T) / nF).ravel(), 10 * X.sum(0), 10 * sum(mul(x, x) for x in X),
                                       [10 * (np.sum(X * X) / N - 1.0), max(0.0, 0.05 - nF / N) * 10]])
            sol = least_squares(res, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
            X = sol.x[:8 * N].reshape(N, 8); F = F_assoc(X, 1); rg, D = releq_g2(X, F)
            row = dict(N=N, trial=trial, nfev=int(sol.nfev), releq_g2=rg, speed=float(np.linalg.norm(F) / N),
                       closure=closure_resid(X), k=k_of(X), assoc=assoc_max(X),
                       min_norm=float(np.linalg.norm(X, axis=1).min()))
            ok = rg < 1e-8 and row['closure'] < 1e-8 and row['speed'] > 1e-6
            row['solved'] = bool(ok)
            if ok:
                w = np.sort(np.abs(np.linalg.eigvals(D).imag))[::-1]
                row['freqs'] = sorted(set(np.round(w[w > 1e-9], 8).tolist()), reverse=True)
                row['centralizer_dim'] = centralizer_dim(D)
                sig = float(max(np.linalg.norm(F, axis=1) / np.linalg.norm(X, axis=1))); dt = 2 * math.pi / (400 * sig)
                for k, (h, n) in enumerate(((dt, 4000), (dt / 4, 16000))):
                    Y, gs = run_retract(X, h, n); t = np.arange(1, n + 1) * h
                    row['dt' if k == 0 else 'dt_quarter'] = dict(h=h, shape_max=float(gs.max()), growth_dec_per_time=slope(t, gs),
                                                                 closure_end=closure_resid(Y), assoc_end=assoc_max(Y))
                row['X'] = X.tolist()
            out.append(row)
            print('B', json.dumps({k: v for k, v in row.items() if k != 'X'}), flush=True)
    return out


def assoc3(x, y, z):
    return mul(mul(x, y), z) - mul(x, mul(y, z))


def oexp(A):
    t = np.linalg.norm(A)
    return E[0].copy() if t < 1e-300 else math.cos(t) * E[0] + math.sin(t) * A / t


def part_C(rng):
    out = {}
    worst = 0.0
    for _ in range(500):
        oe, of, og, w = [v / np.linalg.norm(v) for v in rng.normal(size=(4, 8))]
        a, b, c = mul(of, oe * CONJ), mul(og, of * CONJ), mul(oe, og * CONJ)
        lhs = mul(c, mul(b, mul(a, w))) - mul(mul(mul(c, b), a), w)
        rhs = -mul(assoc3(c, b, a), w) - assoc3(c, mul(b, a), w) - mul(c, assoc3(b, a, w))
        worst = max(worst, float(np.linalg.norm(lhs - rhs)))
    out['identity_max_err'] = worst
    rows = []
    for eps in np.logspace(-3, -0.5, 11):
        H_, P_, W_ = [], [], []
        for _ in range(40):
            A, B = rng.normal(size=(2, 8)); A[0] = B[0] = 0; A /= np.linalg.norm(A); B /= np.linalg.norm(B)
            oe, of, og = E[0], oexp(eps * A), oexp(eps * B)
            a, b, c = mul(of, oe * CONJ), mul(og, of * CONJ), mul(oe, og * CONJ)
            H = np.array([mul(c, mul(b, mul(a, E[k]))) for k in range(8)]).T
            H_.append(float(np.linalg.norm(H - np.eye(8))))
            P_.append(float(np.linalg.norm(assoc3(c, b, a))))
            W_.append(float(np.sqrt(sum(np.linalg.norm(assoc3(c, mul(b, a), E[k])) ** 2 + np.linalg.norm(mul(c, assoc3(b, a, E[k]))) ** 2 for k in range(8)))))
        rows.append(dict(eps=float(eps), holonomy=float(np.median(H_)), plates_bracketing=float(np.median(P_)), frame_terms=float(np.median(W_))))
    out['small_triangles'] = rows
    le = np.log10([r['eps'] for r in rows])
    for key in ('holonomy', 'plates_bracketing', 'frame_terms'):
        out[f'slope_{key}'] = float(np.polyfit(le[:6], np.log10([r[key] for r in rows[:6]]), 1)[0])
    print('C', json.dumps({k: v for k, v in out.items() if k != 'small_triangles'}), flush=True)
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(SEED)
    res = dict(seed=SEED)
    res['C_holonomy_order'] = part_C(rng)
    res['A_centralizer'] = part_A()
    res['B_closure_releq'] = part_B(rng)
    json.dump(res, open(os.path.join(HERE, 'results', 'assoc_stage3.json'), 'w'), ensure_ascii=False, indent=1)
