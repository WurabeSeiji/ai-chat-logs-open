#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ASSOC（結合子の作用）の第 2 段（2026-10-20）。作用・初期値は octonion_selfconsistent_actions_20261020.py をそのまま読み込む。

(1) 固定点のまわりの線形化：N4・NULL28・ORTHO_AB の ASSOC（κ = ±1）。目的は自己無撞着な軌道の候補を出すこと：
    純虚の固有値の対 ±iω は、それぞれ固定点から出る周期 2π/ω の周期軌道の族を与える（Lyapunov の中心定理、非共鳴のとき）。
    振動数の集合と、比が有理数か（共鳴＝周期的な族の重なり）を記録する。ヤコビアンは中心差分で h = 1e-5 と 1e-6 の二通り
    （実部が h とともに縮めば差分誤差。実部が残る場合は、その方向には周期軌道の族が出ないことだけを記録する）。
(2) g₂ の相対平衡を解く：u_e(X) = F_e(X) = D X_e（D ∈ g₂、全 e 共通）を、初期値 X（各 |x_e| = 1）と D の 14 係数を
    同時に未知数とした最小二乗（scipy least_squares）で解く。残差は ‖F − DX‖/‖F‖。N = 4..8、各 N で乱数の初期点 12 個。
    g₂ の基底は八元数の微分 D(xy) = (Dx)y + x(Dy) の解空間（14 次元）をここで計算する。
    見つかった解は RK4 で積分し、形（グラム行列）が保たれるか、結合子が 0 でないか（八元数のまま）を確かめる。
(3) N5 の長い走行：刻み dt = 2π/(2000σ)、60000 歩、κ = ±1。形の量 G₀₁(t) = ⟨x₀, x₁⟩ のパワースペクトル（離散の線か連続か）と、
    1e-8 だけずらした二軌道の離れ方（指数的なら不規則）。
乱数シード 20261020（本体と同じ値。流れは本体と別にこのプログラムの先頭から）。
"""
import importlib.util
import json
import math
import os

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('sca', os.path.join(HERE, 'octonion_selfconsistent_actions_20261020.py'))
sca = importlib.util.module_from_spec(spec); spec.loader.exec_module(sca)
F_assoc, mul, E, releq_resid, assoc_max, rk4 = sca.F_assoc, sca.mul, sca.E, sca.releq_resid, sca.assoc_max, sca.rk4
SEED = sca.SEED


def g2_basis():
    rows = []
    for i in range(8):
        for j in range(8):
            p = mul(E[i], E[j])
            for c in range(8):
                r = np.zeros(64)
                for m in range(8):
                    r[c * 8 + m] += p[m]
                    r[m * 8 + i] -= mul(E[m], E[j])[c]
                    r[m * 8 + j] -= mul(E[i], E[m])[c]
                rows.append(r)
    u, s, vt = np.linalg.svd(np.array(rows))
    return [v.reshape(8, 8) for v in vt[int(np.sum(s > 1e-9 * s[0])):]]


G2 = g2_basis()


def releq_g2(X, F):
    nF = np.linalg.norm(F)
    if nF < 1e-300:
        return 0.0, None
    A = np.stack([(X @ D.T).ravel() for D in G2], axis=1)
    c, *_ = np.linalg.lstsq(A, F.ravel(), rcond=None)
    return float(np.linalg.norm(F.ravel() - A @ c) / nF), sum(ci * D for ci, D in zip(c, G2))


def jac(F, X, h):
    x = X.ravel(); n = len(x); J = np.zeros((n, n))
    for i in range(n):
        d = np.zeros(n); d[i] = h
        J[:, i] = (F((x + d).reshape(X.shape)) - F((x - d).reshape(X.shape))).ravel() / (2 * h)
    return J


def part1(cfg):
    out = {}
    for n in ('N4', 'NULL28', 'ORTHO_AB'):
        X = cfg[n]
        for kappa in (1, -1):
            F = lambda Y: F_assoc(Y, kappa)
            r = dict(F0=float(np.linalg.norm(F(X))))
            for h in (1e-5, 1e-6):
                ev = np.linalg.eigvals(jac(F, X, h))
                r[f'h{h:.0e}'] = dict(max_Re=float(ev.real.max()), max_absIm=float(np.abs(ev.imag).max()),
                                      n_absRe_gt_1em3=int(np.sum(np.abs(ev.real) > 1e-3)),
                                      distinct_absIm=sorted(set(np.round(np.abs(ev.imag[np.abs(ev.imag) > 1e-6]), 4).tolist()))[:12])
            out[f'{n}_k{kappa:+d}'] = r
            print('lin', n, kappa, json.dumps(r), flush=True)
    return out


def part2(rng):
    out = []
    for N in range(4, 9):
        for trial in range(12):
            x0 = rng.normal(size=(N, 8)); c0 = rng.normal(size=14) * 0.1

            def res(p):
                X = p[:8 * N].reshape(N, 8); X = X / np.linalg.norm(X, axis=1, keepdims=True)
                D = sum(ci * M for ci, M in zip(p[8 * N:], G2))
                F = F_assoc(X, 1)
                nF = max(np.linalg.norm(F), 1e-12)
                return np.concatenate([((F - X @ D.T) / nF).ravel(), [max(0.0, 0.05 - nF / N) * 10]])
            sol = least_squares(res, np.concatenate([x0.ravel(), c0]), xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=3000)
            X = sol.x[:8 * N].reshape(N, 8); X = X / np.linalg.norm(X, axis=1, keepdims=True)
            F = F_assoc(X, 1); rg, D = releq_g2(X, F)
            row = dict(N=N, trial=trial, nfev=int(sol.nfev), speed=float(np.linalg.norm(F) / N), releq_g2=rg,
                       releq_so8=releq_resid(X, F), assoc=assoc_max(X))
            if rg < 1e-8 and row['speed'] > 1e-6:
                w = np.sort(np.abs(np.linalg.eigvals(D).imag))[::-1]
                row['freqs'] = sorted(set(np.round(w[w > 1e-9], 8).tolist()), reverse=True)
                # 積分して確かめる
                sig = float(max(np.linalg.norm(F, axis=1))); dt = 2 * math.pi / (400 * sig); Y = X.copy(); G0 = X @ X.T
                gmax = 0.0
                for _ in range(4000):
                    Y = rk4(lambda Z: F_assoc(Z, 1), Y, dt)
                    gmax = max(gmax, float(np.linalg.norm(Y @ Y.T - G0) / np.linalg.norm(G0)))
                row.update(check_steps=4000, check_dt=dt, check_shape_max=gmax, check_assoc_end=assoc_max(Y),
                           X=X.tolist())
            out.append(row)
            print('g2', json.dumps({k: v for k, v in row.items() if k != 'X'}), flush=True)
    return out


def part3(cfg):
    out = {}
    X0 = cfg['N5']
    for kappa in (1, -1):
        F = lambda Y: F_assoc(Y, kappa)
        sig = float(max(np.linalg.norm(F(X0), axis=1) / np.linalg.norm(X0, axis=1)))
        dt = 2 * math.pi / (2000 * sig); steps = 60000
        X = X0.copy(); Xp = X0.copy(); Xp[0, 1] += 1e-8
        g01, sep, drift = [], [], 0.0
        for t in range(steps):
            X = rk4(F, X, dt); Xp = rk4(F, Xp, dt)
            g01.append(float(X[0] @ X[1]))
            if t % 100 == 0:
                sep.append(float(np.linalg.norm((Xp @ Xp.T) - (X @ X.T))))
        drift = float(np.max(np.abs(np.linalg.norm(X, axis=1) - np.linalg.norm(X0, axis=1))))
        g = np.array(g01) - np.mean(g01); P = np.abs(np.fft.rfft(g * np.hanning(len(g)))) ** 2
        fr = np.fft.rfftfreq(len(g), dt) * 2 * math.pi
        top = np.argsort(P)[::-1][:8]
        Pn = P / P.sum(); ent = float(-(Pn[Pn > 0] * np.log(Pn[Pn > 0])).sum() / math.log(len(Pn)))
        sep = np.array(sep); ls = np.log10(np.maximum(sep, 1e-300))
        # 離れ方：前半と後半の傾き（桁/時間）
        tt = np.arange(len(sep)) * 100 * dt; h = len(sep) // 2
        out[f'k{kappa:+d}'] = dict(sigma=sig, dt=dt, steps=steps, T=steps * dt, norm_drift=drift,
                                   spectral_entropy=ent, top_omegas=[float(fr[i]) for i in top],
                                   top_power_frac=[float(Pn[i]) for i in top],
                                   sep_log10_start=float(ls[1]), sep_log10_mid=float(ls[h]), sep_log10_end=float(ls[-1]),
                                   sep_slope_first_half=float(np.polyfit(tt[1:h], ls[1:h], 1)[0]),
                                   sep_slope_second_half=float(np.polyfit(tt[h:], ls[h:], 1)[0]),
                                   assoc_end=assoc_max(X))
        np.savez_compressed(os.path.join(HERE, 'results', f'assoc_stage2_N5_k{kappa:+d}.npz'), g01=np.array(g01), sep=sep, dt=dt)
        print('N5', kappa, json.dumps(out[f'k{kappa:+d}']), flush=True)
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(SEED)
    cfg = sca.configs(np.random.default_rng(SEED))
    res = dict(seed=SEED, dim_g2=len(G2))
    res['linearization'] = part1(cfg)
    res['g2_relative_equilibria'] = part2(rng)
    res['N5_long'] = part3(cfg)
    json.dump(res, open(os.path.join(HERE, 'results', 'assoc_stage2.json'), 'w'), ensure_ascii=False, indent=1)
