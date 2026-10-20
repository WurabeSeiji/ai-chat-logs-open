#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""高対称な四元数と八元数を初期値にしたインフレーション実験（2026-10-20）。

著者の指摘：高対称でなければ一歩でジャンプして作用が破綻する。極めて高対称だから、揺らぎが幾何級数的に育つ過程が
緩和された形で見える状態がある。複素数の十字（床）でも、作用から見えない状態のわずかな揺らぎが蓄積・拡大した。

作用と配置は既存のプログラムをそのまま読み込む（再実装しない）：
  四元数：../四元数が自己無撞着になる作用と初期値_20261020/quaternion_selfconsistent_actions_20261020.py の F_C1, F_C2, configs
  八元数：../八元数が自己無撞着になる作用と初期値_20261020/assoc_stage2_20261020.py の F_assoc, G2, mul
組：四元数 CROSS8, T24, O48, I120, TWO_C × C1, C2 × κ = ±1。八元数 ORTHO16, NULL28 × ASSOC × κ = ±1。
揺らぎの種：外から足さない。各配置を (a) そのまま、(b) 自己同型でランダムに回したもの（四元数は q ↦ u q ū、u は単位四元数、
  八元数は exp(D)、D ∈ g₂ の乱数）で置く。配置は物理的に同じで、浮動小数の表し方だけが違う（種は丸め誤差だけ）。
刻み：同じ N の一般の配置（乱数、シード固定）での σ_ref = max_e |F_e|/|x_e| から dt = 2π/(200 σ_ref)。
走行：RK4 ＋各歩のあと |x_e| を初期値に戻す射影。四元数 20000 歩、八元数 ORTHO16 20000 歩、NULL28 8000 歩。
測る量：g(t) = ‖G(t) − G₀‖/‖G₀‖（G = XXᵀ、形の変化）。10 歩ごとに記録。速さ ‖F‖/Σ|x|。
  log10 g が −13 から −3 の区間の傾き（桁/時間）＝指数的な育ちの速さ。終端の形の変化、生成する部分代数の次元、
  結合子（八元数）・交換子（四元数）の最大、終端の相対平衡の残差（四元数は so(3)、八元数は g₂）。
乱数シード 20261020。
"""
import importlib.util
import json
import math
import os

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


Q = load('q', os.path.join(HERE, '..', '四元数が自己無撞着になる作用と初期値_20261020', 'quaternion_selfconsistent_actions_20261020.py'))
O = load('o', os.path.join(HERE, '..', '八元数が自己無撞着になる作用と初期値_20261020', 'assoc_stage2_20261020.py'))
SEED = 20261020


def ortho16():
    return np.array([s * O.E[a] for a in range(8) for s in (1, -1)])


def null28():
    E = O.E
    return np.array([s * (E[0] + t * E[a]) / math.sqrt(2) for a in range(1, 8) for t in (1, -1) for s in (1, -1)])


def rk4r(F, X, dt, r0):
    k1 = F(X); k2 = F(X + 0.5 * dt * k1); k3 = F(X + 0.5 * dt * k2); k4 = F(X + dt * k3)
    Y = X + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return Y * (r0 / np.linalg.norm(Y, axis=1))[:, None]


def run(X0, F, sigma_ref, steps, label, extra):
    r0 = np.linalg.norm(X0, axis=1); dt = 2 * math.pi / (200 * sigma_ref); sx = float(r0.sum())
    G0 = X0 @ X0.T; nG = np.linalg.norm(G0); X = X0.copy(); rec_t, rec_g, rec_v = [], [], []
    for t in range(1, steps + 1):
        X = rk4r(F, X, dt, r0)
        if t % 10 == 0 or t == steps:
            rec_t.append(t * dt); rec_g.append(float(np.linalg.norm(X @ X.T - G0) / nG))
            rec_v.append(float(np.linalg.norm(F(X)) / sx))
    tt = np.array(rec_t); g = np.maximum(np.array(rec_g), 1e-300); lg = np.log10(g)
    m = (lg > -13) & (lg < -3)
    slope = float(np.polyfit(tt[m], lg[m], 1)[0]) if m.sum() >= 10 else None
    i3 = np.where(lg > -3)[0]
    out = dict(label=label, N=len(X0), dt=dt, steps=steps, T=steps * dt, F0=float(np.linalg.norm(F(X0)) / sx),
               g_start=float(g[0]), g_end=float(g[-1]), g_max=float(g.max()), growth_dec_per_time=slope,
               t_reach_1e_3=float(tt[i3[0]]) if len(i3) else None,
               speed_end=float(rec_v[-1]), **extra(X))
    np.savez_compressed(os.path.join(HERE, 'results', f'series_{label}.npz'), t=tt, g=np.array(rec_g), v=np.array(rec_v))
    print(json.dumps({k: (float('%.4g' % v) if isinstance(v, float) else v) for k, v in out.items()}, ensure_ascii=False), flush=True)
    return out


def q_extra(X):
    F = Q.F_C2(X, 1)
    return dict(subalgebra_dim_end=Q.subalgebra_dim(X), commutator_max_end=Q.commutator_max(X),
                releq_so3_end_C2=Q.releq(X, F, Q.SO3)[0])


def o_extra(X):
    F = O.F_assoc(X, 1)
    return dict(assoc_max_end=O.assoc_max(X), releq_g2_end=O.releq_g2(X, F)[0])


if __name__ == '__main__':
    rng = np.random.default_rng(SEED)
    res = dict(seed=SEED, runs=[])
    qc = Q.configs(np.random.default_rng(SEED))
    for name in ('CROSS8', 'T24', 'O48', 'I120', 'TWO_C'):
        X = qc[name]; N = len(X)
        u = rng.normal(size=4); u /= np.linalg.norm(u)
        Xrot = np.array([Q.mul(Q.mul(u, x), u * Q.CONJ) for x in X])
        Xg = np.random.default_rng(SEED + N).normal(size=(N, 4))
        for aname, Fa in (('C1', Q.F_C1), ('C2', Q.F_C2)):
            sref = float(max(np.linalg.norm(Fa(Xg, 1), axis=1) / np.linalg.norm(Xg, axis=1)))
            for kappa in (1, -1):
                for tag, X0 in (('asis', X), ('rot', Xrot)):
                    res['runs'].append(run(X0, lambda Y: Fa(Y, kappa), sref, 20000, f'{name}_{aname}_k{kappa:+d}_{tag}', q_extra))
    for name, X, steps in (('ORTHO16', ortho16(), 20000), ('NULL28', null28(), 8000)):
        N = len(X)
        D = sum(c * M for c, M in zip(rng.normal(size=14), O.G2)); g = expm(D)
        Xrot = X @ g.T
        Xg = np.random.default_rng(SEED + N).normal(size=(N, 8))
        sref = float(max(np.linalg.norm(O.F_assoc(Xg, 1), axis=1) / np.linalg.norm(Xg, axis=1)))
        for kappa in (1, -1):
            for tag, X0 in (('asis', X), ('rot', Xrot)):
                res['runs'].append(run(X0, lambda Y: O.F_assoc(Y, kappa), sref, steps, f'{name}_ASSOC_k{kappa:+d}_{tag}', o_extra))
    json.dump(res, open(os.path.join(HERE, 'results', 'symmetric_inflation.json'), 'w'), ensure_ascii=False, indent=1)
