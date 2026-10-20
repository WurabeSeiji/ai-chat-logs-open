#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ASSOC の高速版（2026-10-20）。縮約の恒等式（../八元数が自己無撞着になる作用と初期値_20261020/assoc_contraction_identity_20261020.py で確認、誤差 ≤ 1.2×10⁻¹³）
  ASSOC_e = κ|x_e| Σ_f [w_e, w_f, S_w(w_e w_f)]、S_w = Σ_g w_g w_gᵀ、[a,b,c] = (ab)c − a(bc)
で三重和を二重和にし、全ての e を一度に計算する。同じ式の別の書き方なので値は丸めの桁で一致する（下の検証）が、
丸め誤差の入り方は元の実装（assoc_stage2_20261020.py の F_assoc）と違う。
検証：乱数の配置（N = 5, 16, 28）と高対称な配置（ORTHO16、NULL28、それらを G₂ で回したもの）で ‖F_fast − F_assoc‖ と時間。
"""
import importlib.util
import json
import os
import time

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('si', os.path.join(HERE, 'symmetric_inflation_20261020.py'))
si = importlib.util.module_from_spec(spec); spec.loader.exec_module(si)
O = si.O
C = O.sca.C                                            # 本体（octonion_selfconsistent_actions_20261020.py）の構造定数


def m3(x, y):
    return np.einsum('...a,...b,abc->...c', x, y, C)


def F_assoc_fast(X, kappa):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; S = W.T @ W
    P = m3(W[:, None, :], W[None, :, :])                  # w_e w_f        (N,N,8)
    Qm = P @ S.T                                           # S_w(w_e w_f)
    A = m3(P, Qm) - m3(W[:, None, :], m3(W[None, :, :], Qm))
    return kappa * r[:, None] * A.sum(1)


if __name__ == '__main__':
    rng = np.random.default_rng(20261020)
    out = {}
    cases = {f'RAND{N}': rng.normal(size=(N, 8)) for N in (5, 16, 28)}
    cases['ORTHO16'] = si.ortho16(); cases['NULL28'] = si.null28()
    D = sum(c * M for c, M in zip(rng.normal(size=14), O.G2)); g = expm(D)
    cases['ORTHO16_rot'] = cases['ORTHO16'] @ g.T; cases['NULL28_rot'] = cases['NULL28'] @ g.T
    for n, X in cases.items():
        t = time.time(); Fa = O.F_assoc(X, 1); ta = time.time() - t
        t = time.time(); Ff = F_assoc_fast(X, 1); tf = time.time() - t
        out[n] = dict(N=len(X), max_abs_diff=float(np.abs(Fa - Ff).max()), norm_F=float(np.linalg.norm(Fa)),
                      t_orig=ta, t_fast=tf)
        print(n, json.dumps(out[n]), flush=True)
    os.makedirs(os.path.join(HERE, 'results_long'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'results_long', 'assoc_fast_check.json'), 'w'), ensure_ascii=False, indent=1)
