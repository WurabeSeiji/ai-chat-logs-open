#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ASSOC の縮約の恒等式と、ORTHO_AB・N4 が静止する理由（2026-10-20）。
claude.ai の `定式化の突き合わせと走行の直し方_20261020.md` §4 の未解明の問い（ORTHO_AB がなぜ ASSOC で静止するか）への解析と確認。

解析：重み ⟨w_e w_f, w_g⟩ は g について線形なので、g の和は二次のモーメント S_w = Σ_g w_g w_gᵀ（単位ベクトルの S）に縮約できる：
  Σ_g ⟨w_e w_f, w_g⟩ [w_e, w_f, w_g] = [w_e, w_f, S_w (w_e w_f)]
よって ASSOC_e = κ|x_e| Σ_f [w_e, w_f, S_w(w_e w_f)]。
  ・S_w ∝ I なら [w_e, w_f, w_e w_f] = 0（Artin）で消える（2-デザイン、第 1 段の解析）。
  ・より一般に、各 w_e w_f が S_w の固有ベクトル（あるいは S_w(w_e w_f) が w_e, w_f の生成する部分代数に入る）なら消える。
    ORTHO_AB・N4 は要素が基底 ±e_a の上にあり、S_w が対角、w_e w_f も ±e_c なので固有ベクトル → 静止。
確認：(1) 縮約の式と F_assoc の差（N4, N5, NULL28, ORTHO_AB, RAND8, RAND16）。
      (2) 各配置で「全ての w_e w_f が S_w の固有ベクトルか」の残差と max|F_assoc|。
"""
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('sca', os.path.join(HERE, 'octonion_selfconsistent_actions_20261020.py'))
sca = importlib.util.module_from_spec(spec); spec.loader.exec_module(sca)
mul, F_assoc = sca.mul, sca.F_assoc


def contracted(X, kappa=1):
    r = np.linalg.norm(X, axis=1); W = X / r[:, None]; Sw = W.T @ W
    out = np.zeros_like(X)
    for e in range(len(X)):
        P = mul(W[e][None, :], W)                   # w_e w_f
        Q = P @ Sw.T                                 # S_w (w_e w_f)
        A = mul(mul(W[e][None, :], W), Q) - mul(W[e][None, :], mul(W, Q))
        out[e] = A.sum(0) * r[e]
    return kappa * out


def eig_products_resid(X):
    W = X / np.linalg.norm(X, axis=1, keepdims=True); Sw = W.T @ W; m = 0.0
    for e in range(len(W)):
        for f in range(len(W)):
            p = mul(W[e], W[f]); n = p @ p
            if n > 1e-24:
                m = max(m, float(np.linalg.norm(Sw @ p - (p @ Sw @ p) / n * p)))
    return m


if __name__ == '__main__':
    cfg = sca.configs(np.random.default_rng(sca.SEED))
    out = {}
    for n, X in cfg.items():
        F = F_assoc(X, 1)
        out[n] = dict(contracted_minus_F=float(np.linalg.norm(contracted(X) - F)), max_F=float(max(np.linalg.norm(F, axis=1))),
                      products_eigvec_resid=eig_products_resid(X))
        print(n, json.dumps(out[n]), flush=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'assoc_contraction_identity.json'), 'w'), ensure_ascii=False, indent=1)
