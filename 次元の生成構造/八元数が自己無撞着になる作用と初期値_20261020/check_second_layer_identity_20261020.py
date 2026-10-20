#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude.ai の「第二層」の作用（../八元数のインフレーション試行_20261020/作用の定式化と探索空間_20261020.md §3、
A2：u_e = Σ_{f≠e} Re(w_ef) Im(w_ef)、w_ef = o_f ō_e、ȯ_e = u_e o_e）の照合（2026-10-20）。

解析：各項は o_e と o_f の二つだけで作られるので、Artin の定理で結合的。(o_f ō_e) o_e = |o_e|² o_f から
  u_e o_e = Σ_f ⟨o_f, o_e⟩ (|o_e|² o_f − ⟨o_f, o_e⟩ o_e) = |o_e|² S o_e − (o_eᵀ S o_e) o_e、 S = Σ_f o_f o_fᵀ
（f = e の項は 0 なので和に含めてよい）。したがって A2 は SO(8) 共変（G₂ に限らない）で八元数の非結合性を使わず、
**A2 が厳密に 0 ⟺ すべての要素が S の固有ベクトル**。定理 B（鏡映で閉じる）はこの十分条件の一つ。
確認：A2 を定義どおり omul で計算したものと右辺の差、各配置での max|A2| と「固有ベクトル条件の残差」の一致。
配置は claude.ai のプログラムと同じ作り方（ORTHO16、NULL28、正十六胞体 v/√3、24 胞体、E₈ の根、ランダム閉塞 N=9 は
同じ生成手順と乱数シード 3）に、本フォルダの N4・N5・ORTHO_AB・RAND8・RAND16 を加える。
"""
import itertools
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')
exec(open(SRC).read().split("results = {")[0])
sys.path.insert(0, HERE)
import importlib.util  # noqa: E402
spec = importlib.util.spec_from_file_location('sca', os.path.join(HERE, 'octonion_selfconsistent_actions_20261020.py'))
sca = importlib.util.module_from_spec(spec); spec.loader.exec_module(sca)


def Im(x):
    y = np.array(x, float); y[0] = 0.0; return y


def A2_velocity(S):
    N = len(S); V = np.zeros_like(S)
    for e in range(N):
        u = np.zeros(8)
        for f in range(N):
            if f != e:
                w = omul(S[f], oconj(S[e])); u += w[0] * Im(w)
        V[e] = omul(u, S[e])
    return V


def formula(S):
    M = S.T @ S
    return np.array([(o @ o) * (M @ o) - (o @ M @ o) * o for o in S])


def eig_resid(S):
    M = S.T @ S
    return float(max(np.linalg.norm(M @ o - (o @ M @ o) / (o @ o) * o) for o in S))


def claude_cfgs():
    c = {}
    c['ORTHO16'] = np.vstack([np.vstack([E[i], -E[i]]) for i in range(8)])
    nl = []
    for i in range(1, 8):
        for s in (1, -1):
            v = (E[0] + s * E[i]) / np.sqrt(2); nl.append(v); nl.append(-v)
    c['NULL28'] = np.array(nl)
    c16 = np.zeros((8, 8)); sc = 1 / np.sqrt(3); c16[0, 0] = 1; c16[1, 0] = -1
    for i in range(1, 4):
        c16[2 * i, i] = sc; c16[2 * i + 1, i] = -sc
    c['16-cell v/sqrt3'] = c16
    d24 = []
    for i in range(4):
        for s in (1, -1):
            v = np.zeros(8); v[i] = s; d24.append(v)
    for sg in itertools.product([0.5, -0.5], repeat=4):
        v = np.zeros(8); v[:4] = sg; d24.append(v)
    c['24-cell (in H)'] = np.array(d24)
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (1, -1):
                for sj in (1, -1):
                    v = np.zeros(8); v[i] = si; v[j] = sj; roots.append(v)
    for s in itertools.product([0.5, -0.5], repeat=8):
        if np.prod([1 if x > 0 else -1 for x in s]) > 0:
            roots.append(np.array(s))
    c['E8 roots 240'] = np.array(roots) / np.sqrt(2)

    def gen_closure(N, r, k=7):
        one = np.ones(N) / np.sqrt(N); q0 = r.normal(size=N); q0 -= one * (one @ q0)
        P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0) / (q0 @ q0)
        C = np.array([P @ r.normal(size=N) for _ in range(k)]); C *= np.sqrt((q0 @ q0) / np.sum(C * C))
        bs = []
        for i in range(k):
            v = r.normal(size=8); v[0] = 0
            for b in bs:
                v -= (b @ v) * b
            bs.append(v / np.linalg.norm(v))
        return np.outer(q0, E[0]) + C.T @ np.array(bs)
    c['random closure N=9'] = gen_closure(9, np.random.default_rng(3))
    return c


if __name__ == '__main__':
    cfg = claude_cfgs()
    cfg.update(sca.configs(np.random.default_rng(sca.SEED)))   # 本フォルダの初期値（同じ乱数の流れの先頭から）
    out = {}
    for n, S in cfg.items():
        if len(S) > 60:
            V = formula(S); diff = None
        else:
            V = A2_velocity(S); diff = float(np.linalg.norm(V - formula(S)))
        out[n] = dict(N=len(S), A2_minus_formula_abs=diff, max_A2=float(max(np.linalg.norm(v) for v in V)),
                      eigvec_resid=eig_resid(S))
        print(n, json.dumps(out[n]), flush=True)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'second_layer_identity.json'), 'w'), ensure_ascii=False, indent=1)
