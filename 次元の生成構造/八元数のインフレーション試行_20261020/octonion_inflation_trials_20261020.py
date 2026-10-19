#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""高対称な八元数の配置でインフレーションが起こるかの試行錯誤（2026-10-20）。

公理：基本量は八元数（代数は公理、feedback_algebra_is_axiom_not_derived）。答えは未知なので、組み合わせを試す。
代数：八元数の積は 次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py の定義をそのまま読み込む。

相互作用の相手は全員対全員（A = J − I、自分以外すべて）。

生成子の候補（条件：大きさを保つ／八元数の積と左右の掛け算だけ／一つの複素平面 span(1, ℓ) に収まる配置では
過去の複素数の力学 K_ef = Im(z̄_e z_f) にそのまま戻る）：
  L_amp : K_ef = A_ef ⟨ℓ o_e, o_f⟩（左から ℓ）、Q ← exp(Δτ K) Q（番号の側の実直交回転、8 成分に一様）
  L_ph  : 同じで o を o/|o| に置き換えた位相のみ版（stage123 と同じ型）
  R_amp : K_ef = A_ef ⟨o_e ℓ, o_f⟩（右から ℓ）
  R_ph  : その位相のみ版
  ELEM  : 要素ごとの左の掛け算 o_e ← exp(Δτ Σ_f A_ef Im(o_f ō_e)/|o_e|²) o_e（各 |o_e| を保つ。複素数では蔵本型の位相方程式になり、
          過去の力学には一致しない。比較のために入れる）
  ℓ = e₄（八元数フォルダの規約）。
初期配置：
  ORTHO16 : 8 次元の正軸体 {±1, ±e₁, …, ±e₇}（等方性の閉塞：Σx = 0、グラム = 2I）
  NULL28  : ±(1 + n)/√2、n = ±e₁…±e₇（本当の二乗の閉塞 Σo = Σo² = 0、各要素が 8 次元ヌル、o⁸ = 1）
  CROSS4  : 複素数の十字 {±1, ±ℓ}（対照。過去の複素数の系に相当）
測る量：
  relEq    : 相対平衡の残差 ‖Ω + Ωᵀ‖/‖Ω‖（2026-10-20 追加。selfcons だけでは張りの中の形の崩れを見逃した：NULL28 の位相のみ版が
             数ステップで 27 桁跳んだ原因）
  g(t)     : 要素どうしのグラム行列 QQᵀ の初期値からの相対変化（形の変化。全体の回転では変わらない。2026-10-20 追加）
  selfcons : 自己無撞着の残差。番号の側の力学では ‖KQ − QΩ‖/‖KQ‖（Ω = Q⁺KQ、状態が自分の列の張りの中で回るだけか）。
             ELEM では 1 ステップ後の形の変化（全体の左右の回転を除かない素朴な ‖ΔQ‖ ではなく、グラム行列 QQᵀ の変化）
  f(t)     : 初期の列の張り（R^M の部分空間）の外へ出た割合 ‖(I − P₀)Q‖²/‖Q‖²（インフレーションの指標、H⊥/H の一般化）
             ELEM では要素どうしのグラム行列 G = QQᵀ の初期値からの相対変化
  刻み Δτ = 2π/(200σ)、σ = 初期の K の最大固有値の絶対値（ELEM は初期の一歩の回転の最大）。
"""
import json
import math
import os

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '八元数閉塞と関係の分解_20261019', 'octonion_closure_check_20261019.py')
exec(open(SRC).read().split("results = {")[0])   # omul, oconj, E などを読み込む

ELL = E[4]


def Lmat(a):
    return np.array([omul(a, E[j]) for j in range(8)]).T      # x ↦ a x


def Rmat(a):
    return np.array([omul(E[j], a) for j in range(8)]).T      # x ↦ x a


LL, RL = Lmat(ELL), Rmat(ELL)


def configs():
    c = {}
    c['ORTHO16'] = np.array([s * E[a] for a in range(8) for s in (1, -1)])
    c['NULL28'] = np.array([s * (E[0] + t * E[a]) / math.sqrt(2) for a in range(1, 8) for t in (1, -1) for s in (1, -1)])
    c['CROSS4'] = np.array([E[0], -E[0], ELL, -ELL])
    return c


def kernel(Q, kind):
    M = len(Q); A = np.ones((M, M)) - np.eye(M)
    U = Q / np.linalg.norm(Q, axis=1, keepdims=True) if kind.endswith('_ph') else Q
    T = U @ (LL.T if kind.startswith('L') else RL.T)          # 各行 = ℓ u_e（左）または u_e ℓ（右）
    return A * (T @ U.T)                                        # K_ef = ⟨ℓ u_e, u_f⟩


def step_uniform(Q, kind, dt):
    return expm(dt * kernel(Q, kind)) @ Q


def elem_gen(Q):
    M = len(Q); out = []
    for e in range(M):
        g = np.zeros(8)
        for f in range(M):
            if f != e:
                g += omul(Q[f], oconj(Q[e]))
        g[0] = 0.0
        out.append(g / np.dot(Q[e], Q[e]))
    return np.array(out)


def step_elem(Q, dt):
    G = elem_gen(Q)
    return np.array([expm(dt * Lmat(G[e])) @ Q[e] for e in range(len(Q))])


def trial(name, Q0, kind, steps=4000):
    Q0 = Q0.astype(np.float64)
    if kind == 'ELEM':
        G0 = elem_gen(Q0); sigma = float(np.max(np.linalg.norm(G0, axis=1)))
        if sigma < 1e-12:
            return dict(config=name, gen=kind, sigma=sigma, note='生成子がゼロ（全要素で関係の虚部の和が 0）＝ 静止', steps=0)
        dt = 2 * math.pi / (200 * sigma)
        Gram0 = Q0 @ Q0.T; Q = Q0.copy(); fs = []
        for t in range(steps):
            Q = step_elem(Q, dt)
            fs.append(float(np.linalg.norm(Q @ Q.T - Gram0) / np.linalg.norm(Gram0)))
        gs = fs
        Q1 = step_elem(Q0, dt)
        return dict(config=name, gen=kind, sigma=sigma, dt=dt, steps=steps,
                    selfcons=float(np.linalg.norm(Q1 @ Q1.T - Gram0) / np.linalg.norm(Gram0)), f=fs, g=gs)
    K0 = kernel(Q0, kind); sigma = float(np.max(np.abs(np.linalg.eigvals(K0))))
    KQ = K0 @ Q0
    if np.linalg.norm(KQ) < 1e-14:
        return dict(config=name, gen=kind, sigma=sigma, note='KQ = 0（生成子が状態に作用しない）＝ 静止', steps=0)
    Om = np.linalg.pinv(Q0) @ KQ
    selfcons = float(np.linalg.norm(KQ - Q0 @ Om) / np.linalg.norm(KQ))
    relEq = float(np.linalg.norm(Om + Om.T) / np.linalg.norm(Om))   # 相対平衡なら Ω は反対称（0）
    dt = 2 * math.pi / (200 * sigma)
    U, s, _ = np.linalg.svd(Q0, full_matrices=False); r = s > 1e-12 * s[0]; P0 = U[:, r] @ U[:, r].T
    Gram0 = Q0 @ Q0.T
    Q = Q0.copy(); fs = []; gs = []
    for t in range(steps):
        Q = step_uniform(Q, kind, dt)
        Qp = Q - P0 @ Q
        fs.append(float(np.sum(Qp * Qp) / np.sum(Q * Q)))
        gs.append(float(np.linalg.norm(Q @ Q.T - Gram0) / np.linalg.norm(Gram0)))
    return dict(config=name, gen=kind, sigma=sigma, dt=dt, steps=steps, selfcons=selfcons, relEq=relEq,
                rank0=int(r.sum()), f=fs, g=gs)


def summarize(r):
    if not r.get('f'):
        return r
    f = np.maximum(np.array(r['f']), 1e-40); lg = np.log10(f)
    out = {k: v for k, v in r.items() if k not in ('f', 'g')}
    g = np.maximum(np.array(r['g']), 1e-40); lgg = np.log10(g)
    out.update(g_step1=float(g[0]), g_max=float(g.max()), log10_g_at=[float(lgg[i]) for i in (0, 9, 99, 499, 999, 1999, len(lgg) - 1)])
    out.update(f_step1=float(f[0]), f_max=float(f.max()), f_final=float(f[-1]),
               log10_f_at=[float(lg[i]) for i in (0, 9, 99, 499, 999, 1999, len(lg) - 1)])
    i1 = np.where(lg > lg[0] + 3)[0]
    out['steps_to_grow_3_decades'] = int(i1[0]) + 1 if len(i1) else -1
    if len(i1):
        a, b = 0, int(i1[0])
        out['mean_growth_dec_per_step_initial'] = float((lg[b] - lg[a]) / max(b - a, 1))
    return out


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    rows, raw = [], {}
    for cname, Q0 in configs().items():
        for kind in ('L_amp', 'L_ph', 'R_amp', 'R_ph', 'ELEM'):
            r = trial(cname, Q0, kind)
            raw[f'{cname}_{kind}_f'] = r.get('f', []); raw[f'{cname}_{kind}_g'] = r.get('g', [])
            s = summarize(r); rows.append(s)
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in s.items()}, ensure_ascii=False), flush=True)
    json.dump(rows, open(os.path.join(HERE, 'results', 'trials_summary.json'), 'w'), ensure_ascii=False, indent=1)
    np.savez_compressed(os.path.join(HERE, 'results', 'trials_f_series.npz'), **{k: np.array(v) for k, v in raw.items()})
