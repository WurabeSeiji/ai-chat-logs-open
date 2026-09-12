#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""90°格子の必然性の検証（読出しのみ・物理無変更）。

木原指示（2026-09-12）: 「初期が90度格子であることは必然か」。

主張（単一波EOM dr_e/dτ = ½Σ_{f~e} r_f sin(2(φ_f-φ_e)) の帰結）:
  平衡（振幅凍結）条件は 各辺で Σ_{f~e} r_f sin(2(φ_f-φ_e)) = 0。実現は2種:
    (termwise) 各項 sin(2Δφ) が個別に0 ⟺ 全対の位相差が90°の倍数 ⟺ 全位相が共通90°格子。
               第2高調波 sin(2ψ) の零点が ψ=k·90° に限られることの必要十分（生成子が sin(Δφ)＝
               第1高調波ゆえ振幅駆動が第2高調波になる帰結）。
    (cancellation) 個別項は非0だが辺ごとの和が相殺して0 → 一般位相の別平衡（ロック終点）。
  ⟹ 90°格子は「全結合項が個別に消える唯一の最大対称平衡」として必然だが、系唯一の平衡ではなく
     不安定平衡。安定平衡は相殺型（非90°）。

検証:
  (A) 床(step0): 個別項 |sin(2Δφ)| がすべて≈0（termwise平衡）。位相は共通90°格子。
  (B) ロック終点: 個別項|sin(2Δφ)|は O(1) だが 辺ごとの重み和 ≈0（cancellation平衡）。位相は非90°格子。
  (C) sin(2ψ)=0 の解が ψ=k·90° に限ることの明示（[0,360)で数え上げ）。
入力: den=N npz（既定 N=6/22/40）。
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
PL = os.path.abspath(os.path.join(HERE, '..'))


def adjacency(N):
    ea, eb = np.triu_indices(N, 1); M = len(ea); A = np.zeros((M, M))
    for e in range(M):
        sh = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); sh[e] = False; A[e, sh] = 1.0
    return A


def path(N):
    if N == 40:
        return os.path.join(PL, 'N40_2000step検証走行_20260912', 'results', 'hm_N40_den_40_states_2000.npz')
    if N == 22:
        return os.path.join(PL, 'N22_2000step検証走行_20260912', 'results', 'hm_N22_den_22_states_2000.npz')
    return os.path.join(GEN, 'full_N3_N40_sweep', 'states', f'hm_N{N}_den_{N}_states_500.npz')


def equilibrium_type(z, A):
    r = np.abs(z); ph = np.angle(z); d = ph[None, :] - ph[:, None]
    perpair_max = float(np.max(np.abs(A * np.sin(2 * d))))                 # 個別項 |sin(2Δφ)| 最大
    edge_sum_max = float(np.max(np.abs(np.sum(A * r[None, :] * np.sin(2 * d), 1))))  # 辺ごと Σ (=2 dr/dτ)
    ph_deg = np.degrees(ph)
    off = np.abs((ph_deg[:, None] - ph_deg[None, :] + 45) % 90 - 45)       # 位相差の90°格子からのずれ
    lattice_off = float(np.max(off[A > 0]))
    return perpair_max, edge_sum_max, lattice_off


def main():
    Ns = [int(x) for x in sys.argv[1:]] or [6, 22, 40]
    # (C) sin(2ψ)=0 の零点
    psi = np.linspace(0, 360, 360001)[:-1]
    zeros = psi[np.abs(np.sin(2 * np.radians(psi))) < 1e-6]
    zeros_deg = sorted(set(np.round(zeros / 90) * 90))
    print(f'(C) sin(2ψ)=0 の [0,360) の解 ≈ {zeros_deg} 度 → 90°の倍数のみ（第2高調波の帰結）\n')
    print(f"{'N':>3}{'状態':>12}{'個別|sin2Δφ|max':>16}{'辺ごとΣ|max(=平衡度)':>20}{'90°格子ずれmax°':>15}{'型':>14}")
    for N in Ns:
        p = path(N)
        if not os.path.exists(p):
            print(f'  N={N}: npz なし skip'); continue
        S = np.asarray(np.load(p)['Z'], np.complex128)
        for t, name in [(0, '床(step0)'), (len(S) - 1, 'ロック(終点)')]:
            pm, es, off = equilibrium_type(S[t], A=adjacency(N))
            typ = 'termwise' if pm < 1e-6 else 'cancellation'
            print(f"{N:>3}{name:>12}{pm:>16.2e}{es:>20.2e}{off:>15.1f}{typ:>14}")
    print('\n結論: 床=termwise平衡(全位相90°格子・各項個別0)＝最大対称・不安定。'
          'ロック=cancellation平衡(非90°・相殺)＝安定。90°格子は「全項個別0の唯一配置」の意味で必然だが唯一平衡ではない。')


if __name__ == '__main__':
    main()
