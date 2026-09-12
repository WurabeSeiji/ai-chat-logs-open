#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""単一波の運動方程式の導出の数値検証（読出しのみ・物理無変更）。

木原指示（2026-09-12）: 「なぜ1つの波が今回の軌道を取るのか」が物理。剛体回転を除いた τ発展方向の
位相の変位と振幅の変化を説明できる定式化がポイント。

正本 one_step（../../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py）:
  H = H_of(exp(i·angle(z)), A) = A ∘ (conj(u) u),  u=exp(i·angle(z))  → H_ef = A_ef e^{i(φ_f-φ_e)}
  H ← i·Im(H)  →  H_ef = i·A_ef sin(φ_f-φ_e) = i·K_ef,  K 実反対称
  z_{n+1} = exp(-iΔτ·H) z = exp(Δτ·K) z,  Δτ = 2π/den   （実直交回転）

導出（本スクリプトで検証する主張）:
  生成子   K_ef = A_ef sin(φ_f-φ_e)          （位相のみ。振幅は生成子に入らない）
  単一波   dz_e/dτ = (K z)_e = Σ_{f~e} sin(φ_f-φ_e) z_f
  z_e=r_e e^{iφ_e} と分解し RHS=Σ A_ef sin(φ_f-φ_e) r_f e^{i(φ_f-φ_e)} を実虚に分けて:
    振幅  dr_e/dτ      = (1/2) Σ_{f~e} r_f sin(2(φ_f-φ_e))
    位相  r_e dφ_e/dτ  =        Σ_{f~e} r_f sin^2(φ_f-φ_e)     （sin^2≥0 → 位相は単調前進＝回転）
  床（位相が 0/90/180/270 の90°格子）では全ての位相差が90°の倍数 → sin(2Δφ)=0 → dr_e/dτ=0（相対平衡）。

検証（機械精度で一致すれば導出は厳密）:
  (T1) 厳密写像: exp(Δτ K(φ_n)) z_n が正本軌道 z_{n+1} を再現。
  (T2) 生成子分解: 微小サブステップ (exp(εK)z - z)/ε の振幅/位相成分が上記EOMと一致。
  (T3) 床の振幅凍結: 90°格子で dr_e/dτ ≈ 0。
入力: den=N・状態npz。既定 N=6（500step）と N=40（2000step, 本系列で走行）。
"""
import math
import os
import sys

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
PL = os.path.abspath(os.path.join(HERE, '..'))


def adjacency(N):
    ea, eb = np.triu_indices(N, 1); M = len(ea); A = np.zeros((M, M))
    for e in range(M):
        sh = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); sh[e] = False; A[e, sh] = 1.0
    return A


def K_of(z, A):
    ph = np.angle(z); d = ph[None, :] - ph[:, None]
    return A * np.sin(d)                       # K_ef = A_ef sin(φ_f-φ_e), 実反対称


def amp_eom(z, A):
    r = np.abs(z); ph = np.angle(z); d = ph[None, :] - ph[:, None]
    return 0.5 * np.sum(A * r[None, :] * np.sin(2 * d), axis=1)


def phase_eom(z, A):
    r = np.abs(z); ph = np.angle(z); d = ph[None, :] - ph[:, None]
    return np.sum(A * r[None, :] * np.sin(d) ** 2, axis=1)   # = r_e dφ_e/dτ


def path(N):
    if N == 6:
        return os.path.join(GEN, 'full_N3_N40_sweep', 'states', 'hm_N6_den_6_states_500.npz')
    if N == 40:
        return os.path.join(PL, 'N40_2000step検証走行_20260912', 'results', 'hm_N40_den_40_states_2000.npz')
    if N == 22:
        return os.path.join(PL, 'N22_2000step検証走行_20260912', 'results', 'hm_N22_den_22_states_2000.npz')
    return os.path.join(GEN, 'full_N3_N40_sweep', 'states', f'hm_N{N}_den_{N}_states_500.npz')


def run(N):
    A = adjacency(N); dtau = 2 * math.pi / N
    S = np.asarray(np.load(path(N))['Z'], np.complex128); T = len(S)
    # (T1) 厳密写像の再現
    err = [np.max(np.abs(expm(dtau * K_of(S[t], A)) @ S[t] - S[t + 1])) for t in range(0, min(T - 1, 300))]
    print(f'  [T1] N={N}: 厳密写像 exp(Δτ·K)z の軌道再現誤差 最大={max(err):.2e} 中央={np.median(err):.2e}')
    # (T2) 生成子分解を微小サブステップで
    eps = 1e-6; ea = ep = 0.0
    for t in np.linspace(1, T - 2, 8).astype(int):
        z = S[t]; zk = expm(eps * K_of(z, A)) @ z; r = np.abs(z)
        dr = (np.abs(zk) - r) / eps
        dphi = np.angle(zk * np.conj(z)) / eps
        ea = max(ea, np.max(np.abs(dr - amp_eom(z, A))))
        ep = max(ep, np.max(np.abs(r * dphi - phase_eom(z, A))))
    print(f'  [T2] N={N}: 振幅EOM最大誤差={ea:.2e}  位相EOM最大誤差={ep:.2e}')
    # (T3) 床の振幅凍結
    ph0 = sorted(set(np.round(np.degrees(np.angle(S[0])) % 90, 1)))
    print(f'  [T3] N={N}: 床の位相 mod90={ph0}（0付近＝90°格子）  床 |dr/dτ|最大={np.max(np.abs(amp_eom(S[0], A))):.2e}')


def main():
    Ns = [int(x) for x in sys.argv[1:]] or [6, 40]
    print('単一波EOMの検証（機械精度で一致すれば導出は厳密）:')
    for N in Ns:
        if os.path.exists(path(N)):
            run(N)
        else:
            print(f'  N={N}: npz なし skip')


if __name__ == '__main__':
    main()
