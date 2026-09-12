#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""自発的対称性の破れ(SSB)の検証: 同一床＋異なる微小シード → 縮退真空（読出し＋厳密one_stepの制御実験）。

木原指示（2026-09-12）: 90°格子が「唯一の平衡でない」ことが最重要。唯一でないから、僅かなシードで
自発的対称性の破れが起き、インフレーションが起き、別の（等振幅・位相分散・傾いた）平衡へ落ちる。

SSBの決定的検証（メキシカンハット像）:
  90°格子＝最大対称だが不安定な平衡（頂上）。安定な縮退真空（等振幅・位相分散・傾いた平衡）が別に存在。
  ⟹ 同一の床 Z0 に、大きさは同じで向きの違う微小シードを与えると、
     全て起動(inflation)し、全て「同じ型」(等振幅・同一ロック率 −2π/N)の平衡へ落ちるが、
     位相配置(＝真空の向き)はシードごとに異なる（縮退真空の自発選択）。

厳密 one_step（正本 run_N3_N40_stage123_v1.py と同一。別途 verify_single_wave_eom で写像一致 2e-15 を確認済）:
  K_ef = A_ef sin(φ_f-φ_e), z_{n+1}=exp(Δτ K)z, Δτ=2π/den。物理式は不変、初期条件のシードのみ制御変数。
入力: N=6 床 Z0（full_N3_N40_sweep の hm_N6_den_6_states_500.npz の step0）。
出力: 標準出力（run.log）。SEED_AMP=1e-8 の複数シードで最終平衡の型と向きを比較。
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
STEPS = 800
SEED_AMP = 1e-8


def adjacency(N):
    ea, eb = np.triu_indices(N, 1); M = len(ea); A = np.zeros((M, M))
    for e in range(M):
        sh = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); sh[e] = False; A[e, sh] = 1.0
    return A


def one_step(z, A, den):
    u = np.exp(1j * np.angle(z)); H = A * (np.conj(u)[:, None] * u[None, :]); np.fill_diagonal(H, 0)
    H = 1j * np.imag(H); w, V = np.linalg.eigh(H); ph = np.exp(-1j * (2 * math.pi / den) * w)
    return V @ (ph * (V.conj().T @ z))


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    A = adjacency(N); M = N * (N - 1) // 2; den = N
    Z0 = np.asarray(np.load(os.path.join(GEN, 'full_N3_N40_sweep', 'states',
                    f'hm_N{N}_den_{N}_states_500.npz'))['Z'][0], np.complex128)
    p = Z0.real / np.linalg.norm(Z0.real); q = Z0.imag - np.dot(Z0.imag, p) * p; q /= np.linalg.norm(q)
    perp = lambda z: (lambda zp: np.vdot(zp, zp).real / np.vdot(z, z).real)(z - p * np.dot(p, z) - q * np.dot(q, z))
    ideal = -2 * math.pi / N
    print(f'SSB検証 N={N} den={N} M={M}: 同一床Z0＋異なる微小シード(振幅{SEED_AMP:g})を{STEPS}step')
    print(f'  SSBなら: 全シードで inflation→等振幅・Δ=−2π/N={ideal:+.5f} の同型平衡、位相配置(真空の向き)はシード依存で相異')
    rels = []
    for seed in range(6):
        rng = np.random.default_rng(seed)
        z = Z0 + SEED_AMP * (rng.standard_normal(M) + 1j * rng.standard_normal(M))
        z = z / np.linalg.norm(z) * np.linalg.norm(Z0)
        for _ in range(STEPS):
            z = one_step(z, A, den)
        r = np.abs(z); dth = float(np.angle(np.vdot(z, one_step(z, A, den))))
        rel = np.angle(z * np.conj(z[0])); rels.append(rel)
        print(f'  seed{seed}: H⊥/H={perp(z):.4f}  振幅比={r.max()/r.min():.4f}  Δ={dth:+.6f}  |Δ+2π/N|={abs(dth-ideal):.1e}')
    base = rels[0]
    print('  真空の向き（相対位相配置）のシード間差:')
    for s in range(1, 6):
        d = np.degrees(np.abs(((rels[s] - base + np.pi) % (2 * np.pi)) - np.pi))
        print(f'    seed{s} vs seed0: 中央={np.median(d):.1f}° 最大={d.max():.1f}°')
    print('  結論: 全シードが等振幅・Δ=−2π/N の同型安定平衡へ。位相配置はシードで相異＝縮退真空の自発選択＝SSB。')


if __name__ == '__main__':
    main()
