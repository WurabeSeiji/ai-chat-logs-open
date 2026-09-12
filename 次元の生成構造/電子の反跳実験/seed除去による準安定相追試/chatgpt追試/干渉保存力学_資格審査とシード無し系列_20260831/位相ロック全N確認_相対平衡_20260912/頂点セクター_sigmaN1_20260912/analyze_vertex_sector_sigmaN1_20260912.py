#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""追試C: 辺空間の3分解 1⊕(N-1)⊕N(N-3)/2（線グラフ/Johnsonスキーム）と、終端 σ=N-1 との対応を検証
（読出しのみ・物理無変更・新規走行なし。既存 10000step 終端状態を使用）。

木原指示（2026-09-12）: 論文8 の「頂点セクター（N-1次元）と σ=N-1 定理の一致は要検証仮説」を既存データで検定。

検証項目（N=3..40）:
  (S1) 次元検算: 1 + (N-1) + N(N-3)/2 = M。
  (S2) 線グラフ L(K_N) の隣接固有値が {2(N-2)[×1], N-4[×(N-1)], -2[×N(N-3)/2]} と一致（3分解の実在）。
  (S3) 頂点セクター射影: 接続行列 B(N×M, B[v,e]=1 if v∈e) から
       P_row = Bᵀ(BBᵀ)⁻¹B（頂点由来 N 次元＝trivial⊕頂点セクター）、
       P_vert = P_row − (1/M)11ᵀ（trivial を除いた (N-1) 次元頂点セクター）。
  (S4) 終端 iK の最大固有値 σ_max（→ N-1 のはず、実験9-1）と、終端状態 z の頂点セクター内割合
       frac_vert = ‖P_vert Re z‖²+‖P_vert Im z‖²）/‖z‖²、trivial 割合、残差セクター割合。
判定: σ_max≈N-1 かつ z が頂点セクターに集中するなら「σ=N-1 ↔ 頂点セクター」を支持。集中しなければ次元一致は偶合。
正本 run_N3_N40_stage123_v1.py（SHA照合）から adjacency/one_step を import。
"""
import csv
import hashlib
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.abspath(os.path.join(HERE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py'))
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
LONG = os.path.abspath(os.path.join(HERE, '..', '..', 'N3_N40_long10000_20260905', 'results'))

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
adjacency, one_step = g['adjacency'], g['one_step']


def incidence(N):
    ea, eb = np.triu_indices(N, 1); M = len(ea); B = np.zeros((N, M))
    B[ea, np.arange(M)] = 1.0; B[eb, np.arange(M)] = 1.0
    return B


def main():
    rows = []
    for N in range(3, 41):
        M = N * (N - 1) // 2
        A = adjacency(N)
        # (S1)(S2) 静的3分解
        dim_ok = (1 + (N - 1) + N * (N - 3) // 2 == M)
        evA = np.round(np.sort(np.linalg.eigvalsh(A))[::-1], 6)
        # 期待固有値と多重度
        exp = {round(2 * (N - 2), 6): 1, round(N - 4, 6): N - 1, round(-2, 6): N * (N - 3) // 2}
        from collections import Counter
        cnt = Counter(np.round(evA, 4).tolist())
        spec_ok = all(abs(cnt.get(round(k, 4), 0) - v) == 0 for k, v in exp.items())
        # (S3) 頂点セクター射影
        B = incidence(N)
        Prow = B.T @ np.linalg.inv(B @ B.T) @ B
        Ptriv = np.ones((M, M)) / M
        Pvert = Prow - Ptriv
        # (S4) 終端状態と σ_max
        st = os.path.join(LONG, f'hm_N{N}_den_{N}_states_10000.npz')
        if not os.path.exists(st):
            continue
        z = np.asarray(np.load(st)['Z'][-1], np.complex128)
        u = np.exp(1j * np.angle(z)); H = A * (np.conj(u)[:, None] * u[None, :]); np.fill_diagonal(H, 0)
        H = 1j * np.imag(H); w = np.linalg.eigvalsh(H); sigma_max = float(np.max(np.abs(w)))
        nz = np.vdot(z, z).real
        frac_vert = float((np.vdot(Pvert @ z.real, Pvert @ z.real).real + np.vdot(Pvert @ z.imag, Pvert @ z.imag).real) / nz)
        frac_triv = float((np.vdot(Ptriv @ z.real, Ptriv @ z.real).real + np.vdot(Ptriv @ z.imag, Ptriv @ z.imag).real) / nz)
        frac_rest = 1.0 - frac_vert - frac_triv
        rows.append(dict(N=N, M=M, dim_ok=dim_ok, spec_ok=spec_ok, sigma_max=sigma_max, target_Nm1=N - 1,
                         frac_trivial=frac_triv, frac_vertex=frac_vert, frac_rest=frac_rest))
        print(f'N={N:>2} M={M:>4}: 3分解次元={dim_ok} L(K_N)固有値={spec_ok} | σ_max={sigma_max:.4f}(→N-1={N-1}) '
              f'| z成分割合 trivial={frac_triv:.4f} 頂点セクター={frac_vert:.4f} 残余={frac_rest:.4f}')
    with open(os.path.join(HERE, 'results', 'vertex_sector_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    fv = np.array([r['frac_vertex'] for r in rows])
    print(f'\n頂点セクター割合 frac_vertex: 範囲 {fv.min():.4f}–{fv.max():.4f}, 中央 {np.median(fv):.4f}')
    print('判定: σ_max→N-1 は成立。z が頂点セクターに集中すれば「σ=N-1↔頂点セクター」支持、'
          '残余セクターに広がれば次元ラベルN-1の一致は偶合（動力学の固有値と静的既約分解は別物）。')


if __name__ == '__main__':
    main()
