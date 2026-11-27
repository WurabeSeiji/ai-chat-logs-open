#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""90度理論床の振幅・σの代数的閉形式の検証（2026-09-09、木原指示）。

整数 K（成分0,±1）の σ_max 固有モードとして解いた90度床（make_90deg_floor_analytic_v1.py）の
σ² と振幅² が、以下の閉形式に厳密一致することを sympy 有理演算で確認する。

σ²（＝iK の最大固有値² ＝ 整数 K² の固有値の符号反転）:
  偶数N: σ² = N(N−2)
  奇数N: σ² = N²−2N−1

振幅²（規格化 ‖Z‖²=1、値×本数）:
  偶数N（2クラス, 各クラス総パワー 1/2）:
    2/N²                       × (N/2)² 本
    2/(N(N−2))                 × N(N−2)/4 本
  奇数N（3クラス, N=3 は最終クラス0本）:
    2/((N−1)(N+1))             × (N²−1)/4 本
    2(N−1)/((N+1)(N²−2N−1))    × (N²−1)/8 本
    2(N+1)/((N−1)(N²−2N−1))    × (N−1)(N−3)/8 本

振幅そのものは √（有理数）、σ は √（整数）。|z_e|² は厳密に有理数
（整数 K² の固有ベクトル u は有理、z=σu+iKu ゆえ |z_e|²=σ²u_e²+(Ku)_e² が有理）。
"""
import os
from collections import Counter

import numpy as np
import sympy as sp

BASE = os.path.dirname(os.path.abspath(__file__))
DST = os.path.join(BASE, 'parents_90deg_floor_analytic')


def integer_K(labels, N):
    ea, eb = np.triu_indices(N, k=1)
    M = len(ea)
    A = np.zeros((M, M), dtype=int)
    for e in range(M):
        s = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        s[e] = False
        A[e, s] = 1
    K = sp.zeros(M, M)
    for e in range(M):
        for f in range(M):
            if A[e, f]:
                dfe = (int(labels[f]) - int(labels[e])) % 4
                K[e, f] = 1 if dfe == 1 else (-1 if dfe == 3 else 0)
    return K, M


def closed_form(N):
    if N % 2 == 0:
        sig2 = N * (N - 2)
        classes = {sp.Rational(2, N * N): (N // 2) ** 2,
                   sp.Rational(2, N * (N - 2)): N * (N - 2) // 4}
    else:
        sig2 = N * N - 2 * N - 1
        classes = {sp.Rational(2, (N - 1) * (N + 1)): (N * N - 1) // 4,
                   sp.Rational(2 * (N - 1), (N + 1) * sig2): (N * N - 1) // 8,
                   sp.Rational(2 * (N + 1), (N - 1) * sig2): (N - 1) * (N - 3) // 8}
        classes = {k: v for k, v in classes.items() if v > 0}
    return sig2, classes


def main():
    ok_all = True
    for N in range(3, 21):
        d = np.load(os.path.join(DST, f'parent_90deg_floor_N{N:05d}_analytic.npz'))
        labels = np.array(d['family'], dtype=int)
        K, M = integer_K(labels, N)
        sig2_cf, cls_cf = closed_form(N)
        # 厳密: 整数 K² の固有値 -σ² のヌル空間から振幅² を有理で作る
        u = (K * K + sig2_cf * sp.eye(M)).nullspace()[0]
        Ku = K * u
        amp2 = [sp.nsimplify(sig2_cf * u[e] ** 2 + Ku[e] ** 2) for e in range(M)]
        tot = sum(amp2)
        cnt = Counter(sp.nsimplify(a / tot) for a in amp2)
        eig_ok = (K * K * u == -sig2_cf * u)
        cf_ok = (dict(cnt) == cls_cf)
        pow_ok = (sum(k * v for k, v in cnt.items()) == 1)
        ok_all = ok_all and eig_ok and cf_ok and pow_ok
        print(f'N={N:>2} M={M:>3}: σ²={sig2_cf:>3} 固有={eig_ok} 閉形式={cf_ok} 総パワー1={pow_ok}')
    print('\n全N(3..20) 閉形式厳密一致:', ok_all)


if __name__ == '__main__':
    main()
