#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import numpy as np


def edges(N):
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)


def adjacency(N):
    ea, eb = edges(N)
    M = len(ea)
    A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        share[e] = False
        A[e, share] = 1.0
    return A


def K_of(z, A):
    return (A * np.imag(np.conj(z)[:, None] * z[None, :])).astype(np.float64, copy=False)


def one_step(z, A, den):
    K = K_of(z, A)

    # H=iK は exp(dtau*K) を full dense Hermitian eigendecomposition で
    # 直接評価するための数値表現であり、物理的な追加相互作用ではない。
    H = (1j * K).astype(np.complex128, copy=False)

    # 全固有値・全固有ベクトルを毎step計算する。部分固有分解は禁止。
    w, V = np.linalg.eigh(H)

    dtau = np.float64(2.0 * math.pi / den)
    phase = np.exp(-1j * dtau * w)

    return (V @ (phase * (V.conj().T @ z))).astype(np.complex128, copy=False)
