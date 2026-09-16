#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""正本相互作用カーネル（密参照実装）— 忠実コピー。

出典（正本、変更禁止）:
  次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/
  干渉保存力学_資格審査とシード無し系列_20260831/
  make_parent型初期値_自己無撞着構造_20260907/
  run_N3_N40_stage123_v1_CANONICAL_REFERENCE.py
    - edges():     14-15 行目
    - adjacency(): 16-20 行目
    - H_of():      21-22 行目
    - one_step():  23-28 行目

本ファイルは上記 4 関数の bit 同一コピーであり、物理式・数値手順は一切変更しない。
Δτ = 2π/den（den は呼び出し側が指定: 主系列 den=P=L、副系列 den=40）。
"""
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


def H_of(z, A):
    H = A * (np.conj(z)[:, None] * z[None, :])
    np.fill_diagonal(H, 0.0)
    return H.astype(np.complex128, copy=False)


def one_step(z, A, den):
    # 正本コメント: 位相のみ生成子 Ĥ の虚部だけを取る H=i·K（K=sin(Δθ) 実反対称）。
    # exp(-iΔτ·iK)=exp(Δτ·K) の実直交回転となり、Z^T Z（零閉塞）と ‖Z‖ を厳密保存する。
    H = H_of(np.exp(1j * np.angle(z)), A)
    H = (1j * np.imag(H)).astype(np.complex128, copy=False)
    w, V = np.linalg.eigh(H)
    phase = np.exp(-1j * np.float64(2.0 * math.pi / den) * w)
    return (V @ (phase * (V.conj().T @ z))).astype(np.complex128, copy=False)
