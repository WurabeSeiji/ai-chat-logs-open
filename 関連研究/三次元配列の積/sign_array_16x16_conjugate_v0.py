#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
八元数の積表の符号の 8x8 (S) と、その共役 (一列目以外の符号を反転した Sc) から 16x16 の ±1 の配列を組み、
行と列の和が全部同じになる書き方 (行・列の符号の反転) があるかを全数で調べる。
  H = [[S, Sc], [Sc, -S]]
乱数なし。使い方: python3 sign_array_16x16_conjugate_v0.py
"""
import itertools
import numpy as np

S = np.array([[1, 1, 1, 1, 1, 1, 1, 1], [1, -1, 1, -1, 1, -1, -1, 1], [1, -1, -1, 1, 1, 1, -1, -1],
              [1, 1, -1, -1, 1, -1, 1, -1], [1, -1, -1, -1, -1, 1, 1, 1], [1, 1, -1, 1, -1, -1, -1, 1],
              [1, 1, 1, -1, -1, 1, -1, -1], [1, -1, 1, 1, -1, -1, 1, -1]])
Sc = S * np.array([1, -1, -1, -1, -1, -1, -1, -1])
H = np.block([[S, Sc], [Sc, -S]])
I16 = 16 * np.eye(16, dtype=int)
print("H H^T = 16 I か:", np.array_equal(H @ H.T, I16), "  H^T H = 16 I か:", np.array_equal(H.T @ H, I16))

found = []
for s in itertools.product((1, -1), repeat=15):
    s = np.array((1,) + s)                       # 全体の符号は同じものなので先頭を + に固定
    col = H.T @ s
    if np.all(np.abs(col) == 4):
        found.append((s, col))
print("行の符号の反転 2^15 通りのうち、列の和の大きさが全部 4 になるもの:", len(found), "通り")
s, col = found[0]
T = (H * s[:, None]) * np.sign(col)[None, :]     # 行を反転し、列の和が全部 +4 になるよう列も反転
print("一つ目を採る。行の和:", T.sum(axis=1).tolist())
print("             列の和:", T.sum(axis=0).tolist())
print("T T^T = 16 I か:", np.array_equal(T @ T.T, I16), "  T^T T = 16 I か:", np.array_equal(T.T @ T, I16))
print("T は対称 (転置しても同じ) か:", np.array_equal(T, T.T))
print("T =")
for row in T:
    print("  " + " ".join("+1" if v > 0 else "-1" for v in row))
