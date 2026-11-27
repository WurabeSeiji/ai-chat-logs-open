#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
24 について、三つの条件だけで高対称な組を作って確かめる (木原さんの指示、2026-11-27)。
  重心が原点                  : どの列 (成分) も、足すと 0
  等方性とモーメントの釣り合い  : 列どうしが直交し、列の長さが同じ
  振幅が同じ                  : 行 (元) の長さが同じ
行どうしの直交、実数軸、ブロックの構造は条件にしない。成分は ±1。乱数なし。

事実: 行 24 本 × 列 24 本の正方の表では、三つを同時に満たせない。
  どの列も足して 0 なら、24 本の列はすべて「全成分が 1 の列」に直交する。
  互いに直交する列は、24 行の表の中に 24 本までしか取れないので、「全成分が 1 の列」と合わせた 25 本は入らない。
そこで二通りを作る。
  [甲] 行 24 本 × 列 23 本 (24 個の元、成分 23 個)。
  [乙] 行 48 本 × 列 24 本 (48 個の元、成分 24 個)。H の 24 行と、符号を全部反転した 24 行。
もとにする 24x24 は sign_array_24x24_v0.py と同じ作り方 (23 を法とした平方数)。
使い方: python3 balanced_set_24_v0.py
"""
p = 23
squares = {(x * x) % p for x in range(1, p)}
chi = lambda x: 0 if x % p == 0 else (1 if x % p in squares else -1)
n = p + 1
H = [[1] * n] + [[-1] + [chi((j - 1) - (i - 1)) + (1 if i == j else 0) for j in range(1, n)] for i in range(1, n)]


def check(T, name):
    R, C = len(T), len(T[0])
    col_sum = [sum(T[i][j] for i in range(R)) for j in range(C)]
    G = [[sum(T[t][i] * T[t][j] for t in range(R)) for j in range(C)] for i in range(C)]
    col_len = sorted(set(G[i][i] for i in range(C)))
    col_off = sorted(set(G[i][j] for i in range(C) for j in range(C) if i != j))
    row_len = sorted(set(sum(v * v for v in r) for r in T))
    row_off = sorted(set(sum(a * b for a, b in zip(T[i], T[j])) for i in range(R) for j in range(R) if i != j))
    print("[%s] 行 %d 本 × 列 %d 本、成分は ±1" % (name, R, C))
    print("   重心          : 列の和 = %s  -> %s" % (sorted(set(col_sum)), "満たす" if set(col_sum) == {0} else "満たさない"))
    print("   等方性・モーメント: 列の長さの二乗 = %s、違う二列の内積 = %s  -> %s"
          % (col_len, col_off, "満たす" if len(col_len) == 1 and col_off == [0] else "満たさない"))
    print("   振幅          : 行の長さの二乗 = %s  -> %s" % (row_len, "満たす" if len(row_len) == 1 else "満たさない"))
    print("   (条件ではない) 違う二行の内積 = %s" % row_off)


check(H, "参考: もとの 24x24")
# [甲] 一行目の符号を反転して一列目を全部 -1 にそろえ、その列を取り除く
A = [[-v for v in H[0]]] + [r[:] for r in H[1:]]
assert all(r[0] == -1 for r in A)
KOU = [r[1:] for r in A]
check(KOU, "甲")
# [乙] H と -H
OTSU = [r[:] for r in H] + [[-v for v in r] for r in H]
check(OTSU, "乙")
print("\n[甲] の表 (+ は +1、- は -1):")
for r in KOU:
    print("  " + " ".join("+" if v > 0 else "-" for v in r))
