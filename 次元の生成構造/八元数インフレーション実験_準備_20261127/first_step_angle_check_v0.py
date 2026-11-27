#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配列を自分に一回掛けたとき (O × O)、各行を八元数のように o = A(cosθ + n sinθ) と読んだ θ が
「一本だけ 0 度か 180 度、残りは全部 90 度」になるかどうかを、いくつかの配列と、実部と読む列の選び方について調べる。
cosθ = (O×O の i 行・実部の列の成分) / (行の長さ)。この成分は「O の i 行目」と「O の実部の列」の内積に等しい。
乱数なし。厳密な整数で計算する。使い方: python3 first_step_angle_check_v0.py
"""
import itertools
from fractions import Fraction
from balanced_init_pm1_v2 import square_array, square_array_ordered

S = [[1, 1, 1, 1, 1, 1, 1, 1], [1, -1, 1, -1, 1, -1, -1, 1], [1, -1, -1, 1, 1, 1, -1, -1], [1, 1, -1, -1, 1, -1, 1, -1],
     [1, -1, -1, -1, -1, 1, 1, 1], [1, 1, -1, 1, -1, -1, -1, 1], [1, 1, 1, -1, -1, 1, -1, -1], [1, -1, 1, 1, -1, -1, 1, -1]]


def cos_list(K, col):
    """O×O の各行について cosθ (実部を col 列目と読む)。厳密な有理数。"""
    N = len(K)
    return [Fraction(sum(K[i][t] * K[t][col] for t in range(N)), N) for i in range(N)]   # 行の長さは N


def describe(K, col):
    cs = cos_list(K, col)
    n0 = sum(1 for c in cs if abs(c) == 1)
    n90 = sum(1 for c in cs if c == 0)
    other = sorted({str(c) for c in cs if c != 0 and abs(c) != 1})
    ok = (n0 == 1 and n90 == len(K) - 1)
    return ok, "0度/180度 %d 本、90度 %d 本%s" % (n0, n90, "" if not other else "、それ以外の cosθ: " + ", ".join(other))


def report(name, K):
    N = len(K)
    res = [describe(K, c) for c in range(N)]
    good = [c + 1 for c in range(N) if res[c][0]]
    print("%s (N=%d)" % (name, N))
    print("   実部を 1 列目と読む: %s" % res[0][1])
    print("   「一本だけ 0/180 度、残り 90 度」になる実部の列: %s" % ("全部の列" if len(good) == N else (good if good else "なし")))
    bad = [c for c in range(N) if not res[c][0]]
    if bad:
        print("   ならない列の例 (%d 列目): %s" % (bad[0] + 1, res[bad[0]][1]))


report("手作りの S", S)
for N in (8, 12, 16, 24):
    H, how = square_array(N)
    report("汎用初期化 [%s]" % how, H)
K, how = square_array_ordered(8, 1)
report("汎用初期化 [%s]" % how, K)
K, how = square_array_ordered(8, "1,6,4,3,7,2,5,8")
report("汎用初期化 [%s]" % how, K)
# 列の符号を一本だけ反転した S
T = [[-v if j == 3 else v for j, v in enumerate(r)] for r in S]
report("S の 4 列目の符号を反転したもの", T)
# 行の並べ替え全部 (N=8、実部は 1 列目)
H8, _ = square_array(8)
cnt = sum(1 for order in itertools.permutations(range(8)) if describe([H8[i] for i in order], 0)[0])
print("二進の 8x8 の行の並べ替え 40320 通りのうち、実部を 1 列目と読んで「一本だけ 0/180 度、残り 90 度」になる数:", cnt)
