#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N = 8 の二進の配列 H (balanced_init_pm1_v2.py の square_array(8)) について、行の並べ替え 8! = 40,320 通りを全部調べ、
並べ替えた配列を自分に掛ける演算がどんな回転を持つかで分類する。仕様は「思考実験.md」§32。

分類の基準: 並べ替えた配列の固有値は、どれも大きさ √8。位相が 0 度・180 度でない固有値 (回転) の角度が
  何種類あるかを数える (共役の組は一つと数える)。
    0 種類 : 回転なし (自分に掛けると単位の 8 倍)。向きは一点。
    1 種類 : 一つの角度で回る。向きは一本の楕円の上。            <- 最初の手作りの配列 S と同じ型
    2 種類以上: 角度の違う回転が混ざる。一本の楕円にならない。
固有値は浮動小数で求め、角度は 1e-9 度に丸めて数える (分類のための読出し。演算ではない)。乱数なし。
使い方: python3 survey_row_orders_N8_v0.py
"""
import itertools
import collections
import numpy as np
from balanced_init_pm1_v2 import square_array

N = 8
H, _ = square_array(N)
H = np.array(H, dtype=float)


def cycle_type(order):
    seen, lens = set(), []
    for i in range(N):
        if i in seen:
            continue
        L, j = 0, i
        while j not in seen:
            seen.add(j)
            j = order[j]
            L += 1
        lens.append(L)
    return tuple(sorted(lens, reverse=True))


by_kind = collections.Counter()
by_angle = collections.Counter()
by_cycle = collections.defaultdict(collections.Counter)
examples = {}
for order in itertools.permutations(range(N)):
    ev = np.linalg.eigvals(H[list(order)])
    ang = np.degrees(np.angle(ev))
    rot = sorted({round(abs(a), 6) for a in ang if 1e-6 < abs(a) < 180 - 1e-6})
    k = len(rot)
    by_kind[k] += 1
    if k == 1:
        mult = int(sum(1 for a in ang if abs(abs(a) - rot[0]) < 1e-6) // 2)
        key = (rot[0], mult)
        by_angle[key] += 1
        ct = cycle_type(order)
        by_cycle[key][ct] += 1
        examples.setdefault(key, order)

print("N = 8、行の並べ替え %d 通り" % sum(by_kind.values()))
print("回転の角度の種類の数ごとの、並べ替えの数:")
for k in sorted(by_kind):
    print("   %d 種類: %6d 通り%s" % (k, by_kind[k], "   <- 一本の楕円になる型" if k == 1 else ("   <- 回転なし (一点)" if k == 0 else "")))
print("\\n一種類の角度で回る並べ替えの内訳 (角度、その角度で回る平面の数):")
for (a, mult), c in sorted(by_angle.items()):
    ex = ",".join(str(v + 1) for v in examples[(a, mult)])
    cts = "  ".join("%s:%d" % ("+".join(map(str, ct)), n) for ct, n in sorted(by_cycle[(a, mult)].items(), key=lambda kv: -kv[1]))
    print("   角度 %10.6f°  平面 %d 個: %5d 通り   例 %s" % (a, mult, c, ex))
    print("        並べ替えの巡回の形 (長さの組:通り数)  " + cts)
print("\\n参考: 最初の手作りの配列 S は、角度 arccos(-1/√8) = %.6f°、平面 3 個。" % np.degrees(np.arccos(-1 / np.sqrt(8))))
