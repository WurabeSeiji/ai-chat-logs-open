#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
八元数の積表の符号 (8x8、成分は ±1) の 8 行を八元数 8 個として読み、実部を r に置き換えた
x_i = r + v_i (v_i は i 行目の虚数部、成分は ±1) の釣り合いを厳密に調べる。r = 1 と r = √7 を比べる。

積表はケーリー・ディクソンの構成 (a,b)(c,d) = (ac - d* b, da + b c*) で作る (480 通りある書き方の一つ)。
乱数・浮動小数・正規化なし (sympy の厳密計算)。使い方: python3 octonion_sign_table_sqrt7_v0.py
"""
import itertools
import numpy as np
import sympy as sp


def cd(c):
    n = c.shape[0]
    C = np.zeros((2 * n, 2 * n, 2 * n), int)
    conj = np.array([1] + [-1] * (n - 1))
    mul = lambda x, y: np.einsum('i,j,ijk->k', x, y, c)
    E = np.eye(n, dtype=int)
    Z = 0 * E[0]
    for i in range(2 * n):
        for j in range(2 * n):
            a, b = (E[i], Z) if i < n else (Z, E[i - n])
            cc, d = (E[j], Z) if j < n else (Z, E[j - n])
            C[i, j, :n] = mul(a, cc) - mul(conj * d, b)
            C[i, j, n:] = mul(d, a) + mul(b, conj * cc)
    return C


O = cd(cd(cd(np.ones((1, 1, 1), int))))
S = O.sum(axis=2)                      # 積表の符号の 8x8


def mul(x, y):
    return [sp.expand(sum(int(O[i, j, k]) * x[i] * y[j] for i in range(8) for j in range(8) if O[i, j, k])) for k in range(8)]


def add(xs):
    return [sp.expand(sum(x[k] for x in xs)) for k in range(8)]


def power(x, m):
    p = x
    for _ in range(m - 1):
        p = mul(p, x)                  # 同じ数の冪は括弧に依らない
    return p


def assoc_nonzero(a, b, c):
    l, r = mul(mul(a, b), c), mul(a, mul(b, c))
    return any(sp.expand(u - v) != 0 for u, v in zip(l, r))


print("積表の符号の 8x8:")
print(S)
for rname, r in (("1", sp.Integer(1)), ("√7", sp.sqrt(7))):
    X8 = [[r] + [sp.Integer(int(S[i, k])) for k in range(1, 8)] for i in range(8)]
    X16 = X8 + [[-c for c in x] for x in X8]
    print("\n=== 実部 r = %s ===" % rname)
    for name, X in (("8 個  x_i = r + v_i", X8), ("16 個 ±(r + v_i)", X16)):
        print(" [%s]" % name)
        print("   振幅の二乗:", sorted(set(sp.expand(sum(c * c for c in x)) for x in X)))
        G = sorted(set(sp.expand(sum(a * b for a, b in zip(X[i], X[j]))) for i in range(len(X)) for j in range(len(X)) if i != j), key=lambda v: float(v))
        print("   互いの内積:", G)
        for m in (1, 2, 3, 4):
            print("   Σ x^%d = %s" % (m, add([power(x, m) for x in X])))
        tr = list(itertools.combinations(range(len(X)), 3))
        print("   三つ組 %d 通りのうち括弧の差が 0 でない組: %d" % (len(tr), sum(1 for t in tr if assoc_nonzero(X[t[0]], X[t[1]], X[t[2]]))))
    x = X8[1]
    print(" 一個の冪 (x = 2 行目):")
    for m in (2, 4, 8):
        print("   x^%d = %s" % (m, power(x, m)))
