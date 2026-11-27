#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2x2x2 の立方体に、複素数三つ a, b, c (と対角の d) を置いたときの積を、記号のまま厳密に計算する。

升目 8 個の分類 (n=2 では三つとも違う升目は無い):
  対角      (x,x,x)            2 個  -> d
  型 1      (y,x,x), y != x    2 個  -> a   (1 番目の軸だけが違う)
  型 2      (x,y,x), y != x    2 個  -> b   (2 番目の軸だけが違う)
  型 3      (x,x,y), y != x    2 個  -> c   (3 番目の軸だけが違う)

三項の積: Prod(A,B,C)[i,j,k] = Σ_t A[i,t,k] B[i,j,t] C[t,j,k]
二項の積: (A*B)[i,j,s] = Σ_{a+b=s mod 2} Σ_t A[i,t,a] B[t,j,b]
乱数・割り算・正規化なし。使い方: python3 cube222_three_numbers_v0.py
"""
import itertools
import sympy as sp

a, b, c, d = sp.symbols('a b c d')
n = 2


def kind(i, j, k):
    if i == j == k:
        return 'd'
    if j == k:
        return 'a'
    if i == k:
        return 'b'
    return 'c'


SYM = {'d': d, 'a': a, 'b': b, 'c': c}
T = {(i, j, k): SYM[kind(i, j, k)] for i, j, k in itertools.product(range(n), repeat=3)}


def bm(A, B, C):
    return {(i, j, k): sp.expand(sum(A[i, t, k] * B[i, j, t] * C[t, j, k] for t in range(n)))
            for i, j, k in itertools.product(range(n), repeat=3)}


def tprod(A, B):
    R = {}
    for i, j, s in itertools.product(range(n), repeat=3):
        R[i, j, s] = sp.expand(sum(A[i, t, x] * B[t, j, (s - x) % n] for t in range(n) for x in range(n)))
    return R


def report(R, name):
    print(name)
    for kd in 'dabc':
        vals = {sp.factor(R[idx]) for idx in R if kind(*idx) == kd}
        print("   %s の升目 (2 個): %s%s" % (kd, " , ".join(str(v) for v in vals), "" if len(vals) == 1 else "   <- 二つの升目で値が違う"))


print("置き方:")
for idx in sorted(T):
    print("   T%s = %s" % (idx, T[idx]))
report(bm(T, T, T), "\n[三項の積] Prod(T,T,T)")
report(tprod(T, T), "\n[二項の積] T*T")

print("\n[三項の積] 具体的な値で繰り返す (厳密)")
w = sp.Rational(-1, 2) + sp.sqrt(3) / 2 * sp.I
cases = [("d=0, (a,b,c)=(1,ω,ω²)", {d: 0, a: 1, b: w, c: w**2}),
         ("d=1, (a,b,c)=(1,ω,ω²)", {d: 1, a: 1, b: w, c: w**2}),
         ("d=1, (a,b,c)=(i,-1,-i)", {d: 1, a: sp.I, b: -1, c: -sp.I}),
         ("d=0, (a,b,c)=(i,-1,-i)", {d: 0, a: sp.I, b: -1, c: -sp.I})]
for name, sub in cases:
    v = {s: sp.nsimplify(sp.simplify(sp.sympify(x))) for s, x in sub.items()}
    print("  " + name)
    for step in range(1, 5):
        dn = sp.simplify(sp.expand(v[d]**3 + v[a] * v[b] * v[c]))
        an = sp.simplify(sp.expand(v[a]**2 * v[d] + v[a] * v[b] * v[c]))
        bn = sp.simplify(sp.expand(v[b]**2 * v[d] + v[a] * v[b] * v[c]))
        cn = sp.simplify(sp.expand(v[c]**2 * v[d] + v[a] * v[b] * v[c]))
        v = {d: dn, a: an, b: bn, c: cn}
        print("     %d 回目: d=%s  a=%s  b=%s  c=%s" % (step, dn, an, bn, cn))
