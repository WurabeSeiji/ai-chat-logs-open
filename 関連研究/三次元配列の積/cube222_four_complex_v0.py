#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2x2x2 の立方体に、複素数四つ z_m = a_m + i b_m をそのまま置く (木原さんの置き方、2026-11-27)。

  (1,1,1) = a1   (1,1,2) = b1        1・2 番目の番号 (i,j) が四つの複素数の番号、
  (1,2,1) = a2   (1,2,2) = b2        3 番目の番号 k が 1 なら実部 a、2 なら虚部 b。
  (2,1,1) = a3   (2,1,2) = b3
  (2,2,1) = a4   (2,2,2) = b4

三つを掛ける積 (Bhattacharya-Mesner): Prod(A,B,C)[i,j,k] = Σ_t A[i,t,k] B[i,j,t] C[t,j,k]
二つを掛ける積 (t 積, Kilmer-Martin):   (A*B)[i,j,s] = Σ_{x+y=s mod 2} Σ_t A[i,t,x] B[t,j,y]   (番号は 0 始まりで計算)

このプログラムは一回分の積を記号のまま出すだけ。繰り返しはしない。乱数・割り算・正規化なし。
使い方: python3 cube222_four_complex_v0.py
"""
import itertools
import sympy as sp

a = sp.symbols('a1:5', real=True)
b = sp.symbols('b1:5', real=True)
IDX = list(itertools.product((1, 2), repeat=3))
num = {(1, 1): 0, (1, 2): 1, (2, 1): 2, (2, 2): 3}


def place(av, bv):
    return {(i, j, k): (av[num[i, j]] if k == 1 else bv[num[i, j]]) for i, j, k in IDX}


def bm(A, B, C):
    return {(i, j, k): sp.expand(sum(A[i, t, k] * B[i, j, t] * C[t, j, k] for t in (1, 2))) for i, j, k in IDX}


def tprod(A, B):
    R = {}
    for i, j, s in IDX:
        R[i, j, s] = sp.expand(sum(A[i, t, x] * B[t, j, ((s - 1) - (x - 1)) % 2 + 1] for t in (1, 2) for x in (1, 2)))
    return R


def show(R, title):
    print(title)
    for m, (i, j) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)), 1):
        print("   複素数 %d (升目 (%d,%d,·)):  実部 = %s   虚部 = %s" % (m, i, j, sp.factor(R[i, j, 1]), sp.factor(R[i, j, 2])))


T = place(a, b)
print("置き方:")
for idx in IDX:
    print("   T%s = %s" % (idx, T[idx]))
show(bm(T, T, T), "\n[三つを掛ける積] Prod(T,T,T) の一回分 (記号のまま)")
show(tprod(T, T), "\n[二つを掛ける積] T*T の一回分 (記号のまま)")

print("\n参考: 複素数の 2x2 行列 Z = [[z1,z2],[z3,z4]] のふつうの積 Z·Z")
z = [a[m] + sp.I * b[m] for m in range(4)]
Z = sp.Matrix([[z[0], z[1]], [z[2], z[3]]])
ZZ = (Z * Z).applyfunc(sp.expand)
for m, (i, j) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1)), 1):
    print("   複素数 %d:  実部 = %s   虚部 = %s" % (m, sp.factor(sp.re(ZZ[i, j])), sp.factor(sp.im(ZZ[i, j]))))

print("\n十字を番号順に置く: (z1,z2,z3,z4) = (1, i, -1, -i)  ->  a = (1,0,-1,0), b = (0,1,0,-1)")
C = place((1, 0, -1, 0), (0, 1, 0, -1))
for idx in IDX:
    print("   T%s = %s" % (idx, C[idx]))
for name, R in (("三つを掛ける積 Prod(T,T,T)", bm(C, C, C)), ("二つを掛ける積 T*T", tprod(C, C))):
    print("  " + name)
    for m, (i, j) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)), 1):
        print("     複素数 %d: %s + (%s) i" % (m, R[i, j, 1], R[i, j, 2]))
