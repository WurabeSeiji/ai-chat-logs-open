#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三つの八元数 A, B, C の積の、順序と括弧の全部の組み合わせ (12 通り) の検算。v0。
仕様は「思考実験.md」§12。積の表は 11/19・11/20 の論文、11/22 のプログラムと同じ。
使い方: python3 triple_product_orderings_v0.py
"""
import itertools
import numpy as np

rng = np.random.default_rng(20261127)
TRI = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
T = np.zeros((8, 8, 8))
T[0, 0, 0] = 1
for i in range(1, 8):
    T[0, i, i] = 1; T[i, 0, i] = 1; T[i, i, 0] = -1
for a, b, c in TRI:
    for (i, j, k), s in [((a, b, c), 1), ((b, c, a), 1), ((c, a, b), 1), ((b, a, c), -1), ((a, c, b), -1), ((c, b, a), -1)]:
        T[i, j, k] = s
mul = lambda x, y: np.einsum('ijk,i,j->k', T, x, y)
conj = lambda x: x * np.r_[1, -np.ones(7)]
L = lambda x, y, z: mul(mul(x, y), z)      # (xy)z
R = lambda x, y, z: mul(x, mul(y, z))      # x(yz)
n = np.linalg.norm
sgn = lambda p: np.linalg.det(np.eye(3)[list(p)])

A, B, C = (rng.normal(size=8) for _ in range(3))
V = [A, B, C]
perms = list(itertools.permutations(range(3)))

print("=== (1) 12 通りの積 (6 つの順序 × 2 つの括弧) ===")
prods = {}
for p in perms:
    x, y, z = (V[i] for i in p)
    prods[p] = (L(x, y, z), R(x, y, z))
    print("  順序 %s: |(xy)z - x(yz)| = %.4f   Re の差 = %.1e" % ("".join("ABC"[i] for i in p), n(prods[p][0] - prods[p][1]), abs(prods[p][0][0] - prods[p][1][0])))
assoc = L(A, B, C) - R(A, B, C)
print("  結合子 [A,B,C] = (AB)C - A(BC): 実部 %.1e  大きさ %.4f" % (assoc[0], n(assoc)))
print("  順序を入れ替えた結合子は 符号×[A,B,C]: 最大のずれ %.1e" % max(n(prods[p][0] - prods[p][1] - sgn(p) * assoc) for p in perms))

print("\n=== (2) 12 通りを同じ重みで全部足す ===")
SL = sum(prods[p][0] for p in perms); SR = sum(prods[p][1] for p in perms)
print("  Σ(xy)z と Σx(yz) の差: %.1e   (括弧の違いは、和では消える)" % n(SL - SR))
cube = lambda X: mul(mul(X, X), X)
pol = cube(A + B + C) - cube(A + B) - cube(B + C) - cube(A + C) + cube(A) + cube(B) + cube(C)
print("  Σ(xy)z と、(A+B+C)^3 の交差項 (A,B,C を一つずつ含む項) の差: %.1e" % n(SL - pol))

print("\n=== (3) N 個を無名のまま全部足す。 P = Σ O ===")
N = 7
O = [rng.normal(size=8) for _ in range(N)]
P = sum(O)
S3 = sum(L(O[i], O[j], O[k]) for i in range(N) for j in range(N) for k in range(N))
print("  Σ_ijk (O_i O_j) O_k と P^3 の差: %.1e" % n(S3 - cube(P)))
S2 = sum(mul(O[i], O[j]) for i in range(N) for j in range(N))
print("  Σ_ij O_i O_j と P^2 の差: %.1e" % n(S2 - mul(P, P)))
Pm = P / n(P)
print("  P^3 は P の複素平面 (1 と P の虚部の向き) の中にある: 外への成分 %.1e" % n(cube(P)[1:] - (cube(P)[1:] @ Pm[1:] / (Pm[1:] @ Pm[1:])) * Pm[1:]))

print("\n=== (4) N 個に書き戻す。自分 O_i を一つ、残り二つは全員の和 P ===")
x = O[0]
cands = {"(O P)P": L(x, P, P), "O(P P)": R(x, P, P), "(P O)P": L(P, x, P), "P(O P)": R(P, x, P), "(P P)O": L(P, P, x), "P(P O)": R(P, P, x)}
print("  (O P)P と O(P P) の差 %.1e | (P O)P と P(O P) の差 %.1e | (P P)O と P(P O) の差 %.1e"
      % (n(cands["(O P)P"] - cands["O(P P)"]), n(cands["(P O)P"] - cands["P(O P)"]), n(cands["(P P)O"] - cands["P(P O)"])))
print("  残る三つは互いに違う: |O P^2 - P^2 O| = %.3f  |O P^2 - P O P| = %.3f  |P^2 O - P O P| = %.3f"
      % (n(cands["O(P P)"] - cands["(P P)O"]), n(cands["O(P P)"] - cands["(P O)P"]), n(cands["(P P)O"] - cands["(P O)P"])))
h = 1e-6
deriv = (cube(P + h * x) - cube(P - h * x)) / (2 * h)
print("  O P^2 + P O P + P^2 O と、P^3 を O の向きに微分したものの差: %.1e" % n(cands["O(P P)"] + cands["(P O)P"] + cands["(P P)O"] - deriv))
print("  Σ_i (O_i P^2 + P O_i P + P^2 O_i) と 3 P^3 の差: %.1e" % n(sum(R(o, P, P) + L(P, o, P) + L(P, P, o) for o in O) - 3 * cube(P)))

print("\n=== (5) 相手の和が二種類あるとき (重みが違う和 P1, P2) ===")
w1 = rng.uniform(size=N); w2 = rng.uniform(size=N)
P1 = sum(w * o for w, o in zip(w1, O)); P2 = sum(w * o for w, o in zip(w2, O))
print("  |(O P1)P2 - O(P1 P2)| = %.4f   (括弧の違いが残る。結合子が生きる)" % n(L(x, P1, P2) - R(x, P1, P2)))
print("  同じ和なら |(O P1)P1 - O(P1 P1)| = %.1e" % n(L(x, P1, P1) - R(x, P1, P1)))

print("\n=== (6) 三つの積を、面積 conj(P)P で一度だけ割る ===")
PP = mul(conj(P), P)
print("  conj(P)P は実数: 虚部 %.1e  値 %.4f = |P|^2 %.4f" % (n(PP[1:]), PP[0], P @ P))
for name in ("O(P P)", "(P O)P", "(P P)O"):
    y = cands[name] / PP[0]
    print("  %s / (conj(P)P): |結果|/|O| = %.12f" % (name, n(y) / n(x)))
print("  全員に施したとき、Σ|O_i'|^2 / Σ|O_i|^2 = %.12f (右)" % (sum(n(R(o, P, P) / PP[0]) ** 2 for o in O) / sum(n(o) ** 2 for o in O)))
