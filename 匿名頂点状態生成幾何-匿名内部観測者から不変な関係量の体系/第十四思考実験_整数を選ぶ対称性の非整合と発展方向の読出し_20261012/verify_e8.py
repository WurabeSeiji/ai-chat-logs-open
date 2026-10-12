#!/usr/bin/env python3
"""Checks for Thought Experiment 14, Sections 10.3 and 7.6 (E8), 2026-10-12.
1. Build the 240 roots of E8 in R^8: 112 of type (+-1,+-1,0^6) and 128 of type (+-1/2)^8 with an even number of minus signs.
2. Every root has squared length 2; every inner product is an integer (the lattice is even and integral).
3. The Gram matrix of eight simple roots has determinant 1 (unimodular).
4. dim E8 = 8 + 240 = 248, Coxeter number = 240 / 8 = 30.
5. 2 cos(pi/5) equals the golden ratio.
6. (Section 7.6) The Coxeter element c of E8 satisfies Phi_30(c) = 0, so its eigenvalues are
   exp(2 pi i m / 30) with m in {1,7,11,13,17,19,23,29} (the exponents of E8).
7. (Section 7.6) Integer relations: 23 + 8 = 31 = h + 1, 124 = 4 * 31, 248 = 8 * 31,
   23 pi/124 + 8 pi/124 = pi/4, and the breakdown of the 137 four-dimensional cells by coordinate type.
Standard library only (exact arithmetic with fractions).
"""
from fractions import Fraction as F
from itertools import combinations, product
import math

roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product((1, -1), repeat=2):
        v = [F(0)] * 8
        v[i], v[j] = F(si), F(sj)
        roots.append(tuple(v))
n_int = len(roots)
for signs in product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(tuple(F(s, 2) for s in signs))
n_half = len(roots) - n_int

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))

assert (n_int, n_half, len(roots)) == (112, 128, 240)
assert all(dot(r, r) == 2 for r in roots)
assert all(dot(u, v).denominator == 1 for u in roots for v in roots)

h = F(1, 2)
e = lambda k: tuple(F(1) if t == k else F(0) for t in range(8))
sub = lambda u, v: tuple(a - b for a, b in zip(u, v))
add = lambda u, v: tuple(a + b for a, b in zip(u, v))
simple = [
    (h, -h, -h, -h, -h, -h, -h, h),
    add(e(0), e(1)),
    sub(e(1), e(0)), sub(e(2), e(1)), sub(e(3), e(2)),
    sub(e(4), e(3)), sub(e(5), e(4)), sub(e(6), e(5)),
]
assert all(s in roots for s in simple)
G = [[dot(a, b) for b in simple] for a in simple]

def det(M):
    M = [row[:] for row in M]
    n, d = len(M), F(1)
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return d

assert det(G) == 1
golden = (1 + math.sqrt(5)) / 2
assert abs(2 * math.cos(math.pi / 5) - golden) < 1e-15

print('roots: 112 (integer type) + 128 (half-integer type) =', len(roots))
print('all squared lengths 2, all inner products integer: even integral lattice')
print('det(Gram of simple roots) =', det(G), '(unimodular)')
print('dim E8 = 8 + 240 =', 8 + len(roots), '; Coxeter number = 240/8 =', len(roots) // 8)
print('2cos(pi/5) =', format(2 * math.cos(math.pi / 5), '.15f'), '= golden ratio')

# ---- Section 7.6 ----
def reflect(r):
    rr = dot(r, r)
    return [[(F(1) if i == j else F(0)) - 2 * r[i] * r[j] / rr for j in range(8)] for i in range(8)]

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(8)) for j in range(8)] for i in range(8)]

Id = [[F(1) if i == j else F(0) for j in range(8)] for i in range(8)]
C = Id
for r in simple:
    C = matmul(C, reflect(r))
P = [Id]
for k in range(8):
    P.append(matmul(P[-1], C))
# Phi_30(x) = x^8 + x^7 - x^5 - x^4 - x^3 + x + 1
coef = {8: 1, 7: 1, 5: -1, 4: -1, 3: -1, 1: 1, 0: 1}
Phi = [[sum(c * P[k][i][j] for k, c in coef.items()) for j in range(8)] for i in range(8)]
assert all(x == 0 for row in Phi for x in row)
exponents = [m for m in range(1, 30) if math.gcd(m, 30) == 1]
assert exponents == [1, 7, 11, 13, 17, 19, 23, 29]
print('Coxeter element: Phi_30(c) = 0; exponents of E8 =', exponents)

assert 23 + 8 == 31 == len(roots) // 8 + 1 and 124 == 4 * 31 and 248 == 8 * 31
assert F(23, 124) + F(8, 124) == F(1, 4)
print('23 + 8 = 31 = h + 1; 124 = 4*31; 248 = 8*31; 23pi/124 + 8pi/124 = pi/4')

from collections import Counter
cells = Counter()
for c in product(range(-3, 4), repeat=4):
    if sum((abs(v) + 0.5) ** 2 for v in c) <= 9:
        cells[tuple(sorted(map(abs, c), reverse=True))] += 1
assert sum(cells.values()) == 137
assert cells[(0, 0, 0, 0)] == 1 and cells[(1, 0, 0, 0)] == 8
outer = {k: v for k, v in cells.items() if max(k) > 0 and k != (1, 0, 0, 0)}
assert sum(outer.values()) == 128
print('137 = 1 + 8 + 128; outer 128 by type:', dict(sorted(outer.items())))
