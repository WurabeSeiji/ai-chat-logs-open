#!/usr/bin/env python3
"""Checks for Thought Experiment 14, Section 11.1 (2026-10-12).
Sum x^2 = 0 and U^n = I as consequences, not assumptions.
1. Integer + orthogonal: every integer matrix that preserves the Euclidean norm (n = 2, 3)
   has finite order (the set is a finite group of signed permutations).
2. Multiplication + norm: for a Cayley-Dickson algebra A of dimension k and integer elements p, q,
   t = |p|^2 + |q|^2, z = |p|^2 - |q|^2, w = 2 conj(p) q  (k components)
   satisfy t^2 - z^2 - |w|^2 = 4(|p|^2 |q|^2 - |conj(p) q|^2).
   This vanishes identically for k = 1, 2, 4, 8 (R, C, H, O: the norm is multiplicative),
   giving null vectors in dimensions 3, 4, 6, 10; it fails for k = 16 (sedenions).
3. Euclid: Z = a + ib, (Re Z^2)^2 + (Im Z^2)^2 = |Z|^4, i.e. (a^2-b^2, 2ab, a^2+b^2) is null.
Standard library only, exact integer arithmetic.
"""
from itertools import product
import random

# 1. integer orthogonal matrices have finite order
def matmul(A, B):
    n = len(A)
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)) for i in range(n))

for n in (2, 3):
    Id = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    group = []
    for entries in product((-1, 0, 1), repeat=n * n):
        M = tuple(tuple(entries[i * n:(i + 1) * n]) for i in range(n))
        MT = tuple(tuple(M[j][i] for j in range(n)) for i in range(n))
        if matmul(M, MT) == Id:
            group.append(M)
    orders = set()
    for M in group:
        P, k = M, 1
        while P != Id:
            P, k = matmul(P, M), k + 1
        orders.add(k)
    print(f'n={n}: {len(group)} integer orthogonal matrices, all of finite order, orders {sorted(orders)}')

# 2. Cayley-Dickson construction
def conj(x):
    if len(x) == 1:
        return x
    h = len(x) // 2
    return conj(x[:h]) + [-v for v in x[h:]]

def mul(x, y):
    if len(x) == 1:
        return [x[0] * y[0]]
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    first = [u - v for u, v in zip(mul(a, c), mul(conj(d), b))]
    second = [u + v for u, v in zip(mul(d, a), mul(b, conj(c)))]
    return first + second

def norm2(x):
    return sum(v * v for v in x)

random.seed(0)
for k in (1, 2, 4, 8, 16):
    worst = 0
    for _ in range(300):
        p = [random.randint(-9, 9) for _ in range(k)]
        q = [random.randint(-9, 9) for _ in range(k)]
        t = norm2(p) + norm2(q)
        z = norm2(p) - norm2(q)
        w = [2 * v for v in mul(conj(p), q)]
        worst = max(worst, abs(t * t - z * z - norm2(w)))
    status = 'null identically' if worst == 0 else 'NOT null'
    print(f'k={k:2d} (dimension {k + 2:2d}): max |t^2 - z^2 - |w|^2| = {worst} -> {status}')
    assert (worst == 0) == (k <= 8)

# 3. Euclid's formula
for a, b in product(range(-12, 13), repeat=2):
    assert (a * a - b * b) ** 2 + (2 * a * b) ** 2 == (a * a + b * b) ** 2
print('Euclid: (a^2-b^2)^2 + (2ab)^2 = (a^2+b^2)^2 for all |a|,|b| <= 12; e.g. (a,b)=(2,1) -> 3,4,5')
