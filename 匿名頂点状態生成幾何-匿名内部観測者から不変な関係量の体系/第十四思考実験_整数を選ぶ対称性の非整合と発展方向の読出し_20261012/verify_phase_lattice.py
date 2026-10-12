#!/usr/bin/env python3
"""Checks for Thought Experiment 14, Sections 3.3(c), 10.1 and 10.1.1 (2026-10-12).
1. Finite orders of all integer 2x2 matrices with det = +-1 and |entries| <= 4.
2. Orders of the Z4 map (a,b)->(-b,a) and the Z6 map (a,b)->(a-b,a).
3. psi(m) <= d (crystallographic restriction in dimension d) for d = 2, 3, 4.
Standard library only.
"""
from itertools import product

I = ((1, 0), (0, 1))

def mul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)) for i in range(2))

def order(M, limit=24):
    P = I
    for k in range(1, limit + 1):
        P = mul(P, M)
        if P == I:
            return k
    return None

orders = set()
for a, b, c, d in product(range(-4, 5), repeat=4):
    if a * d - b * c in (1, -1):
        o = order(((a, b), (c, d)))
        if o:
            orders.add(o)
assert orders == {1, 2, 3, 4, 6}
print('finite orders in GL(2,Z), |entries|<=4:', sorted(orders))

Z4 = ((0, -1), (1, 0))   # (a,b) -> (-b, a)
Z6 = ((1, -1), (1, 0))   # (a,b) -> (a-b, a)
assert order(Z4) == 4 and order(Z6) == 6 and order(mul(Z6, Z6)) == 3
print('Z4 (a,b)->(-b,a): trace', Z4[0][0] + Z4[1][1], 'order', order(Z4))
print('Z6 (a,b)->(a-b,a): trace', Z6[0][0] + Z6[1][1], 'order', order(Z6), '; squared: order', order(mul(Z6, Z6)))

def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r

def psi(m):
    if m == 1:
        return 0
    s, n, p = 0, m, 2
    while n > 1:
        if n % p == 0:
            q = 1
            while n % p == 0:
                n //= p
                q *= p
            s += phi(q)
        p += 1
    if m % 4 == 2:
        s -= 1
    return s

for d in (2, 3, 4):
    allowed = [m for m in range(1, 61) if psi(m) <= d]
    print(f'dimension {d}: allowed rotation orders', allowed)
assert [m for m in range(1, 61) if psi(m) <= 2] == [1, 2, 3, 4, 6]
assert [m for m in range(1, 61) if psi(m) <= 4] == [1, 2, 3, 4, 5, 6, 8, 10, 12]
